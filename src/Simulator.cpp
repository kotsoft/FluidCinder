//
//  Simulator.cpp
//  FluidCinder
//
//  Created by Grant Kot on 3/29/12.
//  Copyright (c) 2012 Grant Kot. All rights reserved.
//
#define numMaterials 4
#include <vector>
#include <math.h>
using namespace std;

struct Material {
	float mass, restDensity, stiffness, bulkViscosity, surfaceTension, kElastic, maxDeformation, meltRate, viscosity, damping, friction, stickiness, smoothing, gravity;
	int materialIndex;
	
	Material() : mass(1), restDensity(3), stiffness(1), bulkViscosity(1), surfaceTension(0), kElastic(0), maxDeformation(0), meltRate(0), viscosity(0), damping(0), friction(0), stickiness(0), smoothing(0), gravity(.05) {};
};

struct Particle {
	Material* mat;
	float x, y, u, v, T00, T01, T11;
	int cx, cy, gi;
	float px[3];
	float py[3];
	float gx[3];
	float gy[3];
	
	Particle(Material* mat) : mat(mat), x(0), y(0), u(0), v(0), T00(0), T01(0), T11(0), cx(0), cy(0), gi(0) {
		memset(px, 0, 12*sizeof(float));
	}
	
	Particle(Material* mat, float x, float y) : mat(mat), x(x), y(y), u(0), v(0), T00(0), T01(0), T11(0), cx(0), cy(0), gi(0) {
		memset(px, 0, 12*sizeof(float));
	}
	
	Particle(Material* mat, float x, float y, float u, float v) : mat(mat), x(x), y(y), u(u), v(v), T00(0), T01(0), T11(0), cx(0), cy(0), gi(0) {
		memset(px, 0, 12*sizeof(float));
	}
};

struct Node {
	float mass, particleDensity, gx, gy, u, v, ax, ay;
	float cgx[numMaterials];
	float cgy[numMaterials];
	bool active;
	Node() : active(false) {}
};

class Simulator {
	int gSizeX, gSizeY, gSizeY_3;
	Node* grid;
	vector<Node*> active;
	Material materials[numMaterials];
	float uscip(float p00, float x00, float y00, float p01, float x01, float y01, float p10, float x10, float y10, float p11, float x11, float y11, float u, float v)
	{
		float dx = x00 - x01;
		float dy = y00 - y10;
		float a = p01 - p00;
		float b = p11 - p10 - a;
		float c = p10 - p00;
		float d = y11 - y01;
		return ((((d - 2 * b - dy) * u - 2 * a + y00 + y01) * v +
				 ((3 * b + 2 * dy - d) * u + 3 * a - 2 * y00 - y01)) * v +
				((((2 * c - x00 - x10) * u + (3 * b + 2 * dx + x10 - x11)) * u - b - dy - dx) * u + y00)) * v +
		(((x11 - 2 * (p11 - p01 + c) + x10 + x00 + x01) * u +
		  (3 * c - 2 * x00 - x10)) * u +
		 x00) * u + p00;
	}
public:
	vector<Particle*> particles;
	Simulator() {
		materials[0].materialIndex = 0;
		materials[1].materialIndex = 1;
		materials[2].materialIndex = 2;
		materials[3].materialIndex = 3;
	}
	void initializeGrid(int sizeX, int sizeY) {
		gSizeX = sizeX;
		gSizeY = sizeY;
		gSizeY_3 = sizeY - 3;
		grid = new Node[gSizeX*gSizeY];
		for (int i = 0; i < gSizeX*gSizeY; i++) {
			grid[i] = Node();
		}
	}
	void addParticles() {
		for (int i = 0; i < 200; i++) {
			for (int j = 0; j < 1000; j++) {
				particles.push_back(new Particle(&materials[0], i*.3 +5.5, j*.3 + 5.5));
			}
		}
	}
	void update() {
		// Reset grid nodes
		int nActive = active.size();
		#pragma omp parallel for
		for (int i = 0; i < nActive; i++) {
			active[i]->active = false;
		}
		active.clear();
		
		// Add particle mass, velocity and density gradient to grid
		int nParticles = particles.size();
		for (int pi = 0; pi < nParticles; pi++) {
			Particle &p = *particles[pi];
			float m =  p.mat->mass;
			float mu = m * p.u;
			float mv = m * p.v;
			int mi = p.mat->materialIndex;
			float *px = p.px;
			float *gx = p.gx;
			float *py = p.py;
			float *gy = p.gy;
			Node* n = &grid[p.gi];
			for (int i = 0; i < 3; i++, n += gSizeY_3) {
				float pxi = px[i];
				float gxi = gx[i];
				for (int j = 0; j < 3; j++, n++) {
					float pyj = py[j];
					float gyj = gy[j];
					float phi = pxi * pyj;
					if (n->active) {
						n->mass += phi * m;
						n->particleDensity += phi;
						n->u += phi * mu;
						n->v += phi * mv;
						n->cgx[mi] += gxi * pyj;
						n->cgy[mi] += pxi * gyj;
					} else {
						n->active = true;
						active.push_back(n);
						n->mass = phi * m;
						n->particleDensity = phi;
						n->u = phi * mu;
						n->v = phi * mv;
						memset(n->cgx, 0, 2 * numMaterials * sizeof(float));
						n->cgx[mi] = gxi * pyj;
						n->cgy[mi] = pxi * gyj;
					}
				}
			}
		}
		
		nActive = active.size();
		
		// Update node velocities and density gradient
		#pragma omp parallel for
		for (int i = 0; i < nActive; i++) {
			Node& n = *active[i];
			n.ax = n.ay = 0;
			n.gx = 0;
			n.gy = 0;
			if (n.mass > 0) {
				n.u /= n.mass;
				n.v /= n.mass;
				for (int j = 0; j < numMaterials; j++) {
					n.gx += n.cgx[j];
					n.gy += n.cgy[j];
				}
				for (int j = 0; j < numMaterials; j++) {
					n.cgx[j] -= n.gx - n.cgx[j];
					n.cgy[j] -= n.gy - n.cgy[j];
				}
			}
		}
		
		// Calculate pressure and add forces to grid
		#pragma omp parallel for
		for (int pi = 0; pi < nParticles; pi++) {
			Particle& p = *particles[pi];
			
			float fx = 0, fy = 0, dudx = 0, dudy = 0, dvdx = 0, dvdy = 0, sx = 0, sy = 0;
			Node* n = &grid[p.gi];
			float *px = p.px;
			float *gx = p.gx;
			float *py = p.py;
			float *gy = p.gy;
			for (int i = 0; i < 3; i++, n += gSizeY_3) {
				float pxi = px[i];
				float gxi = gx[i];
				for (int j = 0; j < 3; j++, n++) {
					float pyj = py[j];
					float gyj = gy[j];
					float phi = pxi * pyj;
				}
			}
			
			int cx = (int)p.x;
			int cy = (int)p.y;
			int gi = cx * gSizeY + cy;
			
			Node& n1 = grid[gi];
			Node& n2 = grid[gi+1];
			Node& n3 = grid[gi+gSizeY];
			Node& n4 = grid[gi+gSizeY+1];
			float density = uscip(n1.particleDensity, n1.gx, n1.gy, n2.particleDensity, n2.gx, n2.gy, n3.particleDensity, n3.gx, n3.gy, n4.particleDensity, n4.gx, n4.gy, p.x - cx, p.y - cy);
			
			float pressure = p.mat->stiffness / p.mat->restDensity * (density - p.mat->restDensity);
			if (pressure > 2) {
				pressure = 2;
			}
			
			if (p.x < 4) {
				fx += (4 - p.x);
			} else if (p.x > gSizeX - 5) {
				fx += (gSizeX - 5 - p.x);
			}
			if (p.y < 4) {
				fy += (4 - p.y);
			} else if (p.y > gSizeY - 5) {
				fy += (gSizeY - 5 - p.y);
			}
			
			n = &grid[p.gi];
			for (int i = 0; i < 3; i++, n += gSizeY_3) {
				float pxi = px[i];
				float gxi = gx[i];
				for (int j = 0; j < 3; j++, n++) {
					float pyj = py[j];
					float gyj = gy[j];
					float phi = pxi * pyj;
					
					float gxm = gxi * pyj;
					float gym = pxi * gyj;
					n->ax += -gxm * pressure + fx*phi;
					n->ay += -gym * pressure + fy*phi;
				}
			}
		}
		
		// Update acceleration of nodes
		#pragma omp parallel for
		for (int i = 0; i < nActive; i++) {
			Node& n = *active[i];
			n.u = 0;
			n.v = 0;
			if (n.mass > 0) {
				n.ax /= n.mass;
				n.ay /= n.mass;
			}
		}
		
		#pragma omp parallel for
		for (int pi = 0; pi < nParticles; pi++) {
			Particle& p = *particles[pi];
			
			// Update particle velocities
			Node* n = &grid[p.gi];
			float *px = p.px;
			float *py = p.py; 
			for (int i = 0; i < 3; i++, n += gSizeY_3) {
				float pxi = px[i];
				for (int j = 0; j < 3; j++, n++) {
					float pyj = py[j];
					float phi = pxi * pyj;
					p.u += phi * n->ax;
					p.v += phi * n->ay;
				}
			}
			
			p.v += .03;
			
			float m =  p.mat->mass;
			float mu = m * p.u;
			float mv = m * p.v;
			
			// Add particle velocities back to the grid
			n = &grid[p.gi];
			for (int i = 0; i < 3; i++, n += gSizeY_3) {
				float pxi = px[i];
				for (int j = 0; j < 3; j++, n++) {
					float pyj = py[j];
					float phi = pxi * pyj;
					n->u += phi * mu;
					n->v += phi * mv;
				}
			}
		}
		
		// Update node velocities
		#pragma omp parallel for
		for (int i = 0; i < nActive; i++) {
			Node& n = *active[i];
			if (n.mass > 0) {
				n.u /= n.mass;
				n.v /= n.mass;
			}
		}
		
		// Advect particles
		#pragma omp parallel for
		for (int pi = 0; pi < nParticles; pi++) {
			Particle& p = *particles[pi];
			
			float gu = 0, gv = 0;
			Node* n = &grid[p.gi];
			float *px = p.px;
			float *py = p.py;
			float* gx = p.gx;
			float* gy = p.gy;
			for (int i = 0; i < 3; i++, n += gSizeY_3) {
				float pxi = px[i];
				for (int j = 0; j < 3; j++, n++) {
					float pyj = py[j];
					float phi = pxi * pyj;
					gu += phi * n->u;
					gv += phi * n->v;
				}
			}
			
			p.x += gu;
			p.y += gv;
			
			p.u += .1*(gu-p.u);
			p.v += .1*(gv-p.v);
			
			// Update grid cell index and kernel weights
			float kx = p.x;
			if (kx < 1) {
				kx = 1 + .01*rand()/RAND_MAX;
			} else if (kx > gSizeX - 2) {
				kx = gSizeX - 2 - .01*rand()/RAND_MAX;
			}
			float ky = p.y;
			if (ky < 1) {
				ky = 1 + .01*rand()/RAND_MAX;
			} else if (ky > gSizeY - 2) {
				ky = gSizeY - 2 - .01*rand()/RAND_MAX;
			}
			
			p.x = kx;
			p.y = ky;
			
			
			int cx = p.cx = (int)(kx - .5f);
			int cy = p.cy = (int)(ky - .5f);
			p.gi = cx * gSizeY + cy;
			
			float x = cx - kx;
			float y = cy - ky;
			
			// Quadratic interpolation kernel - Don't change these constants
			px[0] = .5f * x * x + 1.5f * x + 1.125f;
			gx[0] = x + 1.5f;
			x++;
			px[1] = -x * x + .75f;
			gx[1] = -2 * x;
			x++;
			px[2] = .5f * x * x - 1.5f * x + 1.125f;
			gx[2] = x - 1.5f;
			
			py[0] = .5f * y * y + 1.5f * y + 1.125f;
			gy[0] = y + 1.5f;
			y++;
			py[1] = -y * y + .75f;
			gy[1] = -2 * y;
			y++;
			py[2] = .5f * y * y - 1.5f * y + 1.125f;
			gy[2] = y - 1.5f;
		}
	}
};