//
//  Simulator.cpp
//  FluidCinder
//
//  Created by Grant Kot on 3/29/12.
//  Copyright (c) 2012 Grant Kot. All rights reserved.
//
#define numMaterials 4
#include <vector>
using namespace std;

struct Material {
	float mass, restDensity, stiffness, bulkViscosity, surfaceTension, kElastic, maxDeformation, meltRate, viscosity, damping, friction, stickiness, smoothing, gravity;
	int materialIndex;
};

struct Particle {
	Material* mat;
	float x, y, u, v, T00, T01, T11;
	int cx, cy, gi;
	float px[3];
	float py[3];
	float gx[3];
	float gy[3];
	
	Particle() : x(0), y(0), u(0), v(0), T00(0), T01(0), T11(0), cx(0), cy(0), gi(0) {
		
	}
};

struct Node {
	float mass, particleDensity, gx, gy, u, v, ax, ay;
	float cgx[numMaterials];
	float cgy[numMaterials];
	bool active;
};

class Simulator {
	int gSizeX, gSizeY;
	Node* grid;
	vector<Node*> active;
	vector<Particle*> particles;
public:
	Simulator() {
		
	}
	void initializeGrid(int sizeX, int sizeY) {
		gSizeX = sizeX;
		gSizeY = sizeY;
		free(grid);
		grid = (Node*)malloc(gSizeX * gSizeY * sizeof(Node));
	}
	void update() {
		// Reset grid nodes
		int nActive = active.size();
		for (int i = 0; i < nActive; i++) {
			active[i]->active = false;
		}
		active.clear();
		
		// Add particle mass, velocity and density gradient to grid
		int nParticles = particles.size();
		for (int pi = 0; pi < nParticles; pi++) {
			Particle* p = particles[pi];
			float m =  p->mat->mass;
			float mu = m * p->u;
			float mv = m * p->v;
			int mi = p->mat->materialIndex;
			Node* n = &grid[p->gi];
			float *px = p->px;
			float *gx = p->gx;
			for (int i = 0; i < 3; i++, n += gSizeY) {
				float pxi = *(px++);
				float gxi = *(gx++);
				float *py = p->py;
				float *gy = p->gy;
				for (int j = 0; j < 3; j++, n++) {
					float pyj = *(py++);
					float gyj = *(gy++);
					float phi = pxi * pyj * m;
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
						memset(n->cgx, 0, numMaterials*sizeof(float));
						memset(n->cgy, 0, numMaterials*sizeof(float));
						n->cgx[mi] = gxi * pyj;
						n->cgy[mi] = pxi * gyj;
					}
				}
			}
		}
		
		for (int i = 0; i < nActive; i++) {
			Node* n = active[i];
		}
		
		for (int pi = 0; pi < nParticles; pi++) {
			Particle* p = particles[pi];
			
			// Update grid cell index and kernel weights
			int cx = p->cx = (int)(p->x - .5f);
			int cy = p->cy = (int)(p->y - .5f);
			p->gi = cx * gSizeY + cy;
			
			float x = cx - p->x;
			float y = cy - p->y;
			float* px = p->px;
			float* py = p->py;
			float* gx = p->gx;
			float* gy = p->gy;
			
			// Quadratic interpolation kernel - Don't change these constants!
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