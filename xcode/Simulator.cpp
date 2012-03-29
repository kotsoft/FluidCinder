//
//  Simulator.cpp
//  FluidCinder
//
//  Created by Grant Kot on 3/29/12.
//  Copyright (c) 2012 Grant Kot. All rights reserved.
//
#include <vector>
using namespace std;

struct Material {
	float mass, restDensity, stiffness, bulkViscosity, surfaceTension, kElastic, maxDeformation, meltRate, viscosity, damping, friction, stickiness, smoothing, gravity;
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
		grid = (Node*)malloc(gSizeX*gSizeY*sizeof(Node));
	}
	void update() {
		int nActive = active.size();
		for (int i = 0; i < nActive; i++) {
			active[i]->active = false;
		}
		active.clear();
		
		int nParticles = particles.size();
		for (int pi = 0; pi < nParticles; pi++) {
			Particle* p = particles[pi];
			float m =  p->mat->mass;
			float mu = m * p->u;
			float mv = m * p->v;
			
			Node* gi = &grid[p->gi];
			for (int i = 0; i < 3; i++, gi += gSizeY) {
				for (int j = 0; j < 3; j++, gi++) {
					
				}
			}
		}
	}
};