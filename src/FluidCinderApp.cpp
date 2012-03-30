#include "cinder/app/AppBasic.h"
#include "cinder/gl/gl.h"
#include "Simulator.cpp"
#include <vector>

using namespace ci;
using namespace ci::app;
using namespace std;

class FluidCinderApp : public AppBasic {
	Simulator s;
  public:
	void setup();
	void mouseDown( MouseEvent event );	
	void update();
	void draw();
};

void FluidCinderApp::setup()
{
	s.initializeGrid(30,30);
	s.addParticles();
}

void FluidCinderApp::mouseDown( MouseEvent event )
{
}

void FluidCinderApp::update()
{
	s.update();
}

void FluidCinderApp::draw()
{
	// clear out the window with black
	gl::clear( Color( 0, 0, 0 ) );
	gl::color(1,1,1);
	gl::begin(GL_POINTS);
	vector<Particle*>& particles = s.particles;
	int nParticles = particles.size();
	for (int i = 0; i < nParticles; i++) {
		gl::vertex(particles[i]->x*10, particles[i]->y*10);
	}
	gl::end();
}


CINDER_APP_BASIC( FluidCinderApp, RendererGl )
