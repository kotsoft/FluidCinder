#include "cinder/app/AppBasic.h"
#include "cinder/gl/gl.h"
#include "Simulator.cpp"
#include <vector>
#include <omp.h>

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
	void prepareSettings(Settings *settings);
};

void FluidCinderApp::setup()
{
	s.initializeGrid(400,400);
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
		gl::vertex(particles[i]->x*3, particles[i]->y*3);
	}
	gl::end();
}

void FluidCinderApp::prepareSettings( Settings *settings ) {
	settings->setWindowSize( 1200, 1200 );
    settings->setFrameRate( 60.0f );
}


CINDER_APP_BASIC( FluidCinderApp, RendererGl )
