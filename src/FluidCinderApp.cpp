#include "cinder/app/AppBasic.h"
#include "cinder/gl/gl.h"
#include "Simulator.cpp"
#include <vector>

#include "cinder/gl/Vbo.h"
//#include <omp.h>

using namespace ci;
using namespace ci::app;
using namespace std;

class FluidCinderApp : public AppBasic {
	Simulator s;
	int n;
	GLfloat* vertices;
  public:
	void setup();
	void mouseDown( MouseEvent event );	
	void update();
	void draw();
	void prepareSettings(Settings *settings);
};

void FluidCinderApp::setup()
{
	s.initializeGrid(800,400);
	s.addParticles();
	n = s.particles.size();
	gl::VboMesh::Layout layout;
	layout.setDynamicPositions();
	
	vertices = new GLfloat[n*4];
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
	gl::color(.1,.5,1);
	vector<Particle>& particles = s.particles;
	int nParticles = particles.size();
	float* vi = vertices;
	for (int i = 0; i < nParticles; i++) {
		Particle& p = particles[i];
		*(vi++) = p.x*2;
		*(vi++) = p.y*2;
		*(vi++) = (p.x-p.gu)*2;
		*(vi++) = (p.y-p.gv)*2;
	}
	
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(2, GL_FLOAT, 0, vertices);
	
	glDrawArrays(GL_LINES, 0, n*2);
	
	glDisableClientState(GL_VERTEX_ARRAY);
}

void FluidCinderApp::prepareSettings( Settings *settings ) {
	settings->setWindowSize( 1600, 800 );
}


CINDER_APP_BASIC( FluidCinderApp, RendererGl )
