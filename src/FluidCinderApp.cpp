#include "cinder/app/AppBasic.h"
#include "cinder/gl/gl.h"

using namespace ci;
using namespace ci::app;
using namespace std;

class FluidCinderApp : public AppBasic {
  public:
	void setup();
	void mouseDown( MouseEvent event );	
	void update();
	void draw();
};

void FluidCinderApp::setup()
{
}

void FluidCinderApp::mouseDown( MouseEvent event )
{
}

void FluidCinderApp::update()
{
}

void FluidCinderApp::draw()
{
	// clear out the window with black
	gl::clear( Color( 0, 0, 0 ) ); 
}


CINDER_APP_BASIC( FluidCinderApp, RendererGl )
