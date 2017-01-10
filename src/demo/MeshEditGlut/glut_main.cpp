//
//  glut_main.cpp
//  GEL
//
//  Created by J. Andreas Bærentzen on 04/10/13.
//
//
/*
 *  MeshEdit is a small application which allows you to load and edit a mesh.
 *  The mesh will be stored in GEL's half edge based Manifold data structure.
 *  A number of editing operations are supported. Most of these are accessible from the
 *  console that pops up when you hit 'esc'.
 *
 *  Created by J. Andreas Bærentzen on 15/08/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <string>

#include <GEL/GL/glew.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <GEL/GLGraphics/MeshEditor.h>
#include <GEL/Util/ArgExtracter.h>


using namespace std;
using namespace HMesh;
using namespace Geometry;
using namespace GLGraphics;
using namespace CGLA;
using namespace Util;

MeshEditor me;

void reshape(int W, int H)
{
    me.reshape(W,H);
}

Registry::variable<string> display_render_mode("normal");
Registry::variable<int> display_smooth_shading(true);
Registry::variable<float> display_gamma(2.2);

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
    
    me.display();
    
    glutSwapBuffers();
}

void animate()
{
    //usleep( (int)1e4 );
    me.try_spinning_ball();
    glutPostRedisplay();
}


void mouse(int button, int state, int x, int y)
{
    Vec2i pos(x,y);
    if (state==GLUT_DOWN)
    {
        if (button==GLUT_LEFT_BUTTON && glutGetModifiers() == 0)
           me.grab_ball(ROTATE_ACTION,pos);
        else if (button==GLUT_MIDDLE_BUTTON || glutGetModifiers() == GLUT_ACTIVE_CTRL)
            me.grab_ball(ZOOM_ACTION,pos);
        else if (button==GLUT_RIGHT_BUTTON || glutGetModifiers() == GLUT_ACTIVE_ALT)
            me.grab_ball(PAN_ACTION,pos);
    }
    else if (state==GLUT_UP)
        me.release_ball();
}

void motion(int x, int y) {
    Vec2i pos(x,y);
    me.roll_ball(Vec2i(x,y));
}


void keyboard_spec(int key, int x, int y)
{
    switch (key) {
        case GLUT_KEY_UP:
            me.key_up();
            break;
        case GLUT_KEY_DOWN:
            me.key_down();
            break;
        case GLUT_KEY_LEFT:
            me.key_left();
            break;
        case GLUT_KEY_RIGHT:
            me.key_right();
            break;
        case GLUT_KEY_HOME:
            me.key_home();
            break;
        case GLUT_KEY_END:
            me.key_end();
            break;
            
        default:
            break;
    }
    glutPostRedisplay();
}


void keyboard(unsigned char key, int x, int y)
{
    me.keyparse(key);
    glutPostRedisplay();
}

void init_glut(int argc, char** argv)
{
    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH|GLUT_ALPHA);
    glutInitWindowSize(800, 600);
    glutInit(&argc, argv);
    glutCreateWindow("MeshEdit");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(keyboard_spec);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutIdleFunc(animate);
}
void init_gl()
{
    glewInit();
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
	
    // Set the value of a uniform
    //glUniform2f(glGetUniformLocation(prog_P0,"WIN_SCALE"), win_size_x/2.0, win_size_y/2.0);
	
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glClearColor(1,1,1, 0.f);
    glColor4f(1.0f, 1.0f, 1.0f, 0.f);
    float material[4] = {1,1,1,1};
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material);
    glEnable(GL_DEPTH_TEST);
    
}

int main(int argc, char** argv)
{
    ArgExtracter ae(argc, argv);
	
    init_glut(argc, argv);
    init_gl();
    
    me.init();
	
    glutMainLoop();
    return 0;
}





