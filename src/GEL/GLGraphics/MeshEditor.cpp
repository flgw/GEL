//
//  MeshEditor.cpp
//  GEL
//
//  Created by J. Andreas Bærentzen on 09/10/13.
//
//
#include <thread>
#if !defined (WIN32)
#include <stdarg.h>
#include <unistd.h>
#endif

#include "MeshEditor.h"
#include "consoleFunctions.h"
#include "VisObj.h"

#include <functional>
#include <vector>
#include <queue>
#include <algorithm>
#include <regex>

#include <GL/glew.h>

#include "glsl_shader.h"
#include "ShadowBuffer.h"

#include <CGLA/CGLA.h>
#include <HMesh/HMesh.h>
#include <Util/Timer.h>

using namespace std;
using namespace CGLA;
using namespace HMesh;
using namespace Util;

namespace GLGraphics {

    
    void MeshEditor::register_console_function(const std::string& name,
                                               const std::function<void(MeshEditor*, const std::vector<std::string>&)>& con_fun,
                                               const std::string& help_txt)
    {
        std::function<void (const std::vector<std::string>&)> f = std::bind(con_fun, this, placeholders::_1);
        theConsole->reg_cmdN(name, f, help_txt);
    }
    
    void MeshEditor::keyparse(unsigned short key){
        //toggle console with ESC
        if (key == 27)
        {
            console_visible = !console_visible;
            return;
        }
        
        if (console_visible)
        {
            static Timer tim;
            if(key==13)
                tim.start();
            theConsole->keyboard(key);
            if(key == 13)
            {
                active_visobj().post_create_display_list();
                double t = tim.get_secs();
                printf("%f seconds",t);
            }
            return;
        }
        else {

            switch(key) {
                case '\033':
                    console_visible = false;
                    break;
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                    active = key - '1';
                    break;
                case ' ':
                    active_visobj().clear_vertex_selection();
                    active_visobj().clear_face_selection();
                    active_visobj().clear_halfedge_selection();
                    post_create_display_list();
                    break;
                case 'f':
                    display_smooth_shading = !display_smooth_shading;
                    post_create_display_list();
                    break;
                case 'w':
                case 'n':
                case 'i':
                case 'r':
                case 't':
                case 'g':
                case 'a':
                case 'c':
                case 's':
                case 'd':
                case 'C':
                case 'M':
                case 'G':
                    display_render_mode = hotkey_to_string[key];
                    post_create_display_list();
                    break;
                case 'R':
                    theConsole->key_up();
                    theConsole->keyboard(13);
                    post_create_display_list();
                    break;
            }
            
        }
        
    }
    
    void MeshEditor::printf(const char* format, ...)
    {
        //format text
        char buffer[1024];
        va_list args;
        va_start(args, format);
        vsprintf(buffer, format, args);
        va_end(args);
        theConsole->print(buffer);
    }

    void MeshEditor::key_up(){
        if(console_visible)
            theConsole->key_up();
        else
        {
            int w[4];
            glGetIntegerv(GL_VIEWPORT, w);
            active = 0;
            active_view_control().reshape(w[2],w[3]);
        }
    }
    void MeshEditor::key_down(){
        if(console_visible)
            theConsole->key_down();
        else
        {
            int w[4];
            glGetIntegerv(GL_VIEWPORT, w);
            active = NO_MESHES-1;
            active_view_control().reshape(w[2],w[3]);
        }
    }
    void MeshEditor::key_left(){
        if(console_visible)
            theConsole->key_left();
        else
        {
            int w[4];
            glGetIntegerv(GL_VIEWPORT, w);
            active = max(0, active-1);
            active_view_control().reshape(w[2],w[3]);
        }
    }
    void MeshEditor::key_right(){
        if(console_visible)
            theConsole->key_right();
        else
        {
            int w[4];
            glGetIntegerv(GL_VIEWPORT, w);
            active = min(NO_MESHES-1, active+1);
            active_view_control().reshape(w[2],w[3]);
        }
    }
    void MeshEditor::key_home(){
        theConsole->key_home();
    }
    void MeshEditor::key_end(){
        theConsole->key_end();
    }
   
    bool MeshEditor::listen_commands() {
        static Timer tim;
        tim.start();
        if(theConsole->listen_commands()) {
            double t = tim.get_secs();
            printf("%f seconds",t);
            post_create_display_list();
            return true;
        }
        return false;
    }
    
    void MeshEditor::grab_ball(TrackBallAction action, const CGLA::Vec2i& pos){
        active_view_control().grab_ball(action, pos);
    }
    void MeshEditor::roll_ball(const CGLA::Vec2i& pos){
        active_view_control().roll_ball(pos);
    }
    void MeshEditor::release_ball(){
        active_view_control().release_ball();
    }
    bool MeshEditor::try_spinning_ball(){
        return active_view_control().try_spin();
    }
    
    void MeshEditor::save_ball() {
        ofstream ofs("trackball.bin", ios_base::binary);
        active_view_control().save(ofs);
        
    }
    void MeshEditor::load_ball() {
        ifstream ifs("trackball.bin", ios_base::binary);
        active_view_control().load(ifs);
    
    }

    
    bool MeshEditor::grab_mesh(const CGLA::Vec2i& pos)
    {
        if(depth_pick(pos[0], pos[1], depth))
        {
            dragging = true;
            mouse_x = pos[0];
            mouse_y = pos[1];
            Vec3d p0 = screen2world(mouse_x, mouse_y, depth);
            Manifold& m = active_mesh();
            active_visobj().save_old();
            Vec3d c;
            float r;
            bsphere(m, c, r);
            for(auto vid : m.vertices())
            {
                double l = sqr_length(p0-m.pos(vid));
                weight_vector[vid] = exp(-l/(brush_size*r*r));
            }
            return true;
        }
        return false;
    }
    
    bool MeshEditor::drag_mesh(const CGLA::Vec2i& pos)
    {
        auto deform_mesh = [&](Manifold& m, VertexAttributeVector<float>& wv, Vec3d& v)
        {
            const Manifold& m_old = active_visobj().mesh_old();
            VertexAttributeVector<Vec3d> new_pos;
            for(auto vid : m_old.vertices())
                new_pos[vid] = m_old.pos(vid) + weight_vector[vid] * v;
            m.positions_attribute_vector() = new_pos;
        };
        
        
        auto inflate_mesh = [&](Manifold& m, VertexAttributeVector<float>& wv, const Vec3d& p0)
        {
            Vec3d c;
            float r;
            bsphere(m, c, r);
            VertexAttributeVector<Vec3d> new_pos;
            for(auto vid : m.vertices()) {
                Vec3d p = m.pos(vid);
                double l = sqr_length(p0-p);
                double wgt = exp(-l/(brush_size*r*r));
                new_pos[vid] = m.pos(vid) +
                wgt * (0.25 * laplacian(m, vid) + (r*0.0025)*normal(m, vid));
            }
            m.positions_attribute_vector() = new_pos;
        };
        
        
        auto smooth_mesh = [&](Manifold& m, VertexAttributeVector<float>& wv, const Vec3d& p0)
        {
            Vec3d c;
            float r;
            bsphere(m, c, r);
            VertexAttributeVector<Vec3d> new_pos;
            for(auto vid : m.vertices()) {
                Vec3d p = m.pos(vid);
                double l = sqr_length(p0-p);
                double wgt = exp(-l/(brush_size*r*r));
                new_pos[vid] = m.pos(vid) +
                wgt * 0.5 * laplacian(m, vid);
            }
            m.positions_attribute_vector() = new_pos;
        };

        auto paint_mesh = [&](Manifold& m, VertexAttributeVector<float>& wv, const Vec3d& p0)
        {
            Vec3d c;
            float r;
            bsphere(m, c, r);
            auto& col_map = active_visobj().get_color_field_attrib_vector();
            
            /// This is inelegant, but we need to know if the damn thing is initialized.
            if(std::isnan(col_map[*m.vertices().begin()][0])) {
                cout << "col_map.size " << col_map.size() << endl;
                for(auto vid: m.vertices())
                    col_map[vid] = Vec3d(0);
            }
            
            double support_radius = r * brush_size;
            for(auto vid : m.vertices()) {
                Vec3d p = m.pos(vid);
                double l = length(p0-p);
                if(l<support_radius) {
                    double t = l/support_radius;
                    double wgt = 1.0 - (3*t*t - 2 * t*t*t);
                    col_map[vid] += 0.1*wgt*Vec3d(paint_color);
                }
            }
        };

        
        if(dragging)
        {
            Vec3d p0 = screen2world(mouse_x, mouse_y, depth);
            Vec3d p1 = screen2world(pos[0], pos[1], depth);
            Vec3d v = p1-p0;
            Manifold& m = active_mesh();
            auto vset = active_visobj().get_vertex_selection();
            auto hset = active_visobj().get_halfedge_selection();
            auto fset = active_visobj().get_face_selection();
            
            if(!vset.empty())
                for(auto vid : vset)
                    m.pos(vid) = active_visobj().mesh_old().pos(vid) + v;
            else if(!hset.empty())
                for(auto h : hset) {
                    Walker w = m.walker(h);
                    auto vid0 = w.vertex();
                    auto vid1 = w.opp().vertex();
                    m.pos(vid0) = active_visobj().mesh_old().pos(vid0) + v;
                    m.pos(vid1) = active_visobj().mesh_old().pos(vid1) + v;
                }
            else if(!fset.empty())
                for(auto f: fset)
                {
					circulate_face_ccw(m, f, std::function<void(VertexID)>( [&](VertexID vid){
                        m.pos(vid) = active_visobj().mesh_old().pos(vid) + v;
                    }) );
                }
            else {
                if(string(brush_type) == "smooth")
                    smooth_mesh(m, weight_vector, p1);
                else if(string(brush_type) == "inflate")
                    inflate_mesh(m, weight_vector, p1);
                else if(string(brush_type) == "deform")
                    deform_mesh(m, weight_vector, v);
                else if(string(brush_type) == "paint") {
                    float d;
                    if(depth_pick(pos[0], pos[1], d))
                        paint_mesh(m, weight_vector, screen2world(pos[0], pos[1], d));
                }
            }
            
            post_create_display_list();
            return true;
        }
        return false;
    }
    
    void MeshEditor::release_mesh()
    {
        dragging = false;
    }
    
    
    
    
//    mutex parallel_work;
    void MeshEditor::display(int scale){

        // Clear screen.
        glClearColor(1, 1, 1, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        //        parallel_work.lock();
        //        active_visobj().post_create_display_list();
        
        
        // Display object.
        active_visobj().display(display_render_mode, *theConsole, display_smooth_shading, display_gamma);
        
        // Draw console
        if(console_visible)
        {
            glUseProgram(0);
            theConsole->display(scale);
        }
        SocketConnector::execute_command(theConsole.get());
        
        // Static variable controlling whether we render shadow at all.
        static Registry::variable<int> shadow_enable(0);
        shadow_enable.reg(*theConsole, "display.shadow.enable", "");
        if(shadow_enable) {
            // Static variables that control the shadow display created and linked to console
            static Registry::variable<float> zenith(1.571);
            static Registry::variable<float> azimuth(0);
            static Registry::variable<float> shadow_alpha(0.3);
            azimuth.reg(*theConsole, "display.shadow.azimuth", "");
            zenith.reg(*theConsole, "display.shadow.zenith", "");
            shadow_alpha.reg(*theConsole, "display.shadow.shadow_alpha", "");

            
            // Shadow buffer (really a frame buffer object)
            static ShadowBuffer sb;
            
            // String containing fragment program for shadow rendering
            static string shadow_shdr_fp = "#version 120\n"
            "    uniform sampler2DShadow shadow_map;\n"
            "    uniform mat4 Mat;\n"
            "   uniform float shadow_alpha;\n"
            "    varying vec4 ep;\n"
            "    void main()\n"
            "    {\n"
            "        vec4 light_pos =  Mat * ep;\n"
            "        light_pos.z = max(0.001,min(0.999,light_pos.z-0.003));\n"
            "        if(light_pos.x <=0 || light_pos.x >=1 || light_pos.y <=0 || light_pos.y>=1) gl_FragColor = vec4(1);"
            "   else      gl_FragColor= vec4(1.0-shadow_alpha)+shadow_alpha*0.25*(shadow2D(shadow_map, light_pos.xyz+vec3(0.001,-0.001,0))+shadow2D(shadow_map, light_pos.xyz+vec3(0.001,0.001,0))+shadow2D(shadow_map, light_pos.xyz-vec3(0.001,0.001,0))+shadow2D(shadow_map, light_pos.xyz+vec3(-0.001,0.001,0)));\n"
            "    }\n";
            
            // Shader program for shadow rendering is compiled and linked.
            static GLuint prog = 0;
            if(!prog)
            {
                GLuint vp = create_glsl_shader(GL_VERTEX_SHADER,"#version 120\n varying vec4 ep; void main(){ep = gl_Vertex; gl_Position = ftransform();}");
                GLuint fp = create_glsl_shader(GL_FRAGMENT_SHADER, shadow_shdr_fp);
                prog = glCreateProgram();
                
                glAttachShader(prog, vp);
                glAttachShader(prog, fp);
                glLinkProgram(prog);
            }
            
            
            // Setup OpenGL state for lighting - used when rendering 3D object casting shadow
            Vec4f lpos = Vec4f(cos(azimuth)*cos(zenith),sin(zenith),sin(azimuth)*cos(zenith),0);
            glLightfv(GL_LIGHT0, GL_POSITION, lpos.get());
            Vec4f mamb(.8,.8,.8,0);
            Vec4f lamb(.4,.4,.5,0);
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mamb.get());
            glLightfv(GL_LIGHT0, GL_AMBIENT, lamb.get());

            // Get old viewport dimensions. Setup rendering to FBO and set viewport
            // dimensions for rendering shadows.
            GLint viewp[4];
            glGetIntegerv(GL_VIEWPORT, viewp);
            sb.enable();
            glViewport(0, 0, 1024, 1024);
            
            // Setup object transformations using old school GL.
            Mat4x4f m;
            float r = active_visobj().get_bsphere_radius();
            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glLoadIdentity();
            glMatrixMode(GL_PROJECTION);
            glPushMatrix();
            glLoadIdentity();
            glOrtho(-2*r, 2*r, -2*r, 2*r, -2*r, 2*r);
            gluLookAt(0.1*r*cos(azimuth)*cos(zenith),0.1*r*sin(zenith),0.1*r*sin(azimuth)*cos(zenith), 0,0,0, 0,1,0);
            
            // Copy the transformation matrix to user code.
            glGetFloatv(GL_PROJECTION_MATRIX, m.get());
            
            // Draw the object in light space.
            draw(active_visobj().mesh());
            
            // Restore transformation matrices.
            glPopMatrix();
            glMatrixMode(GL_MODELVIEW);
            glPopMatrix();
            
            // Restore the usual on screen framebuffer.
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
            glDrawBuffer(GL_BACK);
            glViewport(0, 0, viewp[2], viewp[3]);
            
            // Use shadow: bind the shadow texture to unit 0.
            sb.bind_textures(0);
            
            // Use the shadow rendering shader program.
            glUseProgram(prog);
            active_visobj().view_control().set_gl_modelview();
            glUniform1i(glGetUniformLocation(prog, "shadow_map"), 0);
            glUniform1f(glGetUniformLocation(prog, "shadow_alpha"), shadow_alpha);
            m = translation_Mat4x4f(Vec3f(0.5)) * scaling_Mat4x4f(Vec3f(0.5)) * transpose(m);
            glUniformMatrix4fv(glGetUniformLocation(prog, "Mat"), 1, 1, m.get());
            
            
            // Setup blending such that the shadow is alpha blended with model.
            glDepthFunc(GL_LEQUAL);
            glEnable(GL_BLEND);
            glBlendFunc(GL_ZERO, GL_SRC_COLOR);

            // Draw ground plane for shadows.
            Vec3d p0, p7;
            bbox(active_visobj().mesh(), p0, p7);
            glBegin(GL_QUADS);
            glVertex3f(-100*r, p0[1],-100*r);
            glVertex3f(-100*r, p0[1],100*r);
            glVertex3f(100*r, p0[1],100*r);
            glVertex3f(100*r, p0[1],-100*r);
            glEnd();
            
            // Draw model again ... just to add shadow.
            draw(active_visobj().mesh());
            
            // Disable blending and shader program.
            glDisable(GL_BLEND);
            glUseProgram(0);
        }
        //        parallel_work.unlock();
        //        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    
    void MeshEditor::reshape(int w, int h) {
        for(VisObj& v : vo)
            v.view_control().reshape(w, h);
    }
    
    Vec2i MeshEditor::shape() {
        return vo[active].view_control().shape();
    }

    void MeshEditor::init() {
        glewInit();
        
        GLint vp[4];
        glGetIntegerv(GL_VIEWPORT, vp);
        for(VisObj& vis_obj : vo)
            vis_obj.view_control().reshape(vp[2], vp[3]);
        
//        glEnable(GL_CULL_FACE);
//        glCullFace(GL_BACK);
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
        
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glClearColor(1,1,0, 0.f);
        glColor4f(1.0f, 1.0f, 1.0f, 0.f);
        float material[4] = {1,1,1,1};
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material);
        glEnable(GL_DEPTH_TEST);
        
        register_console_function("quit", console_quit,"");
        register_console_function("exit", console_quit,"");
        register_console_function("bye", console_quit,"");
        
        register_console_function("merge_with", console_merge, "");
        register_console_function("clear_mesh", console_clear, "");

        register_console_function("help.controls", console_list_controls, "");
        register_console_function("simplify", console_simplify,"");
        register_console_function("ridge_lines", console_ridge_lines,"");
        
        register_console_function("smooth.laplacian", console_laplacian_smooth,"");
        register_console_function("smooth.taubin", console_taubin_smooth,"");
        register_console_function("smooth.fuzzy_vector_median_anisotropic", console_fvm_anisotropic_smooth ,"");
        register_console_function("smooth.bilateral_anisotropic", console_bilateral_anisotropic_smooth ,"");
        
        register_console_function("optimize.valency", console_optimize_valency,"");
        register_console_function("optimize.minimize_dihedral_angles", console_minimize_dihedral,"");
        register_console_function("optimize.minimize_curvature", console_minimize_curvature,"");
        register_console_function("optimize.maximize_min_angle", console_maximize_min_angle,"");
        register_console_function("load_mesh", console_reload,"");
        register_console_function("add_mesh", console_add_mesh,"");
        
        register_console_function("cleanup.stitch", console_stitch,"");
        register_console_function("cleanup.remove_duplicates", console_remove_duplicates,"");
        register_console_function("cleanup.remove_val2", console_remove_val2, "");
        register_console_function("cleanup.flip_orientation", console_flip_orientation,"");
        register_console_function("cleanup.remove_caps", console_remove_caps,"");
        register_console_function("cleanup.remove_needles", console_remove_needles,"");
        register_console_function("cleanup.close_holes", console_close_holes,"");
        register_console_function("cleanup.compact", console_compact,"");

        register_console_function("triangulate", console_triangulate,"");
        register_console_function("dual", console_dual,"");
        register_console_function("refine.split_edges", console_refine_edges,"");
        register_console_function("refine.split_faces", console_refine_faces,"");
        
        register_console_function("subdivide.catmull_clark", console_cc_subdivide,"");
        register_console_function("subdivide.rootcc", console_root_cc_subdivide,"");
        register_console_function("subdivide.loop", console_loop_subdivide,"");
        register_console_function("subdivide.root3", console_root3_subdivide,"");
        register_console_function("subdivide.doo_sabin", console_doosabin_subdivide,"");
        register_console_function("subdivide.butterfly", console_butterfly_subdivide,"");
        register_console_function("save_mesh", console_save,"");
        register_console_function("noise.perturb_vertices", console_vertex_noise,"");
        register_console_function("noise.perturb_vertices_perpendicular", console_perpendicular_vertex_noise,"");
        register_console_function("noise.perturb_topology", console_noisy_flips,"");
        
        
        register_console_function("align_with", console_align,"");
        register_console_function("undo", console_undo,"");
        
        register_console_function("validity", console_valid,"");
        register_console_function("info", console_info,"");
        register_console_function("info.all_meshes", console_info_all,"");
        
                
        register_console_function("display.refit_trackball", console_refit_trackball, "Resets trackball");
        register_console_function("display.save_trackball", console_save_trackball, "Saves trackball to disk");
        register_console_function("display.load_trackball", console_load_trackball, "Load trackball to disk");
        
        register_console_function("transform.scale", console_scale, "Scale mesh");
        register_console_function("transform.rotate", console_rotate, "Rotate mesh");
        register_console_function("transform.translate", console_translate, "Translate mesh");
        
        register_console_function("edit.selected.merge_1_ring", console_merge_1_ring, "Merge 1-ring of selected vertices");
        
        register_console_function("edit.selected.split_edge", console_split_edge , "");
        register_console_function("edit.selected.collapse_edge", console_collapse_edge, "");
        register_console_function("edit.selected.flip_edge", console_flip_edge, "");
        register_console_function("edit.selected.dissolve_edge", console_dissolve_edge, "");
        register_console_function("edit.selected.stellate_face", console_stellate_face, "");
        register_console_function("edit.selected.split_face", console_split_face, "");
        register_console_function("edit.selected.delete", console_delete, "");
        register_console_function("edit.selected.triangulate_face", console_triangulate_face, "");
        register_console_function("edit.selected.bridge_faces", console_bridge_faces, "");

        register_console_function("test", console_test, "Test some shit");
        
        selection_mode.reg(theConsole, "selection.mode", "The selection mode. 0 = vertex, 1 = halfedge, 2 = face");
        active.reg(theConsole, "active_mesh", "The active mesh");
        display_render_mode.reg(theConsole, "display.render_mode", "Display render mode");
        brush_size.reg(theConsole, "brush.size", "Size of brush used for editing");
        brush_type.reg(theConsole, "brush.type", "Smooth, deform, inflate, or paint");
        paint_color.reg(theConsole, "brush.color", "The color of the brush used in paint mode");
        
        display_smooth_shading.reg(theConsole, "display.smooth_shading", "1 for smooth shading 0 for flat");
        display_gamma.reg(theConsole, "display.gamma", "The gamma setting for the display");

        theConsole->print("Welcome to MeshEdit");
        theConsole->newline();
    }
    
    bool MeshEditor::add_file(const std::string& str)
    {   
        while (active_mesh().no_vertices()>0 && active<NO_MESHES)
            active  = active + 1;
        if(active == NO_MESHES)
            active = 0;
        if(active_visobj().reload(str)) {
            active_visobj().post_create_display_list();
            return true;
        }
        return false;
    }

    bool MeshEditor::reload_active_from_file(const std::string& str)
    {
        if(active_visobj().reload(str)) {
            active_visobj().post_create_display_list();
            return true;
        }
        return false;
    }
    
    bool MeshEditor::add_to_active_from_file(const std::string& str)
    {
        if(active_visobj().add_mesh(str)) {
            active_visobj().post_create_display_list();
            return  true;
        }
        return false;
    }

    
}
