//
// Created by flgw on 1/11/17.
//

#ifndef GEL_CONSOLEFUNCTIONS_H
#define GEL_CONSOLEFUNCTIONS_H

#include <CGLA/CGLA.h>
#include <HMesh/HMesh.h>
#include <Util/Timer.h>

namespace GLGraphics {

    bool wantshelp(const std::vector<std::string> &args) {
        if (args.size() == 0)
            return false;

        std::string str = args[0];

        if (str == "help" || str == "HELP" || str == "Help" || str == "?")
            return true;

        return false;
    }

    static std::map<unsigned short, std::string> hotkey_to_string = {
            {'w', "wireframe render mode"},
            {'n', "normal render mode"},
            {'i', "isophotes render mode"},
            {'r', "reflection render mode"},
            {'t', "toon render mode"},
            {'g', "glazed render mode"},
            {'a', "ambient occlusion render mode"},
            {'c', "color field render mode"},
            {'s', "scalar field render mode"},
            {'d', "debug render mode"},
            {'C', "curvature lines render mode"},
            {'M', "mean curvature render mode"},
            {'G', "gaussian curvature render mode"},
//            {27,  "toggle console"},
            {'1', "switch to mesh 1"},
            {'2', "switch to mesh 2"},
            {'3', "switch to mesh 3"},
            {'4', "switch to mesh 4"},
            {'5', "switch to mesh 5"},
            {'6', "switch to mesh 6"},
            {'7', "switch to mesh 7"},
            {'8', "switch to mesh 8"},
            {'9', "switch to mesh 9"},
            {' ', "clear selections"},
            {'f', "toggle smooth/flat shading"},
            {'R', "repeat last command"}
    };

    void console_list_controls(MeshEditor *me, const std::vector<std::string> &args) {
        me->printf("");
        me->printf("== MeshEdit Controls ==");
        me->printf("");
        me->printf("-- Mouse movement --");
        me->printf("Left button down                           : rotate");
        me->printf("Left button down + SHIFT                   : edit");
        me->printf("Middle button down (or right button + ALT) : zoom");
        me->printf("Right button down                          : pan");
        me->printf("-- Mouse button press --");
        me->printf("Left mouse click + CTRL                 : select");
        me->printf("-- Hotkeys (note some are upper case) --");
        me->printf("  ESC : toggles console");
        me->printf("  <-  : switch to mesh with lower number ");
        me->printf("  ->  : switch to mesh with higher number ");
        for (auto &h: hotkey_to_string) {
            me->printf("  '%c' : %s", h.first, h.second.c_str());
        }
    }

    /// Function that aligns two meshes.
    void console_align(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: align_with <src>");
            me->printf("This function aligns current mesh with src");
            me->printf("In practice the GLViewController of src is copied to current mesh.");
            me->printf("The argument is mandatory and must correspond to a valide mesh entry.");
            me->printf("Note that results might be unexpected if the meshes are not on the same scale");
        }

        int dest = 0;

        if (args.size() > 0) {
            int src = console_arg(args, 0, 1);
            if (src < 0 || src >= me->get_no_meshes()) {
                me->printf("src mesh out of range");
                return;
            }
            me->align(src, me->get_active_no());
        } else
            me->printf("You must enter the mesh number that you want to align with");
    }

    void console_merge(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: merge_with <src>");
            me->printf("merges src into current mesh");
            return;
        }
        int src = console_arg(args, 0, 2);
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        HMesh::Manifold &m_src = me->get_mesh(src);
        m.merge(m_src);
        return;
    }

    void console_clear(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: clear");
            me->printf("clears current mesh");
            return;
        }
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        m.clear();
        return;
    }

    void console_quit(MeshEditor *me, const std::vector<std::string> &args) {
        exit(0);
    }

    void console_ridge_lines(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: ridge_lines");
            return;
        }

        me->save_active_mesh();

        HMesh::Manifold &mani = me->active_mesh();

        HMesh::VertexAttributeVector<CGLA::Mat3x3d> curvature_tensors;
        HMesh::VertexAttributeVector<CGLA::Vec3d> min_curv_direction;
        HMesh::VertexAttributeVector<CGLA::Vec3d> max_curv_direction;
        HMesh::VertexAttributeVector<CGLA::Vec2d> curvature;

        curvature_paraboloids(mani,
                              min_curv_direction,
                              max_curv_direction,
                              curvature);

        for (auto vid : mani.vertices()) {
            CGLA::Vec3d max_curv_dir = normalize(max_curv_direction[vid]);
            CGLA::Vec3d min_curv_dir = normalize(min_curv_direction[vid]);
            double vid_min_pc = curvature[vid][0];
            double vid_max_pc = curvature[vid][1];
            bool ridge = true;
            bool ravine = true;
            HMesh::Walker w = mani.walker(vid);
            CGLA::Vec3d r(0);
            for (; !w.full_circle(); w = w.circulate_vertex_ccw()) {
                CGLA::Vec3d e = (mani.pos(w.vertex()) - mani.pos(vid));

                if (abs(dot(min_curv_dir, e)) > abs(dot(max_curv_dir, e))) {
                    if (curvature[w.vertex()][0] < vid_min_pc + 20)
                        ravine = false;

                } else {
                    if (curvature[w.vertex()][1] > vid_max_pc - 20)
                        ridge = false;
                }
            }
            DebugRenderer::vertex_colors[vid] = CGLA::Vec3f(ridge, ravine, 0.0);
        }
        for (auto fid : mani.faces())
            DebugRenderer::face_colors[fid] = CGLA::Vec3f(.3, .3, .6);
        for (auto hid : mani.halfedges()) {

            HMesh::Walker w = mani.walker(hid);
            CGLA::Vec3f c0 = DebugRenderer::vertex_colors[w.opp().vertex()];
            CGLA::Vec3f c1 = DebugRenderer::vertex_colors[w.vertex()];

            DebugRenderer::edge_colors[hid] = (c0 == c1) ? c0 : CGLA::Vec3f(0.1, 0.1, 0.3);

        }
    }

    void transform_mesh(HMesh::Manifold &mani, const CGLA::Mat4x4d &m) {
        for (HMesh::VertexIDIterator vid = mani.vertices_begin(); vid != mani.vertices_end(); ++vid)
            mani.pos(*vid) = m.mul_3D_point(mani.pos(*vid));
    }

    void console_scale(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: transform.scale sx sy sz");
            me->printf("Note: If only sx is provided, uniform scaling is applied");
            return;
        }

        CGLA::Vec3d s;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            double scale;
            a0 >> scale;
            s = CGLA::Vec3d(scale);
        }
        if (args.size() > 1) {
            std::istringstream a0(args[1]);
            a0 >> s[1];
        }
        if (args.size() > 2) {
            std::istringstream a0(args[2]);
            a0 >> s[2];
        }

        me->save_active_mesh();
        transform_mesh(me->active_mesh(), scaling_Mat4x4d(s));
    }

    void console_rotate(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: transform.rotate axis_x axis_y axis_z angle");
            return;
        }

        CGLA::Vec3d a;
        double angle;

        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> a[0];
        }
        if (args.size() > 1) {
            std::istringstream a0(args[1]);
            a0 >> a[1];
        }
        if (args.size() > 2) {
            std::istringstream a0(args[2]);
            a0 >> a[2];
        }
        if (args.size() > 3) {
            std::istringstream a0(args[3]);
            a0 >> angle;
        }

        me->save_active_mesh();
        CGLA::Quatd q;
        q.make_rot(angle, a);
        transform_mesh(me->active_mesh(), q.get_Mat4x4d());
    }


    void console_translate(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: transform.translate tx ty tz");
            me->printf("Note: recenters if no arguments are provided.");
            return;
        }

        CGLA::Vec3d t;

        if (args.size() == 0) {
            float rad;
            bsphere(me->active_mesh(), t, rad);
        } else {
            if (args.size() > 0) {
                std::istringstream a0(args[0]);
                a0 >> t[0];
            }
            if (args.size() > 1) {
                std::istringstream a0(args[1]);
                a0 >> t[1];
            }
            if (args.size() > 2) {
                std::istringstream a0(args[2]);
                a0 >> t[2];
            }
        }

        me->save_active_mesh();
        transform_mesh(me->active_mesh(), translation_Mat4x4d(-t));
    }


    void console_merge_1_ring(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: edit.selected.merge_1_ring");
            return;
        }
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        auto sel = me->get_vertex_selection();
        for (auto v: sel)
            if (m.in_use(v))
                m.merge_one_ring(v);
    }


    void console_flip_edge(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: edit.selected.flip_edge");
            return;
        }
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        auto sel = me->get_halfedge_selection();
        for (auto h: sel)
            if (m.in_use(h) && precond_flip_edge(m, h))
                m.flip_edge(h);
    }

    void console_collapse_edge(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: edit.selected.collapse_edge");
            return;
        }
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        auto sel = me->get_halfedge_selection();
        for (auto h: sel)
            if (m.in_use(h) && precond_collapse_edge(m, h))
                m.collapse_edge(h, true);
    }

    void console_dissolve_edge(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: edit.selected.dissolve_edge");
            return;
        }
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        auto sel = me->get_halfedge_selection();
        for (auto h: sel)
            if (m.in_use(h))
                m.merge_faces(m.walker(h).face(), h);
    }

    void console_split_edge(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: edit.selected.split_edge");
            return;
        }
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        auto sel = me->get_halfedge_selection();
        for (auto h: sel)
            if (m.in_use(h))
                m.split_edge(h);
    }

    void console_stellate_face(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: edit.selected.stellate_face");
            return;
        }
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        auto sel = me->get_face_selection();
        for (auto f: sel)
            if (m.in_use(f))
                m.split_face_by_vertex(f);
    }


    void console_triangulate_face(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: edit.selected.triangulate_face");
            return;
        }
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        auto sel = me->get_face_selection();
        for (auto f: sel)
            if (m.in_use(f))
                triangulate_face_by_edge_split(m, f);
    }

    void console_bridge_faces(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: edit.selected.bridge_faces");
            return;
        }
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        auto sel = me->get_face_selection();
        if (sel.size() != 2) {
            me->printf("You must select exactly two faces");
            return;
        }

        int i = 0;
        HMesh::FaceID f0, f1;
        for (auto f: sel) {
            int n = no_edges(m, f);
            if (i == 0) {
                f0 = f;
                i += n;
            } else {
                i -= n;
                f1 = f;
            }
        }
        if (i != 0) {
            me->printf("The selected faces must have same number of edges");
            return;
        }
        std::vector<HMesh::VertexID> loop0;
        circulate_face_ccw(m, f0, std::function<void(HMesh::VertexID)>([&](HMesh::VertexID v) {
            loop0.push_back(v);
        }));

        std::vector<HMesh::VertexID> loop1;
        circulate_face_ccw(m, f1, std::function<void(HMesh::VertexID)>([&](HMesh::VertexID v) {
            loop1.push_back(v);
        }));

        std::vector<std::pair<HMesh::VertexID, HMesh::VertexID> > connections;

        size_t L0 = loop0.size();
        size_t L1 = loop1.size();

        assert(L0 == L1);

        size_t L = L0;

        float min_len = FLT_MAX;
        int j_off_min_len = -1;
        for (int j_off = 0; j_off < L; ++j_off) {
            float len = 0;
            for (int i = 0; i < L; ++i)
                len += sqr_length(m.pos(loop0[i]) - m.pos(loop1[(L + j_off - i) % L]));
            if (len < min_len) {
                j_off_min_len = j_off;
                min_len = len;
            }
        }
        for (int i = 0; i < L; ++i)
            connections.push_back(
                    std::pair<HMesh::VertexID, HMesh::VertexID>(loop0[i], loop1[(L + j_off_min_len - i) % L]));
        // Merge the two one rings producing two faces.

        // Bridge the just created faces.
        std::vector<HMesh::HalfEdgeID> newhalfedges = m.bridge_faces(f0, f1, connections);
    }


    void console_delete(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: edit.selected.merge_1_ring");
            return;
        }
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        auto vsel = me->get_vertex_selection();
        auto hsel = me->get_halfedge_selection();
        auto fsel = me->get_face_selection();
        for (auto v: vsel)
            if (m.in_use(v))
                m.remove_vertex(v);
        for (auto h: hsel)
            if (m.in_use(h))
                m.remove_edge(h);

        for (auto f: fsel)
            if (m.in_use(f))
                m.remove_face(f);
    }

    void console_split_face(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: edit.selected.split_face");
            return;
        }
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        auto vsel = me->get_vertex_selection();
        auto fsel = me->get_face_selection();
        for (auto f: fsel)
            if (m.in_use(f)) {
                HMesh::VertexID v0 = HMesh::InvalidVertexID, v1 = HMesh::InvalidVertexID;
                circulate_face_ccw(m, f, std::function<void(HMesh::VertexID)>([&](HMesh::VertexID v) {
                    if (vsel.count(v)) {
                        if (v0 == HMesh::InvalidVertexID)
                            v0 = v;
                        else
                            v1 = v;

                    }

                }));
                m.split_face_by_edge(f, v0, v1);
            }
    }


    void console_refit_trackball(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("display.refit_trackball");
            return;
        }
        me->refit();
    }


    void console_test(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: test");
            return;
        }

        me->save_active_mesh();
        me->active_mesh().slit_edges(me->get_vertex_selection());
    }

    void console_save(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: save <name->x3d|name->obj> ");

            return;
        }
        const std::string &file_name = args[0];
        if (args.size() == 1) {
            if (file_name.substr(file_name.length() - 4, file_name.length()) == ".obj") {
                obj_save(file_name, me->active_mesh());

                return;
            } else if (file_name.substr(file_name.length() - 4, file_name.length()) == ".off") {
                off_save(file_name, me->active_mesh());

                return;
            } else if (file_name.substr(file_name.length() - 4, file_name.length()) == ".x3d") {
                x3d_save(file_name, me->active_mesh());

                return;
            }
            me->printf("unknown format");
            return;
        }
        me->printf("usage: save <name->x3d|name->obj> ");
    }


    void console_refine_edges(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: refine.split_edges <length>");
            me->printf("splits edges longer than <length>; default is 0.5 times average length");
            return;
        }

        me->save_active_mesh();

        float thresh = 0.5f;

        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> thresh;
        }

        float avg_length = average_edge_length(me->active_mesh());

        refine_edges(me->active_mesh(), thresh * avg_length);

        return;

    }

    void console_refine_faces(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: refine.split_faces ");
            me->printf("usage:  Takes no arguments. Inserts a vertex at the centre of each face.");

            return;
        }
        me->save_active_mesh();

        triangulate_by_vertex_face_split(me->active_mesh());

        return;

    }

    void console_cc_subdivide(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: subdivide.catmull_clark ");
            me->printf("Does one step of Catmull-Clark subdivision");

            return;
        }
        me->save_active_mesh();

        cc_split(me->active_mesh(), me->active_mesh());
        cc_smooth(me->active_mesh());

        return;
    }

    void console_root_cc_subdivide(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: subdivide.catmull_clark ");
            me->printf("Does one step of Catmull-Clark subdivision");

            return;
        }
        me->save_active_mesh();

        rootCC_subdivide(me->active_mesh(), me->active_mesh());
        return;
    }


    void console_loop_subdivide(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: subdivide.loop");
            me->printf("Does one step of Loop subdivision");

            return;
        }
        me->save_active_mesh();

        loop_split(me->active_mesh(), me->active_mesh());
        loop_smooth(me->active_mesh());

        return;
    }

    void console_stitch(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: cleanup.stitch <rad>");
            me->printf("Stitches faces");

            return;
        }
        double r = 0.001;

        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> r;
        }

        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        CGLA::Vec3d c;
        float rad;
        bsphere(m, c, rad);
        stitch_mesh(me->active_mesh(), r * rad);
        return;
    }

    void console_remove_duplicates(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: cleanup.remove_duplicates <rad>");
            me->printf("Removes duplicate vertices and incident faces");

            return;
        }
        double r = 0.001;

        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> r;
        }

        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        CGLA::Vec3d c;
        float rad;
        bsphere(m, c, rad);
        remove_duplicates(me->active_mesh(), r * rad);
        return;
    }

    void console_compact(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: cleanup.compact");
            me->printf("Removes unreferenced vertices");

            return;
        }
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        m.cleanup();
        return;
    }


    void console_remove_val2(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: cleanup.remove_val2");
            me->printf("Removes valence 2 vertices");

            return;
        }
        me->save_active_mesh();
        HMesh::Manifold &m = me->active_mesh();
        remove_valence_two_vertices(m);
        return;
    }


    void console_root3_subdivide(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: subdivide.root3");
            me->printf("Does one step of sqrt(3) subdivision");

            return;
        }
        me->save_active_mesh();

        root3_subdivide(me->active_mesh(), me->active_mesh());

        return;
    }


    void console_doosabin_subdivide(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: subdivide.doo_sabin ");
            me->printf("Does one step of Doo-Sabin Subdivision");

            return;
        }
        me->save_active_mesh();

        cc_split(me->active_mesh(), me->active_mesh());
        dual(me->active_mesh());

        return;
    }

    void console_butterfly_subdivide(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: subdivide.butterfly ");
            me->printf("Does one step of Modified Butterfly Subdivision");

            return;
        }
        me->save_active_mesh();

        butterfly_subdivide(me->active_mesh(), me->active_mesh());

        return;
    }

    void console_dual(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: dual ");
            me->printf("Produces the dual by converting each face to a vertex placed at the barycenter.");
            return;
        }
        me->save_active_mesh();

        dual(me->active_mesh());

        return;
    }


    void console_minimize_curvature(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: optimize.minimize_curvature <anneal>");
            me->printf("Flip edges to minimize mean curvature.");
            me->printf("If anneal is true, simulated annealing (slow) is used rather than a greedy scheme");
            return;
        }
        me->save_active_mesh();

        bool anneal = false;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> anneal;
        }

        minimize_curvature(me->active_mesh(), anneal);
        me->post_create_display_list();
        return;
    }

    void console_minimize_dihedral(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: optimize.minimize_dihedral <iter> <anneal> <use_alpha> <gamma> ");
            me->printf("Flip edges to minimize dihedral angles.");
            me->printf("Iter is the max number of iterations. anneal tells us whether to use ");
            me->printf("simulated annealing and not greedy optimization. use_alpha (default=true) ");
            me->printf("means to use angle and not cosine of anglegamma (default=4) is the power ");
            me->printf("to which we raise the dihedral angle");
            return;
        }
        me->save_active_mesh();

        int iter = 1000;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> iter;
        }

        bool anneal = false;
        if (args.size() > 1) {
            std::istringstream a0(args[1]);
            a0 >> anneal;
        }

        bool use_alpha = true;
        if (args.size() > 2) {
            std::istringstream a0(args[2]);
            a0 >> use_alpha;
        }

        float gamma = 4.0f;
        if (args.size() > 3) {
            std::istringstream a0(args[3]);
            a0 >> gamma;
        }


        minimize_dihedral_angle(me->active_mesh(), iter, anneal, use_alpha, gamma);
        return;
    }

    void console_maximize_min_angle(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: optimize.maximize_min_angle <thresh> <anneal>");
            me->printf("Flip edges to maximize min angle - to make mesh more Delaunay.");
            me->printf("If the dot product of the normals between adjacent faces < thresh");
            me->printf("no flip will be made. anneal selects simulated annealing rather ");
            me->printf("nthan greedy optimization.");
            return;
        }
        me->save_active_mesh();

        float thresh = 0.0f;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> thresh;
        }
        bool anneal = false;
        if (args.size() > 1) {
            std::istringstream a0(args[1]);
            a0 >> anneal;
        }
        maximize_min_angle(me->active_mesh(), thresh, anneal);
        return;
    }


    void console_optimize_valency(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: optimize.valency <anneal> ");
            me->printf(
                    "Optimizes valency for triangle meshes. Anneal selects simulated annealing rather than greedy optim.");
            return;
        }
        me->save_active_mesh();

        bool anneal = false;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> anneal;
        }
        optimize_valency(me->active_mesh(), anneal);
        return;
    }


    void console_close_holes(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: cleanup.close_holes");
            me->printf("This function closes holes. It simply follows the loop of halfvectors which");
            me->printf("enclose the hole and add a face to which they all point.");
            return;
        }
        me->save_active_mesh();

        close_holes(me->active_mesh());
        return;
    }

    void console_reload(MeshEditor *me, const std::vector<std::string> &args) {
        std::string file_name = console_arg(args, 0, me->active_visobj().get_file_name());
        if (wantshelp(args)) {
            me->printf("usage:  load <file>");
            me->printf("(Re)loads the current file if no argument is given, but");
            me->printf("if an argument is given, then that becomes the current file");
            return;
        }
        me->save_active_mesh();

        if (me->reload_active_from_file(file_name))
            me->printf("Loaded %s", file_name.c_str());
        else
            me->printf("failed to load: %s", file_name.c_str());
    }


    void console_add_mesh(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage:  add_mesh <file>");
            me->printf("Loads the file but without clearing the mesh. Thus, the loaded mesh is added to the");
            me->printf("current model.");
            return;
        }
        me->save_active_mesh();

        if (!me->add_to_active_from_file(args.size() > 0 ? args[0] : ""))
            me->printf("failed to load");

        return;
    }

    void console_valid(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage:  validity");
            me->printf("Tests validity of Manifold");
            return;
        }
        if (valid(me->active_mesh()))
            me->printf("Mesh is valid");
        else
            me->printf("Mesh is invalid - check console output");
        return;
    }


    const CGLA::Vec3f &get_color(int i) {
        static CGLA::Vec3f ctable[100000];
        static bool was_here;
        CGLA::gel_srand(0);
        if (!was_here) {
            was_here = true;
            ctable[0] = CGLA::Vec3f(0);
            for (int j = 1; j < 100000; ++j)
                ctable[j] = CGLA::Vec3f(0.3) +
                            0.7 * normalize(CGLA::Vec3f(CGLA::gel_rand(), CGLA::gel_rand(), CGLA::gel_rand()));
            ctable[3] = CGLA::Vec3f(1, 0, 0);
            ctable[4] = CGLA::Vec3f(0, 1, 0);
            ctable[5] = CGLA::Vec3f(0, 0, 1);
            ctable[6] = CGLA::Vec3f(1, 0, 1);
        }
        return ctable[i % 100000];
    }

    void console_info_all(MeshEditor *me, const std::vector<std::string> &args) {
        CGLA::Vec3d p0_all(FLT_MAX), p7_all(-FLT_MAX);
        for (int i = 0; i < me->get_no_meshes(); ++i) {
            CGLA::Vec3d p0, p7;
            HMesh::Manifold &m = me->get_mesh(i);
            if (m.no_faces() > 0) {
                bbox(m, p0, p7);
                p0_all = v_min(p0, p0_all);
                p7_all = v_max(p7, p7_all);
            }
        }

        std::stringstream bbox_corners;
        bbox_corners << p0_all << " - " << p7_all << std::endl;
        me->printf("Bounding box corners : %s", bbox_corners.str().c_str());

    }

    void console_info(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage:  info");
            me->printf("Provides information about mesh.");
            return;
        }
        CGLA::Vec3d p0, p7;
        bbox(me->active_mesh(), p0, p7);
        std::stringstream bbox_corners;
        bbox_corners << p0 << " - " << p7 << std::endl;
        me->printf("Bounding box corners : %s", bbox_corners.str().c_str());
        std::map<int, int> val_hist;

        HMesh::Manifold &m = me->active_mesh();

        double avg_len = 0;
        for (HMesh::HalfEdgeID h: m.halfedges()) {
            DebugRenderer::edge_colors[h] = CGLA::Vec3f(0.3);
            avg_len += length(m, h);
        }
        avg_len /= m.no_halfedges();
        me->printf("Avg. edge length: %f", avg_len);
        for (HMesh::VertexID v: m.vertices()) {
            int val = valency(m, v);
            DebugRenderer::vertex_colors[v] = get_color(val);
            ++val_hist[val];

            if (val != 4)
                circulate_vertex_ccw(m, v,
                                     static_cast<std::function<void(HMesh::HalfEdgeID)>>([&](HMesh::HalfEdgeID h) {
                                         HMesh::Walker w = m.walker(h);
                                         DebugRenderer::edge_colors[h] = CGLA::Vec3f(1);
                                         DebugRenderer::edge_colors[w.opp().halfedge()] = CGLA::Vec3f(1);
                                         while (valency(m, w.vertex()) == 4) {
                                             w = w.next().opp().next();
                                             DebugRenderer::edge_colors[w.halfedge()] = CGLA::Vec3f(1);
                                             DebugRenderer::edge_colors[w.opp().halfedge()] = CGLA::Vec3f(1);
                                         }
                                     }));
        }
        std::map<int, int> ngon_hist;
        for (HMesh::FaceID f: m.faces()) {
            int ne = no_edges(m, f);
            ++ngon_hist[ne];
            DebugRenderer::face_colors[f] = 0.7 * get_color(ne);
        }

        me->printf("Valency histogram");
        for (std::map<int, int>::iterator iter = val_hist.begin(); iter != val_hist.end(); ++iter) {
            me->printf("%d, %d", iter->first, iter->second);
        }

        me->printf("Ngon histogram");
        for (std::map<int, int>::iterator iter = ngon_hist.begin(); iter != ngon_hist.end(); ++iter) {
            me->printf("%d, %d", iter->first, iter->second);
        }


        me->printf("Mesh contains %d faces", me->active_mesh().no_faces());
        me->printf("Mesh contains %d halfedges", me->active_mesh().no_halfedges());
        me->printf("Mesh contains %d vertices", me->active_mesh().no_vertices());
        return;
    }


    void console_simplify(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: simplify <fraction> ");
            me->printf("Performs Garland Heckbert (quadric based) mesh simplification.");
            me->printf("The only argument is the fraction of vertices to keep.");
            return;
        }
        me->save_active_mesh();

        float keep_fraction;
        if (args.size() == 0) {
            me->printf("you must specify fraction of vertices to keep");
            return;
        }
        std::istringstream a0(args[0]);
        a0 >> keep_fraction;

        bool optimal_positions = true;
        if (args.size() == 2) {
            std::istringstream a1(args[1]);
            a1 >> optimal_positions;
        }

        CGLA::Vec3d p0, p7;
        HMesh::Manifold m = me->active_mesh();
        bbox(m, p0, p7);
        CGLA::Vec3d d = p7 - p0;
        float s = d.max_coord();
        CGLA::Vec3d pcentre = (p7 + p0) / 2.0;
        for (HMesh::VertexID v: m.vertices()) {
            m.pos(v) -= pcentre;
            m.pos(v) /= s;
        }
        std::cout << "Timing the Garland Heckbert (quadric based) mesh simplication..." << std::endl;
        Util::Timer timer;
        timer.start();

        //simplify
        quadric_simplify(me->active_mesh(), keep_fraction, 0.0001f, optimal_positions);

        std::cout << "Simplification complete, process time: " << timer.get_secs() << " seconds" << std::endl;

        for (HMesh::VertexID v: m.vertices()) {
            m.pos(v) *= s;
            m.pos(v) += pcentre;
        }
        return;
    }

    void console_vertex_noise(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: noise.perturb_vertices <amplitude>");
            me->printf("adds a random vector to each vertex. A random vector in the unit cube is generated and");
            me->printf("to ensure an isotropic distribution, vectors outside the unit ball are discarded.");
            me->printf("The vector is multiplied by the average edge length and then by the amplitude specified.");
            me->printf("If no amplitude is specified, the default (0.5) is used.");
            return;
        }
        me->save_active_mesh();

        float avg_length = average_edge_length(me->active_mesh());

        float noise_amplitude = 0.5f;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> noise_amplitude;
        }

        CGLA::gel_srand(0);
        for (HMesh::VertexIDIterator vi = me->active_mesh().vertices_begin();
             vi != me->active_mesh().vertices_end(); ++vi) {
            CGLA::Vec3d v;
            do {
                v = CGLA::Vec3d(CGLA::gel_rand(), CGLA::gel_rand(), CGLA::gel_rand());
                v /= (float) (CGLA::GEL_RAND_MAX);
                v -= CGLA::Vec3d(0.5);
                v *= 2.0;
            } while (sqr_length(v) > 1.0);

            v *= noise_amplitude;
            v *= avg_length;
            me->active_mesh().pos(*vi) += v;
        }
        return;
    }

    void console_perpendicular_vertex_noise(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: noise.perturb_vertices_perpendicular <amplitude>");
            me->printf("adds the normal times a random scalar times amplitude times");
            me->printf("times average edge length to the vertex. (default amplitude=0.5)");
            return;
        }
        me->save_active_mesh();

        float avg_length = average_edge_length(me->active_mesh());

        float noise_amplitude = 0.5;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> noise_amplitude;
        }

        HMesh::VertexAttributeVector<CGLA::Vec3d> normals;
        for (HMesh::VertexIDIterator vi = me->active_mesh().vertices_begin();
             vi != me->active_mesh().vertices_end(); ++vi)
            normals[*vi] = normal(me->active_mesh(), *vi);

        CGLA::gel_srand(0);
        for (HMesh::VertexIDIterator vi = me->active_mesh().vertices_begin();
             vi != me->active_mesh().vertices_end(); ++vi) {
            float rval = 0.5 - CGLA::gel_rand() / float(CGLA::GEL_RAND_MAX);
            me->active_mesh().pos(*vi) += normals[*vi] * rval * noise_amplitude * avg_length * 2.0;
        }
        return;
    }

    void console_noisy_flips(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage:  noise.perturb_topology <iter>");
            me->printf("Perform random flips. iter (default=1) is the number of iterations.");
            me->printf("mostly for making nasty synthetic test cases.");
            return;
        }
        me->save_active_mesh();

        int iter = 1;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> iter;
        }

        randomize_mesh(me->active_mesh(), iter);
        return;
    }

    void console_laplacian_smooth(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage:  smooth.laplacian <weight> <iter>");
            me->printf("Perform Laplacian smoothing. weight is the scaling factor for the Laplacian.");
            me->printf("default weight = 1.0. Default number of iterations = 1");
            return;
        }
        me->save_active_mesh();

        float t = 1.0;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> t;
        }
        int iter = 1;
        if (args.size() > 1) {
            std::istringstream a0(args[1]);
            a0 >> iter;
        }
        Util::Timer tim;
        tim.start();
        /// Simple laplacian smoothing with an optional weight.
        laplacian_smooth(me->active_mesh(), t, iter);
        std::cout << "It took " << tim.get_secs();
        return;
    }


    void console_taubin_smooth(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage:  smooth.taubin <iter>");
            me->printf("Perform Taubin smoothing. iter (default=1) is the number of iterations.");
            return;
        }
        me->save_active_mesh();

        int iter = 1;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> iter;
        }
        /// Taubin smoothing is similar to laplacian smoothing but reduces shrinkage
        taubin_smooth(me->active_mesh(), iter);

        return;
    }

    void console_fvm_anisotropic_smooth(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: smooth.fuzzy_vector_median <iter>");
            me->printf(
                    "Smooth normals using fuzzy vector median smoothing. iter (default=1) is the number of iterations");
            me->printf("This function does a very good job of preserving sharp edges.");
            return;
        }
        me->save_active_mesh();

        int iter = 1;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> iter;
        }
        // Fuzzy vector median smoothing is effective when it comes to preserving sharp edges.
        anisotropic_smooth(me->active_mesh(), iter, HMesh::FVM_NORMAL_SMOOTH);

        return;
    }

    void console_bilateral_anisotropic_smooth(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: smooth.fuzzy_vector_median <iter>");
            me->printf(
                    "Smooth normals using fuzzy vector median smoothing. iter (default=1) is the number of iterations");
            me->printf("This function does a very good job of preserving sharp edges.");
            return;
        }
        me->save_active_mesh();

        int iter = 1;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> iter;
        }

        anisotropic_smooth(me->active_mesh(), iter, HMesh::BILATERAL_NORMAL_SMOOTH);

        return;
    }

    void console_triangulate(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage:  triangulate");
            me->printf("This function triangulates all non triangular faces of the mesh.");
            me->printf("you may want to call it after hole closing. For a polygon it simply connects");
            me->printf("the two closest vertices in a recursive manner until only triangles remain");
            return;
        }
        me->save_active_mesh();

        shortest_edge_triangulate(me->active_mesh());
        me->active_mesh().cleanup();
        valid(me->active_mesh());
        return;
    }


    void console_remove_caps(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage:  cleanup.remove_caps thresh");
            me->printf("Remove caps (triangles with one very big angle). The thresh argument is the fraction of PI to");
            me->printf("use as threshold for big angle. Default is 0.85. Caps are removed by flipping.");
            return;
        }
        me->save_active_mesh();

        float t = 0.85f;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> t;
        }
        remove_caps(me->active_mesh(), static_cast<float>(M_PI) * t);
        me->active_mesh().cleanup();

        return;
    }

    void console_remove_needles(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: cleanup.remove_needles <thresh>");
            me->printf("Removes very short edges by collapse. thresh is multiplied by the average edge length");
            me->printf("to get the length shorter than which we collapse. Default = 0.1");
            return;
        }
        me->save_active_mesh();

        float thresh = 0.1f;
        if (args.size() > 0) {
            std::istringstream a0(args[0]);
            a0 >> thresh;
        }
        float avg_length = average_edge_length(me->active_mesh());
        remove_needles(me->active_mesh(), thresh * avg_length);
        me->active_mesh().cleanup();

        return;
    }

    void console_flip_orientation(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage: cleanup.flip_orientation");
            me->printf("reorients all faces - flipping normal direction");
            return;
        }
        me->save_active_mesh();
        flip_orientation(me->active_mesh());
        return;
    }


    void console_undo(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage:  undo");
            me->printf("This function undoes one operation. Repeated undo does nothing");
            return;
        }
        me->restore_active_mesh();
        return;
    }

    void console_save_trackball(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage:  display.save_trackball");
            me->printf("This function saves the trackball to disk");
            return;
        }
        me->save_ball();
        return;
    }

    void console_load_trackball(MeshEditor *me, const std::vector<std::string> &args) {
        if (wantshelp(args)) {
            me->printf("usage:  display.load_trackball");
            me->printf("This function loads the trackball from disk");
            return;
        }
        me->load_ball();
        return;
    }

}

#endif //GEL_CONSOLEFUNCTIONS_H
