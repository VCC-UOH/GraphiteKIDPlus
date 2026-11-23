// -----------------------------------------------------------------------
// CRITICAL: ALL Geogram headers BEFORE any OGF header.
// -----------------------------------------------------------------------
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/basic/attributes.h>

#include <OGF/CurveKID/commands/mesh_grob_curvekid_commands.h>
#include <OGF/CurveKID/algo/mesh_geometry.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <limits>
#include <random>
#include <numeric>

namespace OGF {

// Local helper — mirrors the static vertex_point in mesh_geometry_updated.cpp
static inline GEO::vec3 cmd_vertex_point(const GEO::Mesh& mesh, GEO::index_t v) {
    const double* p = mesh.vertices.point_ptr(v);
    return GEO::vec3(p[0], p[1], p[2]);
}

// Normalize file path: replace backslashes with forward slashes.
// Windows paths with backslashes cause escape sequence errors in Lua/GOM.
static inline std::string normalize_path(std::string path) {
    for (char& c : path) if (c == '\\') c = '/';
    return path;
}

// Save mesh as .obj only (OBJ is universally readable; .geogram is Graphite-only).
// base_path may include any extension — it is stripped first.
static inline void save_mesh_both(const GEO::Mesh& mesh, const std::string& base_path) {
    std::string base = base_path;
    auto dot = base.rfind('.');
    if (dot != std::string::npos) {
        std::string ext = base.substr(dot);
        if (ext == ".geogram" || ext == ".obj")
            base = base.substr(0, dot);
    }
    GEO::mesh_save(mesh, base + ".obj");
    GEO::Logger::out("CurveKID") << "  Saved: " << base << ".obj" << std::endl;
}

    MeshGrobCurveKIDCommands::MeshGrobCurveKIDCommands() {
    }

    MeshGrobCurveKIDCommands::~MeshGrobCurveKIDCommands() {
    }

    // -------------------------------------------------------------------
    // inspect_mesh  (unchanged — preserved from working build)
    // -------------------------------------------------------------------
    void MeshGrobCurveKIDCommands::inspect_mesh() {
        MeshGrob* grob = mesh_grob();
        if (!grob) {
            GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl;
            return;
        }
        const GEO::Mesh& mesh = *grob;
        CurveKIDAlgo::MeshStats stats;
        if (!CurveKIDAlgo::compute_mesh_stats(mesh, stats)) {
            GEO::Logger::err("CurveKID")
                << "Mesh has no triangular facets." << std::endl;
            return;
        }
        GEO::Logger::out("CurveKID") << "=== Mesh Statistics ===" << std::endl;
        GEO::Logger::out("CurveKID") << "Vertices      : " << stats.nb_vertices << std::endl;
        GEO::Logger::out("CurveKID") << "Facets        : " << stats.nb_facets   << std::endl;
        GEO::Logger::out("CurveKID") << "Triangles     : " << stats.nb_triangles << std::endl;
        GEO::Logger::out("CurveKID") << "Non-tri facets: " << stats.nb_non_tri   << std::endl;
        GEO::Logger::out("CurveKID") << "Unique edges  : " << stats.nb_edges     << std::endl;
        GEO::Logger::out("CurveKID") << "Border edges  : " << stats.nb_border_edges << std::endl;
        GEO::Logger::out("CurveKID")
            << "Edge length   : min=" << stats.min_edge_length
            << " max=" << stats.max_edge_length
            << " mean=" << stats.mean_edge_length << std::endl;
        GEO::Logger::out("CurveKID")
            << "Face area     : min=" << stats.min_face_area
            << " max=" << stats.max_face_area
            << " mean=" << stats.mean_face_area << std::endl;
        GEO::Logger::out("CurveKID")
            << "Total area    : " << stats.total_surface_area << std::endl;
    }

    // -------------------------------------------------------------------
    // estimate_curvatures  (Stage 1 only — log summary, no visualisation)
    // -------------------------------------------------------------------
    void MeshGrobCurveKIDCommands::estimate_curvatures() {
        MeshGrob* grob = mesh_grob();
        if (!grob) {
            GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl;
            return;
        }
        const GEO::Mesh& mesh = *grob;

        GEO::Logger::out("CurveKID")
            << "Stage 1: computing vertex curvatures..." << std::endl;

        std::vector<CurveKIDAlgo::VertexCurvature> vcurv;
        CurveKIDAlgo::compute_vertex_curvatures(mesh, vcurv);

        std::vector<CurveKIDAlgo::FaceCurvature> fcurv;
        CurveKIDAlgo::compute_face_curvatures(mesh, vcurv, fcurv);

        // Log summary statistics
        double min_K= 1e30, max_K=-1e30, sum_K=0;
        double min_H= 1e30, max_H=-1e30, sum_H=0;
        int cnt = 0;
        for (auto& fc : fcurv) {
            min_K = std::min(min_K, fc.gaussian);
            max_K = std::max(max_K, fc.gaussian);
            sum_K += fc.gaussian;
            min_H = std::min(min_H, fc.mean);
            max_H = std::max(max_H, fc.mean);
            sum_H += fc.mean;
            ++cnt;
        }
        if (cnt > 0) {
            GEO::Logger::out("CurveKID")
                << "Gaussian curvature: min=" << min_K
                << " max=" << max_K
                << " mean=" << sum_K/cnt << std::endl;
            GEO::Logger::out("CurveKID")
                << "Mean curvature    : min=" << min_H
                << " max=" << max_H
                << " mean=" << sum_H/cnt << std::endl;
        }
        GEO::Logger::out("CurveKID") << "Stage 1 done." << std::endl;
    }

    // -------------------------------------------------------------------
    // run_curvekid  (Stages 1-5)
    // -------------------------------------------------------------------
    void MeshGrobCurveKIDCommands::run_curvekid(int nb_patches, int K) {
        MeshGrob* grob = mesh_grob();
        if (!grob) {
            GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl;
            return;
        }
        GEO::Mesh& mesh = *grob;

        // --- Stage 1: Curvature estimation ---
        GEO::Logger::out("CurveKID")
            << "[1/5] Estimating curvatures..." << std::endl;

        std::vector<CurveKIDAlgo::VertexCurvature> vcurv;
        CurveKIDAlgo::compute_vertex_curvatures(mesh, vcurv);

        std::vector<CurveKIDAlgo::FaceCurvature> fcurv;
        CurveKIDAlgo::compute_face_curvatures(mesh, vcurv, fcurv);

        // --- Stage 2: Connected patch growing ---
        GEO::Logger::out("CurveKID")
            << "[2/5] Growing ~" << nb_patches
            << " connected patches..." << std::endl;

        std::vector<GEO::index_t> patch_id;
        GEO::index_t actual_patches = CurveKIDAlgo::grow_patches(
            mesh, fcurv,
            (GEO::index_t)nb_patches,
            patch_id);

        GEO::Logger::out("CurveKID")
            << "      Created " << actual_patches << " patches." << std::endl;

        // --- Stage 3: Sphere fitting per patch ---
        GEO::Logger::out("CurveKID")
            << "[3/5] Fitting spheres to patches..." << std::endl;

        std::vector<CurveKIDAlgo::PatchDescriptor> descs;
        CurveKIDAlgo::compute_patch_descriptors(
            mesh, patch_id, fcurv, actual_patches, descs);

        // Log sphere radius range
        double rmin=1e30, rmax=-1e30, rsum=0; int rcnt=0;
        for (auto& d : descs) {
            if (d.nb_faces == 0) continue;
            rmin = std::min(rmin, d.sphere_radius);
            rmax = std::max(rmax, d.sphere_radius);
            rsum += d.sphere_radius;
            ++rcnt;
        }
        if (rcnt > 0)
            GEO::Logger::out("CurveKID")
                << "      Sphere radius: min=" << rmin
                << " max=" << rmax
                << " mean=" << rsum/rcnt << std::endl;

        // --- Stage 4: K-means clustering ---
        GEO::Logger::out("CurveKID")
            << "[4/5] K-means clustering (K=" << K << ")..." << std::endl;

        double inertia = CurveKIDAlgo::cluster_patches_kmeans(
            descs, (GEO::index_t)K);

        GEO::Logger::out("CurveKID")
            << "      K-means inertia = " << inertia << std::endl;

        // Count faces per template
        std::vector<int> tmpl_count(K, 0);
        for (GEO::index_t f = 0; f < mesh.facets.nb(); ++f) {
            if (mesh.facets.nb_vertices(f) != 3) continue;
            GEO::index_t pid = patch_id[f];
            if (pid < actual_patches &&
                descs[pid].template_id != GEO::index_t(-1))
            {
                int tid = (int)descs[pid].template_id;
                if (tid < K) ++tmpl_count[tid];
            }
        }
        for (int t = 0; t < K; ++t)
            GEO::Logger::out("CurveKID")
                << "      Template " << t
                << ": " << tmpl_count[t] << " faces" << std::endl;

        // --- Stage 5: Write facet attributes ---
        GEO::Logger::out("CurveKID")
            << "[5/5] Writing facet attributes..." << std::endl;

        // Graphite reads Attribute<double> on facets for colour painting.
        // Create three attributes: patch_id, template_id, sphere_residual.
        GEO::Attribute<double> attr_patch(
            mesh.facets.attributes(), "patch_id");
        GEO::Attribute<double> attr_tmpl(
            mesh.facets.attributes(), "template_id");
        GEO::Attribute<double> attr_err(
            mesh.facets.attributes(), "sphere_residual");

        for (GEO::index_t f = 0; f < mesh.facets.nb(); ++f) {
            GEO::index_t pid = patch_id[f];

            attr_patch[f] = (double)pid;

            if (pid < actual_patches) {
                GEO::index_t tid = descs[pid].template_id;
                attr_tmpl[f] = (tid == GEO::index_t(-1))
                               ? -1.0 : (double)tid;
                attr_err[f]  = descs[pid].sphere_residual;
            } else {
                attr_tmpl[f] = -1.0;
                attr_err[f]  =  0.0;
            }
        }

        // Tell Graphite the mesh changed so it redraws
        grob->update();

        GEO::Logger::out("CurveKID")
            << "Done. Attributes written: patch_id, template_id, "
               "sphere_residual." << std::endl;
        GEO::Logger::out("CurveKID")
            << "To visualise: Mesh > Attributes > facets > "
               "template_id (or patch_id / sphere_residual)." << std::endl;
    }

    // -------------------------------------------------------------------
    // validate_curvekid  — 6-question diagnostic (run run_curvekid first)
    // -------------------------------------------------------------------
    void MeshGrobCurveKIDCommands::validate_curvekid(int nb_patches, int K) {
        MeshGrob* grob = mesh_grob();
        if (!grob) {
            GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl;
            return;
        }
        GEO::Mesh& mesh = *grob;

        // Check that run_curvekid has been called (attributes must exist)
        if (!mesh.facets.attributes().is_defined("patch_id") ||
            !mesh.facets.attributes().is_defined("template_id")) {
            GEO::Logger::err("CurveKID")
                << "Attributes not found. Run run_curvekid first." << std::endl;
            return;
        }

        GEO::Attribute<double> attr_patch(mesh.facets.attributes(), "patch_id");
        GEO::Attribute<double> attr_tmpl(mesh.facets.attributes(), "template_id");
        GEO::Attribute<double> attr_err(mesh.facets.attributes(), "sphere_residual");

        GEO::index_t nb_faces = mesh.facets.nb();

        GEO::Logger::out("CurveKID")
            << "========== CurveKID Validation Report ==========" << std::endl;

        // ---- Q1: Count non-empty patches and templates, invalid IDs ----
        std::vector<bool> patch_seen(nb_patches, false);
        std::vector<bool> tmpl_seen(K, false);
        int invalid_patch = 0, invalid_tmpl = 0;
        for (GEO::index_t f = 0; f < nb_faces; ++f) {
            int pid = (int)attr_patch[f];
            int tid = (int)attr_tmpl[f];
            if (pid < 0 || pid >= nb_patches) { ++invalid_patch; continue; }
            patch_seen[pid] = true;
            if (tid < 0 || tid >= K) { ++invalid_tmpl; }
            else tmpl_seen[tid] = true;
        }
        int nonempty_patches = 0, nonempty_tmpls = 0;
        for (int i = 0; i < nb_patches; ++i) if (patch_seen[i]) ++nonempty_patches;
        for (int i = 0; i < K;          ++i) if (tmpl_seen[i])  ++nonempty_tmpls;

        GEO::Logger::out("CurveKID") << "[Q1] Non-empty patches  : "
            << nonempty_patches << " / " << nb_patches << std::endl;
        GEO::Logger::out("CurveKID") << "[Q1] Non-empty templates: "
            << nonempty_tmpls << " / " << K << std::endl;
        GEO::Logger::out("CurveKID") << "[Q1] Invalid patch IDs  : "
            << invalid_patch << std::endl;
        GEO::Logger::out("CurveKID") << "[Q1] Invalid template IDs: "
            << invalid_tmpl << std::endl;

        // ---- Q2: BFS connectivity check per patch ----
        // Build face adjacency via shared edges
        // face_adj[f] = list of adjacent face indices
        std::vector<std::vector<GEO::index_t>> face_adj(nb_faces);
        {
            std::map<std::pair<GEO::index_t,GEO::index_t>, GEO::index_t> edge_face;
            for (GEO::index_t f = 0; f < nb_faces; ++f) {
                GEO::index_t n = mesh.facets.nb_vertices(f);
                for (GEO::index_t lv = 0; lv < n; ++lv) {
                    GEO::index_t va = mesh.facets.vertex(f, lv);
                    GEO::index_t vb = mesh.facets.vertex(f, (lv+1) % n);
                    if (va > vb) std::swap(va, vb);
                    auto key = std::make_pair(va, vb);
                    auto it = edge_face.find(key);
                    if (it != edge_face.end()) {
                        face_adj[f].push_back(it->second);
                        face_adj[it->second].push_back(f);
                        edge_face.erase(it);
                    } else {
                        edge_face[key] = f;
                    }
                }
            }
        }

        int disconnected_patches = 0;
        for (int pid = 0; pid < nb_patches; ++pid) {
            if (!patch_seen[pid]) continue;
            // collect faces of this patch
            std::vector<GEO::index_t> pfaces;
            for (GEO::index_t f = 0; f < nb_faces; ++f)
                if ((int)attr_patch[f] == pid) pfaces.push_back(f);
            if (pfaces.empty()) continue;

            // BFS
            std::unordered_map<GEO::index_t,bool> in_patch;
            for (auto f : pfaces) in_patch[f] = false;
            std::vector<GEO::index_t> queue;
            queue.push_back(pfaces[0]);
            in_patch[pfaces[0]] = true;
            int visited = 1;
            for (int qi = 0; qi < (int)queue.size(); ++qi) {
                GEO::index_t cur = queue[qi];
                for (GEO::index_t nb : face_adj[cur]) {
                    auto it = in_patch.find(nb);
                    if (it != in_patch.end() && !it->second) {
                        it->second = true;
                        ++visited;
                        queue.push_back(nb);
                    }
                }
            }
            if (visited < (int)pfaces.size()) ++disconnected_patches;
        }
        GEO::Logger::out("CurveKID") << "[Q2] Disconnected patches: "
            << disconnected_patches << " (should be 0)" << std::endl;

        // ---- Q3: Histogram of patches and faces per template ----
        std::vector<int> tmpl_patch_count(K, 0);
        std::vector<int> tmpl_face_count(K, 0);
        // per patch, find its template
        std::vector<int> patch_tmpl(nb_patches, -1);
        for (GEO::index_t f = 0; f < nb_faces; ++f) {
            int pid = (int)attr_patch[f];
            int tid = (int)attr_tmpl[f];
            if (pid >= 0 && pid < nb_patches && tid >= 0 && tid < K)
                patch_tmpl[pid] = tid;
        }
        for (int pid = 0; pid < nb_patches; ++pid) {
            int tid = patch_tmpl[pid];
            if (tid >= 0 && tid < K) ++tmpl_patch_count[tid];
        }
        for (GEO::index_t f = 0; f < nb_faces; ++f) {
            int tid = (int)attr_tmpl[f];
            if (tid >= 0 && tid < K) ++tmpl_face_count[tid];
        }
        int max_faces = *std::max_element(tmpl_face_count.begin(), tmpl_face_count.end());
        GEO::Logger::out("CurveKID") << "[Q3] Template histogram (patches | faces | bar):" << std::endl;
        for (int t = 0; t < K; ++t) {
            int bar_len = (max_faces > 0) ? (tmpl_face_count[t] * 30 / max_faces) : 0;
            std::string bar(bar_len, '#');
            GEO::Logger::out("CurveKID")
                << "  T" << t
                << " patches=" << tmpl_patch_count[t]
                << " faces=" << tmpl_face_count[t]
                << " |" << bar << std::endl;
        }

        // ---- Q4: Sphere residuals normalised by bbox diagonal ----
        // Compute bbox diagonal
        double xmin=1e30,xmax=-1e30,ymin=1e30,ymax=-1e30,zmin=1e30,zmax=-1e30;
        for (GEO::index_t v = 0; v < mesh.vertices.nb(); ++v) {
            GEO::vec3 p(mesh.vertices.point_ptr(v));
            xmin=std::min(xmin,p.x); xmax=std::max(xmax,p.x);
            ymin=std::min(ymin,p.y); ymax=std::max(ymax,p.y);
            zmin=std::min(zmin,p.z); zmax=std::max(zmax,p.z);
        }
        double diag = std::sqrt((xmax-xmin)*(xmax-xmin)+
                                (ymax-ymin)*(ymax-ymin)+
                                (zmax-zmin)*(zmax-zmin));

        double sum_res=0, sum_res2=0, max_res=-1e30;
        int res_cnt = 0;
        int degen_patches = 0;
        for (GEO::index_t f = 0; f < nb_faces; ++f) {
            double r = attr_err[f];
            sum_res  += r;
            sum_res2 += r * r;
            max_res   = std::max(max_res, r);
            ++res_cnt;
        }
        double mean_res = (res_cnt>0) ? sum_res/res_cnt : 0;
        double rms_res  = (res_cnt>0) ? std::sqrt(sum_res2/res_cnt) : 0;

        GEO::Logger::out("CurveKID")
            << "[Q4] Bbox diagonal: " << diag << std::endl;
        GEO::Logger::out("CurveKID")
            << "[Q4] Sphere residual (absolute) mean=" << mean_res
            << " RMS=" << rms_res << " max=" << max_res << std::endl;
        GEO::Logger::out("CurveKID")
            << "[Q4] Sphere residual (normalised) mean=" << mean_res/diag
            << " RMS=" << rms_res/diag << " max=" << max_res/diag << std::endl;

        // ---- Q5: Dominant template analysis ----
        int dom_t = (int)(std::max_element(tmpl_face_count.begin(),
                                           tmpl_face_count.end())
                          - tmpl_face_count.begin());
        GEO::Logger::out("CurveKID")
            << "[Q5] Dominant template: T" << dom_t
            << " (" << tmpl_face_count[dom_t] << " faces, "
            << (100.0*tmpl_face_count[dom_t]/nb_faces) << "%)" << std::endl;

        // Find patches in dominant template and collect patch_ids for per-patch stats
        // We need to re-run curvature to get patch descriptors... we don't have them
        // stored, so we just report face-level residuals for the dominant template.
        double dom_res_sum=0; int dom_cnt=0;
        for (GEO::index_t f = 0; f < nb_faces; ++f) {
            if ((int)attr_tmpl[f] == dom_t) {
                dom_res_sum += attr_err[f];
                ++dom_cnt;
            }
        }
        GEO::Logger::out("CurveKID")
            << "[Q5] Dominant template mean sphere_residual: "
            << (dom_cnt>0 ? dom_res_sum/dom_cnt : 0.0) << std::endl;
        GEO::Logger::out("CurveKID")
            << "[Q5] Note: dominant template may reflect near-flat regions "
               "(large-radius spheres ~ flat patches) clustering together." << std::endl;

        // ---- Q6: Verify template_id is consistent within each patch ----
        int inconsistent = 0;
        std::vector<int> patch_first_tmpl(nb_patches, -2);
        for (GEO::index_t f = 0; f < nb_faces; ++f) {
            int pid = (int)attr_patch[f];
            int tid = (int)attr_tmpl[f];
            if (pid < 0 || pid >= nb_patches) continue;
            if (patch_first_tmpl[pid] == -2) {
                patch_first_tmpl[pid] = tid;
            } else if (patch_first_tmpl[pid] != tid) {
                ++inconsistent;
                patch_first_tmpl[pid] = -3; // mark as known-inconsistent
            }
        }
        GEO::Logger::out("CurveKID")
            << "[Q6] Patches with inconsistent template_id: "
            << inconsistent << " (should be 0 — classification is per-patch)" << std::endl;

        GEO::Logger::out("CurveKID")
            << "========== Validation complete ==========" << std::endl;
    }

    // =========================================================================
    // Pottmann-inspired Stage 1: central spheres
    // =========================================================================

    void MeshGrobCurveKIDCommands::inspect_central_spheres(double H_eps)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) {
            GEO::Logger::err("CurveKID") << "No mesh loaded." << std::endl;
            return;
        }
        GEO::Mesh& mesh = *grob;

        const GEO::index_t nv = mesh.vertices.nb();
        const GEO::index_t nf = mesh.facets.nb();
        if (nv == 0 || nf == 0) {
            GEO::Logger::err("CurveKID")
                << "Mesh has no vertices or faces." << std::endl;
            return;
        }

        GEO::Logger::out("CurveKID")
            << "=== Pottmann Stage 1: Central spheres "
            << "(H_eps=" << H_eps << ") ===" << std::endl;

        // --- Stage 1a: vertex curvatures (reuse existing pipeline) ---
        std::vector<VertexCurvature> vcurv;
        compute_vertex_curvatures(mesh, vcurv);

        // --- Stage 1b: central spheres ---
        std::vector<CentralSphere> cspheres;
        compute_central_spheres(mesh, vcurv, cspheres, H_eps);

        // --- Write vertex attributes for Graphite visualisation ---
        // Attribute: central_sphere_radius (double)
        // Attribute: central_sphere_valid  (int, 0=fallback, 1=valid)
        GEO::Attribute<double> attr_radius(
            mesh.vertices.attributes(), "central_sphere_radius");
        GEO::Attribute<int>    attr_valid(
            mesh.vertices.attributes(), "central_sphere_valid");

        int n_valid    = 0;
        int n_fallback = 0;
        double r_min_valid = std::numeric_limits<double>::max();
        double r_max_valid = 0.0;
        double r_sum_valid = 0.0;
        double r_max_fallback = 0.0;
        double r_sum_fallback = 0.0;

        for (GEO::index_t v = 0; v < nv; ++v) {
            attr_radius[v] = cspheres[v].radius;
            attr_valid[v]  = cspheres[v].valid ? 1 : 0;

            if (cspheres[v].valid) {
                ++n_valid;
                r_min_valid = std::min(r_min_valid, cspheres[v].radius);
                r_max_valid = std::max(r_max_valid, cspheres[v].radius);
                r_sum_valid += cspheres[v].radius;
            } else {
                ++n_fallback;
                r_max_fallback = std::max(r_max_fallback, cspheres[v].radius);
                r_sum_fallback += cspheres[v].radius;
            }
        }

        // --- Summary ---
        GEO::Logger::out("CurveKID")
            << "  Vertices total  : " << nv << std::endl;
        GEO::Logger::out("CurveKID")
            << "  Valid (|H|>eps) : " << n_valid
            << "  ("
            << (100.0 * n_valid / (nv > 0 ? nv : 1)) << "%)" << std::endl;
        GEO::Logger::out("CurveKID")
            << "  Fallback (|H|~0): " << n_fallback
            << "  ("
            << (100.0 * n_fallback / (nv > 0 ? nv : 1)) << "%)" << std::endl;

        if (n_valid > 0) {
            GEO::Logger::out("CurveKID")
                << "  Valid radius — min: " << r_min_valid
                << "  max: " << r_max_valid
                << "  mean: " << (r_sum_valid / n_valid) << std::endl;
        }
        if (n_fallback > 0) {
            GEO::Logger::out("CurveKID")
                << "  Fallback radius — max: " << r_max_fallback
                << "  mean: " << (r_sum_fallback / n_fallback) << std::endl;
        }

        GEO::Logger::out("CurveKID")
            << "  Attributes written: central_sphere_radius, central_sphere_valid"
            << std::endl;
        GEO::Logger::out("CurveKID")
            << "  -> In Graphite: Attributes > vertex > central_sphere_radius"
            << std::endl;
        GEO::Logger::out("CurveKID")
            << "=== Stage 1 complete ===" << std::endl;

        grob->update();
    }

    // =========================================================================
    // render_curved_patches
    // =========================================================================

    void MeshGrobCurveKIDCommands::render_curved_patches(
        int arc_subdivisions, int nb_patches, int K)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh." << std::endl; return; }
        GEO::Mesh& mesh = *grob;

        // --- Check required attributes exist ---
        if (!mesh.facets.attributes().is_defined("patch_id") ||
            !mesh.facets.attributes().is_defined("template_id")) {
            GEO::Logger::err("CurveKID")
                << "Attributes not found. Run run_curvekid first." << std::endl;
            return;
        }

        GEO::Attribute<int> attr_patch(mesh.facets.attributes(),   "patch_id");
        GEO::Attribute<int> attr_tmpl (mesh.facets.attributes(),   "template_id");

        const GEO::index_t nf = mesh.facets.nb();
        const int N = std::max(2, arc_subdivisions);

        // --- Rebuild patch sphere descriptors by re-running pipeline ---
        GEO::Logger::out("CurveKID") << "Re-computing patch descriptors..." << std::endl;
        std::vector<VertexCurvature> vcurv;
        compute_vertex_curvatures(mesh, vcurv);
        std::vector<FaceCurvature> fcurv;
        compute_face_curvatures(mesh, vcurv, fcurv);

        // Read patch_id from attribute (already computed by run_curvekid)
        std::vector<GEO::index_t> patch_id(nf);
        for (GEO::index_t f = 0; f < nf; ++f)
            patch_id[f] = (GEO::index_t)attr_patch[f];

        std::vector<PatchDescriptor> descs;
        compute_patch_descriptors(mesh, patch_id, fcurv, (GEO::index_t)nb_patches, descs);

        // --- Compute sphere_quality face attribute ---
        // quality = 1 - clamp(residual / (patch_diameter * 0.1), 0, 1)
        // patch_diameter estimated as sqrt(total_area) * 2
        GEO::Attribute<double> attr_quality(mesh.facets.attributes(), "sphere_quality");
        // per-patch diameter estimate
        std::vector<double> patch_diam(nb_patches, 1.0);
        for (int p = 0; p < nb_patches; ++p) {
            if (p < (int)descs.size() && descs[p].nb_faces > 0)
                patch_diam[p] = std::sqrt(descs[p].total_area) * 2.0;
        }
        for (GEO::index_t f = 0; f < nf; ++f) {
            int pid = attr_patch[f];
            if (pid < 0 || pid >= nb_patches) { attr_quality[f] = 0.0; continue; }
            double diam = patch_diam[pid];
            double res  = (pid < (int)descs.size()) ? descs[pid].sphere_residual : 0.0;
            double q = (diam > 1e-10) ? (1.0 - std::min(1.0, res / (diam * 0.1))) : 0.0;
            attr_quality[f] = std::max(0.0, q);
        }
        GEO::Logger::out("CurveKID")
            << "sphere_quality written as face attribute (1=good sphere fit, 0=poor)."
            << std::endl;

        // --- Find boundary edges (facet corners where adjacent patch differs) ---
        // Use facet adjacency: for each facet f, corner lc, adjacent facet is
        // mesh.facets.adjacent(f, lc)
        // Edge v0-v1 is a boundary edge if adj patch_id differs

        // Build arc polylines into a dedicated MeshGrob.
        // We create it as a plain GEO::Mesh here (not a scene object),
        // then save it as a .geogram file the user can load separately,
        // AND also attempt to inject it into the scene if the API permits.
        GEO::Mesh arcs_mesh;
        GEO::Mesh& arcs = arcs_mesh;

        // Avoid duplicate edges (each boundary edge is shared by 2 faces)
        // Key: (min_v, max_v) → bool processed
        std::unordered_map<GEO::index_t,
            std::unordered_map<GEO::index_t, bool>> edge_done;

        int arc_count = 0;
        int poor_fit  = 0; // patches where sphere_residual/diam > 0.1

        for (GEO::index_t f = 0; f < nf; ++f) {
            int pid_f = attr_patch[f];
            if (pid_f < 0 || pid_f >= nb_patches) continue;
            const PatchDescriptor& desc = descs[pid_f];
            if (desc.sphere_radius < 1e-10) continue;

            GEO::index_t nvc = mesh.facets.nb_vertices(f);
            for (GEO::index_t lc = 0; lc < nvc; ++lc) {
                GEO::index_t adj = mesh.facets.adjacent(f, lc);
                bool is_boundary = (adj == GEO::NO_FACET) ||
                                   ((int)attr_patch[adj] != pid_f);
                if (!is_boundary) continue;

                GEO::index_t v0 = mesh.facets.vertex(f, lc);
                GEO::index_t v1 = mesh.facets.vertex(f, (lc + 1) % nvc);
                GEO::index_t lo = std::min(v0, v1);
                GEO::index_t hi = std::max(v0, v1);
                if (edge_done[lo][hi]) continue;
                edge_done[lo][hi] = true;

                // Great-circle arc on the patch sphere
                const double* a = mesh.vertices.point_ptr(v0);
                const double* b = mesh.vertices.point_ptr(v1);
                GEO::vec3 C = desc.sphere_center;
                double    R = desc.sphere_radius;

                GEO::vec3 pA(a[0]-C.x, a[1]-C.y, a[2]-C.z);
                GEO::vec3 pB(b[0]-C.x, b[1]-C.y, b[2]-C.z);
                double lenA = GEO::length(pA);
                double lenB = GEO::length(pB);
                if (lenA < 1e-12 || lenB < 1e-12) continue;
                pA = pA * (R / lenA);
                pB = pB * (R / lenB);

                // Add N+1 vertices for the arc
                GEO::index_t base_v = arcs.vertices.nb();
                arcs.vertices.create_vertices(N + 1);
                for (int s = 0; s <= N; ++s) {
                    double t = double(s) / double(N);
                    // Slerp between pA and pB on the sphere
                    GEO::vec3 interp = pA * (1.0 - t) + pB * t;
                    double len = GEO::length(interp);
                    if (len < 1e-12) interp = pA;
                    else interp = interp * (R / len);
                    GEO::vec3 world = interp + C;
                    double* pt = arcs.vertices.point_ptr(base_v + s);
                    pt[0] = world.x; pt[1] = world.y; pt[2] = world.z;
                }
                // Add N edge facets
                for (int s = 0; s < N; ++s) {
                    arcs.edges.create_edge(base_v + s, base_v + s + 1);
                }
                ++arc_count;
            }
        }

        // Count poor-fit patches
        for (int p = 0; p < nb_patches && p < (int)descs.size(); ++p) {
            double diam = patch_diam[p];
            if (diam > 1e-10 && descs[p].sphere_residual / (diam * 0.1) > 1.0)
                ++poor_fit;
        }

        // Save arc mesh using the mesh's scene name (e.g. "kitten_arcs.geogram").
        // File lands in Graphite's current working directory.
        std::string arc_path = std::string(grob->name()) + "_arcs";
        save_mesh_both(arcs, arc_path);

        grob->update();

        GEO::Logger::out("CurveKID")
            << "Boundary arcs: " << arc_count
            << " edges, " << N << " subdivisions each." << std::endl;
        GEO::Logger::out("CurveKID")
            << "Poor sphere fit (residual > 10% diam): "
            << poor_fit << " / " << nb_patches << " patches." << std::endl;
        GEO::Logger::out("CurveKID")
            << "-> Visualise: colour mesh by 'sphere_quality' attribute (face)." << std::endl;
        GEO::Logger::out("CurveKID")
            << "-> Arc mesh saved to: " << arc_path << std::endl;
        GEO::Logger::out("CurveKID")
            << "-> In Graphite: File > Open > CurveKID_arcs.geogram to overlay arcs." << std::endl;
    }

    // =========================================================================
    // save_curvekid  — write .ckid (JSON) segmentation file
    // =========================================================================

    void MeshGrobCurveKIDCommands::save_curvekid(std::string ckid_filename)
    {
        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh." << std::endl; return; }
        GEO::Mesh& mesh = *grob;

        if (!mesh.facets.attributes().is_defined("patch_id") ||
            !mesh.facets.attributes().is_defined("template_id")) {
            GEO::Logger::err("CurveKID")
                << "No segmentation found. Run run_curvekid first." << std::endl;
            return;
        }

        // Bind required attributes — these were written as int by run_curvekid
        GEO::Attribute<int> attr_patch, attr_tmpl;
        attr_patch.bind(mesh.facets.attributes(), "patch_id");
        attr_tmpl .bind(mesh.facets.attributes(), "template_id");
        if (!attr_patch.is_bound() || !attr_tmpl.is_bound()) {
            GEO::Logger::err("CurveKID")
                << "Could not bind patch_id/template_id — wrong type? Run run_curvekid again." << std::endl;
            return;
        }

        const GEO::index_t nf = mesh.facets.nb();

        // sphere_quality is optional (only exists if render_curved_patches was run)
        bool has_quality = mesh.facets.attributes().is_defined("sphere_quality");

        // sphere_cx/cy/cz/r are optional too
        bool has_sphere_cx = mesh.facets.attributes().is_defined("sphere_cx");

        // Work on a local copy so we can modify the path
        std::string filename = ckid_filename;
        // If empty, use mesh scene name + .ckid in current working directory
        if (filename.empty()) {
            filename = std::string(grob->name()) + ".ckid";
        }

        std::ofstream out(filename);
        if (!out) {
            GEO::Logger::err("CurveKID")
                << "Cannot open file for writing: " << filename << std::endl;
            return;
        }

        out << std::fixed << std::setprecision(8);
        out << "{\n";
        out << "  \"curvekid_version\": 1,\n";
        out << "  \"mesh_name\": \"" << grob->name() << "\",\n";
        out << "  \"nb_facets\": " << nf << ",\n";
        out << "  \"nb_vertices\": " << mesh.vertices.nb() << ",\n";

        // Bind optional attributes only if they exist
        // Use scoped pointers to avoid Geogram's assert on unbound access
        GEO::Attribute<double> attr_quality;
        if (has_quality)
            attr_quality.bind(mesh.facets.attributes(), "sphere_quality");

        GEO::Attribute<double> s_cx, s_cy, s_cz, s_r;
        if (has_sphere_cx) {
            s_cx.bind(mesh.facets.attributes(), "sphere_cx");
            s_cy.bind(mesh.facets.attributes(), "sphere_cy");
            s_cz.bind(mesh.facets.attributes(), "sphere_cz");
            s_r .bind(mesh.facets.attributes(), "sphere_r");
        }

        // Face data array
        out << "  \"faces\": [\n";
        for (GEO::index_t f = 0; f < nf; ++f) {
            out << "    {\"f\":" << f
                << ",\"patch\":" << attr_patch[f]
                << ",\"tmpl\":"  << attr_tmpl[f];
            if (has_quality && attr_quality.is_bound())
                out << ",\"q\":" << attr_quality[f];
            if (has_sphere_cx && s_cx.is_bound())
                out << ",\"cx\":" << s_cx[f]
                    << ",\"cy\":" << s_cy[f]
                    << ",\"cz\":" << s_cz[f]
                    << ",\"r\":"  << s_r[f];
            out << "}";
            if (f + 1 < nf) out << ",";
            out << "\n";
        }
        out << "  ]\n";
        out << "}\n";
        out.close();

        GEO::Logger::out("CurveKID")
            << "Segmentation saved to: " << filename
            << "  (" << nf << " faces)" << std::endl;
    }

    // =========================================================================
    // load_curvekid  — read .ckid (JSON) and restore face attributes
    // =========================================================================

    void MeshGrobCurveKIDCommands::load_curvekid(std::string ckid_filename)
    {
        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh." << std::endl; return; }
        GEO::Mesh& mesh = *grob;

        // If empty, use mesh scene name + .ckid in current working directory
        std::string filepath = ckid_filename;
        if (filepath.empty()) {
            filepath = std::string(grob->name()) + ".ckid";
            GEO::Logger::out("CurveKID")
                << "Auto-derived path: " << filepath << std::endl;
        }

        std::ifstream in(filepath);
        if (!in) {
            GEO::Logger::err("CurveKID")
                << "Cannot open file: " << filepath << std::endl;
            return;
        }

        // Simple hand-rolled JSON parser — no dependency on a JSON library.
        // Reads the flat "faces" array produced by save_curvekid().
        std::string content(
            (std::istreambuf_iterator<char>(in)),
             std::istreambuf_iterator<char>());
        in.close();

        const GEO::index_t nf = mesh.facets.nb();

        GEO::Attribute<int>    attr_patch(mesh.facets.attributes(), "patch_id");
        GEO::Attribute<int>    attr_tmpl (mesh.facets.attributes(), "template_id");
        GEO::Attribute<double> attr_quality(mesh.facets.attributes(), "sphere_quality");

        // Parse: find each {"f":<n>,"patch":<p>,"tmpl":<t>,...}
        int loaded = 0;
        size_t pos = 0;
        while (pos < content.size()) {
            // Find next { after "faces" array starts
            size_t obj_start = content.find("{\"f\":", pos);
            if (obj_start == std::string::npos) break;
            size_t obj_end   = content.find('}', obj_start);
            if (obj_end   == std::string::npos) break;
            std::string obj = content.substr(obj_start, obj_end - obj_start + 1);
            pos = obj_end + 1;

            // Extract field value by key
            auto get_int = [&](const std::string& key, int def) -> int {
                std::string search = "\"" + key + "\":";
                size_t p = obj.find(search);
                if (p == std::string::npos) return def;
                p += search.size();
                return std::stoi(obj.substr(p));
            };
            auto get_dbl = [&](const std::string& key, double def) -> double {
                std::string search = "\"" + key + "\":";
                size_t p = obj.find(search);
                if (p == std::string::npos) return def;
                p += search.size();
                return std::stod(obj.substr(p));
            };

            int f     = get_int("f",     -1);
            int patch = get_int("patch", -1);
            int tmpl  = get_int("tmpl",  -1);
            double q  = get_dbl("q",     -1.0);

            if (f < 0 || (GEO::index_t)f >= nf) continue;
            attr_patch[f] = patch;
            attr_tmpl[f]  = tmpl;
            if (q >= 0.0) attr_quality[f] = q;
            ++loaded;
        }

        grob->update();
        GEO::Logger::out("CurveKID")
            << "Loaded segmentation from: " << filepath << std::endl;
        GEO::Logger::out("CurveKID")
            << "  " << loaded << " faces restored (mesh has " << nf << ")." << std::endl;
        if (loaded != (int)nf)
            GEO::Logger::warn("CurveKID")
                << "  Face count mismatch — file may be for a different mesh." << std::endl;
        GEO::Logger::out("CurveKID")
            << "  Attributes restored: patch_id, template_id, sphere_quality." << std::endl;
        GEO::Logger::out("CurveKID")
            << "  -> Colour by 'patch_id' or 'template_id' in Attributes panel." << std::endl;
    }

    // =========================================================================
    // build_curved_mesh — Stages 1+2+3: coarse remesh → sphere fit → tessellate
    // =========================================================================

    void MeshGrobCurveKIDCommands::build_curved_mesh(int target_faces, int subdiv)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh." << std::endl; return; }

        GEO::Logger::out("CurveKID") << "=== build_curved_mesh ===" << std::endl;
        GEO::Logger::out("CurveKID")
            << "  target_faces=" << target_faces
            << "  subdiv=" << subdiv << std::endl;

        // --- Stage 0: per-vertex curvatures (needed for adaptive sizing) ---
        std::vector<VertexCurvature> vcurv;
        compute_vertex_curvatures(*grob, vcurv);
        GEO::Logger::out("CurveKID")
            << "  Curvatures computed for " << vcurv.size() << " vertices." << std::endl;

        // --- Stage 1: coarse curvature-adaptive remesh ---
        GEO::Mesh coarse_mesh;
        GEO::index_t nf_coarse = build_coarse_mesh(
            *grob, vcurv, coarse_mesh,
            GEO::index_t(target_faces), /*H_eps=*/0.01);
        GEO::Logger::out("CurveKID")
            << "  Coarse mesh: " << nf_coarse << " faces, "
            << coarse_mesh.vertices.nb() << " vertices." << std::endl;

        if (nf_coarse == 0) {
            GEO::Logger::err("CurveKID")
                << "  Coarse mesh is empty — remeshing failed." << std::endl;
            return;
        }

        // --- Stage 2: fit one sphere per coarse triangle ---
        std::vector<CurvedFace> cfaces;
        fit_curved_faces(*grob, coarse_mesh, cfaces);
        {
            int valid_cnt = 0;
            for (auto& cf : cfaces) if (cf.valid) ++valid_cnt;
            GEO::Logger::out("CurveKID")
                << "  Sphere fitting: " << valid_cnt << "/" << cfaces.size()
                << " faces valid." << std::endl;
            double mean_res = 0.0;
            for (auto& cf : cfaces) if (cf.valid) mean_res += cf.sphere_residual;
            if (valid_cnt > 0) mean_res /= valid_cnt;
            GEO::Logger::out("CurveKID")
                << "  Mean sphere residual: " << mean_res << std::endl;
        }

        // --- Stage 3: tessellate onto spheres for visualization ---
        GEO::Mesh viz_mesh;
        tessellate_curved_faces(coarse_mesh, cfaces, viz_mesh, subdiv);
        GEO::Logger::out("CurveKID")
            << "  Tessellated viz mesh: " << viz_mesh.facets.nb()
            << " faces, " << viz_mesh.vertices.nb() << " vertices." << std::endl;

        // --- Stage 4 (light): approximation error ---
        std::vector<double> errors;
        compute_approximation_error(*grob, cfaces, errors);
        {
            double mean_e = 0.0, max_e = 0.0;
            for (double e : errors) {
                mean_e += e;
                max_e = std::max(max_e, e);
            }
            if (!errors.empty()) mean_e /= errors.size();
            GEO::Logger::out("CurveKID")
                << "  Approx error (vertex-to-nearest-sphere): "
                << "mean=" << mean_e << "  max=" << max_e << std::endl;
        }

        // --- Store the coarse-mesh sphere params as face attributes ---
        {
            GEO::Attribute<double> attr_r(
                viz_mesh.facets.attributes(), "sphere_radius");
            GEO::Attribute<double> attr_res(
                viz_mesh.facets.attributes(), "sphere_residual");
            GEO::Attribute<int>    attr_id(
                viz_mesh.facets.attributes(), "curved_face_id");
            // attr_id is already written by tessellate_curved_faces;
            // fill radius+residual from the parent CurvedFace
            for (GEO::index_t f = 0; f < viz_mesh.facets.nb(); ++f) {
                if (!attr_id.is_bound()) break;
                int cf_id = attr_id[f];
                if (cf_id >= 0 && cf_id < (int)cfaces.size()) {
                    if (attr_r.is_bound())   attr_r[f]   = cfaces[cf_id].sphere_radius;
                    if (attr_res.is_bound()) attr_res[f] = cfaces[cf_id].sphere_residual;
                }
            }
        }

        // --- Create / replace "CurveKID_CurvedMesh" in the scene ---
        // We write to a geogram file and log the path.
        save_mesh_both(viz_mesh,   std::string(grob->name()) + "_curved");
        save_mesh_both(coarse_mesh, std::string(grob->name()) + "_coarse");

        GEO::Logger::out("CurveKID") << "=== Done ===" << std::endl;
    }

    // =========================================================================
    // build_coarse_curved_mesh — PRIMARY paper experiment
    // M (original) -> N (curved) with E_curved(N) vs E_planar(N) comparison
    // =========================================================================

    void MeshGrobCurveKIDCommands::build_coarse_curved_mesh(int N, int subdiv, double alpha)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh." << std::endl; return; }

        const GEO::Mesh& original = *grob;
        const GEO::index_t M = original.facets.nb();

        // N must be strictly less than M — otherwise the coarse mesh IS the
        // original and E_planar = 0 (trivially, no approximation error).
        if ((GEO::index_t)N >= M) {
            GEO::Logger::err("CurveKID")
                << "  N (" << N << ") must be less than M (" << M << "). "
                << "Try N = " << (M / 5) << " for a 5x reduction." << std::endl;
            return;
        }

        GEO::Logger::out("CurveKID") << std::string(60, '=') << std::endl;
        GEO::Logger::out("CurveKID") << "  CurveKID: M -> N curved face reduction" << std::endl;
        GEO::Logger::out("CurveKID") << "  Original mesh M = " << M << " faces" << std::endl;
        GEO::Logger::out("CurveKID") << "  Target   curved N = " << N << " faces" << std::endl;
        GEO::Logger::out("CurveKID") << std::string(60, '=') << std::endl;

        // --- Stage 0: per-vertex curvatures (adaptive remesh sizing) ---
        std::vector<VertexCurvature> vcurv;
        compute_vertex_curvatures(original, vcurv);

        // --- Stage 1: coarse remesh -> N flat triangles ---
        GEO::Mesh coarse_flat;
        GEO::index_t nf_actual = build_coarse_mesh(
            original, vcurv, coarse_flat, GEO::index_t(N), 0.01);

        GEO::Logger::out("CurveKID")
            << "  Coarse flat mesh: " << nf_actual << " faces  ("
            << coarse_flat.vertices.nb() << " vertices)" << std::endl;

        if (nf_actual == 0) {
            GEO::Logger::err("CurveKID") << "  Coarse remesh failed." << std::endl;
            return;
        }

        // --- Topology validation ---
        MeshTopologyReport topo = validate_mesh_topology(coarse_flat, true);
        if (!topo.valid) {
            // Save flat mesh — user can repair in Graphite then use fit_curved_to_coarse
            std::string flat_base =
                std::string(grob->name()) + "_flat_N" + std::to_string(nf_actual);
            save_mesh_both(coarse_flat, flat_base);
            GEO::Logger::err("CurveKID")
                << "  Coarse mesh has connectivity issues (Components="
                << topo.connected_components << " BoundaryEdges=" << topo.boundary_edges << ")." << std::endl;
            GEO::Logger::out("CurveKID")
                << "  FIX: Load it in Graphite -> Surface->Repair->Repair mesh -> Save -> fit_curved_to_coarse" << std::endl;
            GEO::Logger::out("CurveKID")
                << "  OR: Use Surface->Remesh->Remesh isotropic for a clean coarse mesh." << std::endl;
            return;
        }

        // --- E_planar(N): how well do N flat triangles approximate original ---
        std::vector<double> err_planar;
        compute_planar_approximation_error(original, coarse_flat, err_planar);
        double mean_planar = 0.0, max_planar = 0.0, rms_planar = 0.0;
        for (double e : err_planar) {
            mean_planar += e;
            max_planar   = std::max(max_planar, e);
            rms_planar  += e * e;
        }
        if (!err_planar.empty()) {
            mean_planar /= err_planar.size();
            rms_planar   = std::sqrt(rms_planar / err_planar.size());
        }

        // --- Stage 2: fit one sphere per coarse triangle ---
        // Run BOTH methods for comparison (ChatGPT Direction A/B/C experiments).
        // Method B = original unconstrained Pottmann Eq.7
        // Method C = constrained corner-interpolating fit (Experiment C)
        std::vector<CurvedFace> cfaces_B, cfaces_C;
        fit_curved_faces(original, coarse_flat, cfaces_B,
                         CurveKIDAlgo::SphereFittingMode::UNCONSTRAINED_7_SAMPLE);
        fit_curved_faces(original, coarse_flat, cfaces_C,
                         CurveKIDAlgo::SphereFittingMode::CONSTRAINED_CORNER_INTERPOLATION);

        // Log method comparison: corner residual, sphere RMS, radius distribution
        {
            auto log_method = [&](const char* tag,
                                  const std::vector<CurveKIDAlgo::CurvedFace>& cfs,
                                  const GEO::Mesh& cmesh) {
                int valid = 0; double sum_res = 0, sum_R = 0, sum_corner = 0;
                double min_R = 1e18, max_R = 0;
                for (GEO::index_t f = 0; f < cmesh.facets.nb(); ++f) {
                    const auto& cf = cfs[f];
                    if (!cf.valid) continue;
                    ++valid;
                    sum_res    += cf.sphere_residual;
                    sum_R      += cf.sphere_radius;
                    min_R       = std::min(min_R, cf.sphere_radius);
                    max_R       = std::max(max_R, cf.sphere_radius);
                    // Corner incidence residual
                    for (int k = 0; k < 3; ++k) {
                        GEO::index_t vi = cmesh.facets.vertex(f, k);
                        const double* pp = cmesh.vertices.point_ptr(vi);
                        GEO::vec3 vp(pp[0], pp[1], pp[2]);
                        double e = std::fabs(GEO::length(vp - cf.sphere_center)
                                             - cf.sphere_radius);
                        sum_corner += e;
                    }
                }
                if (valid == 0) return;
                GEO::Logger::out("CurveKID")
                    << "  [" << tag << "] valid=" << valid
                    << "  mean_res=" << sum_res/valid
                    << "  mean_R=" << sum_R/valid
                    << "  R_range=[" << min_R << "," << max_R << "]"
                    << "  mean_corner_res=" << sum_corner/(3*valid) << std::endl;
            };
            GEO::Logger::out("CurveKID") << "  --- Fitting method comparison ---" << std::endl;
            log_method("Method-B unconstrained", cfaces_B, coarse_flat);
            log_method("Method-C constrained",   cfaces_C, coarse_flat);
        }

        // Use Method B (original unconstrained) as the active fitting.
        // Method C is logged for comparison but not used for tessellation.
        std::vector<CurvedFace> cfaces = cfaces_B;

        // --- Stage 2b: joint optimisation (Pottmann 2023 Sec 4) ---
        // Only run for closed manifolds (no boundary edges).
        {
            MeshTopologyReport tr = validate_mesh_topology(coarse_flat, false);
            if (tr.boundary_edges == 0) {
                GEO::Logger::out("CurveKID") << "  Running joint sphere optimisation..." << std::endl;
                optimise_joint_spheres(coarse_flat, cfaces, 30);
            } else {
                GEO::Logger::out("CurveKID")
                    << "  Skipping joint opt (open surface, "
                    << tr.boundary_edges << " boundary edges)" << std::endl;
            }
        }

        int valid_cnt = 0, flat_cnt = 0;
        double sum_res = 0.0;
        for (const CurvedFace& cf : cfaces) {
            if (cf.valid) { ++valid_cnt; sum_res += cf.sphere_residual; }
            else          { ++flat_cnt; }
        }
        double mean_res = (valid_cnt > 0) ? sum_res / valid_cnt : 0.0;
        GEO::Logger::out("CurveKID")
            << "  Face partition: " << valid_cnt << " curved + "
            << flat_cnt << " flat = " << cfaces.size() << " total" << std::endl;
        GEO::Logger::out("CurveKID")
            << "  Mean sphere fit residual (curved faces): " << mean_res << std::endl;

        // --- Stage 3: shared arc resolution ---
        std::vector<GEO::index_t> edge_winner;
        resolve_shared_arcs(coarse_flat, cfaces, edge_winner);

        // --- Stage 4: tessellate for visualization ---
        // Mode: BLENDED_SHARED_EDGE_DEBUG (ChatGPT Direction 016, Option 2)
        // Each interior shared edge uses e(t) = 0.5*slerp_A(t) + 0.5*slerp_B(t)
        // → zero seam between adjacent panels.
        // NOTE: output is a blended approximation, NOT exact spherical panels.
        GEO::Mesh viz_mesh;
        tessellate_curved_faces_v2(coarse_flat, cfaces, edge_winner, viz_mesh, subdiv, alpha,
            CurveKIDAlgo::TessellationMode::BLENDED_SHARED_EDGE_DEBUG);
        GEO::Logger::out("CurveKID")
            << "  Viz mesh: " << viz_mesh.facets.nb()
            << " render triangles  (" << subdiv << "x subdiv per curved face)"
            << "  [BLENDED_SHARED_EDGE_DEBUG]"
            << std::endl;

        // --- E_curved(N): how well do N curved spherical faces approximate original ---
        std::vector<double> err_curved;
        compute_approximation_error(original, cfaces, err_curved);
        double mean_curved = 0.0, max_curved = 0.0, rms_curved = 0.0;
        for (double e : err_curved) {
            mean_curved += e;
            max_curved   = std::max(max_curved, e);
            rms_curved  += e * e;
        }
        if (!err_curved.empty()) {
            mean_curved /= err_curved.size();
            rms_curved   = std::sqrt(rms_curved / err_curved.size());
        }

        // --- Summary ---
        GEO::Logger::out("CurveKID") << std::string(60, '-') << std::endl;
        GEO::Logger::out("CurveKID")
            << "  RESULTS: M=" << M << "  N=" << nf_actual << std::endl;
        GEO::Logger::out("CurveKID")
            << "  E_planar(N)  mean=" << mean_planar
            << "  rms=" << rms_planar
            << "  max=" << max_planar << std::endl;
        GEO::Logger::out("CurveKID")
            << "  E_curved(N)  mean=" << mean_curved
            << "  rms=" << rms_curved
            << "  max=" << max_curved << std::endl;
        double improve_mean = (mean_planar > 1e-14) ?
            (1.0 - mean_curved / mean_planar) * 100.0 : 0.0;
        double improve_rms  = (rms_planar > 1e-14) ?
            (1.0 - rms_curved  / rms_planar)  * 100.0 : 0.0;
        GEO::Logger::out("CurveKID")
            << "  Improvement: mean=" << improve_mean << "%"
            << "  rms=" << improve_rms << "%" << std::endl;
        GEO::Logger::out("CurveKID") << std::string(60, '=') << std::endl;

        // --- Save both meshes in .geogram AND .obj ---
        std::string base = std::string(grob->name()) + "_N" + std::to_string(nf_actual);
        save_mesh_both(viz_mesh,    base + "_curved");
        save_mesh_both(coarse_flat, base + "_flat");
    }

    // =========================================================================
    // =========================================================================
    // fit_curved_to_coarse — fit spheres to a pre-built CVT coarse mesh
    // =========================================================================
    void MeshGrobCurveKIDCommands::fit_curved_to_coarse(
        std::string coarse_mesh_name, int subdiv, double alpha)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl; return; }
        const GEO::Mesh& original = *grob;

        // Load the coarse mesh from a file path
        coarse_mesh_name = normalize_path(coarse_mesh_name);
        if (coarse_mesh_name.empty()) {
            GEO::Logger::err("CurveKID")
                << "Specify the path to the coarse mesh file "
                << "(e.g. C:/path/to/sphere_remesh.geogram)." << std::endl;
            return;
        }
        GEO::Mesh coarse_flat;
        if (!GEO::mesh_load(coarse_mesh_name, coarse_flat)) {
            GEO::Logger::err("CurveKID")
                << "Could not load coarse mesh from '" << coarse_mesh_name << "'." << std::endl;
            return;
        }

        GEO::Logger::out("CurveKID") << std::string(60, '=') << std::endl;
        GEO::Logger::out("CurveKID") << "  fit_curved_to_coarse" << std::endl;
        GEO::Logger::out("CurveKID") << "  Original: " << original.facets.nb() << " faces" << std::endl;
        GEO::Logger::out("CurveKID") << "  Coarse  : " << coarse_flat.facets.nb()
                                     << " faces, " << coarse_flat.vertices.nb() << " vertices" << std::endl;

        // --- Topology validation ---
        MeshTopologyReport topo = validate_mesh_topology(coarse_flat, true);
        if (!topo.valid) {
            GEO::Logger::err("CurveKID")
                << "  Coarse mesh failed topology check. "
                << "Use Surface->Remesh->Remesh isotropic to get a valid coarse mesh." << std::endl;
            return;
        }

        // --- E_planar(N) ---
        std::vector<double> err_planar;
        compute_planar_approximation_error(original, coarse_flat, err_planar);
        double mean_planar=0, rms_planar=0, max_planar=0;
        for (double e : err_planar) {
            mean_planar += e; max_planar = std::max(max_planar, e); rms_planar += e*e;
        }
        if (!err_planar.empty()) {
            mean_planar /= err_planar.size();
            rms_planar = std::sqrt(rms_planar / err_planar.size());
        }

        // --- Fit spheres (Pottmann Eq.7) ---
        std::vector<CurvedFace> cfaces;
        fit_curved_faces(original, coarse_flat, cfaces);

        // --- Joint optimisation — closed surfaces only ---
        {
            MeshTopologyReport tr2 = validate_mesh_topology(coarse_flat, false);
            if (tr2.boundary_edges == 0) {
                GEO::Logger::out("CurveKID") << "  Running joint sphere optimisation..." << std::endl;
                optimise_joint_spheres(coarse_flat, cfaces, 30);
            } else {
                GEO::Logger::out("CurveKID")
                    << "  Skipping joint opt (open surface)" << std::endl;
            }
        }

        int curved_cnt = 0;
        for (auto& cf : cfaces) if (cf.valid) ++curved_cnt;
        GEO::Logger::out("CurveKID") << "  Curved faces: " << curved_cnt
                                     << " / " << cfaces.size() << std::endl;

        // --- Shared arcs + tessellation (BLENDED_SHARED_EDGE_DEBUG) ---
        std::vector<GEO::index_t> edge_winner;
        resolve_shared_arcs(coarse_flat, cfaces, edge_winner);
        GEO::Mesh viz_mesh;
        tessellate_curved_faces_v2(coarse_flat, cfaces, edge_winner, viz_mesh, subdiv, alpha,
            CurveKIDAlgo::TessellationMode::BLENDED_SHARED_EDGE_DEBUG);

        // --- E_curved(N) ---
        std::vector<double> err_curved;
        compute_approximation_error(original, cfaces, err_curved);
        double mean_curved=0, rms_curved=0, max_curved=0;
        for (double e : err_curved) {
            mean_curved += e; max_curved = std::max(max_curved, e); rms_curved += e*e;
        }
        if (!err_curved.empty()) {
            mean_curved /= err_curved.size();
            rms_curved = std::sqrt(rms_curved / err_curved.size());
        }

        double improve_mean = (mean_planar>1e-14) ? (1.0 - mean_curved/mean_planar)*100.0 : 0.0;
        double improve_rms  = (rms_planar >1e-14) ? (1.0 - rms_curved /rms_planar )*100.0 : 0.0;

        GEO::Logger::out("CurveKID") << std::string(60, '-') << std::endl;
        GEO::Logger::out("CurveKID") << "  E_planar : mean=" << mean_planar << " rms=" << rms_planar << " max=" << max_planar << std::endl;
        GEO::Logger::out("CurveKID") << "  E_curved : mean=" << mean_curved << " rms=" << rms_curved << " max=" << max_curved << std::endl;
        GEO::Logger::out("CurveKID") << "  Improvement: mean=" << improve_mean << "% rms=" << improve_rms << "%" << std::endl;
        GEO::Logger::out("CurveKID") << std::string(60, '=') << std::endl;

        save_mesh_both(viz_mesh, coarse_mesh_name + "_curved");
    }

    // =========================================================================
    // classify_curved_faces — K-set template classification (2nd contribution)
    // =========================================================================
    void MeshGrobCurveKIDCommands::classify_curved_faces(
        std::string curved_mesh_path, std::string flat_mesh_path, int K)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl; return; }
        const GEO::Mesh& original = *grob;

        // Load curved viz mesh (has curved_face_id attribute)
        curved_mesh_path = normalize_path(curved_mesh_path);
        if (curved_mesh_path.empty()) {
            GEO::Logger::err("CurveKID") << "Specify path to curved viz mesh." << std::endl;
            return;
        }
        GEO::Mesh curved_mesh;
        if (!GEO::mesh_load(curved_mesh_path, curved_mesh)) {
            GEO::Logger::err("CurveKID") << "Could not load: " << curved_mesh_path << std::endl;
            return;
        }

        // Re-fit curved faces to get CurvedFace descriptors
        // We need the flat coarse mesh — use provided path or auto-derive from name
        flat_mesh_path = normalize_path(flat_mesh_path);
        std::string flat_path;
        if (!flat_mesh_path.empty()) {
            flat_path = flat_mesh_path;
        } else {
            // Auto-derive: replace "_curved_" with "_flat_" in curved mesh path
            flat_path = curved_mesh_path;
            size_t pos = flat_path.find("_curved_");
            if (pos != std::string::npos) flat_path.replace(pos, 8, "_flat_");
            else {
                // Try replacing "_curved." with "_flat."
                pos = flat_path.find("_curved.");
                if (pos != std::string::npos) flat_path.replace(pos, 8, "_flat.");
                else {
                    GEO::Logger::err("CurveKID")
                        << "Cannot auto-find flat mesh. Provide flat_mesh_path explicitly." << std::endl;
                    return;
                }
            }
        }
        GEO::Mesh coarse_flat;
        if (!GEO::mesh_load(flat_path, coarse_flat)) {
            GEO::Logger::err("CurveKID") << "Could not load flat mesh: " << flat_path << std::endl;
            return;
        }
        GEO::Logger::out("CurveKID") << "  Flat mesh: " << flat_path << std::endl;

        const GEO::index_t N = coarse_flat.facets.nb();
        GEO::Logger::out("CurveKID") << std::string(60,'=') << std::endl;
        GEO::Logger::out("CurveKID") << "  STAGE 2: K-set classification" << std::endl;
        GEO::Logger::out("CurveKID") << "  N=" << N << " curved faces -> K=" << K << " templates" << std::endl;

        // Recompute sphere fits
        std::vector<CurvedFace> cfaces;
        fit_curved_faces(original, coarse_flat, cfaces);

        int valid_cnt = 0;
        for (auto& cf : cfaces) if (cf.valid) ++valid_cnt;
        GEO::Logger::out("CurveKID") << "  Curved faces: " << valid_cnt << " / " << N << std::endl;

        // Build pose-invariant descriptors matching SphericalFace::desc_* fields.
        // 6D feature vector per face:
        //   [0] log(R+1)         — curvature magnitude (log-compressed)
        //   [1] arc₀ (smallest)  } sorted great-circle arc lengths
        //   [2] arc₁ (middle)    } SSS congruence: 3 arcs determine spherical
        //   [3] arc₂ (largest)   } triangle up to reflection
        //   [4] solid_angle      — subtended area on unit sphere
        //   [5] mean_H_local     — mean curvature of original surface neighbourhood
        // This matches the CurveKID SphericalFace descriptor design.
        int nf = (int)N;
        const int FDIM = 6;
        std::vector<double> feat((size_t)nf * FDIM, 0.0);
        for (int i = 0; i < nf; ++i) {
            const CurvedFace& cf = cfaces[(GEO::index_t)i];
            if (!cf.valid) continue;
            GEO::vec3 C = cf.sphere_center;
            double R = cf.sphere_radius;
            GEO::vec3 p[3];
            for (int k=0;k<3;++k) p[k] = cmd_vertex_point(coarse_flat, cf.v[k]);
            // Arc lengths (great-circle): R * acos(dot(unit(p-C), unit(q-C)))
            double arcs[3] = {0,0,0};
            for (int e=0;e<3;++e) {
                GEO::vec3 a=p[e]-C, b=p[(e+1)%3]-C;
                double la=GEO::length(a), lb=GEO::length(b);
                if (la>1e-14 && lb>1e-14) {
                    double cs=GEO::dot(a,b)/(la*lb);
                    cs=std::max(-1.0,std::min(1.0,cs));
                    arcs[e]=R*std::acos(cs);
                }
            }
            std::sort(arcs,arcs+3); // sort → pose invariant
            // Solid angle (Van Oosterom & Strackee)
            GEO::vec3 va=p[0]-C, vb=p[1]-C, vc=p[2]-C;
            double la=GEO::length(va),lb=GEO::length(vb),lc=GEO::length(vc);
            double solid_angle=0.0;
            if (la>1e-14&&lb>1e-14&&lc>1e-14) {
                double num=std::fabs(GEO::dot(va,GEO::cross(vb,vc)));
                double den=la*lb*lc+GEO::dot(va,vb)*lc+GEO::dot(vb,vc)*la+GEO::dot(vc,va)*lb;
                solid_angle=2.0*std::atan2(num,den);
                if (solid_angle<0) solid_angle+=4.0*3.14159265358979323846;
            }
            feat[(size_t)i*FDIM+0]=std::log1p(R);
            feat[(size_t)i*FDIM+1]=arcs[0];
            feat[(size_t)i*FDIM+2]=arcs[1];
            feat[(size_t)i*FDIM+3]=arcs[2];
            feat[(size_t)i*FDIM+4]=solid_angle;
            feat[(size_t)i*FDIM+5]=cf.mean_H_local;
        }
        // Compute mean and std per feature
        // Normalise each feature by its std across all faces
        std::vector<double> meanF(FDIM,0.0), stdF(FDIM,1.0);
        for (int fi=0;fi<FDIM;++fi) {
            for (int i=0;i<nf;++i) meanF[fi]+=feat[(size_t)i*FDIM+fi];
            meanF[fi]/=nf;
            for (int i=0;i<nf;++i){ double d2=feat[(size_t)i*FDIM+fi]-meanF[fi]; stdF[fi]+=d2*d2; }
            stdF[fi]=std::sqrt(stdF[fi]/nf);
            if (stdF[fi]<1e-10) stdF[fi]=1.0;
        }
        for (int i=0;i<nf;++i)
            for (int fi=0;fi<FDIM;++fi)
                feat[(size_t)i*FDIM+fi]=(feat[(size_t)i*FDIM+fi]-meanF[fi])/stdF[fi];

        // K-means++ init + Lloyd iterations
        std::mt19937 rng(42);
        std::vector<int> centres;
        centres.push_back((int)(rng()%(unsigned)nf));
        std::vector<double> dists(nf,1e30);
        auto fdist=[&](int a,int b){
            double s=0;
            for(int fi=0;fi<FDIM;++fi){double d2=feat[(size_t)a*FDIM+fi]-feat[(size_t)b*FDIM+fi];s+=d2*d2;}
            return s;
        };
        while ((int)centres.size()<std::min(K,nf)) {
            int last=centres.back();
            double sum=0;
            for (int i=0;i<nf;++i){ double d2=fdist(i,last); if(d2<dists[i])dists[i]=d2; sum+=dists[i]; }
            double r=std::uniform_real_distribution<>(0,sum)(rng);
            double acc=0; int ch=0;
            for (int i=0;i<nf;++i){acc+=dists[i];if(acc>=r){ch=i;break;}}
            centres.push_back(ch);
        }
        int Kact=(int)centres.size();
        std::vector<double> ctr((size_t)Kact*FDIM);
        for (int ki=0;ki<Kact;++ki)
            for (int fi=0;fi<FDIM;++fi) ctr[(size_t)ki*FDIM+fi]=feat[(size_t)centres[ki]*FDIM+fi];

        std::vector<int> assign(nf,0);
        for (int iter=0;iter<50;++iter) {
            bool changed=false;
            for (int i=0;i<nf;++i) {
                double best=1e30; int bk=0;
                for (int ki=0;ki<Kact;++ki){
                    double s=0;
                    for(int fi=0;fi<FDIM;++fi){double d2=feat[(size_t)i*FDIM+fi]-ctr[(size_t)ki*FDIM+fi];s+=d2*d2;}
                    if(s<best){best=s;bk=ki;}
                }
                if(assign[i]!=bk){changed=true;assign[i]=bk;}
            }
            if(!changed) break;
            std::fill(ctr.begin(),ctr.end(),0.0);
            std::vector<int> cnt(Kact,0);
            for (int i=0;i<nf;++i){
                int ki=assign[i]; cnt[ki]++;
                for(int fi=0;fi<FDIM;++fi) ctr[(size_t)ki*FDIM+fi]+=feat[(size_t)i*FDIM+fi];
            }
            for (int ki=0;ki<Kact;++ki)
                if(cnt[ki]>0)
                    for(int fi=0;fi<FDIM;++fi) ctr[(size_t)ki*FDIM+fi]/=cnt[ki];
        }

        // Log class sizes
        std::vector<int> class_size(Kact,0);
        for (int i=0;i<nf;++i) class_size[assign[i]]++;
        GEO::Logger::out("CurveKID") << "  K=" << Kact << " template classes:" << std::endl;
        for (int ki=0;ki<Kact;++ki)
            GEO::Logger::out("CurveKID") << "    Class " << ki << ": " << class_size[ki] << " faces" << std::endl;

        // Write template_id attribute to curved viz mesh (one value per coarse face,
        // broadcast to all sub-triangles of that face via curved_face_id attribute)
        GEO::Attribute<int> attr_face_id;
        bool has_fid = curved_mesh.facets.attributes().is_defined("curved_face_id");
        if (has_fid)
            attr_face_id.bind(curved_mesh.facets.attributes(), "curved_face_id");

        GEO::Attribute<int> attr_tmpl(curved_mesh.facets.attributes(), "template_id");
        for (GEO::index_t tf=0; tf<curved_mesh.facets.nb(); ++tf) {
            int fid = has_fid ? attr_face_id[tf] : -1;
            if (fid>=0 && fid<nf) attr_tmpl[tf]=assign[fid];
            else attr_tmpl[tf]=-1;
        }

        // Save classified mesh in .geogram and .obj
        std::string out_base = curved_mesh_path.substr(0, curved_mesh_path.rfind('.')) + "_classified";
        save_mesh_both(curved_mesh, out_base);
        GEO::Logger::out("CurveKID") << "  Paint by 'template_id' attribute to see K=" << Kact << " classes." << std::endl;
        GEO::Logger::out("CurveKID") << std::string(60,'=') << std::endl;
    }

    // =========================================================================
    // build_sphere_ground_truth — perfect watertight curved sphere
    // =========================================================================
    // =========================================================================
    // K=1 ICOSAHEDRAL BASELINE (ChatGPT Direction 017)
    // =========================================================================
    void MeshGrobCurveKIDCommands::build_icosahedral_baseline(int subdiv, double radius)
    {
        using namespace CurveKIDAlgo;

        GEO::Logger::out("CurveKID") << "=== K=1 Icosahedral Baseline ===" << std::endl;
        GEO::Logger::out("CurveKID") << "  radius=" << radius << "  subdiv=" << subdiv << std::endl;

        // ---- Step 1: Build icosahedron vertices ----
        // 12 vertices of a regular icosahedron inscribed in unit sphere.
        // Golden ratio phi = (1 + sqrt(5)) / 2
        const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
        const double norm = std::sqrt(1.0 + phi*phi);
        // Normalised so all vertices lie on the unit sphere
        const double a = 1.0 / norm;
        const double b = phi / norm;

        std::vector<GEO::vec3> verts = {
            { 0,  a,  b}, { 0, -a,  b}, { 0,  a, -b}, { 0, -a, -b},
            { a,  b,  0}, {-a,  b,  0}, { a, -b,  0}, {-a, -b,  0},
            { b,  0,  a}, {-b,  0,  a}, { b,  0, -a}, {-b,  0, -a}
        };
        // Scale to desired radius
        for (auto& v : verts) v = v * radius;

        // 20 triangular faces of icosahedron (vertex indices, CCW outward)
        // Using flat int array to avoid MSVC brace-initialisation issues
        const int face_data[20][3] = {
            {0,8,1},{0,1,9},{0,9,5},{0,5,4},{0,4,8},
            {3,2,11},{3,11,7},{3,7,6},{3,6,10},{3,10,2},
            {1,8,6},{1,6,7},{1,7,9},{9,7,11},{9,11,5},
            {5,11,2},{5,2,4},{4,2,10},{4,10,8},{8,10,6}
        };

        GEO::Logger::out("CurveKID") << "  Icosahedron: 12 vertices, 20 faces, 30 edges" << std::endl;

        // ---- Step 2: Build coarse mesh ----
        GEO::Mesh coarse; coarse.clear();
        coarse.vertices.create_vertices(12);
        for (int i = 0; i < 12; ++i) {
            double* p = coarse.vertices.point_ptr(i);
            p[0]=verts[i].x; p[1]=verts[i].y; p[2]=verts[i].z;
        }
        for (int fi_ = 0; fi_ < 20; ++fi_) {
            GEO::index_t fi = coarse.facets.create_triangle(
                face_data[fi_][0], face_data[fi_][1], face_data[fi_][2]);
            (void)fi;
        }
        coarse.facets.connect();

        // ---- Step 3: Verify all vertices lie on sphere ----
        double max_v_err = 0.0;
        for (int i = 0; i < 12; ++i) {
            double r = GEO::length(verts[i]);
            max_v_err = std::max(max_v_err, std::fabs(r - radius));
        }
        GEO::Logger::out("CurveKID")
            << "  Vertex-on-sphere error: " << max_v_err
            << "  (should be < 1e-14)" << std::endl;

        // ---- Step 4: Assign ONE sphere to ALL 20 faces ----
        // The supporting sphere for every face IS the enclosing sphere.
        // This guarantees K=1 (all panels congruent) and C0 compatibility
        // (all shared edges are great-circle arcs on the same sphere).
        std::vector<CurvedFace> cfaces(20);
        for (int f = 0; f < 20; ++f) {
            CurvedFace& cf = cfaces[f];
            cf.face_id = f;
            cf.v[0] = face_data[f][0];
            cf.v[1] = face_data[f][1];
            cf.v[2] = face_data[f][2];
            cf.sphere_center   = GEO::vec3(0,0,0);
            cf.sphere_radius   = radius;
            cf.sphere_residual = 0.0;
            cf.valid = true;
        }

        // ---- Step 5: Verify K=1 rigid congruence ----
        // All faces are congruent iff all edge lengths are equal.
        // Icosahedron: all edges have the same length.
        double edge_len = GEO::length(verts[face_data[0][1]] - verts[face_data[0][0]]);
        double max_edge_dev = 0.0;
        for (int fi_ = 0; fi_ < 20; ++fi_)
            for (int e = 0; e < 3; ++e) {
                double l = GEO::length(verts[face_data[fi_][(e+1)%3]] - verts[face_data[fi_][e]]);
                max_edge_dev = std::max(max_edge_dev, std::fabs(l - edge_len));
            }
        GEO::Logger::out("CurveKID")
            << "  Edge length: " << edge_len
            << "  max deviation: " << max_edge_dev
            << "  → K=1 rigid congruence: "
            << (max_edge_dev < 1e-10 ? "CONFIRMED" : "FAILED") << std::endl;

        // ---- Step 6: Verify C0 boundary compatibility ----
        // All edges are great-circle arcs on the SAME sphere → C0 exact.
        // Check: for each shared edge (va,vb), both endpoints lie on sphere.
        int n_interior = 0; double max_ep_err = 0.0;
        for (GEO::index_t f = 0; f < coarse.facets.nb(); ++f) {
            for (GEO::index_t le = 0; le < 3; ++le) {
                GEO::index_t adj = coarse.facets.adjacent(f, le);
                if (adj == GEO::NO_FACET || adj <= f) continue;
                ++n_interior;
                GEO::index_t va = coarse.facets.vertex(f, le);
                GEO::index_t vb = coarse.facets.vertex(f, (le+1)%3);
                const double* pa = coarse.vertices.point_ptr(va);
                const double* pb = coarse.vertices.point_ptr(vb);
                double ea = std::fabs(std::sqrt(pa[0]*pa[0]+pa[1]*pa[1]+pa[2]*pa[2]) - radius);
                double eb = std::fabs(std::sqrt(pb[0]*pb[0]+pb[1]*pb[1]+pb[2]*pb[2]) - radius);
                max_ep_err = std::max(max_ep_err, std::max(ea,eb));
            }
        }
        GEO::Logger::out("CurveKID")
            << "  Interior edges: " << n_interior
            << "  max endpoint-on-sphere error: " << max_ep_err
            << "  → C0 compatible: "
            << (max_ep_err < 1e-10 ? "CONFIRMED" : "FAILED") << std::endl;

        // ---- Step 7: Tessellate ----
        std::vector<GEO::index_t> edge_winner;
        resolve_shared_arcs(coarse, cfaces, edge_winner);
        GEO::Mesh viz_mesh;
        tessellate_curved_faces_v2(coarse, cfaces, edge_winner, viz_mesh, subdiv, 1.0,
            TessellationMode::BLENDED_SHARED_EDGE_DEBUG);

        GEO::Logger::out("CurveKID")
            << "  Tessellated: " << viz_mesh.facets.nb()
            << " triangles (" << subdiv << "x subdiv)" << std::endl;

        // ---- Step 8: Fidelity vs dense sphere ----
        // For each tessellated vertex, measure distance to sphere surface.
        double sum_e2 = 0.0; double max_e = 0.0;
        int nv_viz = viz_mesh.vertices.nb();
        for (int v = 0; v < nv_viz; ++v) {
            const double* p = viz_mesh.vertices.point_ptr(v);
            double r = std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
            double e = std::fabs(r - radius);
            sum_e2 += e*e;
            max_e = std::max(max_e, e);
        }
        double rms_e = std::sqrt(sum_e2 / nv_viz);
        GEO::Logger::out("CurveKID")
            << "  Fidelity vs sphere: rms=" << rms_e
            << "  max=" << max_e << std::endl;

        // ---- Step 9: Self-intersection check (simple: check no face has zero area) ----
        int degenerate = 0;
        for (GEO::index_t f = 0; f < viz_mesh.facets.nb(); ++f) {
            GEO::index_t v0=viz_mesh.facets.vertex(f,0),
                         v1=viz_mesh.facets.vertex(f,1),
                         v2=viz_mesh.facets.vertex(f,2);
            const double* p0=viz_mesh.vertices.point_ptr(v0);
            const double* p1=viz_mesh.vertices.point_ptr(v1);
            const double* p2=viz_mesh.vertices.point_ptr(v2);
            GEO::vec3 e1(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2]);
            GEO::vec3 e2(p2[0]-p0[0],p2[1]-p0[1],p2[2]-p0[2]);
            double area = GEO::length(GEO::cross(e1,e2)) * 0.5;
            if (area < 1e-14) ++degenerate;
        }
        GEO::Logger::out("CurveKID")
            << "  Degenerate triangles: " << degenerate
            << "  (self-intersection proxy)" << std::endl;

        // ---- Step 10: Save and display ----
        std::string out_name = "icosahedral_K1_r" + std::to_string((int)radius)
                             + "_s" + std::to_string(subdiv);
        save_mesh_both(viz_mesh, out_name);

        // Add to scene
        MeshGrob* out_grob = MeshGrob::find_or_create(scene_graph(), out_name);
        out_grob->copy(viz_mesh);
        out_grob->update();

        GEO::Logger::out("CurveKID")
            << "=== Icosahedral baseline complete ===" << std::endl
            << "  K=1 congruence: " << (max_edge_dev < 1e-10 ? "YES" : "NO") << std::endl
            << "  C0 compatible:  " << (max_ep_err < 1e-10 ? "YES" : "NO") << std::endl
            << "  Fidelity rms:   " << rms_e << std::endl
            << "  Degenerate:     " << degenerate << std::endl;
    }

    void MeshGrobCurveKIDCommands::build_sphere_ground_truth(
        std::string coarse_mesh_path, int subdiv)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl; return; }
        const GEO::Mesh& original = *grob;

        coarse_mesh_path = normalize_path(coarse_mesh_path);
        if (coarse_mesh_path.empty()) {
            GEO::Logger::err("CurveKID") << "Provide path to coarse mesh." << std::endl; return;
        }
        GEO::Mesh coarse;
        if (!GEO::mesh_load(coarse_mesh_path, coarse)) {
            GEO::Logger::err("CurveKID") << "Cannot load: " << coarse_mesh_path << std::endl; return;
        }
        if (subdiv < 1) subdiv = 1;

        GEO::Logger::out("CurveKID") << std::string(60,'=') << std::endl;
        GEO::Logger::out("CurveKID") << "  GROUND TRUTH: single-sphere curved mesh" << std::endl;
        GEO::Logger::out("CurveKID") << "  Coarse: " << coarse.facets.nb() << " faces, "
                                     << coarse.vertices.nb() << " vertices" << std::endl;

        // --- Step 1: fit ONE sphere to the entire original mesh ---
        const GEO::index_t nv_orig = original.vertices.nb();
        double ATA[4][4]={}, ATb[4]={};
        for (GEO::index_t v=0; v<nv_orig; ++v) {
            const double* p = original.vertices.point_ptr(v);
            double row[4]={2*p[0],2*p[1],2*p[2],1.0};
            double rhs=p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
            for (int i=0;i<4;++i){ ATb[i]+=row[i]*rhs; for(int j=0;j<4;++j) ATA[i][j]+=row[i]*row[j]; }
        }
        double M[4][5];
        for (int i=0;i<4;++i){ for(int j=0;j<4;++j) M[i][j]=ATA[i][j]; M[i][4]=ATb[i]; }
        for (int col=0;col<4;++col) {
            int piv=col;
            for (int r=col+1;r<4;++r) if(std::fabs(M[r][col])>std::fabs(M[piv][col])) piv=r;
            std::swap(M[col],M[piv]);
            double inv=1.0/M[col][col];
            for(int j=col;j<=4;++j) M[col][j]*=inv;
            for(int r=0;r<4;++r){ if(r==col) continue; double f=M[r][col]; for(int j=col;j<=4;++j) M[r][j]-=f*M[col][j]; }
        }
        GEO::vec3 C(M[0][4],M[1][4],M[2][4]);
        double R2=M[3][4]+C.x*C.x+C.y*C.y+C.z*C.z;
        double R=std::sqrt(std::max(0.0,R2));
        GEO::Logger::out("CurveKID") << "  Sphere: C=(" << C.x << "," << C.y << "," << C.z
                                     << ")  R=" << R << std::endl;

        // --- Step 2: project all coarse vertices onto the sphere ---
        const GEO::index_t nv_c = coarse.vertices.nb();
        std::vector<GEO::vec3> proj_verts(nv_c);
        for (GEO::index_t v=0; v<nv_c; ++v) {
            const double* p = coarse.vertices.point_ptr(v);
            GEO::vec3 pt(p[0],p[1],p[2]);
            GEO::vec3 d = pt - C;
            double len = GEO::length(d);
            proj_verts[v] = (len > 1e-14) ? C + d*(R/len) : pt;
        }

        // --- Step 3: build watertight tessellation with shared edges ---
        // For each coarse edge (va,vb), generate (subdiv-1) interior points
        // using slerp. Store by canonical edge key (min,max) with direction flag.
        // Each coarse face then looks up its 3 edges' points and builds sub-triangles.
        GEO::Mesh viz;
        // Add projected coarse vertices first (indices 0..nv_c-1)
        viz.vertices.create_vertices(nv_c);
        for (GEO::index_t v=0; v<nv_c; ++v) {
            double* vp = viz.vertices.point_ptr(v);
            vp[0]=proj_verts[v].x; vp[1]=proj_verts[v].y; vp[2]=proj_verts[v].z;
        }

        // Edge point table: edge (a,b) with a<b → vector of (subdiv-1) point indices
        // stored from a→b direction
        std::map<std::pair<GEO::index_t,GEO::index_t>, std::vector<GEO::index_t>> edge_pts;

        auto get_edge_pts = [&](GEO::index_t va, GEO::index_t vb)
            -> std::vector<GEO::index_t>
        {
            // Returns interior points from va to vb (not including endpoints)
            bool flipped = (va > vb);
            GEO::index_t a = flipped ? vb : va;
            GEO::index_t b = flipped ? va : vb;
            auto key = std::make_pair(a, b);
            auto it = edge_pts.find(key);
            if (it == edge_pts.end()) {
                // Generate interior points via slerp on sphere
                GEO::vec3 pa = proj_verts[a] - C;
                GEO::vec3 pb = proj_verts[b] - C;
                double la = GEO::length(pa), lb = GEO::length(pb);
                std::vector<GEO::index_t> ids;
                for (int s = 1; s < subdiv; ++s) {
                    double t = double(s) / double(subdiv);
                    // Slerp
                    GEO::vec3 pt;
                    if (la > 1e-14 && lb > 1e-14) {
                        GEO::vec3 ua = pa*(1.0/la), ub = pb*(1.0/lb);
                        double dot = std::max(-1.0, std::min(1.0, GEO::dot(ua, ub)));
                        double theta = std::acos(dot);
                        if (std::fabs(theta) < 1e-10) {
                            pt = C + ua*R;
                        } else {
                            GEO::vec3 interp = ua*(std::sin((1-t)*theta)/std::sin(theta))
                                             + ub*(std::sin(t*theta)    /std::sin(theta));
                            pt = C + interp * R;
                        }
                    } else {
                        pt = C + pa*(1-t) + pb*t; // fallback linear
                    }
                    GEO::index_t idx = viz.vertices.create_vertices(1);
                    double* vp = viz.vertices.point_ptr(idx);
                    vp[0]=pt.x; vp[1]=pt.y; vp[2]=pt.z;
                    ids.push_back(idx);
                }
                edge_pts[key] = ids;
                it = edge_pts.find(key);
            }
            std::vector<GEO::index_t> result = it->second;
            if (flipped) std::reverse(result.begin(), result.end());
            return result;
        };

        // For each coarse face: build (subdiv x subdiv) sub-triangles
        const GEO::index_t nf_c = coarse.facets.nb();
        for (GEO::index_t f = 0; f < nf_c; ++f) {
            if (coarse.facets.nb_vertices(f) != 3) continue;
            GEO::index_t cv[3];
            for (int k=0;k<3;++k) cv[k] = coarse.facets.vertex(f,k);

            // Build (subdiv+1) x (subdiv+1) grid of vertex indices
            // Row i, col j: barycentric (i/subdiv towards v2, j/subdiv towards v1)
            // Grid: row 0 = edge v0→v1, col 0 = edge v0→v2, etc.
            // Use triangle barycentric rows:
            //   grid[i][j] where i+j <= subdiv
            //   i=0,j=0 → v0;  i=subdiv,j=0 → v2;  i=0,j=subdiv → v1
            int S = subdiv;
            // Store grid as flat array indexed by (i*(S+1)+j) but only i+j<=S
            // Map (i,j) to a vertex index
            std::vector<std::vector<GEO::index_t>> grid(S+1, std::vector<GEO::index_t>(S+1, GEO::NO_VERTEX));

            // Corners
            grid[0][0] = cv[0];
            grid[0][S] = cv[1];
            grid[S][0] = cv[2];

            // Edge v0→v1 (i=0, j=1..S-1)
            {
                auto ep = get_edge_pts(cv[0], cv[1]);
                for (int j=1;j<S;++j) grid[0][j] = ep[j-1];
            }
            // Edge v0→v2 (j=0, i=1..S-1)
            {
                auto ep = get_edge_pts(cv[0], cv[2]);
                for (int i=1;i<S;++i) grid[i][0] = ep[i-1];
            }
            // Edge v1→v2 (i+j=S, i=1..S-1)
            {
                auto ep = get_edge_pts(cv[1], cv[2]);
                // ep[0] is closest to cv[1], ep[S-2] closest to cv[2]
                for (int i=1;i<S;++i) grid[i][S-i] = ep[i-1];
            }

            // Interior points via slerp barycentric on sphere
            for (int i=1;i<S;++i) {
                for (int j=1;j<S-i;++j) {
                    double bv0 = 1.0 - double(i+j)/S;
                    double bv1 = double(j)/S;
                    double bv2 = double(i)/S;
                    // Spherical barycentric: weighted slerp
                    // Simple: normalise linear combination and project onto sphere
                    GEO::vec3 p = proj_verts[cv[0]]*bv0
                                + proj_verts[cv[1]]*bv1
                                + proj_verts[cv[2]]*bv2;
                    GEO::vec3 d = p - C;
                    double len = GEO::length(d);
                    GEO::vec3 pt = (len>1e-14) ? C + d*(R/len) : p;
                    GEO::index_t idx = viz.vertices.create_vertices(1);
                    double* vp = viz.vertices.point_ptr(idx);
                    vp[0]=pt.x; vp[1]=pt.y; vp[2]=pt.z;
                    grid[i][j] = idx;
                }
            }

            // Build sub-triangles from grid.
            // Grid uses coordinates (i,j) where i=row from v0→v2, j=col from v0→v1.
            // Valid points satisfy i+j <= S.
            // Each quad (i,j),(i,j+1),(i+1,j),(i+1,j+1) gives:
            //   upward triangle:   (i,j), (i,j+1), (i+1,j)       — always valid if i+j+1<=S
            //   downward triangle: (i,j+1),(i+1,j+1),(i+1,j)     — valid only if i+j+2<=S-1 (not on hypotenuse)
            for (int i=0; i<S; ++i) {
                for (int j=0; j<S-i; ++j) {
                    // Upward triangle: always present
                    GEO::index_t a=grid[i][j], b=grid[i][j+1], c2=grid[i+1][j];
                    if (a!=GEO::NO_VERTEX && b!=GEO::NO_VERTEX && c2!=GEO::NO_VERTEX)
                        viz.facets.create_triangle(a, b, c2);
                    // Downward triangle: valid when (i+1)+(j+1) <= S, i.e. i+j <= S-2
                    if (i+j+2 <= S) {
                        GEO::index_t d2=grid[i+1][j+1];
                        if (b!=GEO::NO_VERTEX && d2!=GEO::NO_VERTEX && c2!=GEO::NO_VERTEX)
                            viz.facets.create_triangle(b, d2, c2);
                    }
                }
            }
        }
        viz.facets.connect();

        // --- Step 4: Validate ---
        GEO::Logger::out("CurveKID") << "  --- Validation ---" << std::endl;
        MeshTopologyReport topo = validate_mesh_topology(viz, true);

        // Radial error: all vertices should be at distance R from C
        double max_radial_err = 0.0, sum_radial_err = 0.0;
        GEO::index_t nv_viz = viz.vertices.nb();
        for (GEO::index_t v=0; v<nv_viz; ++v) {
            const double* p = viz.vertices.point_ptr(v);
            double dist = std::sqrt((p[0]-C.x)*(p[0]-C.x)+(p[1]-C.y)*(p[1]-C.y)+(p[2]-C.z)*(p[2]-C.z));
            double err = std::fabs(dist - R);
            sum_radial_err += err; max_radial_err = std::max(max_radial_err, err);
        }
        GEO::Logger::out("CurveKID")
            << "  Radial error: mean=" << sum_radial_err/nv_viz
            << "  max=" << max_radial_err
            << "  (should be ~0)" << std::endl;
        GEO::Logger::out("CurveKID")
            << "  Viz mesh: " << viz.facets.nb() << " triangles, "
            << viz.vertices.nb() << " vertices" << std::endl;

        // Save in .geogram and .obj
        std::string out_base = normalize_path(coarse_mesh_path);
        size_t dot2 = out_base.rfind('.');
        if (dot2 != std::string::npos) out_base = out_base.substr(0, dot2);
        out_base += "_ground_truth";
        save_mesh_both(viz, out_base);
        GEO::Logger::out("CurveKID") << std::string(60,'=') << std::endl;
    }

    // build_curved_mesh_synthetic — circumsphere-per-face, shared arcs
    // =========================================================================

    void MeshGrobCurveKIDCommands::build_curved_mesh_synthetic(
        int target_faces, double t_factor, int subdiv)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh." << std::endl; return; }

        GEO::Logger::out("CurveKID") << "=== build_curved_mesh_synthetic ===" << std::endl;
        GEO::Logger::out("CurveKID")
            << "  target_faces=" << target_faces
            << "  t_factor="     << t_factor
            << "  subdiv="       << subdiv << std::endl;

        // --- Stage 0: per-vertex curvatures (needed for adaptive sizing) ---
        std::vector<VertexCurvature> vcurv;
        compute_vertex_curvatures(*grob, vcurv);
        GEO::Logger::out("CurveKID")
            << "  Curvatures computed for " << vcurv.size() << " vertices." << std::endl;

        // --- Stage 1: coarse curvature-adaptive remesh ---
        GEO::Mesh coarse_mesh;
        GEO::index_t nf_coarse = build_coarse_mesh(
            *grob, vcurv, coarse_mesh,
            GEO::index_t(target_faces), /*H_eps=*/0.01);
        GEO::Logger::out("CurveKID")
            << "  Coarse mesh: " << nf_coarse << " faces, "
            << coarse_mesh.vertices.nb() << " vertices." << std::endl;

        if (nf_coarse == 0) {
            GEO::Logger::err("CurveKID")
                << "  Coarse mesh is empty — remeshing failed." << std::endl;
            return;
        }

        // --- Stage 2 (SYNTHETIC): circumsphere per coarse face ---
        std::vector<CurvedFace> cfaces;
        fit_curved_faces_synthetic(coarse_mesh, cfaces, t_factor);
        {
            int valid_cnt = 0;
            double max_res = 0.0;
            for (const CurvedFace& cf : cfaces) {
                if (cf.valid) {
                    ++valid_cnt;
                    max_res = std::max(max_res, cf.sphere_residual);
                }
            }
            GEO::Logger::out("CurveKID")
                << "  Circumsphere fitting: " << valid_cnt
                << "/" << cfaces.size() << " faces valid." << std::endl;
            GEO::Logger::out("CurveKID")
                << "  Max corner-to-sphere error (should be ~0): "
                << max_res << std::endl;
        }

        // --- Stage 3: resolve shared arcs (pick winning sphere per edge) ---
        std::vector<GEO::index_t> edge_winner;
        resolve_shared_arcs(coarse_mesh, cfaces, edge_winner);
        GEO::Logger::out("CurveKID")
            << "  Shared arcs resolved." << std::endl;

        // --- Stage 4: tessellate with shared-arc awareness ---
        GEO::Mesh viz_mesh;
        tessellate_curved_faces_v2(coarse_mesh, cfaces, edge_winner, viz_mesh, subdiv, 1.0);
        GEO::Logger::out("CurveKID")
            << "  Tessellated viz mesh: " << viz_mesh.facets.nb()
            << " faces, " << viz_mesh.vertices.nb() << " vertices." << std::endl;

        // --- Store sphere radius as face attribute for coloring ---
        {
            GEO::Attribute<double> attr_r(
                viz_mesh.facets.attributes(), "sphere_radius");
            GEO::Attribute<int> attr_id(
                viz_mesh.facets.attributes(), "curved_face_id");
            for (GEO::index_t f = 0; f < viz_mesh.facets.nb(); ++f) {
                if (!attr_id.is_bound()) break;
                int cf_id = attr_id[f];
                if (cf_id >= 0 && cf_id < (int)cfaces.size()) {
                    if (attr_r.is_bound())
                        attr_r[f] = cfaces[cf_id].sphere_radius;
                }
            }
        }

        // --- Save both meshes in .geogram and .obj ---
        save_mesh_both(viz_mesh,   std::string(grob->name()) + "_synthetic");
        save_mesh_both(coarse_mesh, std::string(grob->name()) + "_coarse_synth");
        GEO::Logger::out("CurveKID") << "=== Done ===" << std::endl;
    }

    // =========================================================================
    // build_curved_edges — render great-circle arcs only, no interior fill
    // =========================================================================
    void MeshGrobCurveKIDCommands::build_curved_edges(
        std::string coarse_mesh_path, int arc_pts)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl; return; }
        const GEO::Mesh& original = *grob;

        coarse_mesh_path = normalize_path(coarse_mesh_path);
        if (coarse_mesh_path.empty()) {
            GEO::Logger::err("CurveKID") << "Provide path to coarse mesh." << std::endl; return;
        }
        GEO::Mesh coarse_flat;
        if (!GEO::mesh_load(coarse_mesh_path, coarse_flat)) {
            GEO::Logger::err("CurveKID") << "Cannot load: " << coarse_mesh_path << std::endl; return;
        }
        if (arc_pts < 2) arc_pts = 2;

        // Fit spheres
        std::vector<CurvedFace> cfaces;
        fit_curved_faces(original, coarse_flat, cfaces);
        optimise_joint_spheres(coarse_flat, cfaces, 30);

        // Build edge arc mesh: for each coarse edge, generate arc_pts points
        // via slerp on the owner face's sphere, stored as line segments.
        // Output is a mesh of edges (polylines), saved as .obj (edges) and .geogram.
        GEO::Mesh arc_mesh;

        // Helper: slerp on sphere
        auto slerp_sphere = [](GEO::vec3 C, double R,
                                GEO::vec3 pa, GEO::vec3 pb, double t) -> GEO::vec3 {
            GEO::vec3 da = pa-C, db = pb-C;
            double la = GEO::length(da), lb = GEO::length(db);
            if (la<1e-14||lb<1e-14) return pa*(1-t)+pb*t;
            GEO::vec3 ua=da*(1/la), ub=db*(1/lb);
            double dot=std::max(-1.0,std::min(1.0,GEO::dot(ua,ub)));
            double theta=std::acos(dot);
            if (std::fabs(theta)<1e-10) return C+ua*R;
            GEO::vec3 interp=ua*(std::sin((1-t)*theta)/std::sin(theta))
                            +ub*(std::sin(t*theta)/std::sin(theta));
            return C+interp*R;
        };

        // Build canonical edge set (each edge generated once)
        std::map<std::pair<GEO::index_t,GEO::index_t>, GEO::index_t> edge_owner;
        for (GEO::index_t f=0; f<(GEO::index_t)cfaces.size(); ++f) {
            const CurvedFace& cf = cfaces[f];
            if (!cf.valid) continue;
            for (int e=0; e<3; ++e) {
                GEO::index_t va=cf.v[e], vb=cf.v[(e+1)%3];
                GEO::index_t lo=std::min(va,vb), hi=std::max(va,vb);
                auto key=std::make_pair(lo,hi);
                if (edge_owner.find(key)==edge_owner.end())
                    edge_owner[key]=f;
            }
        }

        for (auto& [key, owner_f] : edge_owner) {
            const CurvedFace& cf = cfaces[owner_f];
            GEO::vec3 pa = cmd_vertex_point(coarse_flat, key.first);
            GEO::vec3 pb = cmd_vertex_point(coarse_flat, key.second);

            // Generate arc_pts points along this arc
            GEO::index_t base = arc_mesh.vertices.create_vertices(arc_pts);
            for (int s=0; s<arc_pts; ++s) {
                double t = double(s)/double(arc_pts-1);
                GEO::vec3 pt;
                if (cf.valid)
                    pt = slerp_sphere(cf.sphere_center, cf.sphere_radius, pa, pb, t);
                else
                    pt = pa*(1-t)+pb*t;
                double* vp = arc_mesh.vertices.point_ptr(base+s);
                vp[0]=pt.x; vp[1]=pt.y; vp[2]=pt.z;
            }
            // Create line segments connecting consecutive arc points
            for (int s=0; s<arc_pts-1; ++s) {
                arc_mesh.edges.create_edge(base+s, base+s+1);
            }
        }

        // Save as .obj and .geogram
        std::string base_name = coarse_mesh_path;
        // strip extension
        auto dot = base_name.rfind('.');
        if (dot != std::string::npos) base_name = base_name.substr(0, dot);
        base_name += "_arcs";

        GEO::Logger::out("CurveKID")
            << "  Curved edges: " << edge_owner.size() << " arcs, "
            << arc_pts << " pts/arc" << std::endl;
        save_mesh_both(arc_mesh, base_name);
    }

    // =========================================================================
    // torus_experiment — four-variant quantitative comparison
    // A: planar coarse  B: independent sphere  C: joint opt  D: error-selected
    // =========================================================================
    void MeshGrobCurveKIDCommands::torus_experiment(
        std::string coarse_mesh_path, int subdiv)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl; return; }
        const GEO::Mesh& original = *grob;

        coarse_mesh_path = normalize_path(coarse_mesh_path);
        if (coarse_mesh_path.empty()) {
            GEO::Logger::err("CurveKID") << "Provide path to coarse mesh." << std::endl; return;
        }
        GEO::Mesh coarse_base;
        if (!GEO::mesh_load(coarse_mesh_path, coarse_base)) {
            GEO::Logger::err("CurveKID") << "Cannot load: " << coarse_mesh_path << std::endl; return;
        }

        std::string base = coarse_mesh_path.substr(0, coarse_mesh_path.rfind('.'));

        // Helper: compute E_planar and E_curved stats
        auto compute_stats = [&](const GEO::Mesh& coarse, const std::vector<CurvedFace>& cfs,
                                  double& mean_p, double& mean_s, int& n_curved, int& n_flat) {
            std::vector<double> ep, es;
            compute_planar_approximation_error(original, coarse, ep);
            // Curved error: for each original vertex, distance to nearest face sphere
            const GEO::index_t nv = original.vertices.nb();
            es.resize(nv, 0.0);
            for (GEO::index_t v = 0; v < nv; ++v) {
                GEO::vec3 pt; { const double* p = original.vertices.point_ptr(v); pt = GEO::vec3(p[0],p[1],p[2]); }
                double best = std::numeric_limits<double>::max();
                for (auto& cf : cfs) {
                    if (!cf.valid) continue;
                    double d = std::fabs(GEO::length(pt - cf.sphere_center) - cf.sphere_radius);
                    best = std::min(best, d*d);
                }
                es[v] = (best < 1e30) ? best : ep[v]*ep[v];
            }
            mean_p = 0; for (auto e : ep) mean_p += e*e; mean_p = std::sqrt(mean_p/ep.size());
            mean_s = 0; for (auto e : es) mean_s += e;   mean_s = std::sqrt(mean_s/es.size());
            n_curved = 0; n_flat = 0;
            for (auto& cf : cfs) { if (cf.valid) ++n_curved; else ++n_flat; }
        };

        GEO::Logger::out("CurveKID") << std::string(60,'=') << std::endl;
        GEO::Logger::out("CurveKID") << "  TORUS EXPERIMENT: four variants from same coarse mesh" << std::endl;
        GEO::Logger::out("CurveKID") << "  Coarse: " << coarse_base.facets.nb() << " faces" << std::endl;
        GEO::Logger::out("CurveKID") << std::string(60,'-') << std::endl;

        // --- Variant A: planar (no sphere fitting) ---
        {
            GEO::Logger::out("CurveKID") << "  Variant A: PLANAR (baseline)" << std::endl;
            std::vector<double> ep;
            compute_planar_approximation_error(original, coarse_base, ep);
            double mean_p=0, rms_p=0, max_p=0;
            for (auto e:ep){mean_p+=e;rms_p+=e*e;max_p=std::max(max_p,e);}
            mean_p/=ep.size(); rms_p=std::sqrt(rms_p/ep.size());
            GEO::Logger::out("CurveKID")
                << "    E_plane: mean=" << mean_p << " rms=" << rms_p << " max=" << max_p << std::endl;
            save_mesh_both(coarse_base, base + "_varA");
        }

        // --- Variant B: independent sphere fitting (no joint opt) ---
        {
            GEO::Logger::out("CurveKID") << "  Variant B: INDEPENDENT sphere fitting" << std::endl;
            GEO::Mesh coarse_B; coarse_B.copy(coarse_base, false);
            std::vector<CurvedFace> cfs_B;
            fit_curved_faces(original, coarse_B, cfs_B);
            std::vector<GEO::index_t> ew_B; resolve_shared_arcs(coarse_B, cfs_B, ew_B);
            GEO::Mesh viz_B;
            tessellate_curved_faces_v2(coarse_B, cfs_B, ew_B, viz_B, subdiv, 1.0);
            double mp,ms; int nc,nf;
            compute_stats(coarse_B, cfs_B, mp, ms, nc, nf);
            GEO::Logger::out("CurveKID")
                << "    Curved=" << nc << " Flat=" << nf
                << "  E_plane_rms=" << mp << "  E_sphere_rms=" << ms
                << "  Improvement=" << (mp>1e-14?(1.0-ms/mp)*100.0:0.0) << "%" << std::endl;
            save_mesh_both(viz_B, base + "_varB");
        }

        // --- Variant C: joint optimisation ---
        {
            GEO::Logger::out("CurveKID") << "  Variant C: JOINT optimisation" << std::endl;
            GEO::Mesh coarse_C; coarse_C.copy(coarse_base, false);
            std::vector<CurvedFace> cfs_C;
            fit_curved_faces(original, coarse_C, cfs_C);
            MeshTopologyReport tr = validate_mesh_topology(coarse_C, false);
            if (tr.boundary_edges == 0)
                optimise_joint_spheres(coarse_C, cfs_C, 30);
            else
                GEO::Logger::out("CurveKID") << "    (open surface — joint opt skipped)" << std::endl;
            std::vector<GEO::index_t> ew_C; resolve_shared_arcs(coarse_C, cfs_C, ew_C);
            GEO::Mesh viz_C;
            tessellate_curved_faces_v2(coarse_C, cfs_C, ew_C, viz_C, subdiv, 1.0);
            double mp,ms; int nc,nf;
            compute_stats(coarse_C, cfs_C, mp, ms, nc, nf);
            GEO::Logger::out("CurveKID")
                << "    Curved=" << nc << " Flat=" << nf
                << "  E_plane_rms=" << mp << "  E_sphere_rms=" << ms
                << "  Improvement=" << (mp>1e-14?(1.0-ms/mp)*100.0:0.0) << "%" << std::endl;
            save_mesh_both(viz_C, base + "_varC");
        }

        // --- Variant D: error-selected (sphere only if E_S < E_P) ---
        {
            GEO::Logger::out("CurveKID") << "  Variant D: ERROR-SELECTED (sphere iff E_S < E_P)" << std::endl;
            GEO::Mesh coarse_D; coarse_D.copy(coarse_base, false);
            std::vector<CurvedFace> cfs_D;
            fit_curved_faces(original, coarse_D, cfs_D);
            // Post-filter: for each valid face, compare E_S vs E_P over 7 samples
            int forced_flat = 0;
            for (auto& cf : cfs_D) {
                if (!cf.valid) continue;
                GEO::vec3 p0 = cmd_vertex_point(coarse_D, cf.v[0]);
                GEO::vec3 p1 = cmd_vertex_point(coarse_D, cf.v[1]);
                GEO::vec3 p2 = cmd_vertex_point(coarse_D, cf.v[2]);
                GEO::vec3 cen = (p0+p1+p2)*(1.0/3.0);
                GEO::vec3 nf2 = GEO::cross(p1-p0, p2-p0);
                double nl = GEO::length(nf2);
                if (nl < 1e-14) continue;
                nf2 = nf2 * (1.0/nl);
                // 7 samples
                GEO::vec3 spts[7] = {p0,p1,p2,cen,(p0+p1)*0.5,(p1+p2)*0.5,(p2+p0)*0.5};
                double e_plane=0, e_sphere=0;
                for (auto& q : spts) {
                    double dp = GEO::dot(q-cen, nf2); e_plane  += dp*dp;
                    double ds = GEO::length(q - cf.sphere_center) - cf.sphere_radius; e_sphere += ds*ds;
                }
                // If sphere not better than plane → force planar
                if (e_sphere >= e_plane) { cf.valid = false; ++forced_flat; }
            }
            GEO::Logger::out("CurveKID") << "    Forced flat by E_S>=E_P: " << forced_flat << std::endl;
            std::vector<GEO::index_t> ew_D; resolve_shared_arcs(coarse_D, cfs_D, ew_D);
            GEO::Mesh viz_D;
            tessellate_curved_faces_v2(coarse_D, cfs_D, ew_D, viz_D, subdiv, 1.0);
            double mp,ms; int nc,nf;
            compute_stats(coarse_D, cfs_D, mp, ms, nc, nf);
            GEO::Logger::out("CurveKID")
                << "    Curved=" << nc << " Flat=" << nf
                << "  E_plane_rms=" << mp << "  E_sphere_rms=" << ms
                << "  Improvement=" << (mp>1e-14?(1.0-ms/mp)*100.0:0.0) << "%" << std::endl;
            save_mesh_both(viz_D, base + "_varD");
        }

        GEO::Logger::out("CurveKID") << std::string(60,'=') << std::endl;
        GEO::Logger::out("CurveKID") << "  Saved: " << base << "_varA/B/C/D.obj" << std::endl;
    }

    // =========================================================================
    // alpha_sweep — E_rendered evaluation for alpha in {0,0.2,0.4,0.6,0.8,1.0}
    // ChatGPT Research Direction 008
    // =========================================================================
    void MeshGrobCurveKIDCommands::alpha_sweep(
        std::string coarse_mesh_path, int subdiv)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl; return; }
        const GEO::Mesh& original = *grob;

        coarse_mesh_path = normalize_path(coarse_mesh_path);
        if (coarse_mesh_path.empty()) {
            GEO::Logger::err("CurveKID") << "Provide path to coarse mesh." << std::endl; return;
        }
        GEO::Mesh coarse_base;
        if (!GEO::mesh_load(coarse_mesh_path, coarse_base)) {
            GEO::Logger::err("CurveKID") << "Cannot load: " << coarse_mesh_path << std::endl; return;
        }

        // Fit spheres once — same for all alpha values
        GEO::Mesh coarse_work; coarse_work.copy(coarse_base, false);
        std::vector<CurvedFace> cfaces;
        fit_curved_faces(original, coarse_work, cfaces);
        {
            MeshTopologyReport tr = validate_mesh_topology(coarse_work, false);
            if (tr.boundary_edges == 0) optimise_joint_spheres(coarse_work, cfaces, 30);
        }
        std::vector<GEO::index_t> ew; resolve_shared_arcs(coarse_work, cfaces, ew);

        // Baseline: flat coarse mesh errors
        std::vector<double> ep;
        compute_planar_approximation_error(original, coarse_base, ep);
        double rms_flat=0, max_flat=0;
        for (auto e:ep){rms_flat+=e*e;max_flat=std::max(max_flat,e);}
        rms_flat=std::sqrt(rms_flat/ep.size());

        std::string base = coarse_mesh_path.substr(0, coarse_mesh_path.rfind('.'));

        GEO::Logger::out("CurveKID") << std::string(60,'=') << std::endl;
        GEO::Logger::out("CurveKID") << "  ALPHA SWEEP — E_rendered evaluation" << std::endl;
        GEO::Logger::out("CurveKID") << "  Coarse: " << coarse_base.facets.nb() << " faces  subdiv=" << subdiv << std::endl;
        GEO::Logger::out("CurveKID") << "  Baseline flat: rms_M2V=" << rms_flat << "  max=" << max_flat << std::endl;
        GEO::Logger::out("CurveKID") << std::string(60,'-') << std::endl;
        GEO::Logger::out("CurveKID")
            << "  alpha | rms_M2V  | max_M2V  | rms_V2M  | max_V2M  | rms_SphDev | max_SphDev" << std::endl;

        double alphas[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
        for (double a : alphas) {
            GEO::Mesh viz;
            tessellate_curved_faces_v2(coarse_work, cfaces, ew, viz, subdiv, a);

            double rms_M2V, max_M2V, rms_V2M, max_V2M, rms_sd, max_sd;
            compute_rendered_errors(original, viz, cfaces,
                rms_M2V, max_M2V, rms_V2M, max_V2M, rms_sd, max_sd);

            GEO::Logger::out("CurveKID")
                << "  " << a
                << " | " << rms_M2V  << " | " << max_M2V
                << " | " << rms_V2M  << " | " << max_V2M
                << " | " << rms_sd   << " | " << max_sd << std::endl;

            // Save best (alpha=0.6) and extremes (0.0 and 1.0)
            if (a == 0.0 || a == 0.6 || a == 1.0) {
                std::string tag = "_alpha" + std::to_string((int)(a*10));
                save_mesh_both(viz, base + tag);
            }
        }
        GEO::Logger::out("CurveKID") << std::string(60,'=') << std::endl;
        GEO::Logger::out("CurveKID")
            << "  Definitions:\n"
            << "    E_M2V     = RMS distance from original vertices to rendered mesh triangles\n"
            << "    E_V2M     = RMS distance from rendered vertices to original mesh triangles\n"
            << "    E_SphDev  = RMS distance of rendered vertices from their supporting spheres\n"
            << "    alpha=0   = flat tessellation  |  alpha=1 = full sphere projection" << std::endl;
    }

    // =========================================================================
    // diagnose_panel_geometry — ChatGPT Direction 018
    // =========================================================================
    // =========================================================================
    // construct_reference_curved_patches — Phase A, ChatGPT Architecture 019
    // Dr. KID style: subdivision + projection onto original dense mesh.
    // No sphere fitting. No bulging. Geometrically valid by construction.
    // =========================================================================
    void MeshGrobCurveKIDCommands::construct_reference_curved_patches(
        std::string coarse_mesh_path, int subdiv,
        int smooth_iters, double smooth_lambda)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh." << std::endl; return; }
        const GEO::Mesh& original = *grob;

        if (coarse_mesh_path.empty()) {
            GEO::Logger::err("CurveKID") << "  coarse_mesh_path required." << std::endl; return;
        }
        GEO::Mesh coarse;
        if (!GEO::mesh_load(coarse_mesh_path, coarse)) {
            GEO::Logger::err("CurveKID") << "  Cannot load: " << coarse_mesh_path << std::endl; return;
        }
        if (subdiv < 1) subdiv = 1;

        GEO::Logger::out("CurveKID") << "=== Phase A: Reference Curved Patches ===" << std::endl;
        GEO::Logger::out("CurveKID")
            << "  M_GT=" << original.facets.nb() << " faces"
            << "  M_coarse=" << coarse.facets.nb() << " faces"
            << "  subdiv=" << subdiv
            << "  smooth=" << smooth_iters << "x lambda=" << smooth_lambda << std::endl;

        const GEO::index_t nf_c = coarse.facets.nb();
        const GEO::index_t nv_c = coarse.vertices.nb();

        // ---- Build closest-point lookup for original mesh ----
        // Pre-cache original mesh triangle data for fast brute-force queries.
        // For typical meshes (< 5000 faces) and ~7200 queries this is fast enough.
        const GEO::index_t nf_orig = original.facets.nb();
        struct OrigTri { GEO::vec3 a, b, c; };
        std::vector<OrigTri> orig_tris(nf_orig);
        for (GEO::index_t f = 0; f < nf_orig; ++f) {
            auto vp = [&](GEO::index_t vi) -> GEO::vec3 {
                const double* p = original.vertices.point_ptr(vi);
                return GEO::vec3(p[0],p[1],p[2]);
            };
            orig_tris[f] = {
                vp(original.facets.vertex(f,0)),
                vp(original.facets.vertex(f,1)),
                vp(original.facets.vertex(f,2))
            };
        }
        // Closest point on a triangle (returns nearest point on triangle to pt)
        auto closest_on_tri = [](GEO::vec3 pt,
                                  GEO::vec3 a, GEO::vec3 b, GEO::vec3 c) -> GEO::vec3 {
            GEO::vec3 ab=b-a, ac=c-a, ap=pt-a;
            double d1=GEO::dot(ab,ap), d2=GEO::dot(ac,ap);
            if (d1<=0 && d2<=0) return a;
            GEO::vec3 bp=pt-b;
            double d3=GEO::dot(ab,bp), d4=GEO::dot(ac,bp);
            if (d3>=0 && d4<=d3) return b;
            GEO::vec3 cp=pt-c;
            double d5=GEO::dot(ab,cp), d6=GEO::dot(ac,cp);
            if (d6>=0 && d5<=d6) return c;
            double vc=d1*d4-d3*d2;
            if (vc<=0 && d1>=0 && d3<=0){
                double v2=d1/(d1-d3); return a+ab*v2;}
            double vb=d5*d2-d1*d6;
            if (vb<=0 && d2>=0 && d6<=0){
                double w2=d2/(d2-d6); return a+ac*w2;}
            double va=d3*d6-d5*d4;
            if (va<=0 && (d4-d3)>=0 && (d5-d6)>=0){
                double w2=(d4-d3)/((d4-d3)+(d5-d6));
                return b+(c-b)*w2;}
            double denom=1.0/(va+vb+vc);
            double v3=vb*denom, w3=vc*denom;
            return a+ab*v3+ac*w3;
        };

        // Project a point onto the original dense mesh (brute-force closest point)
        // Defined here (before Step 1) so it can be used to snap coarse corners.
        auto project_to_gt = [&](GEO::vec3 pt) -> GEO::vec3 {
            double best_d2 = 1e30;
            GEO::vec3 best_pt = pt;
            for (auto& tri : orig_tris) {
                GEO::vec3 cp = closest_on_tri(pt, tri.a, tri.b, tri.c);
                GEO::vec3 diff = cp - pt;
                double d2 = GEO::dot(diff,diff);
                if (d2 < best_d2) { best_d2=d2; best_pt=cp; }
            }
            return best_pt;
        };

        // ---- Output mesh ----
        GEO::Mesh out; out.clear();

        // Attributes: parent_face_id per vertex and per facet
        GEO::Attribute<int> vattr_parent; // set after vertices created
        GEO::Attribute<int> fattr_parent; // set after facets created

        // ---- Step 0: Validate coarse mesh topology ----
        {
            coarse.facets.connect();
            int n_be = 0, n_nme = 0;
            for (GEO::index_t f = 0; f < coarse.facets.nb(); ++f)
                for (GEO::index_t le = 0; le < 3; ++le)
                    if (coarse.facets.adjacent(f,le) == GEO::NO_FACET) ++n_be;
            int Vc = coarse.vertices.nb();
            int Fc = coarse.facets.nb();
            int Ec = (3*Fc + n_be) / 2;
            int chi_c = Vc - Ec + Fc;
            GEO::Logger::out("CurveKID")
                << "  Coarse mesh: V=" << Vc << " F=" << Fc << " E=" << Ec
                << " Euler=" << chi_c << " BoundaryEdges=" << n_be << std::endl;
            if (n_be > 0)
                GEO::Logger::out("CurveKID")
                    << "  WARNING: coarse mesh has " << n_be
                    << " boundary edges — open surface." << std::endl;
        }

        // ---- Step 1: Copy coarse vertices, snapping them onto M_GT ----
        // CVT remeshing may leave coarse vertices slightly off the surface.
        // Project each coarse vertex onto the nearest point on M_GT so that
        // corners are exactly on the surface (same as edge/interior vertices).
        // These vertices are FIXED during smoothing (never moved after this).
        out.vertices.create_vertices(nv_c);
        for (GEO::index_t v = 0; v < nv_c; ++v) {
            const double* p = coarse.vertices.point_ptr(v);
            GEO::vec3 pt(p[0],p[1],p[2]);
            GEO::vec3 pt_snapped = project_to_gt(pt);
            double* q = out.vertices.point_ptr(v);
            q[0]=pt_snapped.x; q[1]=pt_snapped.y; q[2]=pt_snapped.z;
        }

        // ---- Step 2: Generate shared edge subdivision vertices ----
        // For coarse edge (lo,hi) with lo<hi: (subdiv-1) interior points
        // stored once, shared by both adjacent faces.
        std::map<std::pair<GEO::index_t,GEO::index_t>,
                 std::vector<GEO::index_t>> edge_verts;

        // Build edge→face map to find owner for each edge
        std::map<std::pair<GEO::index_t,GEO::index_t>, GEO::index_t> edge_to_face_c;
        for (GEO::index_t f = 0; f < nf_c; ++f) {
            for (int le = 0; le < 3; ++le) {
                GEO::index_t va = coarse.facets.vertex(f,le);
                GEO::index_t vb = coarse.facets.vertex(f,(le+1)%3);
                auto key = std::make_pair(std::min(va,vb), std::max(va,vb));
                if (edge_to_face_c.find(key) == edge_to_face_c.end())
                    edge_to_face_c[key] = f;
            }
        }

        // Build set of coarse boundary edges (edges with only one adjacent face).
        // Boundary edge subdivision vertices must NOT be projected onto M_GT —
        // they lie on the open boundary of the mesh and projection would snap
        // them to interior triangles, causing the border to shrink inward.
        // Instead, boundary edges use pure linear interpolation between the
        // (already-snapped) corner endpoints.
        std::set<std::pair<GEO::index_t,GEO::index_t>> coarse_boundary_edges;
        {
            std::map<std::pair<GEO::index_t,GEO::index_t>,int> edge_count;
            for (GEO::index_t f = 0; f < nf_c; ++f) {
                for (int le = 0; le < 3; ++le) {
                    GEO::index_t ea = coarse.facets.vertex(f,le);
                    GEO::index_t eb = coarse.facets.vertex(f,(le+1)%3);
                    auto key = std::make_pair(std::min(ea,eb), std::max(ea,eb));
                    edge_count[key]++;
                }
            }
            for (auto& kv : edge_count)
                if (kv.second == 1) coarse_boundary_edges.insert(kv.first);
        }
        GEO::Logger::out("CurveKID")
            << "  Coarse boundary edges (pinned, no projection): "
            << coarse_boundary_edges.size() << std::endl;

        // For each edge, generate (subdiv-1) interior points by lerp + project.
        // Uses snapped corner positions (from out, not raw coarse) for consistency.
        // Boundary edges: lerp only (no projection) to preserve open boundaries.
        auto get_edge_verts_phase_a = [&](GEO::index_t va, GEO::index_t vb)
            -> std::vector<GEO::index_t>
        {
            bool flipped = (va > vb);
            GEO::index_t lo = std::min(va,vb), hi = std::max(va,vb);
            auto key = std::make_pair(lo,hi);
            auto it = edge_verts.find(key);
            if (it == edge_verts.end()) {
                // Use SNAPPED positions (stored in out) not raw coarse positions
                const double* pa = out.vertices.point_ptr(lo);
                const double* pb = out.vertices.point_ptr(hi);
                GEO::vec3 vpa(pa[0],pa[1],pa[2]), vpb(pb[0],pb[1],pb[2]);
                bool is_boundary = (coarse_boundary_edges.count(key) > 0);
                std::vector<GEO::index_t> ids;
                for (int s = 1; s < subdiv; ++s) {
                    double t = double(s)/double(subdiv);
                    GEO::vec3 pt_lerp = vpa*(1.0-t) + vpb*t;
                    // Boundary edges: keep lerp exactly (no projection).
                    // Interior edges: project onto M_GT surface.
                    GEO::vec3 pt_final = is_boundary ? pt_lerp : project_to_gt(pt_lerp);
                    GEO::index_t nv = out.vertices.create_vertices(1);
                    double* vp = out.vertices.point_ptr(nv);
                    vp[0]=pt_final.x; vp[1]=pt_final.y; vp[2]=pt_final.z;
                    ids.push_back(nv);
                }
                edge_verts[key] = ids;
                it = edge_verts.find(key);
            }
            std::vector<GEO::index_t> result = it->second;
            if (flipped) std::reverse(result.begin(), result.end());
            return result;
        };

        // ---- Step 3: Tessellate each coarse face ----
        // The hypotenuse assignment uses ev[S-1-i] (not ev[i-1]) so that
        // the hypotenuse vertices are stored in FORWARD order from v1 to v2.
        // This ensures consistent child-edge creation with adjacent faces
        // regardless of which slot (axis j=0, axis i=0, or hypotenuse)
        // the shared coarse edge occupies in each face.
        // See tet_fix2.cpp for the proof-of-correctness test.
        std::vector<int> face_parents; // parent face id per output triangle

        int S = subdiv;
        for (GEO::index_t f = 0; f < nf_c; ++f) {
            GEO::index_t v0 = coarse.facets.vertex(f,0);
            GEO::index_t v1 = coarse.facets.vertex(f,1);
            GEO::index_t v2 = coarse.facets.vertex(f,2);

            const double* p0_ = coarse.vertices.point_ptr(v0);
            const double* p1_ = coarse.vertices.point_ptr(v1);
            const double* p2_ = coarse.vertices.point_ptr(v2);
            GEO::vec3 p0(p0_[0],p0_[1],p0_[2]);
            GEO::vec3 p1(p1_[0],p1_[1],p1_[2]);
            GEO::vec3 p2(p2_[0],p2_[1],p2_[2]);

            // Grid: grid[i][j] = output vertex index, i+j <= S
            // (i,j,S-i-j) = barycentric coords * S
            // corner: (0,0)=v0, (S,0)=v1, (0,S)=v2
            std::vector<std::vector<GEO::index_t>> grid(
                S+1, std::vector<GEO::index_t>(S+1, GEO::NO_VERTEX));

            // Corners
            grid[0][0] = v0;
            grid[S][0] = v1;
            grid[0][S] = v2;

            // Edge v0→v1 (j=0, i=1..S-1)
            // get_edge_verts returns vertices in direction v0→v1
            // (reversed if v0 > v1). grid[i][0] corresponds to
            // parameter t=i/S moving from v0 (i=0) toward v1 (i=S).
            { auto ev = get_edge_verts_phase_a(v0,v1);
              // ev[k] = vertex at parameter (k+1)/S from v0 toward v1
              for (int i=1;i<S;++i) grid[i][0] = ev[i-1]; }

            // Edge v0→v2 (i=0, j=1..S-1)
            // grid[0][j] = parameter j/S from v0 toward v2
            { auto ev = get_edge_verts_phase_a(v0,v2);
              for (int j=1;j<S;++j) grid[0][j] = ev[j-1]; }

            // Edge v1→v2 (hypotenuse i+j=S):
            // The hypotenuse path from v1 (at g[S][0]) to v2 (at g[0][S])
            // visits g[S-k][k] for k=1..S-1 in order of increasing distance from v1.
            // = g[S-1][1], g[S-2][2], ..., g[1][S-1]
            // get_edge_verts_phase_a(v1,v2) returns ev[k] = vertex at t=(k+1)/S from v1.
            // So ev[0] is closest to v1, ev[S-2] is closest to v2.
            // Mapping: g[S-k][k] = ev[k-1] for k=1..S-1
            // i.e. g[i][S-i] = ev[S-1-i] for i=1..S-1  (FIX: reversed vs naive i-1)
            { auto ev = get_edge_verts_phase_a(v1,v2);
              for (int i=1;i<S;++i) grid[i][S-i] = ev[S-1-i]; }

            // Interior points: barycentric lerp + project
            for (int i=1;i<S;++i) {
                for (int j=1;j<S-i;++j) {
                    double b0 = 1.0 - double(i+j)/S;
                    double b1 = double(i)/S;
                    double b2 = double(j)/S;
                    GEO::vec3 pt_lerp = p0*b0 + p1*b1 + p2*b2;
                    GEO::vec3 pt_proj = project_to_gt(pt_lerp);
                    GEO::index_t nv = out.vertices.create_vertices(1);
                    double* vp = out.vertices.point_ptr(nv);
                    vp[0]=pt_proj.x; vp[1]=pt_proj.y; vp[2]=pt_proj.z;
                    grid[i][j] = nv;
                }
            }

            // Emit triangles using standard barycentric subdivision pattern.
            // For each cell (i,j) with i+j < S:
            //   Upward triangle:   grid[i][j], grid[i+1][j], grid[i][j+1]
            //   Downward triangle: grid[i+1][j], grid[i+1][j+1], grid[i][j+1]
            //     (only when i+j+2 <= S, i.e. grid[i+1][j+1] exists)
            // This covers all S^2 sub-triangles with no gaps.
            for (int i=0;i<S;++i) {
                for (int j=0;j<S-i;++j) {
                    // Upward triangle (always exists)
                    GEO::index_t a=grid[i][j], b=grid[i+1][j], c=grid[i][j+1];
                    if (a!=GEO::NO_VERTEX && b!=GEO::NO_VERTEX && c!=GEO::NO_VERTEX) {
                        out.facets.create_triangle(a,b,c);
                        face_parents.push_back((int)f);
                    }
                    // Downward triangle (exists when i+j+1 < S)
                    if (i+j+1 < S) {
                        GEO::index_t d=grid[i+1][j], e=grid[i+1][j+1], ff2=grid[i][j+1];
                        if (d!=GEO::NO_VERTEX && e!=GEO::NO_VERTEX && ff2!=GEO::NO_VERTEX) {
                            out.facets.create_triangle(d,e,ff2);
                            face_parents.push_back((int)f);
                        } else if (f < 2) { // debug first 2 faces only
                            GEO::Logger::out("CurveKID")
                                << "  SKIP downward f=" << f
                                << " i=" << i << " j=" << j
                                << " d=" << (d==GEO::NO_VERTEX?"NONE":std::to_string(d))
                                << " e=" << (e==GEO::NO_VERTEX?"NONE":std::to_string(e))
                                << " ff2=" << (ff2==GEO::NO_VERTEX?"NONE":std::to_string(ff2))
                                << std::endl;
                        }
                    }
                }
            }
        }

        // ---- Raw index buffer edge count (before connect) ----
        {
            std::map<std::pair<GEO::index_t,GEO::index_t>,int> raw_edge_count;
            int degenerate_tris = 0;
            int degenerate_edges = 0;
            for (GEO::index_t f = 0; f < out.facets.nb(); ++f) {
                GEO::index_t v0=out.facets.vertex(f,0);
                GEO::index_t v1=out.facets.vertex(f,1);
                GEO::index_t v2=out.facets.vertex(f,2);
                if (v0==v1||v1==v2||v0==v2) { ++degenerate_tris; continue; }
                for (int le=0;le<3;++le) {
                    GEO::index_t va=out.facets.vertex(f,le);
                    GEO::index_t vb=out.facets.vertex(f,(le+1)%3);
                    if (va==vb) { ++degenerate_edges; continue; }
                    raw_edge_count[{std::min(va,vb),std::max(va,vb)}]++;
                }
            }
            int r1=0,r2=0,rgt2=0;
            for (auto& kv:raw_edge_count)
                if(kv.second==1)++r1; else if(kv.second==2)++r2; else ++rgt2;
            GEO::Logger::out("CurveKID")
                << "  Raw index buffer: total_edges=" << raw_edge_count.size()
                << "  1-face=" << r1 << "  2-face=" << r2 << "  >2=" << rgt2
                << "  degenerate_tris=" << degenerate_tris
                << "  degenerate_edges=" << degenerate_edges << std::endl;

            // Print first boundary edge details
            int printed = 0;
            for (GEO::index_t f = 0; f < out.facets.nb() && printed < 3; ++f) {
                GEO::index_t v0=out.facets.vertex(f,0);
                GEO::index_t v1=out.facets.vertex(f,1);
                GEO::index_t v2=out.facets.vertex(f,2);
                if (v0==v1||v1==v2||v0==v2) continue;
                for (int le=0;le<3;++le) {
                    GEO::index_t va=out.facets.vertex(f,le);
                    GEO::index_t vb=out.facets.vertex(f,(le+1)%3);
                    auto key=std::make_pair(std::min(va,vb),std::max(va,vb));
                    if (raw_edge_count[key]==1 && printed < 3) {
                        // Get parent face info
                        int pf = (f < face_parents.size()) ? face_parents[f] : -1;
                        GEO::Logger::out("CurveKID")
                            << "  Boundary edge: f=" << f << " parent=" << pf
                            << " va=" << va << " vb=" << vb
                            << " va<nv_c=" << (va<nv_c) << " vb<nv_c=" << (vb<nv_c)
                            << std::endl;
                        ++printed;
                    }
                }
            }
        }

        // ---- Pre-connect diagnostic: check raw triangle index buffer ----
        // Print one shared coarse edge from both adjacent faces
        {
            bool found = false;
            for (GEO::index_t f = 0; f < nf_c && !found; ++f) {
                for (GEO::index_t le = 0; le < 3 && !found; ++le) {
                    GEO::index_t adj = coarse.facets.adjacent(f, le);
                    if (adj == GEO::NO_FACET || adj <= f) continue;
                    // Found a shared coarse edge between face f and face adj
                    GEO::index_t va = coarse.facets.vertex(f, le);
                    GEO::index_t vb = coarse.facets.vertex(f, (le+1)%3);
                    // Get edge verts from both sides
                    auto ev_f   = get_edge_verts_phase_a(va, vb);
                    // Find the same edge from adj face
                    GEO::index_t va_adj=GEO::NO_VERTEX, vb_adj=GEO::NO_VERTEX;
                    for (GEO::index_t le2=0; le2<3; ++le2) {
                        GEO::index_t pa = coarse.facets.vertex(adj, le2);
                        GEO::index_t pb = coarse.facets.vertex(adj, (le2+1)%3);
                        if ((pa==va&&pb==vb)||(pa==vb&&pb==va)) {
                            va_adj=pa; vb_adj=pb; break;
                        }
                    }
                    auto ev_adj = get_edge_verts_phase_a(va_adj, vb_adj);
                    GEO::Logger::out("CurveKID")
                        << "  DEBUG shared edge coarse(" << va << "," << vb << "):" << std::endl
                        << "    face " << f   << " dir(" << va << "->" << vb << "): [" << va;
                    for (auto v : ev_f)   GEO::Logger::out("CurveKID") << "," << v;
                    GEO::Logger::out("CurveKID") << "," << vb << "]" << std::endl;
                    GEO::Logger::out("CurveKID")
                        << "    face " << adj << " dir(" << va_adj << "->" << vb_adj << "): [" << va_adj;
                    for (auto v : ev_adj) GEO::Logger::out("CurveKID") << "," << v;
                    GEO::Logger::out("CurveKID") << "," << vb_adj << "]" << std::endl;
                    found = true;
                }
            }
        }

        out.facets.connect();

        // ---- Edge incidence diagnostic (ChatGPT Direction 023) ----
        {
            // Build canonical edge → face count map
            std::map<std::pair<GEO::index_t,GEO::index_t>, int> edge_count;
            for (GEO::index_t f = 0; f < out.facets.nb(); ++f) {
                for (int le = 0; le < 3; ++le) {
                    GEO::index_t va = out.facets.vertex(f, le);
                    GEO::index_t vb = out.facets.vertex(f, (le+1)%3);
                    auto key = std::make_pair(std::min(va,vb), std::max(va,vb));
                    edge_count[key]++;
                }
            }
            int n1=0, n2=0, ngt2=0;
            int boundary_in_patch=0, boundary_on_coarse=0;
            for (auto& kv : edge_count) {
                if      (kv.second == 1) ++n1;
                else if (kv.second == 2) ++n2;
                else                     ++ngt2;
            }
            GEO::Logger::out("CurveKID")
                << "  Edge incidence histogram:"
                << "  1-face(boundary)=" << n1
                << "  2-face(manifold)=" << n2
                << "  >2-face(nonmanifold)=" << ngt2 << std::endl;
            // Check: are boundary edges along coarse edges or inside patches?
            // A boundary edge whose both endpoints have index < nv_c is a coarse edge
            for (auto& kv : edge_count) {
                if (kv.second != 1) continue;
                bool va_coarse = (kv.first.first  < nv_c);
                bool vb_coarse = (kv.first.second < nv_c);
                // Edge vertex (index >= nv_c) could be a coarse-edge interior vertex
                // Check if it's in edge_verts
                // For now just categorise by endpoint types
                if (va_coarse || vb_coarse) ++boundary_on_coarse;
                else ++boundary_in_patch;
            }
            GEO::Logger::out("CurveKID")
                << "  Boundary edge location:"
                << "  on/near-coarse-vertex=" << boundary_on_coarse
                << "  interior=" << boundary_in_patch << std::endl;
        }

        // Set parent_face_id attribute
        GEO::Attribute<int> fattr(out.facets.attributes(), "parent_face_id");
        for (GEO::index_t f = 0; f < out.facets.nb(); ++f)
            fattr[f] = (f < face_parents.size()) ? face_parents[f] : -1;

        GEO::Logger::out("CurveKID")
            << "  Tessellated: " << out.vertices.nb() << " vertices, "
            << out.facets.nb() << " triangles" << std::endl;

        // ---- Step 4: Laplacian smoothing (boundary-constrained) ----
        // Build vertex adjacency. Shared-edge vertices and coarse-corner
        // vertices are "boundary" in the sense that they are shared between
        // patches — do NOT move them (preserve shared edges exactly).
        // Only move pure interior vertices of each patch.
        std::set<GEO::index_t> fixed_verts;
        // Fix all coarse vertices (corners) — snapped in Step 1, never move again.
        for (GEO::index_t v = 0; v < nv_c; ++v) fixed_verts.insert(v);
        // Fix all shared-edge subdivision vertices — generated once, shared by
        // adjacent patches; moving them would break seam continuity.
        // (This includes boundary edge verts which are pure-lerp and must stay put.)
        for (auto& kv : edge_verts)
            for (auto vi : kv.second) fixed_verts.insert(vi);

        // Build vertex → neighbour map from output mesh triangles
        int nv_out = out.vertices.nb();
        std::vector<std::vector<GEO::index_t>> v2v(nv_out);
        for (GEO::index_t f = 0; f < out.facets.nb(); ++f) {
            GEO::index_t a=out.facets.vertex(f,0),
                         b=out.facets.vertex(f,1),
                         c=out.facets.vertex(f,2);
            v2v[a].push_back(b); v2v[a].push_back(c);
            v2v[b].push_back(a); v2v[b].push_back(c);
            v2v[c].push_back(a); v2v[c].push_back(b);
        }

        for (int iter = 0; iter < smooth_iters; ++iter) {
            std::vector<GEO::vec3> new_pos(nv_out);
            for (int v = 0; v < nv_out; ++v) {
                const double* pp = out.vertices.point_ptr(v);
                new_pos[v] = GEO::vec3(pp[0],pp[1],pp[2]);
            }
            for (int v = 0; v < nv_out; ++v) {
                if (fixed_verts.count(v)) continue;
                if (v2v[v].empty()) continue;
                // Mean of neighbours
                GEO::vec3 mean(0,0,0);
                std::set<GEO::index_t> uniq(v2v[v].begin(), v2v[v].end());
                for (auto nb : uniq) {
                    const double* pp = out.vertices.point_ptr(nb);
                    mean = mean + GEO::vec3(pp[0],pp[1],pp[2]);
                }
                mean = mean * (1.0/uniq.size());
                const double* pp = out.vertices.point_ptr(v);
                GEO::vec3 cur(pp[0],pp[1],pp[2]);
                GEO::vec3 smoothed = cur*(1.0-smooth_lambda) + mean*smooth_lambda;
                // Reproject onto M_GT
                new_pos[v] = project_to_gt(smoothed);
            }
            for (int v = 0; v < nv_out; ++v) {
                double* pp = out.vertices.point_ptr(v);
                pp[0]=new_pos[v].x; pp[1]=new_pos[v].y; pp[2]=new_pos[v].z;
            }
        }

        // ---- Step 5: Validate ----
        int degenerate=0, inverted=0;
        double sum_e2=0, max_e=0;
        int nf_out = out.facets.nb();
        for (GEO::index_t f = 0; f < (GEO::index_t)nf_out; ++f) {
            GEO::index_t va=out.facets.vertex(f,0),
                         vb=out.facets.vertex(f,1),
                         vc=out.facets.vertex(f,2);
            const double* a_=out.vertices.point_ptr(va);
            const double* b_=out.vertices.point_ptr(vb);
            const double* c_=out.vertices.point_ptr(vc);
            GEO::vec3 a(a_[0],a_[1],a_[2]),b(b_[0],b_[1],b_[2]),c(c_[0],c_[1],c_[2]);
            double area = GEO::length(GEO::cross(b-a,c-a))*0.5;
            if (area < 1e-14) ++degenerate;
        }
        // Fidelity: for each output vertex, distance to M_GT
        for (int v = 0; v < nv_out; ++v) {
            const double* pp = out.vertices.point_ptr(v);
            GEO::vec3 pt(pp[0],pp[1],pp[2]);
            GEO::vec3 nearest = project_to_gt(pt);
            double e = GEO::length(nearest - pt);
            sum_e2 += e*e;
            max_e = std::max(max_e, e);
        }
        double rms_e = std::sqrt(sum_e2/nv_out);

        GEO::Logger::out("CurveKID")
            << "  Fidelity vs M_GT: rms=" << rms_e << "  max=" << max_e << std::endl;
        GEO::Logger::out("CurveKID")
            << "  Degenerate triangles: " << degenerate << std::endl;
        GEO::Logger::out("CurveKID")
            << "  Shared edges: " << edge_verts.size()
            << "  (each generated once, shared by both patches)" << std::endl;

        // ---- Step 6: Mesh validity stats ----
        {
            // Euler characteristic, boundary edges, components
            int n_boundary_e = 0;
            out.facets.connect();
            for (GEO::index_t f = 0; f < out.facets.nb(); ++f)
                for (GEO::index_t le = 0; le < 3; ++le)
                    if (out.facets.adjacent(f,le) == GEO::NO_FACET) ++n_boundary_e;
            int V = out.vertices.nb();
            int F = out.facets.nb();
            int E = 0;
            // E = 3*F/2 + boundary/2 for a manifold mesh
            // For manifold: E = (3F + boundary_e) / 2
            E = (3*F + n_boundary_e) / 2;
            int chi = V - E + F;
            GEO::Logger::out("CurveKID")
                << "  Mesh validity: V=" << V << " E=" << E << " F=" << F
                << " Euler=" << chi << " BoundaryEdges=" << n_boundary_e << std::endl;
        }

        // ---- Step 7: Save as .geogram (preserves attributes) ----
        std::string out_name = "ref_curved_patches_s" + std::to_string(subdiv);
        // Save as geogram to preserve parent_face_id attribute
        std::string geogram_path = out_name + ".geogram";
        GEO::mesh_save(out, geogram_path);
        GEO::Logger::out("CurveKID") << "  Saved: " << geogram_path
            << "  (attributes preserved)" << std::endl;
        // Also save .obj for compatibility
        GEO::mesh_save(out, out_name + ".obj");
        GEO::Logger::out("CurveKID") << "  Saved: " << out_name << ".obj" << std::endl;

        // ---- Step 8: Save authoritative coarse-edge vertex map ----
        // Per ChatGPT Decision 2: store each coarse edge's ordered subdivision
        // vertex IDs explicitly so build_phase_a_edges can read them directly
        // without geometric reconstruction.
        // Format: one line per coarse edge:
        //   lo hi v0 v1 ... v_{S-2}
        // where lo<hi are coarse vertex indices, v0..v_{S-2} are output mesh
        // vertex indices in order from lo toward hi (lo and hi are NOT repeated).
        {
            std::string emap_path = out_name + ".edgemap.txt";
            FILE* fp = fopen(emap_path.c_str(), "w");
            if (fp) {
                fprintf(fp, "# CurveKID Phase A edge map: lo hi v0 v1 ... v_{S-2}\n");
                fprintf(fp, "# S=%d  nv_c=%u  nf_c=%u\n", S, (unsigned)nv_c, (unsigned)nf_c);
                for (auto& kv : edge_verts) {
                    fprintf(fp, "%u %u", (unsigned)kv.first.first,
                                        (unsigned)kv.first.second);
                    for (auto vi : kv.second) fprintf(fp, " %u", (unsigned)vi);
                    fprintf(fp, "\n");
                }
                fclose(fp);
                GEO::Logger::out("CurveKID")
                    << "  Saved edge map: " << emap_path
                    << "  (" << edge_verts.size() << " coarse edges)" << std::endl;
            }
        }

        MeshGrob* out_grob = MeshGrob::find_or_create(scene_graph(), out_name);
        out_grob->copy(out);
        out_grob->update();

        GEO::Logger::out("CurveKID")
            << "=== Phase A complete: " << out_name << " ===" << std::endl
            << "  Paint by 'parent_face_id' for patch colour map." << std::endl
            << "  Load .geogram file to preserve attribute." << std::endl;
    }

    // =========================================================================
    // classify_phase_a_patches — Phase B, ChatGPT Direction 021
    // =========================================================================
    void MeshGrobCurveKIDCommands::classify_phase_a_patches(
        std::string phase_a_mesh_path, std::string coarse_mesh_path, int K)
    {
        using namespace CurveKIDAlgo;

        std::string phase_a_mesh = phase_a_mesh_path;
        if (phase_a_mesh.empty() || coarse_mesh_path.empty()) {
            GEO::Logger::err("CurveKID")
                << "  Both phase_a_mesh and coarse_mesh paths required." << std::endl;
            return;
        }

        // Load Phase A mesh (with parent_face_id attribute)
        GEO::Mesh pa_mesh;
        if (!GEO::mesh_load(phase_a_mesh, pa_mesh)) {
            GEO::Logger::err("CurveKID") << "  Cannot load: " << phase_a_mesh << std::endl;
            return;
        }
        // Load coarse mesh
        GEO::Mesh coarse;
        if (!GEO::mesh_load(coarse_mesh_path, coarse)) {
            GEO::Logger::err("CurveKID") << "  Cannot load: " << coarse_mesh_path << std::endl;
            return;
        }

        const GEO::index_t nf_pa = pa_mesh.facets.nb();
        const GEO::index_t nf_c  = coarse.facets.nb();

        GEO::Logger::out("CurveKID") << "=== Phase B: Classify Phase A Patches ===" << std::endl;
        GEO::Logger::out("CurveKID")
            << "  Phase A mesh: " << nf_pa << " triangles"
            << "  Coarse: " << nf_c << " faces  K=" << K << std::endl;

        // Read parent_face_id attribute — must already exist in the file.
        // GEO::Attribute constructor CREATES the attribute if absent, so we
        // must check is_defined() BEFORE constructing it.
        if (!pa_mesh.facets.attributes().is_defined("parent_face_id")) {
            GEO::Logger::err("CurveKID")
                << "  parent_face_id not found in " << phase_a_mesh
                << ". You must supply the .geogram file (not .obj) — "
                << ".obj does not preserve attributes." << std::endl;
            return;
        }
        GEO::Attribute<int> fattr_parent(pa_mesh.facets.attributes(), "parent_face_id");

        // ---- Compute per-patch geometric descriptor ----
        // For each coarse face f (patch), compute:
        //   d[0..2] = sorted coarse edge lengths (shape descriptor, rigid-invariant)
        //   d[3]    = total patch area (sum of tessellated triangle areas)
        //   d[4]    = mean normal variation (curvature proxy)
        //   d[5]    = aspect ratio (max_edge / min_edge)
        // All descriptors normalised to [0,1] range after computing all patches.

        const int DIM = 6;
        // Use plain struct instead of std::array<double,6> for MSVC compatibility
        struct Desc6 {
            double v[6];
            double& operator[](int i){ return v[i]; }
            double  operator[](int i) const { return v[i]; }
            void fill(double val){ for(int i=0;i<6;++i) v[i]=val; }
        };
        std::vector<Desc6> descriptors(nf_c);
        for (auto& d : descriptors) d.fill(0.0);

        // Coarse edge lengths
        for (GEO::index_t f = 0; f < nf_c; ++f) {
            GEO::vec3 p[3];
            for (int k=0;k<3;++k) {
                const double* pp = coarse.vertices.point_ptr(coarse.facets.vertex(f,k));
                p[k] = GEO::vec3(pp[0],pp[1],pp[2]);
            }
            double e[3] = {
                GEO::length(p[1]-p[0]),
                GEO::length(p[2]-p[1]),
                GEO::length(p[0]-p[2])
            };
            std::sort(e, e+3);
            descriptors[f][0] = e[0];
            descriptors[f][1] = e[1];
            descriptors[f][2] = e[2];
            descriptors[f][3] = 0.0; // area (accumulated below)
            descriptors[f][4] = 0.0; // normal variation
            descriptors[f][5] = (e[0] > 1e-14) ? e[2]/e[0] : 1.0; // aspect ratio
        }

        // Accumulate patch area and normal variation from tessellated triangles
        std::vector<GEO::vec3> patch_mean_normal(nf_c, GEO::vec3(0,0,0));
        for (GEO::index_t f = 0; f < nf_pa; ++f) {
            int cf_id = fattr_parent[f];
            if (cf_id < 0 || cf_id >= (int)nf_c) continue;
            GEO::index_t v0=pa_mesh.facets.vertex(f,0),
                         v1=pa_mesh.facets.vertex(f,1),
                         v2=pa_mesh.facets.vertex(f,2);
            const double* a_=pa_mesh.vertices.point_ptr(v0);
            const double* b_=pa_mesh.vertices.point_ptr(v1);
            const double* c_=pa_mesh.vertices.point_ptr(v2);
            GEO::vec3 a(a_[0],a_[1],a_[2]),b(b_[0],b_[1],b_[2]),c(c_[0],c_[1],c_[2]);
            GEO::vec3 n = GEO::cross(b-a,c-a);
            double area = GEO::length(n)*0.5;
            descriptors[cf_id][3] += area;
            if (area > 1e-14) patch_mean_normal[cf_id] = patch_mean_normal[cf_id] + n*(1.0/GEO::length(n));
        }
        // Normalise mean normals and compute variance
        for (GEO::index_t f = 0; f < nf_c; ++f) {
            double nl = GEO::length(patch_mean_normal[f]);
            if (nl > 1e-14) patch_mean_normal[f] = patch_mean_normal[f]*(1.0/nl);
            descriptors[f][4] = 1.0 - std::min(1.0, nl); // 0=flat, 1=random normals
        }

        // Normalise each descriptor dimension to [0,1]
        for (int d=0; d<DIM; ++d) {
            double mn=1e18, mx=-1e18;
            for (auto& desc : descriptors) { mn=std::min(mn,desc[d]); mx=std::max(mx,desc[d]); }
            double range = mx-mn;
            if (range < 1e-14) continue;
            for (auto& desc : descriptors) desc[d] = (desc[d]-mn)/range;
        }

        // ---- K-means clustering ----
        if (K > (int)nf_c) K = nf_c;
        std::vector<int> assignment(nf_c, 0);
        std::vector<Desc6> centroids(K);
        for (auto& c : centroids) c.fill(0.0);

        // Initialise: spread centroids evenly
        for (int k=0; k<K; ++k) {
            int idx = (int)((double)k/K * nf_c);
            centroids[k] = descriptors[idx];
        }

        for (int iter=0; iter<50; ++iter) {
            // Assign
            bool changed = false;
            for (GEO::index_t f=0; f<nf_c; ++f) {
                double best_d = 1e18; int best_k = 0;
                for (int k=0; k<K; ++k) {
                    double d2=0;
                    for (int d=0;d<DIM;++d) {
                        double diff=descriptors[f][d]-centroids[k][d];
                        d2+=diff*diff;
                    }
                    if (d2<best_d) { best_d=d2; best_k=k; }
                }
                if (assignment[f] != best_k) { assignment[f]=best_k; changed=true; }
            }
            // Update centroids
            std::vector<Desc6> new_c(K);
            std::vector<int> cnt(K,0);
            for (auto& c : new_c) c.fill(0.0);
            for (GEO::index_t f=0; f<nf_c; ++f) {
                int k=assignment[f]; ++cnt[k];
                for (int d=0;d<DIM;++d) new_c[k][d]+=descriptors[f][d];
            }
            for (int k=0; k<K; ++k) if (cnt[k]>0)
                for (int d=0;d<DIM;++d) new_c[k][d]/=cnt[k];
            centroids = new_c;
            if (!changed) break;
        }

        // Log class sizes
        GEO::Logger::out("CurveKID") << "  K=" << K << " template classes:" << std::endl;
        std::vector<int> class_sizes(K,0);
        for (int a : assignment) ++class_sizes[a];
        for (int k=0; k<K; ++k)
            GEO::Logger::out("CurveKID")
                << "    Class " << k << ": " << class_sizes[k] << " patches" << std::endl;

        // ---- Rigid registration error per class ----
        // For each class, compute mean descriptor distance to centroid
        GEO::Logger::out("CurveKID") << "  Rigid registration error (descriptor distance):" << std::endl;
        for (int k=0; k<K; ++k) {
            double sum_d2=0; int cnt=0;
            double max_d=0;
            for (GEO::index_t f=0; f<nf_c; ++f) {
                if (assignment[f]!=k) continue;
                double d2=0;
                for (int d=0;d<DIM;++d) {
                    double diff=descriptors[f][d]-centroids[k][d];
                    d2+=diff*diff;
                }
                sum_d2+=d2; ++cnt;
                max_d=std::max(max_d,std::sqrt(d2));
            }
            double rms=(cnt>0)?std::sqrt(sum_d2/cnt):0;
            GEO::Logger::out("CurveKID")
                << "    Class " << k << ": rms=" << rms << "  max=" << max_d << std::endl;
        }

        // ---- Write template_id (int) and template_color (double [0,1]) ----
        // Per ChatGPT Decision 1: keep template_id as integer class label.
        // Add separate template_color in [0,1] for Graphite colormap painting.
        // Diagnose: log actual min, max, and distinct values in both attributes.
        GEO::Attribute<int>    fattr_tmpl_i(pa_mesh.facets.attributes(), "template_id");
        GEO::Attribute<double> fattr_color (pa_mesh.facets.attributes(), "template_color");
        double color_scale = (K > 1) ? 1.0 / (K - 1) : 1.0;
        int invalid_count = 0;
        for (GEO::index_t f=0; f<nf_pa; ++f) {
            int cf_id = fattr_parent[f];
            int cls = (cf_id>=0 && cf_id<(int)nf_c) ? assignment[cf_id] : -1;
            fattr_tmpl_i[f] = cls;
            // Invalid triangles get 0.0 (not -1.0) so they don't corrupt colormap range.
            // Valid classes map to [0,1] uniformly.
            fattr_color[f]  = (cls >= 0) ? cls * color_scale : 0.0;
            if (cls < 0) ++invalid_count;
        }
        if (invalid_count > 0)
            GEO::Logger::out("CurveKID")
                << "  WARNING: " << invalid_count
                << " triangles have invalid parent_face_id — set to class 0 color. "
                << "This usually means the .geogram file was not reloaded after "
                << "construct_reference_curved_patches." << std::endl;
        // Diagnostic: verify actual attribute values
        {
            int    tid_min=INT_MAX, tid_max=INT_MIN;
            double col_min=1e18,    col_max=-1e18;
            std::set<int> distinct_classes;
            for (GEO::index_t f=0; f<nf_pa; ++f) {
                int t = fattr_tmpl_i[f];
                double c = fattr_color[f];
                tid_min = std::min(tid_min, t);
                tid_max = std::max(tid_max, t);
                col_min = std::min(col_min, c);
                col_max = std::max(col_max, c);
                if (t >= 0) distinct_classes.insert(t);
            }
            GEO::Logger::out("CurveKID")
                << "  DIAGNOSTIC template_id: min=" << tid_min
                << " max=" << tid_max
                << " distinct_classes=" << distinct_classes.size() << std::endl;
            GEO::Logger::out("CurveKID")
                << "  DIAGNOSTIC template_color: min=" << col_min
                << " max=" << col_max << std::endl;
            GEO::Logger::out("CurveKID")
                << "  → In Graphite: paint by 'facets.template_color', "
                << "set colormap min=0 max=1 manually." << std::endl;
        }

        // Save updated mesh with template_id
        std::string out_name = "ref_curved_patches_s_classified";
        GEO::mesh_save(pa_mesh, out_name + ".geogram");
        GEO::mesh_save(pa_mesh, out_name + ".obj");
        GEO::Logger::out("CurveKID")
            << "  Saved: " << out_name << ".geogram  (template_id attribute preserved)" << std::endl
            << "  Paint by 'template_id' for K-class colour map." << std::endl;

        // Update in scene
        MeshGrob* out_grob = MeshGrob::find_or_create(scene_graph(), out_name);
        out_grob->copy(pa_mesh);
        out_grob->update();

        GEO::Logger::out("CurveKID")
            << "=== Phase B classification complete ===" << std::endl
            << "  Paint by 'facets.template_color', set colormap min=0 max=1 manually." << std::endl;
    }

    void MeshGrobCurveKIDCommands::diagnose_panel_geometry(
        std::string coarse_mesh_path, int subdiv)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl; return; }
        const GEO::Mesh& original = *grob;

        // Load coarse mesh
        GEO::Mesh coarse_flat;
        if (coarse_mesh_path.empty()) {
            GEO::Logger::err("CurveKID") << "  coarse_mesh_path required." << std::endl; return;
        }
        if (!GEO::mesh_load(coarse_mesh_path, coarse_flat)) {
            GEO::Logger::err("CurveKID") << "  Cannot load: " << coarse_mesh_path << std::endl; return;
        }

        // Fit spheres (Method B)
        std::vector<CurvedFace> cfaces;
        fit_curved_faces(original, coarse_flat, cfaces,
                         SphereFittingMode::UNCONSTRAINED_7_SAMPLE);

        const GEO::index_t nf = coarse_flat.facets.nb();
        GEO::Logger::out("CurveKID") << "=== Panel Geometry Diagnostic ===" << std::endl;
        GEO::Logger::out("CurveKID") << "  Coarse faces: " << nf << std::endl;

        // Per-panel statistics
        double sum_R = 0, min_R = 1e18, max_R = 0;
        double sum_bulge = 0, max_bulge = 0;
        double sum_signed = 0;
        int wrong_branch = 0, inverted = 0;
        int n_valid = 0;

        // Build face normal for each coarse face (outward)
        std::vector<GEO::vec3> face_normals(nf);
        for (GEO::index_t f = 0; f < nf; ++f) {
            GEO::vec3 p0 = cmd_vertex_point(coarse_flat, coarse_flat.facets.vertex(f,0));
            GEO::vec3 p1 = cmd_vertex_point(coarse_flat, coarse_flat.facets.vertex(f,1));
            GEO::vec3 p2 = cmd_vertex_point(coarse_flat, coarse_flat.facets.vertex(f,2));
            GEO::vec3 n = GEO::cross(p1-p0, p2-p0);
            double nl = GEO::length(n);
            face_normals[f] = (nl > 1e-14) ? n*(1.0/nl) : GEO::vec3(0,0,1);
        }

        // Per-face analysis
        struct FaceStat {
            double R, bulge_max, signed_disp, e_to_orig;
            bool wrong_branch;
            int inverted;
        };
        std::vector<FaceStat> stats(nf);

        for (GEO::index_t f = 0; f < nf; ++f) {
            const CurvedFace& cf = cfaces[f];
            FaceStat& st = stats[f];
            st = {0,0,0,0,false,0};
            if (!cf.valid) continue;
            ++n_valid;

            st.R = cf.sphere_radius;
            sum_R += st.R;
            min_R = std::min(min_R, st.R);
            max_R = std::max(max_R, st.R);

            // Sample interior of this face at subdiv^2 points
            GEO::vec3 p0 = cmd_vertex_point(coarse_flat, cf.v[0]);
            GEO::vec3 p1 = cmd_vertex_point(coarse_flat, cf.v[1]);
            GEO::vec3 p2 = cmd_vertex_point(coarse_flat, cf.v[2]);
            GEO::vec3 C = cf.sphere_center;
            double R = cf.sphere_radius;
            GEO::vec3 fn = face_normals[f];

            // Check sphere branch: C should be on the OPPOSITE side from fn
            // i.e. dot(C - centroid, fn) < 0 means C is below the face
            GEO::vec3 centroid = (p0+p1+p2)*(1.0/3.0);
            double side = GEO::dot(C - centroid, fn);
            // Projection goes AWAY from C, so patch bulges toward +fn if C is below
            // This is correct for a convex surface. Flag if reversed.
            bool branch_ok = true; // will check per sample

            for (int si = 1; si < subdiv; ++si) {
                for (int sj = 1; sj < subdiv-si; ++sj) {
                    double bv0 = 1.0 - double(si+sj)/subdiv;
                    double bv1 = double(sj)/subdiv;
                    double bv2 = double(si)/subdiv;
                    GEO::vec3 pt_flat = p0*bv0 + p1*bv1 + p2*bv2;

                    // Project onto sphere
                    GEO::vec3 d = pt_flat - C;
                    double dlen = GEO::length(d);
                    if (dlen < 1e-14) continue;
                    GEO::vec3 pt_proj = C + d*(R/dlen);

                    // Bulge = displacement from flat position
                    GEO::vec3 disp = pt_proj - pt_flat;
                    double bulge = GEO::length(disp);
                    st.bulge_max = std::max(st.bulge_max, bulge);

                    // Signed displacement along face normal
                    double signed_d = GEO::dot(disp, fn);
                    st.signed_disp = std::max(st.signed_disp,
                                              std::fabs(signed_d));

                    // Branch check: projected point should be on the same
                    // side of the original surface as the flat point.
                    // If signed_d < 0 and face is convex outward, patch
                    // is going inward — wrong branch.
                    if (signed_d < -1e-6) {
                        // Patch is pushing INWARD — potential wrong branch
                        // Only flag if bulge is significant
                        if (bulge > 0.01 * R) branch_ok = false;
                    }
                }
            }

            st.wrong_branch = !branch_ok;
            if (st.wrong_branch) ++wrong_branch;

            sum_bulge += st.bulge_max;
            max_bulge = std::max(max_bulge, st.bulge_max);
            sum_signed += st.signed_disp;
        }

        // Summary log
        double mean_R    = (n_valid > 0) ? sum_R    / n_valid : 0;
        double mean_bulge= (n_valid > 0) ? sum_bulge / n_valid : 0;

        GEO::Logger::out("CurveKID")
            << "  Sphere radius: mean=" << mean_R
            << "  min=" << min_R << "  max=" << max_R << std::endl;
        GEO::Logger::out("CurveKID")
            << "  Interior bulge (displacement from flat): mean=" << mean_bulge
            << "  max=" << max_bulge << std::endl;
        GEO::Logger::out("CurveKID")
            << "  Wrong-branch faces: " << wrong_branch << " / " << n_valid << std::endl;

        // Classify faces by bulge severity
        int severe = 0, moderate = 0, mild = 0;
        for (GEO::index_t f = 0; f < nf; ++f) {
            if (!cfaces[f].valid) continue;
            double b = stats[f].bulge_max;
            double R = stats[f].R;
            if (R < 1e-10) continue;
            double rel = b / R;  // bulge relative to sphere radius
            if      (rel > 0.5)  ++severe;
            else if (rel > 0.1)  ++moderate;
            else                 ++mild;
        }
        GEO::Logger::out("CurveKID")
            << "  Bulge severity (relative to R):"
            << "  severe(>50%R)=" << severe
            << "  moderate(10-50%R)=" << moderate
            << "  mild(<10%R)=" << mild << std::endl;

        // Write colour-mapped error mesh
        // Colour: red=severe bulge, yellow=moderate, green=mild
        GEO::Mesh viz_mesh;
        std::vector<GEO::index_t> ew;
        resolve_shared_arcs(coarse_flat, cfaces, ew);
        tessellate_curved_faces_v2(coarse_flat, cfaces, ew, viz_mesh, subdiv, 1.0,
            TessellationMode::BLENDED_SHARED_EDGE_DEBUG);

        // Add bulge attribute
        GEO::Attribute<double> attr_bulge(viz_mesh.facets.attributes(), "bulge");
        GEO::Attribute<int>    attr_branch(viz_mesh.facets.attributes(), "wrong_branch");
        GEO::Attribute<int>    attr_coarse(viz_mesh.facets.attributes(), "curved_face_id");

        // curved_face_id already set by tessellator; copy bulge per rendered face
        if (attr_coarse.is_bound()) {
            for (GEO::index_t f = 0; f < viz_mesh.facets.nb(); ++f) {
                int cf_id = attr_coarse[f];
                if (cf_id >= 0 && cf_id < (int)nf) {
                    attr_bulge[f]  = stats[cf_id].bulge_max;
                    attr_branch[f] = stats[cf_id].wrong_branch ? 1 : 0;
                }
            }
        }

        std::string out_name = "panel_diag";
        save_mesh_both(viz_mesh, out_name);
        GEO::Logger::out("CurveKID")
            << "  Saved diagnostic mesh: " << out_name << std::endl
            << "  Paint by 'bulge' attribute for error colour map." << std::endl;

        GEO::Logger::out("CurveKID") << "=== Diagnosis complete ===" << std::endl;
    }

    // =========================================================================
    // diagnose_sphere_compatibility — ChatGPT Direction 011
    // =========================================================================
    void MeshGrobCurveKIDCommands::diagnose_sphere_compatibility(
        std::string coarse_mesh_path, double endpoint_tol)
    {
        using namespace CurveKIDAlgo;

        MeshGrob* grob = mesh_grob();
        if (!grob) { GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl; return; }
        const GEO::Mesh& original = *grob;

        coarse_mesh_path = normalize_path(coarse_mesh_path);
        if (coarse_mesh_path.empty()) {
            GEO::Logger::err("CurveKID") << "Provide path to coarse mesh." << std::endl; return;
        }
        GEO::Mesh coarse_flat;
        if (!GEO::mesh_load(coarse_mesh_path, coarse_flat)) {
            GEO::Logger::err("CurveKID") << "Cannot load: " << coarse_mesh_path << std::endl; return;
        }

        // Fit spheres (independent only — no joint opt, to test raw compatibility)
        std::vector<CurvedFace> cfaces;
        fit_curved_faces(original, coarse_flat, cfaces);
        int n_valid = 0;
        for (auto& cf : cfaces) if (cf.valid) ++n_valid;
        GEO::Logger::out("CurveKID")
            << "  Fitted spheres: " << n_valid << " / " << cfaces.size() << std::endl;

        // Run diagnostics
        std::vector<CurveKIDAlgo::EdgeDiagnostic> diags;
        CurveKIDAlgo::CompatibilityReport rep =
            CurveKIDAlgo::diagnose_sphere_compatibility(
                coarse_flat, cfaces, diags, endpoint_tol, true);

        // Write colored visualization: color each face by its worst-edge compatibility
        // 0=compatible(green) 1=incompatible(yellow) 2=disjoint(red) 3=degenerate(gray)
        GEO::Attribute<int> compat_attr(coarse_flat.facets.attributes(), "compat");
        for (GEO::index_t f = 0; f < coarse_flat.facets.nb(); ++f) compat_attr[f] = 0;
        for (auto& d : diags) {
            int col = 0;
            switch (d.compat) {
                case CurveKIDAlgo::EdgeCompatibility::PROPER_COMPATIBLE:      col=0; break;
                case CurveKIDAlgo::EdgeCompatibility::PROPER_INCOMPATIBLE:    col=1; break;
                case CurveKIDAlgo::EdgeCompatibility::DISJOINT_OR_CONTAINED:  col=2; break;
                case CurveKIDAlgo::EdgeCompatibility::TANGENT:                col=2; break;
                case CurveKIDAlgo::EdgeCompatibility::NUMERICALLY_DEGENERATE: col=3; break;
                default: col=3; break;
            }
            if (col > 0) {
                // Colour both adjacent faces
                if (d.face_a < coarse_flat.facets.nb())
                    compat_attr[d.face_a] = std::max(compat_attr[d.face_a], col);
                if (d.face_b != GEO::NO_FACET && d.face_b < coarse_flat.facets.nb())
                    compat_attr[d.face_b] = std::max(compat_attr[d.face_b], col);
            }
        }

        std::string out_base = coarse_mesh_path.substr(0, coarse_mesh_path.rfind('.')) + "_compat";
        save_mesh_both(coarse_flat, out_base);
        GEO::Logger::out("CurveKID")
            << "  Coloured mesh saved: " << out_base << ".obj" << std::endl;
        GEO::Logger::out("CurveKID")
            << "  Paint by 'compat' attribute: 0=compatible 1=incompatible 2=disjoint 3=degenerate" << std::endl;
    }

    // =========================================================================
    // diagnose_template_rigidity — ChatGPT Decision 5
    // =========================================================================
    void MeshGrobCurveKIDCommands::diagnose_template_rigidity(
        std::string phase_a_mesh_path, std::string coarse_mesh_path,
        std::string edgemap_path, int subdiv)
    {
        // Load classified Phase A mesh (needs template_id and parent_face_id)
        GEO::Mesh pa_mesh;
        if (phase_a_mesh_path.empty() || !GEO::mesh_load(phase_a_mesh_path, pa_mesh)) {
            GEO::Logger::err("CurveKID") << "  Cannot load: " << phase_a_mesh_path << std::endl;
            return;
        }
        GEO::Mesh coarse;
        if (coarse_mesh_path.empty() || !GEO::mesh_load(coarse_mesh_path, coarse)) {
            GEO::Logger::err("CurveKID") << "  Cannot load: " << coarse_mesh_path << std::endl;
            return;
        }
        // Check both attributes exist BEFORE constructing (constructor creates if absent)
        if (!pa_mesh.facets.attributes().is_defined("parent_face_id") ||
            !pa_mesh.facets.attributes().is_defined("template_id")) {
            GEO::Logger::err("CurveKID")
                << "  Need both parent_face_id and template_id. "
                << "Load the .geogram output of classify_phase_a_patches "
                << "(not .obj — .obj does not preserve attributes)." << std::endl;
            return;
        }
        GEO::Attribute<int> fattr_parent(pa_mesh.facets.attributes(), "parent_face_id");
        GEO::Attribute<int> fattr_tmpl  (pa_mesh.facets.attributes(), "template_id");

        const int S = subdiv;
        const int nf_c = (int)coarse.facets.nb();
        const int nv_pa = (int)pa_mesh.vertices.nb();

        // Find K from max template_id
        int K = 0;
        for (GEO::index_t f = 0; f < pa_mesh.facets.nb(); ++f)
            K = std::max(K, fattr_tmpl[f] + 1);
        if (K <= 0) {
            GEO::Logger::err("CurveKID") << "  No template_id found." << std::endl;
            return;
        }

        GEO::Logger::out("CurveKID")
            << "=== Template Rigidity Diagnostic  K=" << K
            << "  S=" << subdiv << " ===" << std::endl;

        // ---- Step 1: Collect patch assignment from template_id ----
        // For each coarse patch f, find its class from any triangle with parent_face_id==f
        std::vector<int> patch_class(nf_c, -1);
        for (GEO::index_t f = 0; f < pa_mesh.facets.nb(); ++f) {
            int cf = fattr_parent[f];
            if (cf >= 0 && cf < nf_c && patch_class[cf] < 0)
                patch_class[cf] = fattr_tmpl[f];
        }

        // ---- Step 2: Extract ordered boundary vertex list per coarse patch ----
        // For each coarse patch f, the boundary vertices are the coarse corners
        // plus the edge subdivision vertices, ordered around the 3 edges in sequence:
        //   edge v0→v1 (S-1 interior pts), corner v1,
        //   edge v1→v2 (S-1 interior pts), corner v2,
        //   edge v2→v0 (S-1 interior pts)
        // Total boundary vertices per patch: 3 + 3*(S-1) = 3*S
        // Interior: (S-1)*(S-2)/2 vertices per patch
        //
        // We read the .edgemap.txt to get edge vertex IDs, and derive interior from
        // the parent_face_id vertex attribute.

        // Load edgemap — explicit path preferred; auto-derive from phase_a_mesh_path if empty.
        // The edgemap is saved as ref_curved_patches_s6.edgemap.txt (next to the Phase A geogram).
        // If you supply the classified mesh (ref_curved_patches_s_classified.geogram),
        // set edgemap_path explicitly to ref_curved_patches_s6.edgemap.txt.
        std::string emap_path = edgemap_path;
        if (emap_path.empty()) {
            // Auto-derive: strip extension from phase_a_mesh_path
            emap_path = phase_a_mesh_path;
            auto dot = emap_path.rfind('.');
            if (dot != std::string::npos) emap_path = emap_path.substr(0, dot);
            emap_path += ".edgemap.txt";
        }
        GEO::Logger::out("CurveKID") << "  Looking for edgemap: " << emap_path << std::endl;
        // Map: (lo,hi) → ordered interior vertex IDs from lo to hi
        std::map<std::pair<int,int>, std::vector<int>> emap;
        {
            FILE* fp = fopen(emap_path.c_str(), "r");
            if (!fp) {
                GEO::Logger::err("CurveKID")
                    << "  .edgemap.txt not found: " << emap_path
                    << ". Run construct_reference_curved_patches first." << std::endl;
                return;
            }
            char line[4096];
            while (fgets(line, sizeof(line), fp)) {
                if (line[0] == '#') continue;
                std::vector<int> ids;
                char* p = line;
                while (*p) {
                    while (*p==' '||*p=='\t'||*p=='\r'||*p=='\n') ++p;
                    if (!*p) break;
                    char* end; long val=strtol(p,&end,10);
                    if (end==p) break;
                    ids.push_back((int)val); p=end;
                }
                if (ids.size()<2) continue;
                std::vector<int> interior(ids.begin()+2, ids.end());
                emap[{ids[0],ids[1]}] = interior;
            }
            fclose(fp);
        }

        // Per-coarse-patch: boundary vertex list (ordered, 3*S vertices)
        // and interior vertex list
        auto get_patch_verts = [&](int f_c,
                                   std::vector<int>& bnd,
                                   std::vector<int>& interior)
        {
            bnd.clear(); interior.clear();
            int v[3] = {(int)coarse.facets.vertex(f_c,0),
                        (int)coarse.facets.vertex(f_c,1),
                        (int)coarse.facets.vertex(f_c,2)};
            // Edge 0: v0→v1
            bnd.push_back(v[0]);
            auto key01 = std::make_pair(std::min(v[0],v[1]),std::max(v[0],v[1]));
            auto it01 = emap.find(key01);
            if (it01 != emap.end()) {
                auto ev = it01->second;
                if (v[0] > v[1]) std::reverse(ev.begin(),ev.end());
                for (int vi : ev) bnd.push_back(vi);
            }
            // Edge 1: v1→v2
            bnd.push_back(v[1]);
            auto key12 = std::make_pair(std::min(v[1],v[2]),std::max(v[1],v[2]));
            auto it12 = emap.find(key12);
            if (it12 != emap.end()) {
                auto ev = it12->second;
                if (v[1] > v[2]) std::reverse(ev.begin(),ev.end());
                for (int vi : ev) bnd.push_back(vi);
            }
            // Edge 2: v2→v0
            bnd.push_back(v[2]);
            auto key20 = std::make_pair(std::min(v[2],v[0]),std::max(v[2],v[0]));
            auto it20 = emap.find(key20);
            if (it20 != emap.end()) {
                auto ev = it20->second;
                if (v[2] > v[0]) std::reverse(ev.begin(),ev.end());
                for (int vi : ev) bnd.push_back(vi);
            }
            // Interior: all pa_mesh vertices with parent vertex attr == f_c
            // We collect from pa_mesh: vertices >= nv_coarse that belong to this patch.
            // Use vattr_parent (vertex attribute, not face attribute).
            // If vertex attribute is not present, derive from face attribute:
            // interior vertices are those used only by triangles of this patch.
            // Build from face list:
            std::set<int> bnd_set(bnd.begin(), bnd.end());
            std::set<int> patch_verts;
            for (GEO::index_t f = 0; f < pa_mesh.facets.nb(); ++f) {
                if (fattr_parent[f] != f_c) continue;
                for (int k=0;k<3;++k) patch_verts.insert((int)pa_mesh.facets.vertex(f,k));
            }
            for (int vi : patch_verts)
                if (!bnd_set.count(vi)) interior.push_back(vi);
        };

        // ---- Helper: get 3D point from pa_mesh ----
        auto pt = [&](int vi) -> GEO::vec3 {
            const double* p = pa_mesh.vertices.point_ptr((GEO::index_t)vi);
            return GEO::vec3(p[0],p[1],p[2]);
        };

        // ---- Step 3: For each class, pick representative (first patch in class),
        //              rigidly register all others to it, report RMS. ----
        // Rigid registration via SVD (Kabsch algorithm):
        // Given P[i], Q[i] (N corresponding points):
        //   p_bar = mean(P), q_bar = mean(Q)
        //   P~ = P - p_bar, Q~ = Q - q_bar
        //   H = sum_i P~[i] * Q~[i]^T  (3x3)
        //   SVD: H = U S V^T
        //   R = V U^T  (with det correction for reflections)
        //   t = q_bar - R * p_bar
        // Then RMSE = sqrt(mean ||R*P[i]+t - Q[i]||^2)
        //
        // Stored as flat arrays for MSVC compatibility (no std::array<double,9>).
        struct Mat3 {
            double m[9]; // row-major: m[r*3+c]
            double& at(int r,int c){ return m[r*3+c]; }
            double  at(int r,int c) const { return m[r*3+c]; }
            void zero(){ for(int i=0;i<9;++i) m[i]=0; }
        };

        auto mat_mul = [](const Mat3& A, const Mat3& B) -> Mat3 {
            Mat3 C; C.zero();
            for(int r=0;r<3;++r) for(int c=0;c<3;++c)
                for(int k=0;k<3;++k) C.at(r,c)+=A.at(r,k)*B.at(k,c);
            return C;
        };
        auto mat_T = [](const Mat3& A) -> Mat3 {
            Mat3 B; for(int r=0;r<3;++r) for(int c=0;c<3;++c) B.at(r,c)=A.at(c,r);
            return B;
        };
        auto mat_det = [](const Mat3& A) -> double {
            return A.at(0,0)*(A.at(1,1)*A.at(2,2)-A.at(1,2)*A.at(2,1))
                  -A.at(0,1)*(A.at(1,0)*A.at(2,2)-A.at(1,2)*A.at(2,0))
                  +A.at(0,2)*(A.at(1,0)*A.at(2,1)-A.at(1,1)*A.at(2,0));
        };

        // Simple 3x3 SVD using Jacobi iterations
        // Returns U, S (diagonal as vec3), V such that A = U diag(S) V^T
        auto svd3 = [](Mat3 A, Mat3& U, GEO::vec3& S, Mat3& V) {
            // Golub–Reinsch-style: use symmetric eigendecomposition of A^T A
            // Simple Jacobi for 3x3
            // For now use a Gram-Schmidt + power iteration approach
            // We use the analytical 3x3 SVD from the standard Kabsch implementation:
            // AtA = A^T A, diagonalise by Jacobi
            Mat3 AtA; AtA.zero();
            for(int r=0;r<3;++r) for(int c=0;c<3;++c)
                for(int k=0;k<3;++k) AtA.at(r,c)+=A.at(k,r)*A.at(k,c);
            // Jacobi eigendecomposition of AtA (symmetric)
            Mat3 Vj; Vj.zero();
            for(int i=0;i<3;++i) Vj.at(i,i)=1.0;
            for(int iter=0;iter<64;++iter) {
                // Find largest off-diagonal
                int p2=0,q2=1; double maxval=std::abs(AtA.at(0,1));
                if(std::abs(AtA.at(0,2))>maxval){p2=0;q2=2;maxval=std::abs(AtA.at(0,2));}
                if(std::abs(AtA.at(1,2))>maxval){p2=1;q2=2;maxval=std::abs(AtA.at(1,2));}
                if(maxval<1e-14) break;
                double app=AtA.at(p2,p2),aqq=AtA.at(q2,q2),apq=AtA.at(p2,q2);
                double phi=0.5*atan2(2.0*apq, aqq-app);
                double c2=cos(phi), s2=sin(phi);
                // Givens rotation
                Mat3 G; G.zero(); for(int i=0;i<3;++i) G.at(i,i)=1.0;
                G.at(p2,p2)=c2; G.at(q2,q2)=c2; G.at(p2,q2)=-s2; G.at(q2,p2)=s2;
                // AtA = G^T AtA G
                Mat3 tmp; tmp.zero();
                for(int r=0;r<3;++r) for(int c=0;c<3;++c) for(int k=0;k<3;++k)
                    tmp.at(r,c)+=G.at(k,r)*AtA.at(k,c);
                AtA.zero();
                for(int r=0;r<3;++r) for(int c=0;c<3;++c) for(int k=0;k<3;++k)
                    AtA.at(r,c)+=tmp.at(r,k)*G.at(k,c);
                // Vj = Vj * G
                tmp=Vj; Vj.zero();
                for(int r=0;r<3;++r) for(int c=0;c<3;++c) for(int k=0;k<3;++k)
                    Vj.at(r,c)+=tmp.at(r,k)*G.at(k,c);
            }
            // Singular values = sqrt(eigenvalues of AtA)
            double sv[3];
            for(int i=0;i<3;++i) sv[i]=std::sqrt(std::max(0.0,AtA.at(i,i)));
            S = GEO::vec3(sv[0],sv[1],sv[2]);
            V = Vj;
            // U columns: u_i = A v_i / sigma_i
            U.zero();
            for(int i=0;i<3;++i) {
                GEO::vec3 vi(Vj.at(0,i),Vj.at(1,i),Vj.at(2,i));
                GEO::vec3 Avi(0,0,0);
                for(int r=0;r<3;++r) {
                    double s=0; for(int k=0;k<3;++k) s+=A.at(r,k)*vi[k];
                    Avi[r]=s;
                }
                double len=GEO::length(Avi);
                if(len>1e-14) Avi=Avi*(1.0/len);
                else { // zero singular value: use canonical axis as fallback
                    Avi = (i==0)?GEO::vec3(1,0,0):(i==1)?GEO::vec3(0,1,0):GEO::vec3(0,0,1);
                }
                U.at(0,i)=Avi[0]; U.at(1,i)=Avi[1]; U.at(2,i)=Avi[2];
            }
        };

        // Kabsch rigid registration: register P onto Q (find R,t s.t. R*P+t ≈ Q)
        auto kabsch = [&](const std::vector<GEO::vec3>& P,
                          const std::vector<GEO::vec3>& Q,
                          Mat3& R_out, GEO::vec3& t_out) -> double
        {
            int n = (int)std::min(P.size(), Q.size());
            if (n == 0) { R_out.zero(); for(int i=0;i<3;++i)R_out.at(i,i)=1.0; t_out=GEO::vec3(0,0,0); return 0; }
            GEO::vec3 pbar(0,0,0), qbar(0,0,0);
            for(int i=0;i<n;++i){ pbar=pbar+P[i]; qbar=qbar+Q[i]; }
            pbar=pbar*(1.0/n); qbar=qbar*(1.0/n);
            // Cross-covariance H = sum P~^T Q~  (3x3)
            Mat3 H; H.zero();
            for(int i=0;i<n;++i) {
                GEO::vec3 dp=P[i]-pbar, dq=Q[i]-qbar;
                for(int r=0;r<3;++r) for(int c=0;c<3;++c)
                    H.at(r,c)+=dp[r]*dq[c];
            }
            Mat3 U,V; GEO::vec3 S;
            svd3(H,U,S,V);
            // R = V U^T, with det correction
            Mat3 VUt = mat_mul(V, mat_T(U));
            double d = mat_det(VUt);
            if (d < 0) {
                // flip last column of V
                for(int r=0;r<3;++r) V.at(r,2)*=-1.0;
                VUt = mat_mul(V, mat_T(U));
            }
            R_out = VUt;
            // t = qbar - R * pbar
            GEO::vec3 Rpbar(0,0,0);
            for(int r=0;r<3;++r) for(int c=0;c<3;++c) Rpbar[r]+=R_out.at(r,c)*pbar[c];
            t_out = qbar - Rpbar;
            // RMSE
            double sum2=0;
            for(int i=0;i<n;++i) {
                GEO::vec3 rp(0,0,0);
                for(int r=0;r<3;++r) for(int c=0;c<3;++c) rp[r]+=R_out.at(r,c)*P[i][c];
                rp=rp+t_out;
                sum2+=GEO::dot(rp-Q[i],rp-Q[i]);
            }
            return std::sqrt(sum2/n);
        };

        // ---- Step 4: Per-class registration loop ----
        // Group coarse patches by class
        std::vector<std::vector<int>> class_patches(K);
        for(int f=0;f<nf_c;++f)
            if(patch_class[f]>=0 && patch_class[f]<K)
                class_patches[patch_class[f]].push_back(f);

        // Open CSV output
        std::string csv_path = "template_rigidity.csv";
        FILE* csv = fopen(csv_path.c_str(),"w");
        if(csv) fprintf(csv,"class,patches,bnd_rms,bnd_max,int_rms,int_max,singleton\n");

        double global_bnd_rms_sum=0; int global_count=0;
        double global_worst_bnd=0;
        int singletons=0;

        GEO::Logger::out("CurveKID") << "  Class  Patches  BndRMS      BndMax      IntRMS      IntMax" << std::endl;

        for(int k=0;k<K;++k) {
            auto& members = class_patches[k];
            int np = (int)members.size();
            if(np==0){
                GEO::Logger::out("CurveKID") << "  " << k << "  EMPTY" << std::endl;
                if(csv) fprintf(csv,"%d,0,0,0,0,0,0\n",k);
                continue;
            }
            if(np==1){
                GEO::Logger::out("CurveKID")
                    << "  " << k << "  SINGLETON (patch " << members[0] << ")" << std::endl;
                if(csv) fprintf(csv,"%d,1,0,0,0,0,1\n",k);
                ++singletons;
                continue;
            }
            // Representative = first member
            int rep = members[0];
            std::vector<int> rep_bnd, rep_int;
            get_patch_verts(rep, rep_bnd, rep_int);
            std::vector<GEO::vec3> Q_bnd, Q_int;
            for(int vi: rep_bnd) Q_bnd.push_back(pt(vi));
            for(int vi: rep_int) Q_int.push_back(pt(vi));

            double sum_bnd_rms2=0, max_bnd=0, sum_int_rms2=0, max_int=0;
            int registered=0;
            for(int mi=1;mi<np;++mi) {
                std::vector<int> p_bnd, p_int;
                get_patch_verts(members[mi], p_bnd, p_int);
                std::vector<GEO::vec3> P_bnd, P_int;
                for(int vi: p_bnd) P_bnd.push_back(pt(vi));
                for(int vi: p_int) P_int.push_back(pt(vi));
                // Align boundary counts (take min to handle boundary patches)
                int nb = (int)std::min(P_bnd.size(), Q_bnd.size());
                std::vector<GEO::vec3> Pb(P_bnd.begin(),P_bnd.begin()+nb);
                std::vector<GEO::vec3> Qb(Q_bnd.begin(),Q_bnd.begin()+nb);
                Mat3 R; GEO::vec3 t;
                double brms = kabsch(Pb, Qb, R, t);
                sum_bnd_rms2 += brms*brms; max_bnd=std::max(max_bnd,brms);
                // Apply same R,t to interior
                double sum_i2=0; double mi_max=0;
                int ni=(int)std::min(P_int.size(),Q_int.size());
                for(int i=0;i<ni;++i){
                    GEO::vec3 rp(0,0,0);
                    for(int r2=0;r2<3;++r2) for(int c=0;c<3;++c) rp[r2]+=R.at(r2,c)*P_int[i][c];
                    rp=rp+t;
                    double d2=GEO::dot(rp-Q_int[i],rp-Q_int[i]);
                    sum_i2+=d2; mi_max=std::max(mi_max,std::sqrt(d2));
                }
                double irms=(ni>0)?std::sqrt(sum_i2/ni):0.0;
                sum_int_rms2+=irms*irms; max_int=std::max(max_int,irms);
                ++registered;
            }
            double bnd_rms = (registered>0)?std::sqrt(sum_bnd_rms2/registered):0;
            double int_rms = (registered>0)?std::sqrt(sum_int_rms2/registered):0;
            GEO::Logger::out("CurveKID")
                << "  " << k << "  " << np
                << "  bnd_rms=" << bnd_rms << "  bnd_max=" << max_bnd
                << "  int_rms=" << int_rms << "  int_max=" << max_int << std::endl;
            if(csv) fprintf(csv,"%d,%d,%.6f,%.6f,%.6f,%.6f,0\n",
                            k,np,bnd_rms,max_bnd,int_rms,max_int);
            sum_bnd_rms2 += bnd_rms*bnd_rms; global_count++;
            global_worst_bnd=std::max(global_worst_bnd,max_bnd);
        }

        double mean_bnd_rms=(global_count>0)?std::sqrt(global_bnd_rms_sum/global_count):0;
        GEO::Logger::out("CurveKID")
            << "  === Summary ===" << std::endl
            << "  Singleton classes: " << singletons << " / " << K << std::endl
            << "  Mean boundary RMS (across non-singleton classes): " << mean_bnd_rms << std::endl
            << "  Worst boundary max: " << global_worst_bnd << std::endl
            << "  Saved: " << csv_path << std::endl;

        if(csv) fclose(csv);
        GEO::Logger::out("CurveKID") << "=== diagnose_template_rigidity complete ===" << std::endl;
    }

    // =========================================================================
    // build_phase_a_edges — extract actual Phase A patch boundary polylines
    // Per ChatGPT Decision 2: reads the authoritative .edgemap.txt file saved by
    // construct_reference_curved_patches. Each coarse edge's ordered subdivision
    // vertex IDs are stored explicitly — no geometric reconstruction needed.
    // Falls back to parent_face_id reconstruction if .edgemap.txt is not found.
    // =========================================================================
    void MeshGrobCurveKIDCommands::build_phase_a_edges(
        std::string phase_a_mesh_path, std::string coarse_mesh_path, int subdiv)
    {
        // Load Phase A subdivided mesh
        GEO::Mesh pa_mesh;
        if (phase_a_mesh_path.empty() || !GEO::mesh_load(phase_a_mesh_path, pa_mesh)) {
            GEO::Logger::err("CurveKID")
                << "  build_phase_a_edges: cannot load: " << phase_a_mesh_path << std::endl;
            return;
        }
        const GEO::index_t nv_pa = pa_mesh.vertices.nb();

        GEO::Logger::out("CurveKID")
            << "=== build_phase_a_edges ===" << std::endl
            << "  Phase A: " << nv_pa << " vertices" << std::endl;

        // ---- Try authoritative .edgemap.txt file first ----
        // Derive edgemap path from phase_a_mesh_path by replacing extension
        std::string emap_path = phase_a_mesh_path;
        {
            auto dot = emap_path.rfind('.');
            if (dot != std::string::npos) emap_path = emap_path.substr(0, dot);
            emap_path += ".edgemap.txt";
        }

        // Collect all boundary edge segments to emit
        // Each entry: (va, vb) in output mesh vertex indices
        std::vector<std::pair<GEO::index_t,GEO::index_t>> seg_list;

        FILE* fp = fopen(emap_path.c_str(), "r");
        if (fp) {
            GEO::Logger::out("CurveKID")
                << "  Reading authoritative edge map: " << emap_path << std::endl;
            char line[4096];
            int coarse_edges_read = 0;
            while (fgets(line, sizeof(line), fp)) {
                if (line[0] == '#') continue;
                // Parse: lo hi v0 v1 ... v_{S-2}
                std::vector<GEO::index_t> ids;
                char* p = line;
                while (*p) {
                    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;
                    if (!*p) break;
                    char* end;
                    long val = strtol(p, &end, 10);
                    if (end == p) break;
                    ids.push_back((GEO::index_t)val);
                    p = end;
                }
                if (ids.size() < 2) continue;
                // ids[0]=lo (coarse corner), ids[1]=hi (coarse corner),
                // ids[2..] = subdivision verts in order from lo to hi
                // Full polyline: lo → ids[2] → ... → ids.back() → hi
                GEO::index_t lo = ids[0], hi = ids[1];
                std::vector<GEO::index_t> chain;
                chain.push_back(lo);
                for (int i = 2; i < (int)ids.size(); ++i) chain.push_back(ids[i]);
                chain.push_back(hi);
                for (int i = 0; i+1 < (int)chain.size(); ++i)
                    seg_list.push_back({chain[i], chain[i+1]});
                ++coarse_edges_read;
            }
            fclose(fp);
            GEO::Logger::out("CurveKID")
                << "  Loaded " << coarse_edges_read << " coarse edges → "
                << seg_list.size() << " segment pairs." << std::endl;
        } else {
            // Fallback: reconstruct from parent_face_id boundaries
            GEO::Logger::out("CurveKID")
                << "  .edgemap.txt not found at " << emap_path
                << " — falling back to parent_face_id reconstruction." << std::endl;
            if (!pa_mesh.facets.attributes().is_defined("parent_face_id")) {
                GEO::Logger::err("CurveKID")
                    << "  parent_face_id not found. Load .geogram file." << std::endl;
                return;
            }
            GEO::Attribute<int> fattr_parent(pa_mesh.facets.attributes(), "parent_face_id");
            pa_mesh.facets.connect();
            std::set<std::pair<GEO::index_t,GEO::index_t>> boundary_set;
            for (GEO::index_t f = 0; f < pa_mesh.facets.nb(); ++f) {
                int pid_f = fattr_parent[f];
                for (GEO::index_t le = 0; le < 3; ++le) {
                    GEO::index_t adj = pa_mesh.facets.adjacent(f, le);
                    bool is_bnd = (adj == GEO::NO_FACET) ||
                                  (fattr_parent[adj] != pid_f);
                    if (is_bnd) {
                        GEO::index_t va = pa_mesh.facets.vertex(f, le);
                        GEO::index_t vb = pa_mesh.facets.vertex(f, (le+1)%3);
                        boundary_set.insert({std::min(va,vb), std::max(va,vb)});
                    }
                }
            }
            for (auto& e : boundary_set)
                seg_list.push_back(e);
            GEO::Logger::out("CurveKID")
                << "  Fallback: " << seg_list.size() << " boundary segments." << std::endl;
        }

        // Build output polyline mesh — only include vertices actually used by edges.
        // Copying all nv_pa vertices would leave interior patch vertices as
        // isolated dots in the visualisation.
        // Build a compact remapping: old vertex index → new vertex index.
        std::map<GEO::index_t, GEO::index_t> v_remap;
        for (auto& seg : seg_list) {
            if (!v_remap.count(seg.first))  v_remap[seg.first]  = (GEO::index_t)v_remap.size();
            // second insert must happen after first size() has been recorded
            if (!v_remap.count(seg.second)) v_remap[seg.second] = (GEO::index_t)v_remap.size();
        }
        // Note: the two map inserts above are order-dependent; rebuild cleanly:
        v_remap.clear();
        for (auto& seg : seg_list) {
            if (!v_remap.count(seg.first))  { GEO::index_t id=(GEO::index_t)v_remap.size(); v_remap[seg.first]=id;  }
            if (!v_remap.count(seg.second)) { GEO::index_t id=(GEO::index_t)v_remap.size(); v_remap[seg.second]=id; }
        }

        GEO::Mesh out_edges;
        out_edges.clear();
        out_edges.vertices.create_vertices((GEO::index_t)v_remap.size());
        for (auto& kv : v_remap) {
            const double* p = pa_mesh.vertices.point_ptr(kv.first);
            double* q = out_edges.vertices.point_ptr(kv.second);
            q[0]=p[0]; q[1]=p[1]; q[2]=p[2];
        }
        out_edges.edges.create_edges((GEO::index_t)seg_list.size());
        for (GEO::index_t ei = 0; ei < (GEO::index_t)seg_list.size(); ++ei) {
            out_edges.edges.set_vertex(ei, 0, v_remap[seg_list[ei].first]);
            out_edges.edges.set_vertex(ei, 1, v_remap[seg_list[ei].second]);
        }
        GEO::Logger::out("CurveKID")
            << "  Compact mesh: " << v_remap.size() << " vertices used (of "
            << nv_pa << " total), " << seg_list.size() << " edge segments." << std::endl;

        GEO::Logger::out("CurveKID")
            << "  Output: " << out_edges.edges.nb() << " edge segments" << std::endl;

        std::string out_name = "phase_a_boundary_edges";
        GEO::mesh_save(out_edges, out_name + ".geogram");
        MeshGrob* out_grob = MeshGrob::find_or_create(scene_graph(), out_name);
        out_grob->copy(out_edges);
        out_grob->update();

        GEO::Logger::out("CurveKID")
            << "=== build_phase_a_edges complete ===" << std::endl
            << "  Set edges=visible, surface=hidden to see boundary network." << std::endl;
    }

} // namespace OGF
