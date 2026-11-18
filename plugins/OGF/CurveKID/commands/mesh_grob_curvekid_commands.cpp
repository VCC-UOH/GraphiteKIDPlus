// -----------------------------------------------------------------------
// CRITICAL: ALL Geogram headers BEFORE any OGF header.
// -----------------------------------------------------------------------
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
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

    void MeshGrobCurveKIDCommands::build_coarse_curved_mesh(int N, int subdiv)
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

        // --- Stage 2: fit one sphere per coarse triangle (Pottmann Eq.7) ---
        std::vector<CurvedFace> cfaces;
        fit_curved_faces(original, coarse_flat, cfaces);

        // --- Stage 2b: joint optimisation (Pottmann 2023 Sec 4) ---
        // Moves coarse vertices and sphere params so adjacent spheres
        // intersect in real circles → shared arcs are well-defined.
        GEO::Logger::out("CurveKID") << "  Running joint sphere optimisation..." << std::endl;
        optimise_joint_spheres(coarse_flat, cfaces, 30);

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
        GEO::Mesh viz_mesh;
        tessellate_curved_faces_v2(coarse_flat, cfaces, edge_winner, viz_mesh, subdiv);
        GEO::Logger::out("CurveKID")
            << "  Viz mesh: " << viz_mesh.facets.nb()
            << " render triangles  (" << subdiv << "x subdiv per curved face)"
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
        std::string coarse_mesh_name, int subdiv)
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

        // --- Joint optimisation (Pottmann 2023 Sec 4) ---
        GEO::Logger::out("CurveKID") << "  Running joint sphere optimisation..." << std::endl;
        optimise_joint_spheres(coarse_flat, cfaces, 30);

        int curved_cnt = 0;
        for (auto& cf : cfaces) if (cf.valid) ++curved_cnt;
        GEO::Logger::out("CurveKID") << "  Curved faces: " << curved_cnt
                                     << " / " << cfaces.size() << std::endl;

        // --- Shared arcs + tessellation ---
        std::vector<GEO::index_t> edge_winner;
        resolve_shared_arcs(coarse_flat, cfaces, edge_winner);
        GEO::Mesh viz_mesh;
        tessellate_curved_faces_v2(coarse_flat, cfaces, edge_winner, viz_mesh, subdiv);

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
        tessellate_curved_faces_v2(coarse_mesh, cfaces, edge_winner, viz_mesh, subdiv);
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

} // namespace OGF
