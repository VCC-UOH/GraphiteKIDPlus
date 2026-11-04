// -----------------------------------------------------------------------
// CRITICAL: ALL Geogram headers BEFORE any OGF header.
// -----------------------------------------------------------------------
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>
#include <geogram/basic/attributes.h>

#include <OGF/CurveKID/commands/mesh_grob_curvekid_commands.h>
#include <OGF/CurveKID/algo/mesh_geometry.h>

#include <algorithm>
#include <cmath>

namespace OGF {

    MeshGrobCurveKIDCommands::MeshGrobCurveKIDCommands() {}
    MeshGrobCurveKIDCommands::~MeshGrobCurveKIDCommands() {}

    void MeshGrobCurveKIDCommands::inspect_mesh() {
        MeshGrob* grob = mesh_grob();
        if (!grob) {
            GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl;
            return;
        }
        const GEO::Mesh& mesh = *grob;
        CurveKIDAlgo::MeshStats stats;
        if (!CurveKIDAlgo::compute_mesh_stats(mesh, stats)) {
            GEO::Logger::err("CurveKID") << "Mesh has no triangular facets." << std::endl;
            return;
        }
        GEO::Logger::out("CurveKID") << "=== Mesh Statistics ===" << std::endl;
        GEO::Logger::out("CurveKID") << "Vertices      : " << stats.nb_vertices << std::endl;
        GEO::Logger::out("CurveKID") << "Facets        : " << stats.nb_facets << std::endl;
        GEO::Logger::out("CurveKID") << "Triangles     : " << stats.nb_triangles << std::endl;
        GEO::Logger::out("CurveKID") << "Non-tri facets: " << stats.nb_non_tri << std::endl;
        GEO::Logger::out("CurveKID") << "Unique edges  : " << stats.nb_edges << std::endl;
        GEO::Logger::out("CurveKID") << "Border edges  : " << stats.nb_border_edges << std::endl;
        GEO::Logger::out("CurveKID") << "Edge length   : min=" << stats.min_edge_length
            << " max=" << stats.max_edge_length
            << " mean=" << stats.mean_edge_length << std::endl;
        GEO::Logger::out("CurveKID") << "Face area     : min=" << stats.min_face_area
            << " max=" << stats.max_face_area
            << " mean=" << stats.mean_face_area << std::endl;
        GEO::Logger::out("CurveKID") << "Total area    : " << stats.total_surface_area << std::endl;
    }

    void MeshGrobCurveKIDCommands::estimate_curvatures() {
        MeshGrob* grob = mesh_grob();
        if (!grob) {
            GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl;
            return;
        }
        const GEO::Mesh& mesh = *grob;
        GEO::Logger::out("CurveKID") << "Stage 1: computing vertex curvatures..." << std::endl;

        std::vector<CurveKIDAlgo::VertexCurvature> vcurv;
        CurveKIDAlgo::compute_vertex_curvatures(mesh, vcurv);
        std::vector<CurveKIDAlgo::FaceCurvature> fcurv;
        CurveKIDAlgo::compute_face_curvatures(mesh, vcurv, fcurv);

        double min_K=1e30,max_K=-1e30,sum_K=0;
        double min_H=1e30,max_H=-1e30,sum_H=0;
        int cnt=0;
        for (auto& fc:fcurv) {
            min_K=std::min(min_K,fc.gaussian); max_K=std::max(max_K,fc.gaussian);
            sum_K+=fc.gaussian;
            min_H=std::min(min_H,fc.mean); max_H=std::max(max_H,fc.mean);
            sum_H+=fc.mean; ++cnt;
        }
        if (cnt>0) {
            GEO::Logger::out("CurveKID") << "Gaussian curvature: min=" << min_K
                << " max=" << max_K << " mean=" << sum_K/cnt << std::endl;
            GEO::Logger::out("CurveKID") << "Mean curvature    : min=" << min_H
                << " max=" << max_H << " mean=" << sum_H/cnt << std::endl;
        }
        GEO::Logger::out("CurveKID") << "Stage 1 done." << std::endl;
    }

    void MeshGrobCurveKIDCommands::run_curvekid(int nb_patches, int K) {
        MeshGrob* grob = mesh_grob();
        if (!grob) {
            GEO::Logger::err("CurveKID") << "No mesh selected." << std::endl;
            return;
        }
        GEO::Mesh& mesh = *grob;

        GEO::Logger::out("CurveKID") << "[1/5] Estimating curvatures..." << std::endl;
        std::vector<CurveKIDAlgo::VertexCurvature> vcurv;
        CurveKIDAlgo::compute_vertex_curvatures(mesh, vcurv);
        std::vector<CurveKIDAlgo::FaceCurvature> fcurv;
        CurveKIDAlgo::compute_face_curvatures(mesh, vcurv, fcurv);

        GEO::Logger::out("CurveKID") << "[2/5] Growing ~" << nb_patches
            << " connected patches..." << std::endl;
        std::vector<GEO::index_t> patch_id;
        GEO::index_t actual_patches = CurveKIDAlgo::grow_patches(
            mesh, fcurv, (GEO::index_t)nb_patches, patch_id);
        GEO::Logger::out("CurveKID") << "      Created " << actual_patches
            << " patches." << std::endl;

        GEO::Logger::out("CurveKID") << "[3/5] Fitting spheres to patches..." << std::endl;
        std::vector<CurveKIDAlgo::PatchDescriptor> descs;
        CurveKIDAlgo::compute_patch_descriptors(
            mesh, patch_id, fcurv, actual_patches, descs);

        double rmin=1e30,rmax=-1e30,rsum=0; int rcnt=0;
        for (auto& d:descs) {
            if(d.nb_faces==0) continue;
            rmin=std::min(rmin,d.sphere_radius);
            rmax=std::max(rmax,d.sphere_radius);
            rsum+=d.sphere_radius; ++rcnt;
        }
        if (rcnt>0)
            GEO::Logger::out("CurveKID") << "      Sphere radius: min=" << rmin
                << " max=" << rmax << " mean=" << rsum/rcnt << std::endl;

        GEO::Logger::out("CurveKID") << "[4/5] K-means clustering (K="
            << K << ")..." << std::endl;
        double inertia = CurveKIDAlgo::cluster_patches_kmeans(
            descs, (GEO::index_t)K);
        GEO::Logger::out("CurveKID") << "      Inertia = " << inertia << std::endl;

        std::vector<int> tmpl_count(K,0);
        for (GEO::index_t f=0;f<mesh.facets.nb();++f) {
            if(mesh.facets.nb_vertices(f)!=3) continue;
            GEO::index_t pid=patch_id[f];
            if(pid<actual_patches && descs[pid].template_id!=GEO::index_t(-1)) {
                int tid=(int)descs[pid].template_id;
                if(tid<K) ++tmpl_count[tid];
            }
        }
        for (int t=0;t<K;++t)
            GEO::Logger::out("CurveKID") << "      Template " << t
                << ": " << tmpl_count[t] << " faces" << std::endl;

        GEO::Logger::out("CurveKID") << "[5/5] Writing facet attributes..." << std::endl;
        GEO::Attribute<double> attr_patch(mesh.facets.attributes(), "patch_id");
        GEO::Attribute<double> attr_tmpl(mesh.facets.attributes(),  "template_id");
        GEO::Attribute<double> attr_err(mesh.facets.attributes(),   "sphere_residual");

        for (GEO::index_t f=0;f<mesh.facets.nb();++f) {
            GEO::index_t pid=patch_id[f];
            attr_patch[f]=(double)pid;
            if(pid<actual_patches) {
                GEO::index_t tid=descs[pid].template_id;
                attr_tmpl[f]=(tid==GEO::index_t(-1))?-1.0:(double)tid;
                attr_err[f]=descs[pid].sphere_residual;
            } else {
                attr_tmpl[f]=-1.0; attr_err[f]=0.0;
            }
        }

        grob->update();
        GEO::Logger::out("CurveKID")
            << "Done. To visualise: right-click mesh > Attributes > "
               "facets > template_id" << std::endl;
    }

} // namespace OGF