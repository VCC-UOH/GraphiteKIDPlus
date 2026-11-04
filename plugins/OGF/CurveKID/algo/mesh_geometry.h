#ifndef H_OGF_CurveKID_ALGO_MESH_GEOMETRY_H
#define H_OGF_CurveKID_ALGO_MESH_GEOMETRY_H

// -----------------------------------------------------------------------
// CRITICAL: ALL Geogram headers BEFORE any OGF header.
// -----------------------------------------------------------------------
#include <geogram/mesh/mesh.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/attributes.h>

#include <OGF/CurveKID/common/common.h>

#include <vector>
#include <cstdint>

namespace OGF {
namespace CurveKIDAlgo {

// =======================================================================
// Basic mesh statistics
// =======================================================================

struct CurveKID_API MeshStats {
    GEO::index_t nb_vertices      = 0;
    GEO::index_t nb_facets        = 0;
    GEO::index_t nb_triangles     = 0;
    GEO::index_t nb_non_tri       = 0;
    GEO::index_t nb_edges         = 0;
    GEO::index_t nb_border_edges  = 0;
    double min_edge_length  = 0.0;
    double max_edge_length  = 0.0;
    double mean_edge_length = 0.0;
    double min_face_area    = 0.0;
    double max_face_area    = 0.0;
    double mean_face_area   = 0.0;
    double total_surface_area = 0.0;
};

CurveKID_API bool compute_mesh_stats(
    const GEO::Mesh& mesh, MeshStats& stats);

CurveKID_API void compute_face_normals(
    const GEO::Mesh& mesh, std::vector<GEO::vec3>& normals);

CurveKID_API void compute_face_areas(
    const GEO::Mesh& mesh, std::vector<double>& areas);

struct CurveKID_API EdgeInfo {
    GEO::index_t v0 = 0;
    GEO::index_t v1 = 0;
    double length    = 0.0;
    bool is_boundary = false;
};

CurveKID_API void enumerate_triangle_edges(
    const GEO::Mesh& mesh, std::vector<EdgeInfo>& edges);


// =======================================================================
// STAGE 1 — Per-vertex curvature estimation
// =======================================================================

struct CurveKID_API VertexCurvature {
    double gaussian = 0.0;
    double mean     = 0.0;
    double kappa1   = 0.0;
    double kappa2   = 0.0;
};

CurveKID_API void compute_vertex_curvatures(
    const GEO::Mesh& mesh,
    std::vector<VertexCurvature>& vcurv);

struct CurveKID_API FaceCurvature {
    double gaussian = 0.0;
    double mean     = 0.0;
    double kappa1   = 0.0;
    double kappa2   = 0.0;
};

CurveKID_API void compute_face_curvatures(
    const GEO::Mesh& mesh,
    const std::vector<VertexCurvature>& vcurv,
    std::vector<FaceCurvature>& fcurv);

CurveKID_API bool fit_sphere_to_neighborhood(
    const GEO::Mesh& mesh,
    GEO::index_t f,
    GEO::vec3& center,
    double& radius,
    double& residual);


// =======================================================================
// STAGE 2 — Connected patch initialisation
// =======================================================================

CurveKID_API GEO::index_t grow_patches(
    const GEO::Mesh& mesh,
    const std::vector<FaceCurvature>& fcurv,
    GEO::index_t target_count,
    std::vector<GEO::index_t>& patch_id
);


// =======================================================================
// STAGE 3 — Per-patch sphere fitting and descriptor
// =======================================================================

struct CurveKID_API PatchDescriptor {
    GEO::index_t patch_id   = 0;
    GEO::index_t nb_faces   = 0;
    GEO::vec3   sphere_center;
    double      sphere_radius   = 0.0;
    double      sphere_residual = 0.0;
    double mean_gaussian = 0.0;
    double mean_kappa1   = 0.0;
    double mean_kappa2   = 0.0;
    double total_area    = 0.0;
    double mean_edge_len = 0.0;
    GEO::index_t template_id = GEO::index_t(-1);
};

CurveKID_API bool fit_sphere_to_patch(
    const GEO::Mesh& mesh,
    const std::vector<GEO::index_t>& patch_id,
    GEO::index_t pid,
    PatchDescriptor& desc);

CurveKID_API void compute_patch_descriptors(
    const GEO::Mesh& mesh,
    const std::vector<GEO::index_t>& patch_id,
    const std::vector<FaceCurvature>& fcurv,
    GEO::index_t nb_patches,
    std::vector<PatchDescriptor>& descs);


// =======================================================================
// STAGE 4 — K-means clustering
// =======================================================================

CurveKID_API double patch_descriptor_distance(
    const PatchDescriptor& a,
    const PatchDescriptor& b,
    const double inv_std[4]
);

CurveKID_API double cluster_patches_kmeans(
    std::vector<PatchDescriptor>& descs,
    GEO::index_t K,
    int max_iter = 50
);


} // namespace CurveKIDAlgo
} // namespace OGF

#endif // H_OGF_CurveKID_ALGO_MESH_GEOMETRY_H