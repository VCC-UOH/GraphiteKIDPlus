
#ifndef H_OGF_CurveKID_ALGO_MESH_GEOMETRY_H
#define H_OGF_CurveKID_ALGO_MESH_GEOMETRY_H

#include <OGF/CurveKID/common/common.h>

#include <geogram/mesh/mesh.h>
#include <geogram/basic/geometry.h>

#include <vector>

// IMPORTANT:
// All Geogram headers must be included BEFORE
// opening namespace OGF.

namespace OGF {
namespace CurveKIDAlgo {

// ------------------------------------------------------------
// Basic mesh statistics
// ------------------------------------------------------------

struct CurveKID_API MeshStats {

    GEO::index_t nb_vertices = 0;
    GEO::index_t nb_facets = 0;

    GEO::index_t nb_triangles = 0;
    GEO::index_t nb_non_tri = 0;

    GEO::index_t nb_edges = 0;
    GEO::index_t nb_border_edges = 0;

    double min_edge_length = 0.0;
    double max_edge_length = 0.0;
    double mean_edge_length = 0.0;

    double min_face_area = 0.0;
    double max_face_area = 0.0;
    double mean_face_area = 0.0;

    double total_surface_area = 0.0;
};

// Compute statistics for a mesh.
// Returns false if the mesh has no facets.
CurveKID_API bool compute_mesh_stats(
    const GEO::Mesh& mesh,
    MeshStats& stats
);

// ------------------------------------------------------------
// Per-face geometry
// ------------------------------------------------------------

// Compute unit face normals.
// Non-triangular or degenerate facets receive zero vectors.
CurveKID_API void compute_face_normals(
    const GEO::Mesh& mesh,
    std::vector<GEO::vec3>& normals
);

// Compute triangle face areas.
// Non-triangular facets receive zero area.
CurveKID_API void compute_face_areas(
    const GEO::Mesh& mesh,
    std::vector<double>& areas
);

// ------------------------------------------------------------
// Per-edge geometry
// ------------------------------------------------------------

struct CurveKID_API EdgeInfo {

    GEO::index_t v0 = 0;
    GEO::index_t v1 = 0;

    double length = 0.0;
    bool is_boundary = false;
};

// Enumerate unique edges belonging to triangle facets.
// Does not require mesh.edges to be populated.
CurveKID_API void enumerate_triangle_edges(
    const GEO::Mesh& mesh,
    std::vector<EdgeInfo>& edges
);

} // namespace CurveKIDAlgo
} // namespace OGF

#endif // H_OGF_CurveKID_ALGO_MESH_GEOMETRY_H
