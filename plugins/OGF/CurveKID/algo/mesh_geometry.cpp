 
#include <OGF/CurveKID/algo/mesh_geometry.h>

#include <geogram/basic/geometry.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace OGF {
namespace CurveKIDAlgo {

namespace {

// ------------------------------------------------------------
// Internal helper: vertex position
// ------------------------------------------------------------

::GEO::vec3 vertex_point(
    const ::GEO::Mesh& mesh,
    ::GEO::index_t v
) {
    const double* p = mesh.vertices.point_ptr(v);

    return ::GEO::vec3(
        p[0],
        p[1],
        p[2]
    );
}

// ------------------------------------------------------------
// Triangle geometry
// ------------------------------------------------------------

::GEO::vec3 triangle_normal_raw(
    const ::GEO::vec3& p0,
    const ::GEO::vec3& p1,
    const ::GEO::vec3& p2
) {
    return ::GEO::cross(
        p1 - p0,
        p2 - p0
    );
}

double triangle_area(
    const ::GEO::vec3& p0,
    const ::GEO::vec3& p1,
    const ::GEO::vec3& p2
) {
    return 0.5 * ::GEO::length(
        triangle_normal_raw(p0, p1, p2)
    );
}

// ------------------------------------------------------------
// Undirected edge representation
// ------------------------------------------------------------

using EdgeKey = std::pair<
    ::GEO::index_t,
    ::GEO::index_t
>;

struct EdgeKeyHash {

    std::size_t operator()(
        const EdgeKey& edge
    ) const noexcept {

        const std::size_t h0 =
            std::hash<::GEO::index_t>{}(edge.first);

        const std::size_t h1 =
            std::hash<::GEO::index_t>{}(edge.second);

        return h0 ^ (
            h1 +
            static_cast<std::size_t>(
                0x9e3779b9U
            ) +
            (h0 << 6) +
            (h0 >> 2)
        );
    }
};

EdgeKey make_edge_key(
    ::GEO::index_t a,
    ::GEO::index_t b
) {
    if (a > b) {
        std::swap(a, b);
    }

    return EdgeKey(a, b);
}

} // anonymous namespace

// ------------------------------------------------------------
// Enumerate unique triangle edges
// ------------------------------------------------------------

void enumerate_triangle_edges(
    const ::GEO::Mesh& mesh,
    std::vector<EdgeInfo>& edges
) {
    edges.clear();

    std::unordered_map<
        EdgeKey,
        std::size_t,
        EdgeKeyHash
    > edge_index;

    std::vector<std::size_t> adjacency_count;

    for (
        ::GEO::index_t f = 0;
        f < mesh.facets.nb();
        ++f
    ) {
        if (mesh.facets.nb_vertices(f) != 3) {
            continue;
        }

        ::GEO::index_t v[3];

        for (int lv = 0; lv < 3; ++lv) {
            v[lv] = mesh.facets.vertex(f, lv);
        }

        for (int e = 0; e < 3; ++e) {

            const EdgeKey key = make_edge_key(
                v[e],
                v[(e + 1) % 3]
            );

            auto it = edge_index.find(key);

            if (it == edge_index.end()) {

                EdgeInfo edge;

                edge.v0 = key.first;
                edge.v1 = key.second;

                edge.length = ::GEO::length(
                    vertex_point(mesh, edge.v1) -
                    vertex_point(mesh, edge.v0)
                );

                edge.is_boundary = false;

                const std::size_t index =
                    edges.size();

                edges.push_back(edge);
                adjacency_count.push_back(1);

                edge_index.emplace(key, index);

            } else {

                ++adjacency_count[it->second];

            }
        }
    }

    for (
        std::size_t i = 0;
        i < edges.size();
        ++i
    ) {
        edges[i].is_boundary =
            (adjacency_count[i] == 1);
    }
}

// ------------------------------------------------------------
// Compute mesh statistics
// ------------------------------------------------------------

bool compute_mesh_stats(
    const ::GEO::Mesh& mesh,
    MeshStats& stats
) {
    // Keep caller's result unchanged on failure.
    MeshStats result{};

    result.nb_vertices = mesh.vertices.nb();
    result.nb_facets = mesh.facets.nb();

    if (result.nb_facets == 0) {
        return false;
    }

    std::vector<double> areas;

    compute_face_areas(mesh, areas);

    double sum_area = 0.0;

    double min_area =
        std::numeric_limits<double>::max();

    double max_area = 0.0;

    for (
        ::GEO::index_t f = 0;
        f < mesh.facets.nb();
        ++f
    ) {
        if (mesh.facets.nb_vertices(f) != 3) {

            ++result.nb_non_tri;

            continue;
        }

        ++result.nb_triangles;

        const double area = areas[f];

        sum_area += area;

        min_area = std::min(min_area, area);
        max_area = std::max(max_area, area);
    }

    if (result.nb_triangles == 0) {
        return false;
    }

    result.min_face_area =
        min_area;

    result.max_face_area =
        max_area;

    result.mean_face_area =
        sum_area /
        static_cast<double>(result.nb_triangles);

    result.total_surface_area =
        sum_area;

    // --------------------------------------------------------
    // Edge statistics
    // --------------------------------------------------------

    std::vector<EdgeInfo> edges;

    enumerate_triangle_edges(mesh, edges);

    result.nb_edges =
        static_cast<::GEO::index_t>(edges.size());

    double sum_edge = 0.0;

    double min_edge =
        std::numeric_limits<double>::max();

    double max_edge = 0.0;

    for (const EdgeInfo& edge : edges) {

        sum_edge += edge.length;

        min_edge = std::min(
            min_edge,
            edge.length
        );

        max_edge = std::max(
            max_edge,
            edge.length
        );

        if (edge.is_boundary) {
            ++result.nb_border_edges;
        }
    }

    if (!edges.empty()) {

        result.min_edge_length = min_edge;
        result.max_edge_length = max_edge;

        result.mean_edge_length =
            sum_edge /
            static_cast<double>(edges.size());
    }

    // Commit only a successful result.
    stats = result;

    return true;
}

// ------------------------------------------------------------
// Compute per-face normals
// ------------------------------------------------------------

void compute_face_normals(
    const ::GEO::Mesh& mesh,
    std::vector<::GEO::vec3>& normals
) {
    normals.assign(
        mesh.facets.nb(),
        ::GEO::vec3(0.0, 0.0, 0.0)
    );

    for (
        ::GEO::index_t f = 0;
        f < mesh.facets.nb();
        ++f
    ) {
        if (mesh.facets.nb_vertices(f) != 3) {
            continue;
        }

        const ::GEO::vec3 p0 = vertex_point(
            mesh,
            mesh.facets.vertex(f, 0)
        );

        const ::GEO::vec3 p1 = vertex_point(
            mesh,
            mesh.facets.vertex(f, 1)
        );

        const ::GEO::vec3 p2 = vertex_point(
            mesh,
            mesh.facets.vertex(f, 2)
        );

        const ::GEO::vec3 normal =
            triangle_normal_raw(p0, p1, p2);

        const double magnitude =
            ::GEO::length(normal);

        if (magnitude > 1e-14) {

            normals[f] =
                normal / magnitude;
        }
    }
}

// ------------------------------------------------------------
// Compute per-face areas
// ------------------------------------------------------------

void compute_face_areas(
    const ::GEO::Mesh& mesh,
    std::vector<double>& areas
) {
    areas.assign(
        mesh.facets.nb(),
        0.0
    );

    for (
        ::GEO::index_t f = 0;
        f < mesh.facets.nb();
        ++f
    ) {
        if (mesh.facets.nb_vertices(f) != 3) {
            continue;
        }

        const ::GEO::vec3 p0 = vertex_point(
            mesh,
            mesh.facets.vertex(f, 0)
        );

        const ::GEO::vec3 p1 = vertex_point(
            mesh,
            mesh.facets.vertex(f, 1)
        );

        const ::GEO::vec3 p2 = vertex_point(
            mesh,
            mesh.facets.vertex(f, 2)
        );

        areas[f] = triangle_area(
            p0,
            p1,
            p2
        );
    }
}

} // namespace CurveKIDAlgo
} // namespace OGF
