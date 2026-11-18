 
#pragma once

#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh.h>
#include <geogram/basic/logger.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <utility>
#include <vector>

namespace OGF {
namespace CurveKID {

struct CurvedMeshValidation {
    GEO::index_t vertices = 0;
    GEO::index_t faces = 0;
    GEO::index_t edges = 0;
    GEO::index_t boundary_edges = 0;
    GEO::index_t nonmanifold_edges = 0;
    GEO::index_t degenerate_faces = 0;
    GEO::index_t inverted_faces = 0;
    GEO::index_t connected_components = 0;

    double mean_radial_error = 0.0;
    double rms_radial_error = 0.0;
    double max_radial_error = 0.0;
};

inline CurvedMeshValidation validate_spherical_mesh(
    const GEO::Mesh& mesh,
    const GEO::vec3& sphere_center,
    double sphere_radius,
    double tolerance = 1e-10
) {
    CurvedMeshValidation result;

    result.vertices = mesh.vertices.nb();
    result.faces = mesh.facets.nb();

    if (sphere_radius <= 0.0) {
        GEO::Logger::err("CurveKID")
            << "Invalid reference sphere radius."
            << std::endl;
        return result;
    }

    const double scale = std::max(1.0, sphere_radius);
    const double area_epsilon = tolerance * scale * scale;

    std::map<
        std::pair<GEO::index_t, GEO::index_t>,
        std::vector<GEO::index_t>
    > edge_faces;

    std::vector<std::vector<GEO::index_t>>
        adjacency(result.faces);

    double error_sum = 0.0;
    double error_sq_sum = 0.0;

    for (GEO::index_t v = 0; v < result.vertices; ++v) {
        GEO::vec3 p = mesh.vertices.point(v);
        double e = std::abs(
            GEO::length(p - sphere_center) - sphere_radius
        );

        error_sum += e;
        error_sq_sum += e * e;
        result.max_radial_error =
            std::max(result.max_radial_error, e);
    }

    if (result.vertices > 0) {
        result.mean_radial_error =
            error_sum / double(result.vertices);

        result.rms_radial_error = std::sqrt(
            error_sq_sum / double(result.vertices)
        );
    }

    for (GEO::index_t f = 0; f < result.faces; ++f) {
        GEO::index_t nv = mesh.facets.nb_vertices(f);

        if (nv != 3) {
            GEO::Logger::warn("CurveKID")
                << "Non-triangular face: " << f
                << std::endl;
            continue;
        }

        GEO::index_t ids[3] = {
            mesh.facets.vertex(f, 0),
            mesh.facets.vertex(f, 1),
            mesh.facets.vertex(f, 2)
        };

        GEO::vec3 a = mesh.vertices.point(ids[0]);
        GEO::vec3 b = mesh.vertices.point(ids[1]);
        GEO::vec3 c = mesh.vertices.point(ids[2]);

        GEO::vec3 normal = GEO::cross(b - a, c - a);
        double area2 = GEO::length(normal);

        if (area2 <= area_epsilon) {
            ++result.degenerate_faces;
        } else {
            GEO::vec3 centroid = (a + b + c) / 3.0;
            GEO::vec3 outward = centroid - sphere_center;

            if (GEO::dot(normal, outward) < 0.0) {
                ++result.inverted_faces;
            }
        }

        for (int k = 0; k < 3; ++k) {
            GEO::index_t u = ids[k];
            GEO::index_t v = ids[(k + 1) % 3];

            if (u > v) {
                std::swap(u, v);
            }

            edge_faces[{u, v}].push_back(f);
        }
    }

    result.edges =
        static_cast<GEO::index_t>(edge_faces.size());

    for (const auto& entry : edge_faces) {
        const auto& incident = entry.second;

        if (incident.size() == 1) {
            ++result.boundary_edges;
        } else if (incident.size() > 2) {
            ++result.nonmanifold_edges;
        }

        if (incident.size() >= 2) {
            for (size_t i = 1; i < incident.size(); ++i) {
                GEO::index_t f0 = incident[0];
                GEO::index_t f1 = incident[i];

                adjacency[f0].push_back(f1);
                adjacency[f1].push_back(f0);
            }
        }
    }

    std::vector<bool> visited(result.faces, false);

    for (GEO::index_t start = 0;
         start < result.faces; ++start) {

        if (visited[start]) {
            continue;
        }

        ++result.connected_components;

        std::queue<GEO::index_t> q;
        q.push(start);
        visited[start] = true;

        while (!q.empty()) {
            GEO::index_t f = q.front();
            q.pop();

            for (GEO::index_t next : adjacency[f]) {
                if (!visited[next]) {
                    visited[next] = true;
                    q.push(next);
                }
            }
        }
    }

    GEO::Logger::out("CurveKID")
        << "\n=== Spherical Mesh Validation ===\n"
        << "Vertices: " << result.vertices << "\n"
        << "Faces: " << result.faces << "\n"
        << "Edges: " << result.edges << "\n"
        << "Boundary edges: " << result.boundary_edges << "\n"
        << "Nonmanifold edges: " << result.nonmanifold_edges << "\n"
        << "Degenerate faces: " << result.degenerate_faces << "\n"
        << "Inverted faces: " << result.inverted_faces << "\n"
        << "Connected components: "
        << result.connected_components << "\n"
        << "Mean radial error: "
        << result.mean_radial_error << "\n"
        << "RMS radial error: "
        << result.rms_radial_error << "\n"
        << "Maximum radial error: "
        << result.max_radial_error << "\n"
        << "===============================\n"
        << std::endl;

    return result;
}

} // namespace CurveKID
} // namespace OGF
