#ifndef H_OGF_CurveKID_ALGO_MESH_GEOMETRY_H
#define H_OGF_CurveKID_ALGO_MESH_GEOMETRY_H

// -----------------------------------------------------------------------
// CRITICAL: ALL Geogram headers BEFORE any OGF header.
// Including Geogram after OGF opens namespace OGF, causing the compiler
// to see GEO types as OGF::GEO::* — which breaks everything.
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
// Basic mesh statistics  (unchanged from working build)
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
    double gaussian = 0.0;  // angle-deficit Gaussian curvature  K = 2π - Σθᵢ
    double mean     = 0.0;  // cotangent-formula mean curvature  H
    double kappa1   = 0.0;  // min principal curvature
    double kappa2   = 0.0;  // max principal curvature
};

// Compute per-vertex curvatures for all vertices.
// Uses angle-deficit for K and cotangent formula for H,
// then extracts kappa1 = H - sqrt(H²-K), kappa2 = H + sqrt(H²-K).
// Falls back to kappa1=kappa2=H when discriminant < 0 (numerical noise).
CurveKID_API void compute_vertex_curvatures(
    const GEO::Mesh& mesh,
    std::vector<VertexCurvature>& vcurv);

// Average vertex curvatures to face level (mean over 3 corners).
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

// Fit a sphere to the 1-ring neighbourhood of face f.
// Collects all vertices in adjacent faces and fits a least-squares sphere.
// Returns false if neighbourhood is too small or degenerate.
CurveKID_API bool fit_sphere_to_neighborhood(
    const GEO::Mesh& mesh,
    GEO::index_t f,
    GEO::vec3& center,
    double& radius,
    double& residual);


// =======================================================================
// STAGE 2 — Connected patch initialisation
// =======================================================================

// Grow ~target_count connected patches using region-growing BFS.
// Each face gets a patch_id in [0, actual_count).
// Seeds are placed by furthest-point sampling in curvature space.
CurveKID_API GEO::index_t grow_patches(
    const GEO::Mesh& mesh,
    const std::vector<FaceCurvature>& fcurv,
    GEO::index_t target_count,          // e.g. 100
    std::vector<GEO::index_t>& patch_id // size = nb_facets
);


// =======================================================================
// STAGE 3 — Per-patch sphere fitting and descriptor
// =======================================================================

struct CurveKID_API PatchDescriptor {
    GEO::index_t patch_id   = 0;
    GEO::index_t nb_faces   = 0;

    // Sphere fit over all patch vertices (least-squares)
    GEO::vec3   sphere_center;
    double      sphere_radius   = 0.0;
    double      sphere_residual = 0.0;  // RMS distance from vertices to sphere

    // Curvature summary (mean over patch faces)
    double mean_gaussian = 0.0;
    double mean_kappa1   = 0.0;
    double mean_kappa2   = 0.0;

    // Geometry summary
    double total_area    = 0.0;
    double mean_edge_len = 0.0;

    // Template assigned in Stage 4
    GEO::index_t template_id = GEO::index_t(-1);
};

// Fit a sphere to all vertices belonging to one patch (via patch_id array).
// Uses iterative least-squares: initialise at centroid + mean distance,
// then refine with Gauss-Newton.
CurveKID_API bool fit_sphere_to_patch(
    const GEO::Mesh& mesh,
    const std::vector<GEO::index_t>& patch_id,
    GEO::index_t pid,
    PatchDescriptor& desc);

// Build all patch descriptors at once.
CurveKID_API void compute_patch_descriptors(
    const GEO::Mesh& mesh,
    const std::vector<GEO::index_t>& patch_id,
    const std::vector<FaceCurvature>& fcurv,
    GEO::index_t nb_patches,
    std::vector<PatchDescriptor>& descs);


// =======================================================================
// STAGE 4 — K-means clustering into template families
// =======================================================================

// Distance between two patch descriptors in feature space.
// Feature vector: (sphere_radius, mean_kappa1, mean_kappa2, mean_edge_len)
// Normalised by global std before distance computation.
CurveKID_API double patch_descriptor_distance(
    const PatchDescriptor& a,
    const PatchDescriptor& b,
    const double inv_std[4]   // precomputed 1/std for each feature
);

// Run K-means on patch descriptors.
// Assigns template_id to each descriptor.
// Returns final inertia (sum of squared distances to assigned centroid).
CurveKID_API double cluster_patches_kmeans(
    std::vector<PatchDescriptor>& descs,
    GEO::index_t K,            // number of templates
    int max_iter = 50
);


// =======================================================================
// POTTMANN-INSPIRED INIT — Central spheres + spherical face descriptors
// (experimental alternative to the BFS patch pipeline above;
//  existing code is preserved and unchanged)
// =======================================================================

// Per-vertex mean-curvature (central) sphere.
// centre = v + (1/H)*n,  radius = |1/H|.
// When |H| < H_eps (flat or saddle near H=0), the sphere is flagged
// invalid and a fallback osculating-paraboloid radius is stored instead.
struct CurveKID_API CentralSphere {
    GEO::vec3 center;          // sphere centre (position-dependent)
    double    radius  = 0.0;   // |1/H|  or  fallback radius
    double    H_value = 0.0;   // signed mean curvature at this vertex
    bool      valid   = false; // false when |H| < H_eps
};

// Compute one CentralSphere per vertex.
// H_eps: vertices with |H| < H_eps are treated as near-flat.
// Fallback for near-flat vertices: radius estimated from local 1-ring
// geometry (circumradius of 1-ring, as a proxy for curvature radius).
CurveKID_API void compute_central_spheres(
    const GEO::Mesh&                    mesh,
    const std::vector<VertexCurvature>& vcurv,
    std::vector<CentralSphere>&         cspheres,
    double H_eps = 0.01         // tune: smaller = more near-flat vertices
);

// Per-face curvature-adaptive target edge length.
// g(f) = c / sqrt( max(|H_f|, H_eps) ), clamped to [min_len, max_len].
// H_f = mean H over the three face vertices.
// c is chosen so that the mean g equals the mesh mean edge length.
CurveKID_API void compute_target_edge_lengths(
    const GEO::Mesh&                    mesh,
    const std::vector<VertexCurvature>& vcurv,
    std::vector<double>&                target_len, // one per face
    double H_eps   = 0.01
);

// One curved triangular face in the Pottmann-inspired representation.
// Supporting sphere is fitted using the face's three vertices PLUS
// the three vertex normals (as tangent-plane constraints), giving a
// well-determined system even for near-flat faces.
struct CurveKID_API SphericalFace {
    GEO::index_t v[3]         = {0,0,0}; // vertex indices
    GEO::vec3    sphere_center;
    double       sphere_radius  = 0.0;
    double       fit_residual   = 0.0;   // RMS point-to-sphere distance
    bool         valid          = false;

    // Pose-invariant descriptor for K-means classification
    // (does NOT use sphere_center — position-independent)
    double desc_radius          = 0.0;   // sphere_radius
    double desc_arc[3]          = {0,0,0}; // great-circle arc lengths of 3 edges
    double desc_solid_angle     = 0.0;   // solid angle subtended at sphere centre
    double desc_mean_H          = 0.0;   // mean |H| over 3 vertices
};

// =======================================================================
// CURVED MESH CONSTRUCTION  (Stages 1-4 of CurveKID pipeline)
// =======================================================================

// One curved triangular face in the output CurveKID mesh.
// The supporting sphere is fitted to nearby points of the ORIGINAL surface,
// so the curved triangle approximates the original geometry locally.
// Shared boundary arcs: for each edge, the arc is on the sphere of
// whichever adjacent face has the smaller residual (better fit).
struct CurveKID_API CurvedFace {
    // Coarse mesh vertex indices (into the coarse mesh, not original)
    GEO::index_t v[3]            = {0, 0, 0};

    // Supporting sphere fitted to original-surface points near this face
    GEO::vec3    sphere_center;
    double       sphere_radius   = 0.0;
    double       sphere_residual = 0.0;  // RMS original-point to sphere dist
    bool         valid           = false;

    // Index into coarse mesh facets
    GEO::index_t face_id         = 0;

    // Approximation quality: mean distance of original mesh vertices
    // that fall inside this face's Voronoi cell to its sphere
    double approx_error          = 0.0;

    // Mean |H| of original-mesh neighbourhood — used for classification.
    // 0 = flat region (sphere fit skipped), >0 = curved region.
    double mean_H_local;
};

// Validate mesh topology: check Euler characteristic, boundary edges,
// non-manifold edges, degenerate faces, and connectivity.
// Returns true if the mesh is a valid closed manifold triangulation.
// Logs all failures to GEO::Logger.
struct CurveKID_API MeshTopologyReport {
    int  connected_components = 0;
    int  boundary_edges       = 0;
    int  nonmanifold_edges    = 0;
    int  degenerate_faces     = 0;
    int  euler               = 0;   // V - E + F, should be 2 for a sphere
    bool face_count_even      = false;
    bool valid                = false;
};

CurveKID_API MeshTopologyReport validate_mesh_topology(
    const GEO::Mesh& mesh,
    bool log_results = true
);

// Build a coarse, curvature-adaptive triangular mesh from the input surface.
// Uses Geogram's mesh_remesh_smooth with per-face target edge lengths.
// Output coarse_mesh is a NEW mesh (caller clears it first).
// Returns actual number of coarse faces produced.
CurveKID_API GEO::index_t build_coarse_mesh(
    const GEO::Mesh&                    original,
    const std::vector<VertexCurvature>& vcurv,
    GEO::Mesh&                          coarse_mesh,
    GEO::index_t                        target_faces = 500,
    double                              H_eps        = 0.01
);

// Fit one CurvedFace per coarse triangle.
// For each coarse face, collects original-mesh vertices within a search
// Sphere fitting mode for fit_curved_faces().
//
//   UNCONSTRAINED_7_SAMPLE:
//     Original Pottmann Eq.7 — fits sphere to 7 dense-mesh samples
//     (3 nearest-to-corners + 3 edge-midpoint snaps + centroid snap).
//     Minimises algebraic normal-direction residual. Coarse corners do
//     NOT generally lie on the fitted sphere → nonzero endpoint error
//     → 0% PROPER_COMPATIBLE edges. (Current default.)
//
//   CONSTRAINED_CORNER_INTERPOLATION:
//     Experiment C — same 7 dense-mesh samples but sphere centre is
//     constrained to the circumcenter-lift line C(t)=cc+t*n so that
//     all 3 coarse corners lie exactly on the sphere. The lift t is
//     chosen to minimise Σ(|q_i - C(t)|² - R²(t))² (algebraic
//     distance objective). Guarantees corner_res=0 → shared edge
//     endpoints lie on both adjacent spheres → PROPER_COMPATIBLE.
enum class SphereFittingMode {
    UNCONSTRAINED_7_SAMPLE,          // original method B
    CONSTRAINED_CORNER_INTERPOLATION // new method C
};

// radius (estimated from coarse edge length) and fits a least-squares sphere.
// orig_index: an AABB/kd-tree built from the original mesh (passed in so
// it is built once, not per-face).
CurveKID_API void fit_curved_faces(
    const GEO::Mesh&              original,
    const GEO::Mesh&              coarse_mesh,
    std::vector<CurvedFace>&      cfaces,
    SphereFittingMode             mode = SphereFittingMode::UNCONSTRAINED_7_SAMPLE
);

// Joint sphere optimisation (Pottmann 2023, Sec 4).
// After independent per-face fitting, runs alternating LS to move
// coarse vertices and sphere parameters so adjacent spheres intersect
// in real circles — required for well-defined shared arcs.
// coarse_mesh vertices are updated in-place.
CurveKID_API void optimise_joint_spheres(
    GEO::Mesh&                 coarse_mesh,
    std::vector<CurvedFace>&   cfaces,
    int                        max_iter = 20
);

// Tessellate each CurvedFace onto its supporting sphere and write the
// result into viz_mesh (cleared first).  Each coarse triangle becomes
// subdiv*subdiv sub-triangles whose vertices are slerp-projected onto
// the sphere.  Face attribute "curved_face_id" maps sub-faces back to
// the coarse face index.
CurveKID_API void tessellate_curved_faces(
    const GEO::Mesh&                   coarse_mesh,
    const std::vector<CurvedFace>&     cfaces,
    GEO::Mesh&                         viz_mesh,
    int                                subdiv = 4
);

// Compute per-original-vertex approximation error:
// for each original vertex find the nearest CurvedFace sphere and
// record signed distance.  Writes result into errors (size = nb_vertices).
CurveKID_API void compute_approximation_error(
    const GEO::Mesh&                   original,
    const std::vector<CurvedFace>&     cfaces,
    std::vector<double>&               errors   // one per original vertex
);

// =======================================================================
// CURVED MESH — shared arc resolution
// =======================================================================

// For each coarse edge shared by two faces, records which face's sphere
// "wins" (lower residual) and owns the arc.
// edge_winner[f*3 + local_edge] = index of the winning face for that edge
// (may be f itself or the adjacent face).
CurveKID_API void resolve_shared_arcs(
    const GEO::Mesh&               coarse_mesh,
    const std::vector<CurvedFace>& cfaces,
    std::vector<GEO::index_t>&     edge_winner  // size = nf_coarse * 3
);

// Fit circumsphere through the 3 coarse vertices of face f.
// C = circumcenter + t * face_normal,  R = |C - v0|.
// t=0 gives the circumscribed sphere of the triangle in its own plane
//   (infinite radius for a flat triangle lying in one plane — degenerate).
// A small nonzero t gives a sphere that passes exactly through all 3
// corners and bulges by approximately t out of the triangle plane.
// Returns false if the triangle is degenerate.
CurveKID_API bool fit_circumsphere(
    const GEO::vec3& p0,
    const GEO::vec3& p1,
    const GEO::vec3& p2,
    double           t,          // normal offset: try circumradius * 0.5
    GEO::vec3&       center,
    double&          radius
);

// Synthetic mode: fit one CurvedFace per coarse triangle using circumsphere.
// Exact corners preserved. t_factor controls curvature: sphere centre is at
//   circumcenter + t_factor * circumradius * face_normal
// Shared arcs resolved afterward by resolve_shared_arcs.
CurveKID_API void fit_curved_faces_synthetic(
    const GEO::Mesh&         coarse_mesh,
    std::vector<CurvedFace>& cfaces,
    double                   t_factor = 0.5
);

// Tessellate with shared-arc awareness.
// For each coarse edge, the arc vertices come from edge_winner[f*3+e]'s sphere.
// Interior vertices come from the face's own sphere.
// Tessellation mode for tessellate_curved_faces_v2().
//
//   OWNER_SPHERE (original):
//     Each interior edge is owned by one face (lower index).
//     Edge points generated by slerp on the owner's sphere only.
//     Adjacent face uses the same points → no seam in position,
//     but edge curve lies on only ONE sphere → incompatible panels.
//
//   BLENDED_SHARED_EDGE_DEBUG (Option 2 / ChatGPT Direction 016):
//     For each interior edge shared by faces A and B:
//       e(t) = 0.5 * slerp_A(t) + 0.5 * slerp_B(t)
//     Both faces use identical blended positions → zero seam.
//     Interior points use boundary-constrained interpolation:
//     barycentric blend of the 3 (now blended) edge curves,
//     then projected onto the face sphere.
//     NOT exact spherical panels — for visualization/debug only.
//     Label all output from this mode as "blended approximation".
//
//   EXACT_SPHERE_INTERSECTION (Option 1 — future):
//     For PROPER_COMPATIBLE edges: use intersection circle arc.
//     For incompatible edges: returns diagnostic, no silent fallback.
//     Not yet implemented.
enum class TessellationMode {
    OWNER_SPHERE,               // original method
    BLENDED_SHARED_EDGE_DEBUG,  // Option 2: averaged edge curves
    EXACT_SPHERE_INTERSECTION   // Option 1: intersection circle arcs (future)
};

CurveKID_API void tessellate_curved_faces_v2(
    const GEO::Mesh&                   coarse_mesh,
    const std::vector<CurvedFace>&     cfaces,
    const std::vector<GEO::index_t>&   edge_winner,
    GEO::Mesh&                         viz_mesh,
    int                                subdiv = 4,
    double                             alpha  = 1.0,
    TessellationMode                   tess_mode = TessellationMode::OWNER_SPHERE
);

// Compute per-original-vertex approximation error against the FLAT coarse mesh.
// This is E_planar(N): mean/max distance of original vertices to the nearest
// flat triangle of the coarse mesh.  Used to compare against E_curved(N).
// =======================================================================
// SPHERE COMPATIBILITY DIAGNOSTICS (ChatGPT Direction 011)
// =======================================================================

// Edge compatibility classification
enum class EdgeCompatibility {
    PROPER_COMPATIBLE,       // intersection circle exists AND endpoints lie on it
    PROPER_INCOMPATIBLE,     // intersection circle exists but endpoints deviate
    COINCIDENT,              // both spheres identical
    TANGENT,                 // point intersection only
    DISJOINT_OR_CONTAINED,   // no real intersection circle
    BOUNDARY_EDGE,           // only one adjacent face
    NUMERICALLY_DEGENERATE   // ill-conditioned
};

struct EdgeDiagnostic {
    GEO::index_t va, vb;             // coarse vertex indices
    GEO::index_t face_a, face_b;     // adjacent face indices
    EdgeCompatibility compat;
    double e_endpoint;               // max dist of endpoints to intersection circle
    double e_incidence;              // max ||v - C_f|| - R_f over both faces/verts
    double sphere_dist;              // |C1-C2|
    double intersection_circle_r;    // radius of intersection circle (0 if not PROPER)
};

struct CompatibilityReport {
    int n_interior_edges      = 0;
    int n_boundary_edges      = 0;
    int n_proper_compatible   = 0;
    int n_proper_incompatible = 0;
    int n_coincident          = 0;
    int n_tangent             = 0;
    int n_disjoint            = 0;
    int n_degenerate          = 0;
    double rms_endpoint_err   = 0.0;
    double max_endpoint_err   = 0.0;
    double rms_incidence_err  = 0.0;
    double max_incidence_err  = 0.0;
};

// Diagnose sphere compatibility for all interior coarse edges.
// Does NOT modify any geometry.
CurveKID_API CompatibilityReport diagnose_sphere_compatibility(
    const GEO::Mesh&               coarse_mesh,
    const std::vector<CurvedFace>& cfaces,
    std::vector<EdgeDiagnostic>&   diagnostics,
    double                         endpoint_tol = 1e-4,
    bool                           log_results  = true
);

// Compute E_rendered (ChatGPT Feedback 008):
// E_M2V = RMS distance from original vertices to rendered mesh triangles
// E_V2M = RMS distance from rendered vertices to original mesh triangles
// sphere_dev = RMS distance of rendered vertices from their supporting spheres
CurveKID_API void compute_rendered_errors(
    const GEO::Mesh&               original,
    const GEO::Mesh&               viz_mesh,
    const std::vector<CurvedFace>& cfaces,
    double& rms_M2V, double& max_M2V,
    double& rms_V2M, double& max_V2M,
    double& rms_sphere_dev, double& max_sphere_dev
);

CurveKID_API void compute_planar_approximation_error(
    const GEO::Mesh& original,
    const GEO::Mesh& coarse_flat,
    std::vector<double>& errors
);

// Fit one SphericalFace per mesh triangle.
// Uses: vertex positions + vertex normals (from vcurv) to constrain sphere.
// Falls back to CentralSphere from the most-curved corner when LS fails.
CurveKID_API void initialize_spherical_faces(
    const GEO::Mesh&                    mesh,
    const std::vector<VertexCurvature>& vcurv,
    const std::vector<CentralSphere>&   cspheres,
    std::vector<SphericalFace>&         sfaces
);

} // namespace CurveKIDAlgo
} // namespace OGF

#endif // H_OGF_CurveKID_ALGO_MESH_GEOMETRY_H
