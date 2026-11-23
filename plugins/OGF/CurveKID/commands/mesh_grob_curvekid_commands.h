#ifndef H_OGF_CurveKID_COMMANDS_MESH_GROB_CURVEKID_COMMANDS_H
#define H_OGF_CurveKID_COMMANDS_MESH_GROB_CURVEKID_COMMANDS_H

#include <OGF/CurveKID/common/common.h>
#include <OGF/mesh/commands/mesh_grob_commands.h>

namespace OGF {

    gom_class CurveKID_API MeshGrobCurveKIDCommands
        : public MeshGrobCommands {
    public:
        MeshGrobCurveKIDCommands();
        ~MeshGrobCurveKIDCommands() override;

    gom_slots:
        /**
         * \menu /CurveKID
         * \brief Prints basic mesh statistics (vertices, faces, edges, areas).
         */
        void inspect_mesh();

        /**
         * \menu /CurveKID
         * \brief Stage 1: Estimate per-face curvatures and log summary.
         */
        void estimate_curvatures();

        /**
         * \menu /CurveKID
         * \brief Stages 1-5: Grow patches, fit spheres, cluster into K
         *        templates and paint patch_id / template_id / error
         *        as facet attributes.
         * \param nb_patches Target number of connected patches (~100).
         * \param K          Number of template families for K-means.
         */
        void run_curvekid(
            int nb_patches = 100,
            int K          = 10
        );

        /**
         * \menu /CurveKID
         * \brief Validate CurveKID pipeline: print patch/template counts,
         *        connectivity check, histogram, sphere residuals, and
         *        dominant-template analysis. Run run_curvekid first.
         * \param nb_patches Must match the value used in run_curvekid.
         * \param K          Must match the value used in run_curvekid.
         */
        void validate_curvekid(
            int nb_patches = 100,
            int K          = 10
        );

        /**
         * \menu /CurveKID
         * \brief Pottmann-inspired Stage 1: Compute one mean-curvature
         *        (central) sphere per vertex and write radius + validity
         *        as vertex attributes for visualisation in Graphite.
         *        Near-flat vertices (|H| < H_eps) get a 1-ring-based
         *        fallback radius and are flagged valid=0.
         * \param H_eps Threshold below which |H| is treated as near-zero.
         *              Smaller value = fewer fallbacks (default 0.01).
         */
        void inspect_central_spheres(
            double H_eps = 0.01
        );

        /**
         * \menu /CurveKID
         * \brief Render curved patch boundaries: for each boundary edge
         *        (patch_id changes across it) draw a great-circle arc on
         *        the patch's fitted sphere using arc_subdivisions sample
         *        points. Writes a new MeshGrob named "CurveKID_arcs".
         *        Also writes normalised sphere-residual as face attribute
         *        "sphere_quality" (0=bad fit, 1=perfect sphere).
         *        Run run_curvekid first.
         * \param arc_subdivisions Points per boundary edge arc (default 8).
         * \param nb_patches       Must match run_curvekid (default 100).
         * \param K                Must match run_curvekid (default 10).
         */
        void render_curved_patches(
            int arc_subdivisions = 8,
            int nb_patches       = 100,
            int K                = 10
        );

        /**
         * \menu /CurveKID
         * \brief Save CurveKID segmentation to a .ckid file (JSON).
         *        Stores: mesh name, nb_patches, K, per-face patch_id and
         *        template_id, and per-patch sphere (center + radius).
         *        Load later with load_curvekid() without re-running.
         *        Leave filename empty to auto-derive from mesh path.
         * \param filename Output path, e.g. C:/meshes/kitten.ckid
         *   (empty = auto: same folder as mesh, .ckid extension)
         */
        void save_curvekid(
            std::string filename = ""
        );

        /**
         * \menu /CurveKID
         * \brief Load a previously saved .ckid segmentation file and
         *        restore patch_id, template_id and sphere_quality face
         *        attributes on the currently loaded mesh.
         * \param filename Path to the .ckid file to load.
         */
        void load_curvekid(
            std::string filename = ""
        );

        /**
         * \menu /CurveKID
         * \brief Build a NEW curved triangular mesh that approximates the
         *        input surface geometry using spherical curved faces.
         *        Stage 1: coarse curvature-adaptive remesh.
         *        Stage 2: fit one supporting sphere per coarse triangle.
         *        Stage 3: tessellate each face onto its sphere for viz.
         *        Creates a new MeshGrob "CurveKID_CurvedMesh" in the scene.
         *        Also logs approximation error statistics.
         * \param target_faces Target number of coarse triangular faces.
         * \param subdiv       Subdivision level for tessellated visualization.
         */
        void build_curved_mesh(
            int target_faces = 500,
            int subdiv       = 4
        );

        /**
         * \menu /CurveKID
         * \brief PRIMARY EXPERIMENT: build a coarse curved triangular mesh
         *        with N faces and report E_curved(N) vs E_planar(N).
         *        This is the core paper experiment:
         *          Original dense mesh (M faces)
         *          -> New curved triangular mesh (N faces)   [geometric simplification]
         *          -> K reusable templates                   [template reuse, Stage 2]
         *        This command covers ONLY the first reduction: M -> N.
         *        It compares how accurately N curved faces vs N flat faces
         *        approximate the original surface.
         *        Saves:
         *          <mesh>_curved_N.geogram  — tessellated curved viz mesh
         *          <mesh>_flat_N.geogram    — the flat coarse mesh (baseline)
         *        Logs: original face count, N, E_curved(N), E_planar(N).
         * \param N      Target number of curved triangular faces (geometric faces).
         * \param subdiv Subdivision per curved face for visualization only.
         */
        void build_coarse_curved_mesh(
            int    N      = 200,
            int    subdiv = 6,
            double alpha  = 1.0
        );

        /**
         * \menu /CurveKID
         * \brief GROUND TRUTH VALIDATION: build a perfect curved sphere mesh.
         *        Uses ONE exact sphere (fitted to the whole input mesh) for all
         *        coarse faces. Projects all coarse vertices onto the sphere.
         *        Builds shared-edge subdivision — each edge's vertices are
         *        generated once and shared by both adjacent faces, guaranteeing
         *        a watertight, crack-free result.
         *        Use this to verify tessellation correctness before enabling
         *        per-face sphere fitting. Expects a sphere as input.
         *        Validates: topology (Euler=2), watertightness, radial error.
         * \param coarse_mesh_path Path to a valid coarse mesh (.obj or .geogram).
         * \param subdiv           Subdivision level per coarse edge.
         */
        void build_sphere_ground_truth(
            std::string coarse_mesh_path = "",
            int         subdiv           = 4
        );

        /**
         * \menu /CurveKID
         * \brief K=1 ICOSAHEDRAL BASELINE: build 20 congruent spherical
         *        triangular patches on a unit sphere using icosahedral
         *        connectivity. All 20 faces share the SAME supporting sphere
         *        (the unit sphere) → K=1 trivially. All 30 edges are
         *        great-circle arcs → C0 compatible by construction.
         *        Verifies: rigid congruence (K=1), C0 boundary compatibility,
         *        zero self-intersections, fidelity vs dense reference sphere.
         *        This is the positive control for the full CurveKID objective.
         *        (ChatGPT Direction 017 — positive control experiment)
         * \param subdiv  Subdivision per edge for tessellation (default 6).
         * \param radius  Sphere radius (default 1.0).
         */
        void build_icosahedral_baseline(
            int    subdiv = 6,
            double radius = 1.0
        );

        /**
         * \menu /CurveKID
         * \brief PANEL GEOMETRY DIAGNOSTIC: for each coarse face report
         *        sphere parameters, interior displacement, branch validity,
         *        triangle quality. Writes colour-mapped error mesh.
         *        (ChatGPT Direction 018)
         * \param coarse_mesh_path Path to the flat coarse mesh.
         * \param subdiv           Subdivision for tessellation.
         */
        void diagnose_panel_geometry(
            std::string coarse_mesh_path = "",
            int         subdiv           = 6
        );

        /**
         * \menu /CurveKID
         * \brief PHASE A — DR. KID STYLE: construct curved patches via
         *        subdivision + projection onto the original dense mesh.
         *        Does NOT use independent sphere fitting.
         *        Steps: CVT coarse mesh → barycentric subdivision →
         *        closest-point projection onto M_GT → Laplacian smoothing
         *        → reproject → store parent_face_id per vertex/triangle.
         *        Shared edges generated once and shared by both patches.
         *        Returns a geometrically valid curved mesh that lies on
         *        the original surface with zero bulging artefacts.
         *        (ChatGPT Architecture 019, Phase A)
         * \param coarse_mesh_path  Path to the CVT coarse flat mesh.
         * \param subdiv            Subdivision frequency per coarse edge.
         * \param smooth_iters      Laplacian smoothing iterations (default 3).
         * \param smooth_lambda     Smoothing weight 0..1 (default 0.5).
         */
        void construct_reference_curved_patches(
            std::string coarse_mesh_path = "",
            int         subdiv           = 6,
            int         smooth_iters     = 3,
            double      smooth_lambda    = 0.5
        );

        /**
         * \menu /CurveKID
         * \brief TEMPLATE RIGIDITY DIAGNOSTIC (ChatGPT Decision 5):
         *        For each K-means class, extract patches' ordered boundary and
         *        interior vertices, rigidly register each patch to the class
         *        representative (R + t, NO scaling), then report:
         *          - boundary RMS, interior RMS, max mismatch per class
         *          - class occupancy including singleton classes
         *          - overall reuse quality: mean and worst-class boundary RMS
         *        Writes per-class statistics to log and a CSV file.
         *        Run classify_phase_a_patches first (needs template_id attribute).
         * \param phase_a_mesh_path  Path to classified Phase A .geogram mesh (ref_curved_patches_s_classified.geogram).
         * \param coarse_mesh_path   Path to the coarse flat mesh (.obj).
         * \param edgemap_path       Path to the .edgemap.txt file (ref_curved_patches_s6.edgemap.txt). Leave empty to auto-derive.
         * \param subdiv             Subdivision level used in Phase A (must match).
         */
        void diagnose_template_rigidity(
            std::string phase_a_mesh_path = "",
            std::string coarse_mesh_path  = "",
            std::string edgemap_path      = "",
            int         subdiv            = 6
        );

        /**
         * \menu /CurveKID
         * \brief PHASE A EDGES: Extract the actual curved patch boundary polylines
         *        directly from the Phase A subdivided mesh (ref_curved_patches_sN.geogram).
         *        Unlike build_curved_edges (which uses sphere-based great-circle arcs),
         *        this command reads the shared coarse-edge vertex sequences that are
         *        already embedded in the Phase A mesh and emits them as a polyline mesh.
         *        These are the TRUE curved boundaries of each patch — geodesic curves
         *        on M_GT, not sphere arcs. The result is a polyline mesh named
         *        "phase_a_boundary_edges" showing the patch network on the surface.
         * \param phase_a_mesh_path  Path to the Phase A .geogram mesh.
         * \param coarse_mesh_path   Path to the coarse flat mesh (.obj).
         * \param subdiv             Subdivision level used in Phase A (must match).
         */
        void build_phase_a_edges(
            std::string phase_a_mesh_path = "",
            std::string coarse_mesh_path  = "",
            int         subdiv            = 6
        );

        /**
         * \menu /CurveKID
         * \brief PHASE B: Classify Phase A reference curved patches into K
         *        template classes using geometric descriptors (edge lengths,
         *        area, normal variation). Reads parent_face_id from the Phase A
         *        .geogram mesh, runs K-means on 6D descriptors, writes
         *        template_id attribute. Load the output .geogram and paint by
         *        template_id to see K-class colour map.
         *        (ChatGPT Direction 021 — Phase B)
         * \param phase_a_mesh_path  Path to Phase A .geogram mesh (with parent_face_id).
         * \param coarse_mesh_path   Path to the coarse flat mesh (.obj or .geogram).
         * \param K                  Number of template classes (default 10).
         */
        void classify_phase_a_patches(
            std::string phase_a_mesh_path  = "",
            std::string coarse_mesh_path   = "",
            int         K                  = 10
        );

        /**
         * \menu /CurveKID
         * \brief FIT CURVED FACES onto a pre-built coarse mesh.
         *        Use Graphite's Surface->Remesh->Remesh isotropic first to get
         *        a topologically valid coarse mesh, save it, then call this.
         *        The current object is the ORIGINAL dense mesh.
         *        Loads the coarse mesh from coarse_mesh_path (file path).
         *        Reports E_curved(N) vs E_planar(N) and saves curved viz mesh.
         * \param coarse_mesh_path Full path to the coarse mesh .geogram file.
         * \param subdiv           Subdivision per curved face for visualization.
         */
        void fit_curved_to_coarse(
            std::string coarse_mesh_path = "",
            int         subdiv           = 4,
            double      alpha            = 1.0
        );

        /**
         * \menu /CurveKID
         * \brief STAGE 2 — K-SET CLASSIFICATION: cluster N curved faces into
         *        K reusable template classes using pose-invariant descriptors.
         *        Each face descriptor: (sphere_radius, arc_lengths[3], solid_angle).
         *        Runs K-means clustering and writes template_id as face attribute
         *        on the curved viz mesh for colour-coded visualisation.
         *        Run build_coarse_curved_mesh or fit_curved_to_coarse first.
         *        Loads curved mesh from curved_mesh_path, saves result with
         *        template_id attribute as <curved_mesh>_classified.geogram.
         * \param curved_mesh_path Path to the curved viz mesh .geogram file.
         * \param K                Number of template classes (default 10).
         */
        void classify_curved_faces(
            std::string curved_mesh_path = "",
            std::string flat_mesh_path   = "",
            int         K                = 10
        );

        /**
         * \menu /CurveKID
         * \brief SYNTHETIC mode: build a curved mesh where each coarse
         *        triangle is lifted onto a circumsphere (sphere through
         *        all 3 corner vertices). Corner positions are EXACTLY
         *        preserved. Adjacent faces share the same circular arc on
         *        the lower-residual face's sphere. Use this to verify that
         *        spherical interiors and curved edges render correctly before
         *        moving to the geometry-fitting mode.
         * \param target_faces Target number of coarse triangular faces.
         * \param t_factor     Normal offset as a fraction of circumradius.
         *                     Controls curvature: 0=flat, 1=hemisphere-ish,
         *                     try 0.5 for clearly visible bulge.
         * \param subdiv       Subdivision level for tessellated visualization.
         */
        void build_curved_mesh_synthetic(
            int    target_faces = 200,
            double t_factor     = 0.5,
            int    subdiv       = 6
        );

        /**
         * \menu /CurveKID
         * \brief CURVED EDGES ONLY: render the 3 great-circle arc edges of
         *        each curved face as a polyline mesh (no interior tessellation).
         *        Each edge is approximated by arc_pts points via slerp on the
         *        face's supporting sphere. Use this to visualise the pure
         *        curved-triangle boundary structure without subdivision fill.
         *        Saves <mesh>_arcs_N.obj and .geogram.
         * \param coarse_mesh_path Path to the flat coarse mesh .geogram file.
         * \param arc_pts          Points per arc edge (default 20).
         */
        void build_curved_edges(
            std::string coarse_mesh_path = "",
            int         arc_pts          = 20
        );

        /**
         * \menu /CurveKID
         * \brief TORUS EXPERIMENT: generate four variants A/B/C/D from the
         *        same coarse mesh and report quantitative comparison.
         *        A=planar, B=independent sphere, C=joint opt, D=error-selected.
         *        Saves <coarse>_varA/B/C/D.obj and logs error table.
         * \param coarse_mesh_path Path to the coarse flat mesh (.obj/.geogram).
         * \param subdiv           Subdivision for tessellation.
         */
        void torus_experiment(
            std::string coarse_mesh_path = "",
            int         subdiv           = 6
        );

        /**
         * \menu /CurveKID
         * \brief ALPHA SWEEP: evaluate E_rendered for alpha=0,0.2,0.4,0.6,0.8,1.0
         *        on a fixed coarse mesh. Reports E_M2V, E_V2M, sphere_dev for each.
         *        (ChatGPT Research Direction 008)
         * \param coarse_mesh_path Path to the flat coarse mesh.
         * \param subdiv           Subdivision for tessellation.
         */
        void alpha_sweep(
            std::string coarse_mesh_path = "",
            int         subdiv           = 6
        );

        /**
         * \menu /CurveKID
         * \brief DIAGNOSTIC: test sphere compatibility for all interior coarse
         *        edges. Reports PROPER_COMPATIBLE vs PROPER_INCOMPATIBLE counts,
         *        endpoint-to-circle error, sphere incidence error.
         *        Does NOT modify any geometry.
         *        (ChatGPT Direction 011)
         * \param coarse_mesh_path Path to the flat coarse mesh.
         * \param endpoint_tol     Tolerance for endpoint-on-circle test.
         */
        void diagnose_sphere_compatibility(
            std::string coarse_mesh_path = "",
            double      endpoint_tol     = 1e-4
        );
    };

} // namespace OGF

#endif // H_OGF_CurveKID_COMMANDS_MESH_GROB_CURVEKID_COMMANDS_H
