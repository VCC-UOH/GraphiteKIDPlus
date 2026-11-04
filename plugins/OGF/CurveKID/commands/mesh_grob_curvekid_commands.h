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
         * \brief Prints basic mesh statistics.
         */
        void inspect_mesh();

        /**
         * \menu /CurveKID
         * \brief Stage 1: Estimate per-face curvatures and log summary.
         */
        void estimate_curvatures();

        /**
         * \menu /CurveKID
         * \brief Run all 5 stages: grow patches, fit spheres, cluster
         *        into K templates, write patch_id / template_id /
         *        sphere_residual as facet attributes.
         * \param nb_patches Target number of connected patches (~100).
         * \param K          Number of template families for K-means.
         */
        void run_curvekid(
            int nb_patches = 100,
            int K          = 10
        );
    };

} // namespace OGF

#endif