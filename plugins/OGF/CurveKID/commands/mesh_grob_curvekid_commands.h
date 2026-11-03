
#ifndef H_OGF_CurveKID_COMMANDS_MESH_GROB_CURVEKID_COMMANDS_H
#define H_OGF_CurveKID_COMMANDS_MESH_GROB_CURVEKID_COMMANDS_H

#include <OGF/CurveKID/common/common.h>
#include <OGF/mesh/commands/mesh_grob_commands.h>

namespace OGF {

    gom_class CurveKID_API MeshGrobCurveKIDCommands
        : public MeshGrobCommands{

    public:
        MeshGrobCurveKIDCommands();
        ~MeshGrobCurveKIDCommands() override;

    gom_slots:

        /**
         * \menu /CurveKID
         * \brief Prints statistics of the selected mesh.
         */
        void inspect_mesh();
    };

}

#endif
