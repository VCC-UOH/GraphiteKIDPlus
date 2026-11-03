
#include <OGF/CurveKID/commands/mesh_grob_curvekid_commands.h>

#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>

namespace OGF {

    MeshGrobCurveKIDCommands::MeshGrobCurveKIDCommands() {}

    MeshGrobCurveKIDCommands::~MeshGrobCurveKIDCommands() {}

    void MeshGrobCurveKIDCommands::inspect_mesh() {

        MeshGrob* grob = mesh_grob();

        if (grob == nullptr) {
            GEO::Logger::err("CurveKID")
                << "No mesh selected." << std::endl;
            return;
        }

        GEO::Mesh& mesh = *grob;

        GEO::index_t triangles = 0;

        for (
            GEO::index_t f = 0;
            f < mesh.facets.nb();
            ++f
            ) {
            if (mesh.facets.nb_vertices(f) == 3) {
                ++triangles;
            }
        }

        GEO::Logger::out("CurveKID")
            << "Vertices: "
            << mesh.vertices.nb()
            << std::endl;

        GEO::Logger::out("CurveKID")
            << "Facets: "
            << mesh.facets.nb()
            << std::endl;

        GEO::Logger::out("CurveKID")
            << "Triangles: "
            << triangles
            << std::endl;
    }

}
