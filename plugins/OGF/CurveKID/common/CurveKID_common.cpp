
#include <OGF/CurveKID/common/common.h>

#include <OGF/basic/modules/module.h>
#include <OGF/gom/types/gom_defs.h>
#include <OGF/scene_graph/types/scene_graph_library.h>

#include <OGF/CurveKID/commands/mesh_grob_curvekid_commands.h>

namespace OGF {

    void CurveKID_libinit::initialize() {

        Logger::out("Init")
            << "<CurveKID>" << std::endl;

        gom_package_initialize(CurveKID);

        ogf_register_grob_commands<
            OGF::MeshGrob,
            OGF::MeshGrobCurveKIDCommands
        >();

        Module* module_info = new Module;

        module_info->set_name("CurveKID");
        module_info->set_vendor("OGF");
        module_info->set_version("1.0");
        module_info->set_info(
            "K-set near-isometric curved patch decomposition"
        );

        Module::bind_module("CurveKID", module_info);
    }

    void CurveKID_libinit::terminate() {
        Module::unbind_module("CurveKID");
    }

    CurveKID_libinit::CurveKID_libinit() {
        increment_users();
    }

    CurveKID_libinit::~CurveKID_libinit() {
        decrement_users();
    }

    void CurveKID_libinit::increment_users() {
        count_++;
        if (count_ == 1) {
            initialize();
        }
    }

    void CurveKID_libinit::decrement_users() {
        count_--;
        if (count_ == 0) {
            terminate();
        }
    }

    int CurveKID_libinit::count_ = 0;

}

extern "C" void CurveKID_API OGF_CurveKID_initialize(void);

extern "C" void CurveKID_API OGF_CurveKID_initialize() {
    OGF::CurveKID_libinit::increment_users();
}

extern "C" void CurveKID_API OGF_CurveKID_terminate(void);

extern "C" void CurveKID_API OGF_CurveKID_terminate() {
    OGF::CurveKID_libinit::decrement_users();
}
