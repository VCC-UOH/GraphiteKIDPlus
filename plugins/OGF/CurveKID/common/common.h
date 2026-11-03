
#ifndef H_OGF_CurveKID_COMMON_COMMON_H
#define H_OGF_CurveKID_COMMON_COMMON_H

#include <OGF/basic/common/common.h>

#ifdef CurveKID_EXPORTS
#define CurveKID_API OGF_EXPORT
#else
#define CurveKID_API OGF_IMPORT
#endif

namespace OGF {

    static class CurveKID_API CurveKID_libinit {
    public:
        CurveKID_libinit();
        ~CurveKID_libinit();

        static void increment_users();
        static void decrement_users();

    private:
        static void initialize();
        static void terminate();
        static int count_;
    } CurveKID_libinit_instance;

}

#endif
