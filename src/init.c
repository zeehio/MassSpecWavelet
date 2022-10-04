#include "MassSpecWavelet.h"

static const R_CallMethodDef callMethods[]  = {
    {"c_distance_to_next_greater", (DL_FUNC) &distance_to_next_greater, 1},
    {"c_findLocalMaximum", (DL_FUNC) &findLocalMaximum, 1},
    {NULL, NULL, 0}
};


void R_init_MassSpecWavelet(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}