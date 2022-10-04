#include "MassSpecWavelet.h"

static const R_CallMethodDef callMethods[]  = {
    {"c_findLocalMaximum", (DL_FUNC) &findLocalMaximum, 1},
    {NULL, NULL, 0}
};


void R_init_MassSpecWavelet(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}