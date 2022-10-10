#include "MassSpecWavelet.h"

static const R_CallMethodDef callMethods[]  = {
    {"c_findLocalMaxWinSize", (DL_FUNC) &findLocalMaxWinSize, 2},
    {"c_localMaximumSlidingWindow", (DL_FUNC) &localMaximumSlidingWindow, 2},
    {NULL, NULL, 0}
};


void R_init_MassSpecWavelet(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}