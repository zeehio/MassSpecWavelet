#include "MassSpecWavelet.h"

static inline int sign_d(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

static void findLocalMaximum_impl_d(double *x, R_xlen_t xlength, int *outi) {

#ifndef MASSSPECWAVELET_DEBUG
  #define MASSSPECWAVELET_DEBUG 0
#endif
    int in_peak = 0;
    int peak_starts = 0;
    int peak_ends = 0;
    int peak_center = 0;
    int winsize;
    int j;
    if (xlength == 0) {
        return;
    } else if (xlength == 1) {
        outi[0] = 0;
        return;
    }
    int *stack_prev = (int*) R_alloc(xlength, sizeof(int));
    int *stack_next = (int*) R_alloc(xlength, sizeof(int));
    int stack_prev_size = 0;
    int stack_next_ends = 0;
    int stack_next_starts = 0;
    int prev_diff = -2;
    int curr_diff = -2;
    for (int i=0; i < xlength-1; i++) {
        curr_diff = sign_d((x[i+1] - x[i]) > 0);
        for (j = stack_next_starts; j < stack_next_ends; j++) {
            if (stack_next[j] > i) {
                stack_next_starts = j;
                break;
            }
        }
        if (MASSSPECWAVELET_DEBUG) printf("i: %d (stack_prev_size: %d). ", i, stack_prev_size);
        switch(prev_diff) {
        case -2:
            if (MASSSPECWAVELET_DEBUG) printf("signal starts ");
            switch (curr_diff) {
            case -1:
                if (MASSSPECWAVELET_DEBUG) printf("and decreases\n");
                // signal starts decreasing
                outi[i] = 0;
                stack_prev[stack_prev_size++] = i;
                break;
            case  0:
                if (MASSSPECWAVELET_DEBUG) printf("and holds\n");
                // signal starts constant
                stack_prev[stack_prev_size++] = i;
                outi[i] = 0;
                break;
            case  1:
                if (MASSSPECWAVELET_DEBUG) printf("and increases\n");
                // signal starts increasing
                outi[i] = 0;
                stack_prev[stack_prev_size++] = i;
                break;
            }
            break;
        case -1:
            if (MASSSPECWAVELET_DEBUG) printf("signal decreases ");
            switch (curr_diff) {
            case -1:
                // signal was decreasing and is decreasing. 
                if (MASSSPECWAVELET_DEBUG) printf("and decreases\n");
                in_peak = 0;
                outi[i] = 0;
                stack_prev[stack_prev_size++] = i;
                break;
            case  0:
                // signal was decreasing and is stable
                if (MASSSPECWAVELET_DEBUG) printf("and holds\n");
                in_peak = 0;
                outi[i] = 0;
                stack_prev[stack_prev_size++] = i;
                break;
            case  1:
                // signal was decreasing and is increasing (minimum)
                if (MASSSPECWAVELET_DEBUG) printf("and increases\n");
                in_peak = 0;
                outi[i] = 0;
                stack_prev[stack_prev_size++] = i;
                break;
            }
            break;
        case 0:
            if (MASSSPECWAVELET_DEBUG) printf("signal holds ");
            switch (curr_diff) {
            case -1:
                if (MASSSPECWAVELET_DEBUG) printf("and decreases\n");
                // signal was stable and decreases. If peak, set center and peak ends now.
                if (!in_peak) {
                    outi[i] = 0;
                } else {
                    peak_ends = i;
                    peak_center = (peak_starts+peak_ends)/2;
                    winsize = 1;
                    // to the left:
                    for (j = stack_prev_size-1; j >= 0; j--) {
                        if (x[stack_prev[j]] > x[peak_center]) {
                            if (MASSSPECWAVELET_DEBUG) printf("  winsize_left: %d\n", peak_center - stack_prev[j] - 1);
                            winsize += peak_center - stack_prev[j] - 1;
                            break;
                        }
                    }
                    if (j == -1) {
                        winsize += peak_center; // CHECK: -1?
                    }
                    // to the right:
                    for (j = stack_next_starts; j < stack_next_ends; j++) {
                        if (x[stack_next[j]] > x[peak_center]) {
                            if (MASSSPECWAVELET_DEBUG) printf("  winsize_right: %d (a)\n", stack_next[j] - peak_center-1);
                            winsize += stack_next[j] - peak_center-1;
                            break;
                        }
                    }
                    if (j == stack_next_ends) {
                        // we reached the end of the stack, keep looking
                        int last_checked = (stack_next_ends == stack_next_starts) ? i: stack_next[stack_next_ends-1];
                        for (j = last_checked+1; j<xlength;j++) {
                            // FIXME: Performance: Append to stack_next if convenient so next time the stack_next is larger:
                            // FIXME: Performance: Use this loop to discard points to check
                            if (x[j] > x[peak_center]) {
                                if (MASSSPECWAVELET_DEBUG) printf("  winsize_right: %d (b)\n", j - peak_center-1);
                                winsize += j - peak_center-1;
                                break;
                            }
                        }
                        if (j == xlength) {
                            if (MASSSPECWAVELET_DEBUG) printf("  winsize_right: %ld (c)\n", xlength - peak_center-1);
                            winsize += xlength - peak_center-1;
                        }
                    }
                    // set winsize
                    outi[i] = 0;
                    outi[peak_center] = winsize;
                    stack_prev[stack_prev_size++] = i;
                    in_peak = 0;
                }
                break;
            case  0:
                if (MASSSPECWAVELET_DEBUG) printf("and holds\n");
                outi[i] = 0;
                // signal was stable and is stable. Maybe peak.
                break;
            case  1:
                if (MASSSPECWAVELET_DEBUG) printf("and increases\n");
                // signal was stable and increases. If peak, cancel peak.
                outi[i] = 0;
                in_peak = 0;
                break;
            }
            break;
        case 1:
            if (MASSSPECWAVELET_DEBUG) printf("signal increases and ");
            switch (curr_diff) {
            case -1:
                if (MASSSPECWAVELET_DEBUG) printf("decreases\n");
                // signal increased and decreases
                in_peak = 1; // unneded, but nice
                peak_starts = i;
                peak_ends = i;
                peak_center = i;
                winsize = 1;
                // to the left:
                for (j = stack_prev_size-1; j >= 0; j--) {
                    if (MASSSPECWAVELET_DEBUG) printf("   stack_prev[%d]: %d\n", j, stack_prev[j]);
                    if (x[stack_prev[j]] > x[peak_center]) {
                        if (MASSSPECWAVELET_DEBUG) printf("       winsize_left: %d\n", peak_center - stack_prev[j] - 1);
                        winsize += peak_center - stack_prev[j] - 1;
                        break;
                    }
                }
                // to the right:
                for (j = stack_next_starts; j < stack_next_ends; j++) {
                    if (x[stack_next[j]] > x[peak_center]) {
                        if (MASSSPECWAVELET_DEBUG) printf("       winsize_right: %d (a)\n", stack_next[j] - peak_center-1);
                        winsize += stack_next[j] - peak_center-1;
                        break;
                    }
                }
                if (j == stack_next_ends) {
                    // we reached the end of the stack, keep looking
                    int last_checked = (stack_next_ends == stack_next_starts) ? i: stack_next[stack_next_ends-1];
                    for (j = last_checked+1; j<xlength;j++) {
                        // FIXME: Performance: Append to stack_next if convenient:
                        // FIXME: Performance: Use this loop to discard points to check
                        if (x[j] > x[peak_center]) {
                            if (MASSSPECWAVELET_DEBUG) printf("       winsize_right: %d (b)\n", j - peak_center-1);
                            winsize += j - peak_center-1;
                            break;
                        }
                    }
                    if (j == xlength) {
                        if (MASSSPECWAVELET_DEBUG) printf("       winsize_right: %ld (c)\n", xlength - peak_center-1);
                        winsize += xlength - peak_center-1;  // CHECK: -1?
                    }
                }
                outi[peak_center] = winsize;
                stack_prev[stack_prev_size++] = i;
                in_peak = 0;
                break;
            case  0:
                if (MASSSPECWAVELET_DEBUG) printf("holds\n");
                // signal increased and stabilizes
                peak_starts = i;
                in_peak = 1;
                outi[i] = 0;
                break;
            case  1:
                if (MASSSPECWAVELET_DEBUG) printf("increases\n");
                // signal is increasing
                in_peak = 0;
                outi[i] = 0;
                break;
            }
            break;
        }
        prev_diff = curr_diff;
    }
    if (xlength > 0)
        outi[xlength-1] = 0;
    fflush(stdout);
    return;
}

SEXP findLocalMaximum(SEXP s_x) {
    double *xd;
    int *xi;
    R_xlen_t xlength = Rf_length(s_x);
    int is_int = TYPEOF(s_x) == INTSXP;
    if (is_int) {
        xi = INTEGER(s_x);
    } else {
        if (TYPEOF(s_x) != REALSXP) {
            Rf_error("x must be real or integer");
        }
        xd = REAL(s_x);
    }
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, xlength));
    int *outi = INTEGER(out);
    if (is_int) {
        Rf_error("Not yet implemented, please coerce to double");
    } else {
        findLocalMaximum_impl_d(xd, xlength, outi);
    }
    Rf_unprotect(1);
    return out;
}