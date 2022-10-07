#include "MassSpecWavelet.h"

#ifndef MASSSPECWAVELET_DEBUG
#define MASSSPECWAVELET_DEBUG 0
#endif


static inline int sign_d(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

static inline void append_stack_next(int pos, int *stack_next, int *stack_next_ends) {
    if (*stack_next_ends == 0 || stack_next[*stack_next_ends-1] < pos) {
        stack_next[*stack_next_ends] = pos;
        *stack_next_ends = *stack_next_ends +1; 
    }
    return;
}


static inline void append_stack_prev(int pos, int *stack_prev, int *stack_prev_size) {
    if (*stack_prev_size == 0 || stack_prev[*stack_prev_size-1] < pos) {
        stack_prev[*stack_prev_size++] = pos;
    }
    return;
}

static inline void append_plateau_start(int pos, int *stack_plateau, int *stack_plateau_size) {
    if (*stack_plateau_size == 0 || stack_plateau[*stack_plateau_size-1] < pos) {
        #if (MASSSPECWAVELET_DEBUG > 0) 
            printf("  -> adding %d to stack_plateau begins\n", pos);
        #endif
        if (*stack_plateau_size % 2 == 0) {
            stack_plateau[*stack_plateau_size] = pos;
            *stack_plateau_size  =  *stack_plateau_size + 1;
        } else {
            Rf_error("Internal bug in MassSpecWavelet:::append_plateau_start: Tried adding a plateau start in an end position. Please report it.");
        }
    }
    return;
}

static inline void append_plateau_end(int pos, int *stack_plateau, int *stack_plateau_size) {
    // add if
    if (*stack_plateau_size % 2 == 1 && stack_plateau[*stack_plateau_size-1] < pos) {
        #if (MASSSPECWAVELET_DEBUG > 0) 
            printf("  -> adding %d to stack_plateau ends (a)\n", pos);
        #endif
        stack_plateau[*stack_plateau_size] = pos;
        *stack_plateau_size = *stack_plateau_size + 1;
    }
    return;
}

static inline int get_half_left_window(double *x, int peak_center, int *stack_prev, int stack_prev_starts, int capWinSize) {
    int winleft;
    for (int j = stack_prev_starts; j >= 0; j--) {
        winleft = peak_center - stack_prev[j] - 1;
        if (winleft > capWinSize) {
            return capWinSize;
        }
        if (x[stack_prev[j]] > x[peak_center]) {
            #if (MASSSPECWAVELET_DEBUG)
                printf("  winsize_left: %d (a)\n", winleft);
            #endif
            return winleft;
        }
    }
    #if (MASSSPECWAVELET_DEBUG > 0)
        printf("  winsize_left: %d (c)\n", peak_center);
    #endif
    return peak_center;
}

static inline int get_half_right_window_stack(double *x, int peak_center, int capWinSize, int *stack_next, int stack_next_starts, int stack_next_ends) {
    for (int j = stack_next_starts; j < stack_next_ends; j++) {
        int winsizeright = stack_next[j] - peak_center - 1;
        if (winsizeright >= capWinSize) {
            return capWinSize;
        } else if (x[stack_next[j]] > x[peak_center]) {
            #if (MASSSPECWAVELET_DEBUG > 0)
                printf("  winsize_right: %d (a)\n", winsizeright);
            #endif
            return winsizeright;
        }
    }
    return -1;
}

void remove_begin_from_plateau(int *stack_plateau, int *stack_plateau_size) {
    if (*stack_plateau_size %2 == 1) {
        #if (MASSSPECWAVELET_DEBUG > 0)
            printf("  -> removing %d from stack_plateau (hi)\n",stack_plateau[*stack_plateau_size-1]);
        #endif
        stack_plateau[*stack_plateau_size-1] = 0;
        *stack_plateau_size = *stack_plateau_size - 1;
    }
}


static inline void set_stacks(int j,
           int peek_prev_diff, int peek_curr_diff,
           int *stack_prev, int *stack_prev_size,
           int *stack_next, int *stack_next_ends,
           int *stack_plateau, int *stack_plateau_size,
           int *not_a_peak) {
    switch(peek_prev_diff) {
    case -1: // decreasing
        append_stack_prev(j, stack_prev, stack_prev_size);
        not_a_peak[j] = 1;
        break;
    case 0: // holds
        switch(peek_curr_diff) {
        case -1: // and decreases
            append_stack_prev(j, stack_prev, stack_prev_size);
            append_plateau_end(j, stack_plateau, stack_plateau_size);
            not_a_peak[j] = 0; // if it's the end of a plateau, we set the peak at the plateau center here
            break;
        case 0: // and keeps holding
            not_a_peak[j] = 1; // centers of plateaus are checked at the end
            break;
        case 1: // and increases
            // it's not a peak, cancel the plateau
            remove_begin_from_plateau(stack_plateau, stack_plateau_size);
            not_a_peak[j] = 1;
            break;
        }
        break;
    case 1: // increases
        append_stack_next(j, stack_next, stack_next_ends);
        switch(peek_curr_diff) {
        case -1: // and decreases
            // that's a peak
            append_stack_prev(j, stack_prev, stack_prev_size);
            not_a_peak[j] = 0;
            break;
        case 0: // and holds
            append_plateau_start(j, stack_plateau, stack_plateau_size);
            not_a_peak[j] = 1;
            break;
        case 1: // and increases
            not_a_peak[j] = 1;
            break;
        }
        break;
    }
    
}


int get_winsize(double *x, int xlength,
                int pos, int peak_center,
                int *stack_prev, int *stack_prev_starts, int *stack_prev_size,
                int capWinSize, int *stack_next, int *stack_next_starts, int *stack_next_ends,
                int *stack_plateau, int *stack_plateau_size, int *not_a_peak) {
    int winsize = 1;
    // to the left:
    winsize += get_half_left_window(x, peak_center, stack_prev, *stack_prev_starts, capWinSize-winsize);
    
    // to the right:
    int winsize_right = get_half_right_window_stack(x, peak_center, capWinSize - winsize, stack_next, *stack_next_starts, *stack_next_ends);
    if (winsize_right < 0) {
        // stack_next did not have the end, so we explore the signal:
        int last_checked = (*stack_next_ends == *stack_next_starts) ? pos: stack_next[*stack_next_ends-1];
        #if (MASSSPECWAVELET_DEBUG > 0)
            printf("    reached end of stack_next, last_checked: %d\n", last_checked);
        #endif
        int peek_prev_diff;
        int peek_curr_diff;
        peek_prev_diff = sign_d(x[last_checked+1] - x[last_checked]);
        // we reached the end of the stack, keep looking
        int j;
        for (j = last_checked+1; j<xlength;j++) {
            peek_curr_diff = sign_d(x[j+1] - x[j]);
            #if (MASSSPECWAVELET_DEBUG > 0)
                printf("    peeking %d (%d, %d)\n", j, peek_prev_diff, peek_curr_diff);
            #endif
            set_stacks(j, peek_prev_diff, peek_curr_diff, stack_prev, stack_prev_size, stack_next, stack_next_ends, stack_plateau, stack_plateau_size, not_a_peak);
            winsize_right = j - peak_center - 1;
            if (winsize_right >= capWinSize - winsize) {
                break;
            } else if (x[j] > x[peak_center]) {
                #if (MASSSPECWAVELET_DEBUG > 0)
                    printf("  winsize_right: %d (b)\n", winsize_right);
                #endif
                break;
            }
            peek_prev_diff = peek_curr_diff;
        }
        if (j == xlength) {
            winsize_right = xlength - peak_center - 1;
            #if (MASSSPECWAVELET_DEBUG > 0)
                printf("  winsize_right: %d (c)\n", winsize_right);
            #endif
        }
    }
    winsize += winsize_right;
    return winsize;
}

static void findLocalMaximum_impl_d(double *x, R_xlen_t xlength, int *outi, int capWinSize) {

    int in_plateau = 0;
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
    if (capWinSize == R_NaInt) {
        capWinSize = INT_MAX;
    }
    int *stack_prev = (int*) R_alloc(xlength, sizeof(int));
    int *stack_next = (int*) R_alloc(xlength, sizeof(int));
    int *not_a_peak = (int*) R_alloc(xlength, sizeof(int));
    memset(not_a_peak, 0, xlength*sizeof(int));
    int *stack_plateau = (int*) R_alloc(xlength, sizeof(int)); // start-end-start-end indices of the plateau (start=first point in plateau, end = last point in plateau)
    int stack_plateau_size = 0;
    int stack_plateau_starts = 0;
    int stack_prev_starts = 0;
    int stack_prev_size = 0;
    int stack_next_ends = 0;
    int stack_next_starts = 0;
    int prev_diff = -2;
    int curr_diff = -2;
    
    curr_diff = sign_d(x[1] - x[0]);
    #if (MASSSPECWAVELET_DEBUG > 0)
    printf("signal starts and %s\n", curr_diff == -1 ? "decreases" : (curr_diff == 0 ? "holds" : "increases"));
    #endif
    outi[0] = 0;
    append_stack_prev(0, stack_prev, &stack_prev_size);
    prev_diff = curr_diff;
    
    for (int i=1; i < xlength-1; i++) {
        #if (MASSSPECWAVELET_DEBUG > 0)
            printf("plateau limits:\n");
            for (j = 0; j< stack_plateau_size; j++) {
                printf(" - %s: %d\n", (j%2 == 0) ? "starts": "ends", stack_plateau[j]);
            }
            printf("plateau summary ends:\n");
        #endif

#if (MASSSPECWAVELET_DEBUG > 0)
            printf("stack_next [starts: %d - ends: %d]\n", stack_next_starts, stack_next_ends);
            for (int kk=stack_next_starts; kk<stack_next_ends;++kk) {
                printf("  stack_next[%d] = %d\n", kk, stack_next[kk]);
            }
            printf("stack_next ends\n");
#endif
            
#if (MASSSPECWAVELET_DEBUG > 0) 
            printf("i: %d (stack_prev_size: %d). ", i, stack_prev_size);
#endif
            
        if (not_a_peak[i]) {
            outi[i] = 0;
            prev_diff = sign_d(x[i+1] - x[i]);
            #if (MASSSPECWAVELET_DEBUG > 0)
                printf("Skipping\n");
            #endif
            continue;
        }
        curr_diff = sign_d(x[i+1] - x[i]);
        for (j = stack_next_starts; j < stack_next_ends; j++) {
            if (stack_next[j] > i) {
                stack_next_starts = j;
                break;
            }
        }
        if (j == stack_next_ends) {
            stack_next_starts = stack_next_ends;
        }
        for (j = stack_prev_starts; j < stack_prev_size; j++) {
            if (stack_prev[j] >= i) {
                stack_prev_starts = j-1;
                break;
            }
        }
        if (j == stack_prev_size) {
            stack_prev_starts = stack_prev_size - 1;
        }
        
        for (j = stack_plateau_starts; j < stack_plateau_size-1; j = j+2) {
            if (i >= stack_plateau[j] && i <= stack_plateau[j+1]) {
                in_plateau = 1;
                stack_plateau_starts = j;
                break;
            }
        }
        if (stack_plateau_size % 2 == 1 && j >=stack_plateau_size-1) {
            in_plateau = 1;
            stack_plateau_starts = j;
        }
        if (MASSSPECWAVELET_DEBUG) printf("in_plateau: %d. ",in_plateau);
#if (MASSSPECWAVELET_DEBUG > 0)
        printf("stack_next [starts: %d - ends: %d]\n", stack_next_starts, stack_next_ends);
        for (int kk=stack_next_starts; kk<stack_next_ends;++kk) {
            printf("  stack_next[%d] = %d\n", kk, stack_next[kk]);
        }
        printf("stack_next ends\n");
#endif
        
        
        switch(prev_diff) {
        case -1:
            #if (MASSSPECWAVELET_DEBUG > 0)
                printf("signal decreases and %s\n", curr_diff == -1 ? "decreases" : (curr_diff == 0 ? "holds" : "increases"));
            #endif
            outi[i] = 0;
            append_stack_prev(i, stack_prev, &stack_prev_size);
            break;
        case 0:
            #if (MASSSPECWAVELET_DEBUG > 0)
                printf("signal holds and %s\n", curr_diff == -1 ? "decreases" : (curr_diff == 0 ? "holds" : "increases"));
            #endif
            switch (curr_diff) {
            case -1:
                // signal was stable and decreases. If peak, set center and peak ends now.
                // future points will need to check this one for the border:
                if (stack_prev_size == 0 || stack_prev[stack_prev_size-1] < i) {
                    stack_prev[stack_prev_size++] = i;
                }
                if (!in_plateau) {
                    outi[i] = 0;
                } else {
                    append_plateau_end(i, stack_plateau, &stack_plateau_size);
                    peak_starts = stack_plateau[stack_plateau_starts];
                    peak_ends = i;
                    #if (MASSSPECWAVELET_DEBUG > 0)
                        printf("  peak_starts: %d, peak_ends: %d\n", peak_starts, peak_ends);
                    #endif
                    peak_center = (peak_starts+peak_ends)/2;
                    winsize = get_winsize(x, xlength, i, peak_center,
                                          stack_prev, &stack_prev_starts, &stack_prev_size,
                                          capWinSize, stack_next, &stack_next_starts, &stack_next_ends,
                                          stack_plateau, &stack_plateau_size, not_a_peak);
                    // set winsize
                    outi[i] = 0;
                    outi[peak_center] = winsize < capWinSize ? winsize : capWinSize;
                }
                break;
            case  0:
                outi[i] = 0;
                // signal was stable and is stable. Maybe peak.
                break;
            case  1:
                // signal was stable and increases. If peak, cancel peak.
                outi[i] = 0;
                if (in_plateau && stack_plateau_size > 0) {
                    remove_begin_from_plateau(stack_plateau, &stack_plateau_size);
                }
                break;
            }
            break;
        case 1:
            #if (MASSSPECWAVELET_DEBUG > 0)
                printf("signal increases and %s\n", curr_diff == -1 ? "decreases" : (curr_diff == 0 ? "holds" : "increases"));
            #endif
            switch (curr_diff) {
            case -1:
                append_stack_prev(i, stack_prev, &stack_prev_size);
                // signal increased and decreases
                peak_center = i;
                winsize = get_winsize(x, xlength, i, peak_center,
                                      stack_prev, &stack_prev_starts, &stack_prev_size,
                                      capWinSize, stack_next, &stack_next_starts, &stack_next_ends,
                                      stack_plateau, &stack_plateau_size, not_a_peak);
                // set winsize
                outi[peak_center] = winsize < capWinSize ? winsize : capWinSize;
                break;
            case  0:
                // signal increased and stabilizes
                append_plateau_start(i, stack_plateau, &stack_plateau_size);
                outi[i] = 0;
                break;
            case  1:
                // signal is increasing
                outi[i] = 0;
                break;
            }
            break;
        }
        prev_diff = curr_diff;
    }
    outi[xlength-1] = 0;
    #if (MASSSPECWAVELET_DEBUG > 0)
        printf("final plateau limits:\n");
        for (j = 0; j< stack_plateau_size; j++) {
            printf(" - %s: %d\n", (j%2 == 0) ? "starts": "ends", stack_plateau[j]);
        }
        printf("final plateau summary ends:\n");
        fflush(stdout);
    #endif
    return;
}

SEXP findLocalMaximum(SEXP s_x, SEXP s_capWinSize) {
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
    if (TYPEOF(s_capWinSize) != INTSXP) {
        Rf_error("capWinSize must be an integer");
    }
    int capWinSize = INTEGER(s_capWinSize)[0];
    
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, xlength));
    int *outi = INTEGER(out);
    if (is_int) {
        Rf_error("Not yet implemented, please coerce to double");
    } else {
        findLocalMaximum_impl_d(xd, xlength, outi, capWinSize);
    }
    Rf_unprotect(1);
    return out;
}