/*
 * fastCenter.c
 * copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
 */

#include <R.h>
#include <Rdefines.h>

/**
 * Subtract a given vector `mu` from every row of a given matrix `X`.
 */
void fastRowMinus(double *X, int *rows, int *cols, double *mu) {
    for (; *cols > 0; (*cols)--) {
        for (int r = *rows; r > 0; r--) {
            *X++ -= *mu;
        }
        mu++;
    }
}


/**
 * Add a given vector `mu` to every row of a given matrix `X`.
 */
void fastRowPlus(double *X, int *rows, int *cols, double *mu) {
    for (; *cols > 0; (*cols)--) {
        for (int r = *rows; r > 0; r--) {
            *X++ += *mu;
        }
        mu++;
    }
}
