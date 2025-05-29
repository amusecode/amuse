#include <stdio.h>
#include <stdlib.h>

double grid0;
double grid1[10];
double grid2[10][10];

int set0(double in) {
    grid0=in;
    return 0;
}
int get0(double *out) {
    *out=grid0;
    return 0;
}
int get_grid0_range() {
    return 0;
}

