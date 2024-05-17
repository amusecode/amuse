#include <unistd.h>

float _x[10] = { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};

int do_sleep(int in) {
    sleep(in);
    return 0;
}

int return_error(int * out) {
    *out=123;
    return -1;
}

int get_x(int in, float *x) {
    *x=_x[in];
    return 0;
}

int set_x(int in, float x) {
    _x[in]=x;
    return 0;
}

int dummy() {
    return 0;
}

