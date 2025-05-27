#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

int echo_int(int int_in, int * int_out) {
    *int_out = int_in;
    if(int_in < 0) {
        return -1;
    } else {
        return 0;
    }
}
int echo_long_long_int(long long int int_in, long long int * int_out) {
    *int_out = int_in;
    if(int_in < 0) {
        return -1;
    } else {
        return 0;
    }
}

int echo_double(double in, double * out) {
    *out = in;
    return 0;
}

int echo_float(float in, float * out) {
    *out = in;
    return 0;
}
int echo_string(char * in, char ** out) {
    *out = in;
    return 0;
}

int echo_string_int(int inint, char * in, char ** out) {
    *out = in;
    return 0;
}

int echo_string_two(char * in1, char * in2, char ** out1, char ** out2) {
    *out1 = in1;
    *out2 = in2;
    return 0;
}

int print_string(char * in) {
    fprintf(stdout, "%s\\n", in);
    return 0;
}

int print_error_string(char * in) {
    fprintf(stderr, "%s\\n", in);
    return 0;
}

int echo_strings(char ** inout1, char ** inout2) {
    char * tmp;
    tmp = *inout1;
    *inout1 = *inout2;
    *inout2 = tmp;
    
    return 0;
}


void echo_array(int * in, int * out, int len) {
    int i = 0;
    for(i = 0; i < len; i++) {
        out[i] = in[i];
    }
}

int echo_array_with_result(int * in, int *out, int len) {
    int i = 0;
    for(i = 0; i < len; i++) {
        out[i] = in[i];
    }
    return -1;
}
        
int echo_2_int(int * int_in1, int * int_in2, int * int_out1, int * int_out2, int len) {
    int i = 0;
    for(i = 0; i < len; i++) {
        int_out1[i] = int_in1[i];
        int_out2[i] = int_in2[i];
    }
    
    return len;
}
int echo_3_int(int * i, int * j, int * k, int * l, int * m, int * int_out, int len) {
    int x = 0;
    for(x = 0; x < len; x++) {
        int_out[x] = i[x];
    }    
    return len;
}

int dummy_3_int(int i, int j, int k) {
    return 0;
}

int echo_inout_array_with_result(int * inout, int len) {
    int i = 0;
    for(i = 0; i < len; i++) {
        inout[i] = inout[i] + 10;
    }
    return 11;
}


int echo_logical(bool in, bool * out) {
    *out = in;
    return 0;
}
/*
int echo_string_array(char ** in, char ** out, int len) {
    int x = 0;
    for(x = 0; x < len; x++) {
        out[x] = in[x];
    }    
    return len;
}
*/


int sum_doubles(double in1, double in2, double * out) {
    *out = in1 + in2;
    return 0;
}

