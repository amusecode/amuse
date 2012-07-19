
// printbits:  print significant nonzero bits of an unsigned integer.

#include "stdinc.h"

void printbits(unsigned int i)
{
    int n = 31;
    while (n >= 0) {if (GETBIT(i, n)) break; n--;}
    while (n >= 0) {if (GETBIT(i, n)) cerr << "1"; else cerr << "0"; n--;}
}
