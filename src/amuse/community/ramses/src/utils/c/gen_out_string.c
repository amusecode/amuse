// gen_out_string
//=============================================================================
// Author: Michael Rieder (mr@student.ethz.ch)
//
// outputs a string with output_ directory names inside
// given boundaries
//
// example: ./gen_out_string 1 2 5
// gives this:
// output_00001 output_00003 output_00005 
//=============================================================================

#include <stdlib.h>
#include <string.h>
#include <stdio.h>


// Global Variables
#define NUM_DIGITS 5
#define PRE_STRING "output_"

int calc_item_size(const char *pre_string, int num_digits) {
  int result;

  result = strlen(pre_string) + num_digits;

  return result;
}

int calc_size(int startval, int increment, int endval) {
  int result;

  result = endval - startval;
  result = result / increment;

  if ( result < 0 )
    return 0;

  result = result + 1;
  result = result * (calc_item_size(PRE_STRING,NUM_DIGITS) + 1);

  // for the terminating 0 at the end
  result = result + 1;

  return result;
}
 
int main (int argc, char *argv[])
{
  int startval, increment, endval;
  int str_size;
  char *str_buffer;
  int str_pos;
  int i;

	// Init
  if (argc != 4) {
    printf( "Usage: %s startval increment endval\n", argv[0] );
    return 1;
  }

  // get values
  startval = atoi(argv[1]);
  increment = atoi(argv[2]);
  endval = atoi(argv[3]);

  if ( increment == 0 ) {
    printf( "Error: increment is 0!\n" );
    return 1;
  }

  // prepare string buffer
  str_size = calc_size(startval, increment, endval);
  str_buffer = malloc(str_size);
  if ( str_buffer == 0 ) {
    printf( "Error: malloc is 0!\n" );
    return 1;
  }

  // fill buffer
  str_pos = 0;
  memset(str_buffer,0,str_size);
  if (increment>0)
    for (i=startval;i<=endval;i+=increment) {
      sprintf(str_buffer+str_pos,"%s%0*i ",PRE_STRING,NUM_DIGITS,i);

      str_pos = str_pos + calc_item_size(PRE_STRING,NUM_DIGITS)+1;
    }
  else
    for (i=startval;i>=endval;i+=increment) {
      sprintf(str_buffer+str_pos,"%s%0*i ",PRE_STRING,NUM_DIGITS,i);

      str_pos = str_pos + calc_item_size(PRE_STRING,NUM_DIGITS)+1;
    }
  // terminate string at the end
  str_buffer[str_size-2] = 0;

	// done
	printf( "%s\n", str_buffer);
	free( str_buffer );
	return 0;
}

