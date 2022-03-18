extern int echo(int input);

/*
 * Interface code
 */
 
int echo_int(int input, int * output){
    *output = echo(input);
    return 0;
}

