/* This is a trival function that just returns whether we were
   compiled in testing mode or not. This is to let external callers of
   the dynamic library know which form of the interface to use. */
int testingMode() {
#ifdef TESTING_MODE
  return(1);
#else
  return(0);
#endif
}
