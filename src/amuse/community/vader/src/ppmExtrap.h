/**************************************************************************/
/* Generic PPM fitting method for non-uniform mesh. Given an input        */
/* array of cell sizes dx and values v defined at cell centers, this      */
/* routine compute the values vL and vR at cell left and right edges. The */
/* method follows Collela & Woodward (1984). Note the indexing is         */
/* that, for edge-centered quantities, index i corresponds the left       */
/* face of cell-centered cell i.                                          */
/**************************************************************************/

void ppmExtrap(const int nx, const double *dx, double *v, double *vL, 
	       double *vR, double *wksp);
