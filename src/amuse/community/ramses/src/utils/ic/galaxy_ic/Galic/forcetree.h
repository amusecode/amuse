
void force_setpointers(double *_softTable);



void force_treeallocate(int maxnodes);  /* usually maxnodes=2*npart is sufficient */

void force_treefree(void);


int force_treebuild(void **pospointer,int Npart,
		     double thetamax,int costresetflag);



void force_treeevaluate(int target);

void force_treeevaluate_potential(int target);


int force_getcost(void);




void force_testforce(double *targetpart);
