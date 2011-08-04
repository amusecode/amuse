#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main(argc,argv)
int argc;
char *argv[];
{
	int i, nobj;
	float r[7];
	FILE *rvfile;

	if( argc == 2 ) {
		if( !(rvfile = fopen(argv[1],"r")) ) {
			fprintf(stderr,"can't open %s\n",argv[1]);
		}
	}
	else {
		fprintf(stderr,"usage: toascii binary_rvfile > ascii_rvfile\n");
		exit(0);
	}

	fseek(rvfile,0,2); nobj = ftell(rvfile)/28; fseek(rvfile,0,0);

	fprintf(stdout,"%d\n",nobj);
	for(i=0; i<nobj; i++) {
		fread(r,1,28,rvfile);
		fprintf(stdout,"% 13.8e % 13.8e % 13.8e % 13.8e % 13.8e % 13.8e % 13.8e\n",
			r[0], r[1], r[2], r[3], r[4], r[5], r[6]);
	}
	fclose(rvfile);
}


