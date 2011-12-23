#include <iostream>
#include "particle.h"
int main(){
	std::vector<particle> ptcl;
	int idum, nbody;
	std::cin >> idum >> nbody;
	ptcl.reserve(nbody);
	for(int i=0; i<nbody; i++){
		particle::vec pos;
		float h;
		int nnb;
		std::cin >> pos >> h >> nnb;
		h *= 2.f;
		ptcl.push_back(particle(i, pos, h));
		ptcl[i].nnb = nnb;
	}
	for(int i=0; i<nbody; i++){
		int nnb=0;
		for(int j=0; j<nbody; j++){
			// if(i==j) continue;
			float h2 = ptcl[i].h * ptcl[i].h;
			if((ptcl[j].pos - ptcl[i].pos).norm2() < h2) nnb++;
		}
		std::cout << i << " "
			      << ptcl[i].nnb << " " 
				  << nnb << std::endl;
	}
	return 0;
}
