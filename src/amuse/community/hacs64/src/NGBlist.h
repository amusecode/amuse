#ifndef __NGBlist_H__
#define __NGBlist_H__

#include <cstdio>
#include "localassert.h"
#include <vector>

struct NGBlist
{
	public:
		enum {NGB_MAX = 512};
	private:
		int n_list;
		int list[NGB_MAX];

	public:
		void print(const int i, FILE *fp = stderr) const
		{
			fprintf(fp, "%6d%6d :", i, n_list);
			for(int k = 0; k < n_list; k++)
				fprintf(fp, " %d", list[k]);
			fprintf(fp, "\n");
			fflush(fp);
		}


		NGBlist() : n_list(0) {}
		~NGBlist() {}	

		int operator[](const int i) const {return list[i];}
		int&      operator[](const int i)       {return list[i];}
		int size()     const {return n_list;};
		void clear() {n_list = 0;}
		void resize(const int n) {n_list = n;}
		void push_back(const int j) {assert(n_list < NGB_MAX); list[n_list++] = j;}
};

#endif // __NGBlist_H__
