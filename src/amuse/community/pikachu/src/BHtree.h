#ifndef BHTREE_H
#define BHTREE_H

#include"Matrix3.h"
//#include"p2m2.h"
#include"particle.h"
#include"force.h"
#include"mpi_interface.h"
#include"global.h"

#define dump_cout(x) cout<<#x<<" = "<<x<<endl;
#define dump_cerr(x) cerr<<#x<<" = "<<x<<endl;

inline int octant_index(const Vector3 &nodepos, 
			const Vector3 &ppos){
    register int i0,i1,i2;
    i0= (nodepos[0]<ppos[0])? 4:0;
    i1= (nodepos[1]<ppos[1])? 2:0;
    i2= (nodepos[2]<ppos[2])? 1:0;
    return i0+i1+i2;
}

// separation between point and point
inline double separation_squared_point_point(const Vector3 &pos0, 
					     const Vector3 &pos1){
    double r2 = (pos0-pos1)*(pos0-pos1);
    return r2;
}

// shortest separation between box and point
inline double shortest_separation_squared_point_box(const Vector3 &pos_center_box, 
						    const Vector3 &len_box, 
						    const Vector3 &pos){
    double r2 = 0.0;
    double adx, ady, adz;
    adx = fabs(pos_center_box[0] - pos[0]) - 0.5*len_box[0];
    ady = fabs(pos_center_box[1] - pos[1]) - 0.5*len_box[1];
    adz = fabs(pos_center_box[2] - pos[2]) - 0.5*len_box[2];
    if (adx>0.0){r2 += adx*adx;}
    if (ady>0.0){r2 += ady*ady;}
    if (adz>0.0){r2 += adz*adz;}
    return r2;
}

// shortest separation between boxes
inline double shortest_separation_squared_box_box(const Vector3 &pos_center_box0, 
						  const Vector3 &len_box0, 
						  const Vector3 &pos_center_box1, 
						  const Vector3 &len_box1){
    double r2=0.0;
    double dx = fabs(pos_center_box0[0] - pos_center_box1[0]) - 0.5*(len_box0[0] + len_box1[0]);
    double dy = fabs(pos_center_box0[1] - pos_center_box1[1]) - 0.5*(len_box0[1] + len_box1[1]);
    double dz = fabs(pos_center_box0[2] - pos_center_box1[2]) - 0.5*(len_box0[2] + len_box1[2]);
    if (dx>0.0){r2 += dx*dx;}
    if (dy>0.0){r2 += dy*dy;}
    if (dz>0.0){r2 += dz*dz;}
    return r2;
}

inline int are_overlapped_point_box(const Vector3 &pos_center_box0, 
				    const Vector3 &len_box0, 
				    const Vector3 &pos){
    Vector3 dx = pos_center_box0 - pos;
    Vector3 len2 = len_box0 * 0.4999999999999;  
    for(int k=0;k<3;k++){
	if( fabs(dx[k]) > len2[k]) return 0; // not over lapped
    }
    return 1; // over lapped
}

inline int are_overlapped_box_box(const Vector3 &pos_center_box0, 
				  const Vector3 &len_box0,
				  const Vector3 &pos_center_box1, 
				  const Vector3 &len_box1){
    Vector3 len = (len_box0 + len_box1) * 0.4999999999999;  
    Vector3 dx = pos_center_box0 - pos_center_box1;
    for(int k=0;k<3;k++){
	if( fabs(dx[k]) > len[k]) return 0; // not overlapped
    }
    return 1; // over lapped
}


class Cell_Tree{
private:

public:
    Cell_Tree *child[8];
    Particle *prt_first;
    int Nprt;
    int isleaf;
    Vector3 cmpos;
    double cmmass;

    int open;
    int level;
    int Nchild;
    Vector3 pos_center;
    Vector3 size;

    double p2mass[3];
    Vector3 p2pos[3];
    Matrix3 quad;

    Cell_Tree(){

        for(int k=0; k<8; k++){
            child[k] = NULL;
        }
        prt_first = NULL;
        Nprt = 0;
        isleaf = 1;
        cmpos = 0.0;
        cmmass = 0.0;
        open = 0; 
	level = 0;
	Nchild = 0;
        pos_center = 0.0;
	size = 0.0;
        for(int k=0; k<3; k++){
            p2mass[k] = 0.0;
            p2pos[k] = 0.0;
        }
        quad = 0.0;
    }

    void clear(){

        for(int i = 0; i<8;i++){
	    child[i] = NULL;
	}
        prt_first = NULL;
        Nprt = 0;
        isleaf = 1;
        cmpos = 0.0;
        cmmass = 0.0;
        open = 0;
	level = 0;
	Nchild = 0;
	pos_center = 0.0;
	size = 0.0;

        for(int k=0; k<3; k++){
            p2mass[k] = 0.0;
            p2pos[k] = 0.0;
        }
        quad = 0.0;
    }

    void dump(ostream& fout = cout){
	fout<<"pos_center="<<pos_center<<endl;
	fout<<"size="<<size<<endl;

	fout<<"Nprt="<<Nprt<<endl;
	fout<<"level="<<level<<endl;
	fout<<"Nchild="<<Nchild<<endl;
	fout<<"isleaf="<<isleaf<<endl;
	fout<<"open="<<open<<endl;
	/*
        for(int k=0; k<3; k++){
	    fout<<"p2mass[k]"<<p2mass[k]<<endl;
	    fout<<"p2pos[k]"<<p2pos[k]<<endl;
        }
	*/
	fout<<endl;  
    }


    void sanity_check(ostream& fout= cout){
        fout<<endl;
        fout<<"cell info"<<endl;
        dump(fout);
        fout<<endl;
        if(isleaf){
            int Ntmp=0;
            for (Particle *prt_tmp = prt_first; prt_tmp != NULL; prt_tmp = prt_tmp->prt_next){
                Ntmp++;
                fout<<"prt_tmp->index="<<prt_tmp->index<<endl;
                fout<<"prt_tmp->pos="<<prt_tmp->pos<<endl;

                int innode = 0;
                for(int k=0; k<3; k++){
		    if(prt_tmp->pos[k] > pos_center[k] - size[k]*0.5 && prt_tmp->pos[k] < pos_center[k]+size[k]*0.5){
                        innode++;
                    }
                }
                if(innode == 3){
                    fout<<"in cell"<<endl;
                }
                else{
                    fout<<"out of cell...."<<endl;
                }
            }
            fout<<endl;
            fout<<"Ntmp="<<Ntmp<<endl;

            fout<<"-----------------------------"<<endl;
            fout<<endl;
            fout<<endl;
        }
        else{
            for(int i=0; i<8; i++){
                if(child[i] != NULL){
                    child[i]->sanity_check(fout);
                }
            }
        }
    }
     

    void clear_recursive(){
        if(isleaf){clear();}
        else{
            for(int i=0; i<8; i++){
                if(child[i] != NULL){ 
                    child[i]->clear_recursive();
                }
            }
        }
    }

    void add_particles_to_essential_tree(Vector3 pos_list_send[], 
					 double mass_list_send[], 
					 int &nlist, 
					 const int &LIST_COMM_MAX,
					 int index_list_send[]);

    void add_to_essential_tree(const Vector3 &poscenter, // center of pos in box_dest
			       const Vector3 &boxsize, // len of box_dest
			       const double &theta2,
			       Vector3 pos_list_send[],
			       double mass_list_send[],
			       int &list_len,
			       const int &LIST_COMM_MAX,
			       const double &r2_ngh,
			       int index_list_send[],
			       const int &quad_flag,
			       Particle BH_dst[],
			       const int &NBH_dst, 
			       const double &rsearch2_FS_BH); 


    void search_neighbour(const Particle &prt_tree_i, 
                          const double &r2_ngh, 
                          Neighbour_List ngh_list[],
                          int &ngh_list_len,
                          const int &ngh_list_max);

    void search_neighbour_index(Particle &prti,
                                const double &r2_ngh, 
                                Neighbour_List_Index ngh_list_idx[],
                                int &ngh_list_len,
                                const int& ngh_list_max,
				const int& flag_min_box);

};

#endif //TREE_H

