
//// rdc_SeBa: read SeBa short dump
////           output requested binary parameters.
////          
//// Options:     none
//-----------------------------------------------------------------------------
//   version 1:  Sept 1998   Simon Portegies Zwart   spz@grape.c.u-tokyo.ac.jp
//                                                   University of Tokyo
//.............................................................................
//   non-local functions: 
//-----------------------------------------------------------------------------

#include "stdinc.h"
#include "node.h"
#include "double_star.h"
#include "main_sequence.h"
//
// SeBa_hist.h
//

#ifndef    _SeBa_HIST
#  define  _SeBa_HIST

// (SPZ+GN:  1 Aug 2000) Dangerous enum names... 
// i.e., current_time = 2 but not lvalue
enum binary_parameter {identity=0, bin_type, mtrtype,
		       current_time,
		       primary_mass, primary_radius,
		       secondary_mass, secondary_radius,
		       semi_major_axis, eccentricity, mass_ratio
                      };

/*-----------------------------------------------------------------------------
 * SeBa_hist  --  a linked list of SeBa histories.
 *-----------------------------------------------------------------------------
 */
class SeBa_hist {
    protected:

    int  number;
    real time;

    binary_type bin_tpe;
    mass_transfer_type mttype;

    real semi;
    real ecc;

    char label_prim[255];
    stellar_type tpe_prim;
    real m_prim;
    real r_prim;
    
    char label_sec[255];
    stellar_type tpe_sec;
    real m_sec;
    real r_sec;

    SeBa_hist * past;
    SeBa_hist * future;

    public:
       SeBa_hist(SeBa_hist* s=NULL) {
	   if (s) {
  	      past=s;
	      past->future=this;
	      future = NULL;
	   }
	   else {
	      past=NULL;
	      future=NULL;
	  }

	   number = 0;
	   time   = 0;
	   bin_tpe = Detached;
	   mass_transfer_type Unknown;
	   strcpy(label_prim, "1a");
	   strcpy(label_sec, "1b");
	   tpe_prim = tpe_sec = Main_Sequence;
	   m_prim=m_sec=r_prim=r_sec=0;
       }
       ~SeBa_hist(){

         if (future!=NULL) {
	     SeBa_hist *tmp = future;
	     future = NULL;
	     delete tmp;
	 }

	 if (past)
	    past->future = NULL;
       }

       SeBa_hist* get_past() {return past;}
       void set_past(SeBa_hist *sb) {past = sb;}
       SeBa_hist* get_future() {return future;}
       SeBa_hist* get_first() {
           if (past!=NULL)
              return past->get_first();
           else
              return this;
       }
       SeBa_hist* get_last() {
           if (future!=NULL)
              return future->get_last();
           else 
              return this;
       }

       void set_future(SeBa_hist* f) {future = f;}
       void set_last(SeBa_hist* f) {
            if(future!=NULL) 
              future->set_last(f);
            else 
              future=f;
       }

       real get_time()              {return time;}
       int get_number()             {return number;}
       int set_number(int n)        {number = n;}

       real set_stellar_radius(bool);
       void move_SeBa_hist_to(SeBa_hist*);
       bool read_SeBa_hist(istream&);

       void put_history(ostream&, bool);
       void put_single_reverse(ostream&);
       void put_first_formed_left(char*, real);

       bool binary_contains(char*, char *, binary_type, mass_transfer_type);
       bool binary_limits(binary_parameter, real, real);
       real get_parameter(binary_parameter);
       binary_type get_binary_type() { return bin_tpe;}
       stellar_type get_primary_type() { return tpe_prim;}
       stellar_type get_secondary_type() { return tpe_sec;}
       char* get_label_prim() {return label_prim;}
       char* get_label_sec() {return label_sec;}

       SeBa_hist* get_SeBa_hist_at_time(real);

    void add_to_SeBa_hist(SeBa_hist *next_hi);
				  
       void put_state(ostream&);
       friend ostream& operator<<(ostream& s, SeBa_hist&);

       bool operator == (SeBa_hist& ha) const;
       bool operator != (SeBa_hist& ha) const;

//       friend SeBa_hist& operator == (SeBa_hist&);
//       friend istream& operator>>(istream& s, cluster_table& table);
};

#define for_all_SeBa_hist(SeBa_hist, base, SeBa_hist_next)                    \
        for (SeBa_hist* SeBa_hist_next = base;                                \
	     SeBa_hist_next != NULL;                            \
	     SeBa_hist_next = SeBa_hist_next->get_future())


SeBa_hist* get_history(SeBa_hist *hi, istream& is);

bool scenarios_identical(SeBa_hist* hi, SeBa_hist* ha);
void put_state(SeBa_hist * hi, ostream & s);

#endif
