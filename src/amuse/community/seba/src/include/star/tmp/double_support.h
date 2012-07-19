
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  double_support.h: derived class for element evolution systems.
 *          functions as derived class for the real elements.
 *.............................................................................
 *    version 1:  Feb 1994   Simon F. Portegies Zwart
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of class double_support
 *
 *.............................................................................
 */
#ifndef  _DOUBLE_SUPPORT
#  define  _DOUBLE_SUPPORT

#include "stdinc.h"
#include "star_state.h"
       
class double_star;
struct double_init;
enum binary_type  {Strong_Encounter=-1, Unknown_Binary_Type=0, 
		   Synchronized, Detached, 
		   Semi_Detached, Contact, Common_Envelope,
		   Double_Spiral_In,
                   Merged, Disrupted, Spiral_In};

const char * type_string(binary_type);
double_star * get_new_binary(double_init&, const int);
binary_type extract_binary_type_string(const char*);

/*-----------------------------------------------------------------------------
 *  double_hist  -- base struct keeps track of double_star history.
 *-----------------------------------------------------------------------------
 */
struct double_hist
    {
    public:
       real binary_age;
       real semi;
       real eccentricity;

       real donor_timescale;
      
    void put_double_hist();
    void clean() {binary_age=semi=eccentricity=donor_timescale=0;}
    };

/*-----------------------------------------------------------------------------
 *  double_init -- class double_star initial conditions.
 *-----------------------------------------------------------------------------
 */
struct double_init
    {
    public:

       real  start_time;
       real  end_time;
       int   n_steps;

       real  mass_prim;
       real  semi;
       real  q;
       real  eccentricity;

//              More public access functions.
    void clean() {end_time=mass_prim=semi=q=eccentricity=0;}
    void read_element();
    void put_element();
    void dump(ostream &);
    void dump(const char*);

    };

/*-----------------------------------------------------------------------------
 *  double_state -- class double_star initial conditions.
 *-----------------------------------------------------------------------------
 */
struct double_state
    {
      int identity;

      real 		time;
      binary_type	type;
      
      real 		semi;
      real		ecc;
      real		velocity;
//      double_init       init;
      
      real total_mass;

      star_state primary;
      star_state secondary;

      void make_double_state(double_star*, star*, star*);
      void init_double_state(double_star*, star*, star*);
      void clean() {
         type=Detached;
         total_mass=0;
      }
};
double_state make_state(double_star*);
void put_state(double_state, ostream& s = cerr);

/*-----------------------------------------------------------------------------
 *  double_profile -- class double_star profiler.
 *-----------------------------------------------------------------------------
 */
struct double_profile
    {

       double_state init;
       double_state final;

       real mdot;

    void init_double_profile(double_star*, star_state&, star_state&);
    void init_double_profile(double_star*);
    void enhance_double_profile(double_star*, star_state&, star_state&);
    void enhance_double_profile(double_star*);

};

void make_profile(int, real, double_profile&, double_init&);
void put_profile(double_profile&);

// Functions to perform a single binev experiement:
//void  binev(double_star*, real, real, int);

//		Independend initialization and constructor functions.
//double_star * triple_star(double_init&, double_init&, int id=1);

void ppperiod(real period, ostream & s = cerr, const char *p ="Porb");
void pptime(real time, ostream & s = cerr, const char *t = "time");


real period_to_semi(real period, real m_prim, real m_sec);
real semi_to_period(real semi, real m_prim, real m_sec);

#endif 		// _DOUBLE_SUPPORT

