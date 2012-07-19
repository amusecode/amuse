//
// general_support.C
//

#include "general_support.h"
#include "double_star.h"
#include "double_support.h"
#include "star_support.h"

void put_state(double_star* b, 
                star_state& prim, 
                star_state& sec) {

        cout << "status at time = " << b->get_current_time() 
             << " id = " << b->get_identity() << endl;
        cout << type_string(b->get_bin_type());

        switch (b->get_bin_type()) {
           case Merged:
              cout << " (";  
              prim.put_star_state();
              cout << ")" << endl;
              break;
           case Disrupted:
              cout << " )";
              prim.put_star_state();
              cout << ", ";
              sec.put_star_state();
              cout << "(" << endl;
              break;
           default:
              cout << " [a = "<<b->get_semi()<<", e = "<<b->get_eccentricity()
                   << ", v = "<<b->get_velocity() << "]";
              cout << " (";
              prim.put_star_state();
              cout << ", ";
              sec.put_star_state();
              cout << ")" << endl;
           }
     }

void put_initial_conditions(double_init& init) {
   
        cout << "\nInitial conditions: \n";
        if (init.start_time>0)
           cout << "start time = " << init.start_time << endl;
        cout << type_string(Detached)
             << " [a = " << init.semi << ", e = " << init.eccentricity
             << "]" 
             << " (M = " << init.mass_prim << ", m = " 
             << init.q*init.mass_prim << ")" << endl; 

     }
