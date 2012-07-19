
// Counters for SeBa:


#ifndef  SEBA_COUNTERS
#  define  SEBA_COUNTERS

typedef unsigned long step_t;

class seba_counters {
  public:

    real cpu_time;

    step_t add_dstar;                         // number of added dstars
    step_t del_dstar;                         // number of deleted dstars

    step_t add_sstar;                         // number of added sstars
    step_t del_sstar;                         // number of deleted dstars

    step_t step_seba;                         // number of seba calls
    step_t step_sstar;                        // number of sstar steps
    step_t step_dstar;                        // number of dstar steps

    step_t detached;                     // number of detached calls
    step_t semi_detached;                // number of semi-detached calls
    step_t contact;                      // number of contact calls

    step_t dynamic;                      // number of dynamc mass steps
    step_t thermal;                      // number of thermal mass steps
    step_t nuclear;                      // number of nuclear mass steps
    step_t aml_driven;                   // number of aml driven mass steps

    step_t sn_in_sstar;                 // number of supernovae in singles
    step_t snwk_in_sstar;               // number of sn with kicks in singles
    step_t sn_in_dstar;                 // number of supernovae in binaries
    step_t snwk_in_dstar;               // number of sn with kicks in binaries

    step_t first_rlof;                        // number of first contacts
    step_t common_envelope;                   // number of common envelope
    step_t spiral_in;                         // number of spiral in
    step_t double_spiral_in;                  // number of double spiral in
    step_t mergers;                           // number of mergers
    step_t aml_mergers;                       // number of mergers due to AML
    step_t gwr_mergers;                       // number of mergers due to GWR

    step_t recursive_overflow;                // number of recursive overflows

  seba_counters () {
      
    cpu_time = 0;

    add_dstar  =
    del_dstar  = 0;

    step_seba    = 
    step_sstar   = 
    step_dstar   = 0;

    detached      =
    semi_detached =	     
    contact       =
    dynamic       =
    thermal       =
    nuclear       =
    aml_driven    = 0;

    sn_in_sstar   = 
    snwk_in_sstar = 
    sn_in_dstar   =
    snwk_in_dstar = 0;

    first_rlof       =
    common_envelope  =
    spiral_in        =
    double_spiral_in =
    mergers          =
    aml_mergers      =
    gwr_mergers      = 0;

    recursive_overflow = 0;
    }
};

void print_counters(seba_counters* sbc, seba_counters* sbc_prev = NULL);
static seba_counters sc_prev;

#endif   //  SEBA_COUNTERS



