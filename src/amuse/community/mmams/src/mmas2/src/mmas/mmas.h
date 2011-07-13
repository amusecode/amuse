#include "usm/usm.h"
#include <vector>

class mmas {
protected:
  real v_inf, r_p;
  usm *model_a;
  usm *model_b;
  usm *product;
  usm *mixed_product;
  
public: 
  real dm_lost;
 
  /* constructor destructor */
  
  mmas(usm &m_a,
       usm &m_b,
       real v0 = 0,
       real r0 = 0) {
    r_p     = r0;
    v_inf   = v0;
    model_a = &m_a;
    model_b = &m_b;
    product = new usm;
    mixed_product = new usm;
  }
  mmas(real v0 = 0,
       real r0 = 0) {
    r_p     = r0;
    v_inf   = v0;
    model_a = new usm;
    model_b = new usm;
    product = new usm;
    mixed_product = new usm;
  }

  ~mmas() {
    delete model_a;
    delete model_b;
    delete product;
    delete mixed_product;
  };
  
  mmas& operator=(const mmas &original) {
      r_p = original.r_p;
      v_inf = original.v_inf;
      *(model_a) = original.get_model_a();
      *(model_b) = original.get_model_b();
      *(product) = original.get_product();
      *(mixed_product) = original.get_mixed_product();
      dm_lost = original.dm_lost;
      return *this;
  }
  
  void set_model_a(usm &ma) {model_a = &ma;}
  void set_model_b(usm &mb) {model_b = &mb;}
  usm& get_model_a() const {return *model_a;}
  usm& get_model_b() const {return *model_b;}

  usm& get_product() const {return *product;}
  usm& get_mixed_product() const {return *mixed_product;}

  /* mmas functions */

  void compute_extra();
  void sort_model(usm&, vector<real>&, vector<real>&, vector<real>&, vector<real>&);

  int shock_heating(real);
  int shock_heating(usm&, real, real);
  int shock_heating_3();
  int shock_heating_3(usm&, real, real, real);
  int shock_heating_4();
  int shock_heating_4(usm&, real, real, real, real);
  
  real get_lagrad(usm&, real);
  int  merge_stars(real, int);
  void merge_stars_consistently(int, int);
  real mass_loss();
  real compute_stellar_energy(usm &model);
  real compute_orbital_energy();

  void mixing();
  void mixing(usm&);
  void mixing_product(int n = 200);
  void smooth_product();

  void solve_hse(usm& model, int mode);

};
