#include "mmas.h"
#include "sort.h"

vector<int> sort_index (vector<real> data) {
  vector<int> index_vector;
  int n = data.size();
  for (int i = 0; i < n; i++) index_vector.push_back(i);

  sort(index_vector.begin(), index_vector.end(), data_compare(data));

  return index_vector;
  
}

void mmas::sort_model(usm &model, 
		      vector<real> &mass,
		      vector<real> &entr,
		      vector<real> &m_mu,
		      vector<real> &id_s) {
  vector<real> entropy;
  vector<real> dm;
  real m_prev = 0;
  for (int i = 0; i < model.get_num_shells(); i++) {
    mass_shell &shell = model.get_shell(i);
    entropy.push_back(shell.entropy);
    dm.push_back(shell.mass - m_prev);
    m_prev = shell.mass;
  }

  vector<int> index = sort_index(entropy);

  mass.clear();
  entr.clear();
  m_mu.clear();
  id_s.clear();

  real m_cur = 0;
  for (unsigned int i = 0; i < index.size(); i++) {
    int j = index[i];
//     j = i;
    mass_shell &shell = model.get_shell(j);
    m_cur += dm[j];
    mass.push_back(m_cur);
    entr.push_back(entropy[j]);
    m_mu.push_back(shell.mean_mu);
    id_s.push_back(j);
  }
  
}

