#include "include/stdinc.h"
#include <vector>
#include <algorithm>

class data_compare : public std::binary_function<int, int, bool> {

   const vector<real> &data;

public:

  data_compare(const vector<real> &d) : data(d) {};
  ~data_compare() {};
  
  bool operator() (int a, int b) const {
    return data[a] < data[b];
  };
  
};

vector<int> sort_index (vector<real>);

