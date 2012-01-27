#ifndef SORT_H
#define SORT_H

inline void Swap(double value[], int index[], int i, int j){
  double tmp;
  int indextmp;
  tmp=value[i];
  value[i]=value[j];
  value[j]=tmp;
  indextmp=index[i];
  index[i]=index[j];
  index[j]=indextmp;
}

inline void Qsort(double value[], int index[], int i, int j){
  int last;
  if(i>=j) return;
  Swap(value,index,i,(i+j)/2);
  last=i;
  for(int k=i+1; k<=j; k++){
    if(value[k]<value[i]){
      Swap(value, index, ++last, k);
    }
  }
  Swap(value,index,i,last);
  Qsort(value,index,i,last-1);
  Qsort(value,index,last+1,j);
}

inline void Swap_index(int index[], const int i, const int j){
  int tmp = index[i];
  index[i]=index[j];
  index[j]=tmp;
}


template<class T> inline void Qsort_index(T *value, int *index, int first, int last){
  T ref_value;
  int ref_index;
  int i, j;
  int tmp;//,tmplast;

//  tmplast = last;		

  ref_index = (first + last)/2;
  ref_value = value[index[ref_index]];

  i = first;
  j = last;

  for( ; ; ){
    while(value[index[i]] < ref_value) i++;
    while(value[index[j]] > ref_value) j--;
    if(i >= j) break;
    tmp = index[i];
    index[i] = index[j];
    index[j] = tmp;
    i++;
    j--;
  }
  if (first < i-1) Qsort_index(value, index, first, i-1);
  if (j+1 < last) Qsort_index(value, index, j+1, last);
}

#endif //SORT_H
