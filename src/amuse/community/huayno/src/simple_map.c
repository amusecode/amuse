#include <stddef.h>
#include <stdlib.h>
#include "simple_map.h"

int init_map(struct simple_map *m, size_t n)
{
  m->max_id=n;
  m->pindex=(size_t *) malloc(m->max_id*sizeof(size_t));
  if(m->pindex==NULL) return -1;
  for(size_t i=0;i<m->max_id;i++) m->pindex[i]=0;
  return 0;
}

int end_map(struct simple_map *m)
{
  m->max_id=0;
  free(m->pindex);
  return 0;
}

int clear_map(struct simple_map *m)
{
  for(size_t i=0;i<m->max_id;i++) m->pindex[i]=0;
  return 0;
}

int compact_map(struct simple_map *m)
{  
  size_t oldmax=m->max_id;
  size_t *oldindex=m->pindex;
  m->max_id=1;
  for(size_t i=0;i<oldmax;i++) if(oldindex[i]>0 && i+1>m->max_id) m->max_id=i+1;
  m->pindex=(size_t *) malloc(m->max_id*sizeof(size_t));
  if(m->pindex==NULL) return -1;
  for(size_t i=0;i<m->max_id;i++) m->pindex[i]=oldindex[i];
  free(oldindex);
  return 0;
}

int map_update(struct simple_map *m, size_t id,size_t p)
{
  if(id>=m->max_id) return -1;
  if(m->pindex[id]==0) return -2;
  m->pindex[id]=p+1;
  return 0;
}

int map_insert(struct simple_map *m, size_t id,size_t p)
{
 while(id>=m->max_id)
 {
   size_t *new;
   new=(size_t *) realloc( m->pindex, 2*m->max_id*sizeof(size_t) );
   if(new == NULL) return -2;
   for(size_t i=m->max_id;i<2*m->max_id;i++) m->pindex[i]=0;
   m->max_id*=2;
   m->pindex=new;
 }
 m->pindex[id]=p+1;
 return 0;
}

int map_delete(struct simple_map *m, size_t id)
{
  if(id>=m->max_id) return -1;
  if(m->pindex[id]==0) return -2;
  m->pindex[id]=0;
}

int map_lookup(struct simple_map *m, size_t id,size_t *p)
{
  if(id>=m->max_id) return -1;
  *p=m->pindex[id];
  if(*p==0 ) return -2;
  *p-=1;
  return 0; 
}
