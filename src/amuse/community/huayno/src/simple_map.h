#ifndef __SIMPLE_MAP_H__
#define __SIMPLE_MAP_H__ 

struct simple_map
{
  size_t max_id;
  size_t *pindex;
};

int init_map(struct simple_map *map,size_t initialSize);
int end_map(struct simple_map *map);
int clear_map(struct simple_map *map);
int compact_map(struct simple_map *map);

int map_lookup(struct simple_map *map,size_t key, size_t *value);
int map_insert(struct simple_map *map,size_t key, size_t value);
int map_update(struct simple_map *map,size_t key, size_t value);
int map_delete(struct simple_map *map,size_t key);

#endif 
