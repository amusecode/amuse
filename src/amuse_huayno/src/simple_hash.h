#ifndef __SIMPLE_HASH_H__
#define __SIMPLE_HASH_H__ 

#include <stdbool.h>

struct cell
{
  size_t key;
  size_t value;
};

struct simple_hash
{
  struct cell* m_cells;
  size_t m_arraySize;
  size_t m_population;
  bool m_zeroUsed;
  struct cell m_zeroCell;
};

int init_hash(struct simple_hash *hash,size_t initialSize);
int end_hash(struct simple_hash *hash);
int clear_hash(struct simple_hash *hash);
int compact_hash(struct simple_hash *hash);

struct cell * hash_cell_insert(struct simple_hash *hash,size_t key);

int hash_lookup(struct simple_hash *hash,size_t key, size_t *value);
int hash_insert(struct simple_hash *hash,size_t key, size_t value);
int hash_update(struct simple_hash *hash,size_t key, size_t value);
int hash_delete(struct simple_hash *hash,size_t key);

int hash_lookups(struct simple_hash *hash,size_t n, size_t *key, size_t *value, int * errors);
int hash_inserts(struct simple_hash *hash,size_t n, size_t *key, size_t *value);
int hash_updates(struct simple_hash *hash,size_t n, size_t *key, size_t *value);
int hash_deletes(struct simple_hash *hash,size_t n, size_t *key);

#endif 
