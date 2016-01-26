#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

//----------------------------------------------
//  HashTable
//
//  Maps pointer-sized integers to pointer-sized integers.
//  Uses open addressing with linear probing.
//  In the m_cells array, key = 0 is reserved to indicate an unused cell.
//  Actual value for key 0 (if any) is stored in m_zeroCell.
//  The hash table automatically doubles in size when it becomes 75% full.
//  The hash table never shrinks in size, even after Clear(), unless you explicitly call Compact().
//
// code from 
//  http://preshing.com/20130107/this-hash-table-is-faster-than-a-judy-array/
// adapted and translated to C
//
//----------------------------------------------

#include "simple_hash.h"

#define FIRST_CELL(h) (hash->m_cells + ((h) & (hash->m_arraySize - 1)))
#define CIRCULAR_NEXT(c) ((c) + 1 != hash->m_cells + hash->m_arraySize ? (c) + 1 : hash->m_cells)
#define CIRCULAR_OFFSET(a, b) ((b) >= (a) ? (b) - (a) : hash->m_arraySize + (b) - (a))

#define integerHash(x) (sizeof(x)>4? integerHash64(x):integerHash32(x))


// not performance critical!
size_t upper_power_of_two(size_t x)
{
   size_t y=1;
   while(y<x)y*=2;
   return y;
}

// from code.google.com/p/smhasher/wiki/MurmurHash3
uint32_t integerHash32(uint32_t h)
{
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;
	return h;
}

// from code.google.com/p/smhasher/wiki/MurmurHash3
uint64_t integerHash64(uint64_t k)
{
	k ^= k >> 33;
	k *= 0xff51afd7ed558ccd;
	k ^= k >> 33;
	k *= 0xc4ceb9fe1a85ec53;
	k ^= k >> 33;
	return k;
}

int init_hash(struct simple_hash *hash, size_t initialSize)
{
  hash->m_arraySize = upper_power_of_two(4*initialSize/3+1);
  if( hash->m_arraySize==0 || (hash->m_arraySize & (hash->m_arraySize - 1)) != 0) return -1;   // array must be a power of 2
  hash->m_cells = (struct cell*) malloc(hash->m_arraySize*sizeof(struct cell));
  if(hash->m_cells==NULL) return -2;
  memset(hash->m_cells, 0, sizeof(struct cell) * hash->m_arraySize);
  hash->m_population = 0;

// Initialize zero cell
  hash->m_zeroUsed = false;
  hash->m_zeroCell.key = 0;
  hash->m_zeroCell.value = 0;
  return 0;
}

int end_hash(struct simple_hash *hash) 
{
  if(hash->m_cells) free(hash->m_cells);
  hash->m_cells=NULL;
  return 0;
}

struct cell* hash_cell_lookup(struct simple_hash *hash, size_t key)
{
    if(key)
    {
        for (struct cell* cell_ = FIRST_CELL(integerHash(key));; cell_ = CIRCULAR_NEXT(cell_))
        {
            if (cell_->key == key)
                return cell_;
            if (!cell_->key)
                return NULL;
        }
    }
    else // Check zero cell
    {
        if (hash->m_zeroUsed) return &(hash->m_zeroCell);
        return NULL;
    }
};

int hash_lookup(struct simple_hash *hash, size_t key, size_t *value) 
{
  struct cell * cell_=hash_cell_lookup(hash,key);
  if(cell_) {
    *value=cell_->value;
    return 0;
  } else 
    return -2;
}

int repopulate_hash(struct simple_hash *hash, size_t desiredSize)
{
    if( !((desiredSize & (desiredSize - 1)) == 0)) return -1;   // Must be a power of 2
    if( !(hash->m_population * 4  <= desiredSize * 3) ) return -2;

    // Get start/end pointers of old array
    struct cell *oldCells = hash->m_cells;
    struct cell *end = hash->m_cells + hash->m_arraySize;
    // Allocate new array
    struct cell *new = (struct cell*) malloc(desiredSize*sizeof(struct cell));
    if(new==NULL) return -3;
    
    hash->m_arraySize = desiredSize;
    hash->m_cells = new;
    memset(hash->m_cells, 0, sizeof(struct cell) * hash->m_arraySize);

    // Iterate through old array
    for (struct cell* c = oldCells; c != end; c++)
    {
        if (c->key)
        {
            // Insert this element into new array
            for (struct cell* cc = FIRST_CELL(integerHash(c->key));; cc = CIRCULAR_NEXT(cc))
            {
                if (!cc->key)
                {
                    // Insert here
                    *cc = *c;
                    break;
                }
            }
        }
    }

    // Delete old array
    free(oldCells);
    return 0;
}

int compact_hash(struct simple_hash *hash)
{
    return repopulate_hash(hash,upper_power_of_two((hash->m_population * 4 + 3) / 3));
}


struct cell* hash_cell_insert(struct simple_hash *hash, size_t key)
{
    if (key)
    {
        for (;;)
        {
            for (struct cell* c = FIRST_CELL(integerHash(key));; c = CIRCULAR_NEXT(c))
            {
                if (c->key == key)
                    return c;  // Found
                if (c->key == 0)
                {
                    // Insert here
                    if ((hash->m_population + 1) * 4 >= hash->m_arraySize * 3)
                    {
                        // Time to resize
                        int err = repopulate_hash(hash,hash->m_arraySize * 2);
                        if(err != 0) return NULL;
                        break;
                    }
                    hash->m_population++;
                    c->key = key;
                    return c;
                }
            }
        }
    }
    else  // Check zero cell
    {
        if (!hash->m_zeroUsed)
        {
            // Insert here
            hash->m_zeroUsed = true;
            if (++(hash->m_population) * 4 >= hash->m_arraySize * 3)
			      {
				// Even though we didn't use a regular slot, let's keep the sizing rules consistent
                int err = repopulate_hash(hash,hash->m_arraySize * 2);
                if(err != 0) return NULL;
            }
        }
        return &(hash->m_zeroCell);
    }
}

int hash_insert(struct simple_hash *hash, size_t key, size_t value) 
{
  struct cell *cell_=hash_cell_insert(hash,key);
  if(cell_) {
    cell_->value=value;
    return 0;
  } else {
    return -1;
  }
}

int hash_cell_delete(struct simple_hash *hash, struct cell* cell_)
{
    if (cell_ != &(hash->m_zeroCell))
    {
        // Delete from regular cells
        if( !(cell_ >= hash->m_cells && cell_ - hash->m_cells < hash->m_arraySize) ) return -1;
        if( !cell_->key) return -2;

        // Remove this cell by shuffling neighboring cells so there are no gaps in anyone's probe chain
        for (struct cell* neighbor = CIRCULAR_NEXT(cell_);; neighbor = CIRCULAR_NEXT(neighbor))
        {
            if (!neighbor->key)
            {
                // There's nobody to swap with. Go ahead and clear this cell, then return
                cell_->key = 0;
                cell_->value = 0;
                hash->m_population--;
                return 0;
            }
            struct cell* ideal = FIRST_CELL(integerHash(neighbor->key));
            if (CIRCULAR_OFFSET(ideal, cell_) < CIRCULAR_OFFSET(ideal, neighbor))
            {
                // Swap with neighbor, then make neighbor the new cell to remove.
                *cell_ = *neighbor;
                cell_ = neighbor;
            }
        }
    }
    else
    {
        // Delete zero cell
        if( !hash->m_zeroUsed) return -3;
        hash->m_zeroUsed = false;
        cell_->value = 0;
        hash->m_population--;
        return 0;
    }
}

int hash_delete(struct simple_hash *hash, size_t key)
{
  struct cell *cell_ = hash_cell_lookup(hash,key);
  if (cell_) {
    int ret = hash_cell_delete(hash,cell_);
    if(ret != 0) return -2;
    if(hash->m_population<hash->m_arraySize/8) return compact_hash(hash);
    return 0;
  } else 
    return -1;  
}

int hash_update(struct simple_hash *hash, size_t key, size_t value)
{
  struct cell *cell_ = hash_cell_lookup(hash,key);
  if (cell_) {
    if (cell_) cell_->value=value;
    return 0;
  } else 
    return -1;  
}

int clear_hash(struct simple_hash *hash)
{
    // (Does not resize the array)
    // Clear regular cells
    memset(hash->m_cells, 0, sizeof(struct cell) * hash->m_arraySize);
    hash->m_population = 0;
    // Clear zero cell
    hash->m_zeroUsed = false;
    hash->m_zeroCell.value = 0;
    return 0;
}

int hash_lookups(struct simple_hash *hash,size_t n, size_t *key, size_t *value, int * errors)
{
  int err = 0;
  for(size_t i=0;i<n;i++) {
    errors[i]=hash_lookup(hash, *(key+i), value+i); 
    if(errors[i]) {err = errors[i];}
  }
  return err;
}
int hash_inserts(struct simple_hash *hash,size_t n, size_t *key, size_t *value)
{
  int err;
  for(size_t i=0;i<n;i++) {
    err=hash_insert(hash, *(key+i), *(value+i)); 
    if(err!=0) break;
  }
  return err;
  
}
int hash_updates(struct simple_hash *hash,size_t n, size_t *key, size_t *value)
{
  int err;
  for(size_t i=0;i<n;i++) {
    err=hash_update(hash, *(key+i), *(value+i)); 
    if(err!=0) break;
  }
  return err;
  
}
int hash_deletes(struct simple_hash *hash,size_t n, size_t *key)
{
  int err=0, errout=0;
  for(size_t i=0;i<n;i++) {
    err=hash_delete(hash, *(key+i)); 
    if(err<errout) errout=err;
  }
  return errout;
  
}
