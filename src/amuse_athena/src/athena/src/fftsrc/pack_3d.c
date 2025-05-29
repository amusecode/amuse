/* parallel pack functions - 1998, 1999

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs
   (505) 845-7873
   sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level directory of the distribution.
*/

#include "pack_3d.h"

#if !defined(PACK_POINTER) && !defined(PACK_MEMCPY)
#define PACK_ARRAY
#endif

/* ------------------------------------------------------------------- */
/* Pack and unpack functions:

   pack routines copy strided values from data into contiguous locs in buf
   unpack routines copy contiguous values from buf into strided locs in data
   different versions of unpack depending on permutation
     and # of values/element
   PACK_ARRAY routines work via array indices (default)
   PACK_POINTER routines work via pointers
   PACK_MEMCPY routines work via pointers and memcpy function
*/
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* pack/unpack with array indices */
/* ------------------------------------------------------------------- */

#ifdef PACK_ARRAY

/* ------------------------------------------------------------------- */
/* pack from data -> buf */

void pack_3d(double *data, double *buf, struct pack_plan_3d *plan)

{
  register int in,out,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  in = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_plane;
    for (mid = 0; mid < nmid; mid++) {
      out = plane + mid*nstride_line;
      for (fast = 0; fast < nfast; fast++)
	buf[in++] = data[out++];
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data */

void unpack_3d(double *buf, double *data, struct pack_plan_3d *plan)

{
  register int in,out,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_plane;
    for (mid = 0; mid < nmid; mid++) {
      in = plane + mid*nstride_line;
      for (fast = 0; fast < nfast; fast++)
	data[in++] = buf[out++];
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, one axis permutation, 1 value/element */

void unpack_3d_permute1_1(double *buf, double *data, struct pack_plan_3d *plan)

{
  register int in,out,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      in = plane + mid;
      for (fast = 0; fast < nfast; fast++, in += nstride_plane)
	data[in] = buf[out++];
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, one axis permutation, 2 values/element */

void unpack_3d_permute1_2(double *buf, double *data, struct pack_plan_3d *plan)

{
  register int in,out,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      in = plane + 2*mid;
      for (fast = 0; fast < nfast; fast++, in += nstride_plane) {
	data[in] = buf[out++];
	data[in+1] = buf[out++];
      }
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, one axis permutation, nqty values/element */

void unpack_3d_permute1_n(double *buf, double *data, struct pack_plan_3d *plan)

{
  register int in,out,iqty,instart,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane,nqty;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      instart = plane + nqty*mid;
      for (fast = 0; fast < nfast; fast++, instart += nstride_plane) {
	in = instart;
	for (iqty = 0; iqty < nqty; iqty++) data[in++] = buf[out++];
      }
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, two axis permutation, 1 value/element */

void unpack_3d_permute2_1(double *buf, double *data, struct pack_plan_3d *plan)

{
  register int in,out,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      in = slow + mid*nstride_plane;
      for (fast = 0; fast < nfast; fast++, in += nstride_line)
	data[in] = buf[out++];
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, two axis permutation, 2 values/element */

void unpack_3d_permute2_2(double *buf, double *data, struct pack_plan_3d *plan)

{
  register int in,out,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      in = 2*slow + mid*nstride_plane;
      for (fast = 0; fast < nfast; fast++, in += nstride_line) {
	data[in] = buf[out++];
	data[in+1] = buf[out++];
      }
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, two axis permutation, nqty values/element */

void unpack_3d_permute2_n(double *buf, double *data, struct pack_plan_3d *plan)

{
  register int in,out,iqty,instart,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,nqty;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      instart = nqty*slow + mid*nstride_plane;
      for (fast = 0; fast < nfast; fast++, instart += nstride_line) {
	in = instart;
	for (iqty = 0; iqty < nqty; iqty++) data[in++] = buf[out++];
      }
    }
  }
}

#endif

/* ------------------------------------------------------------------- */
/* pack/unpack with pointers */
/* ------------------------------------------------------------------- */

#ifdef PACK_POINTER

/* ------------------------------------------------------------------- */
/* pack from data -> buf */

void pack_3d(double *data, double *buf, struct pack_plan_3d *plan)

{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  in = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_plane;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+mid*nstride_line]);
      end = begin + nfast;
      for (out = begin; out < end; out++)
	*(in++) = *out;
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data */

void unpack_3d(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_plane;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+mid*nstride_line]);
      end = begin + nfast;
      for (in = begin; in < end; in++)
	*in = *(out++);
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, one axis permutation, 1 value/element */

void unpack_3d_permute1_1(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+mid]);
      end = begin + nfast*nstride_plane;
      for (in = begin; in < end; in += nstride_plane)
	*in = *(out++);
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, one axis permutation, 2 values/element */

void unpack_3d_permute1_2(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+2*mid]);
      end = begin + nfast*nstride_plane;
      for (in = begin; in < end; in += nstride_plane) {
	*in = *(out++);
	*(in+1) = *(out++);
      }
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, one axis permutation, nqty values/element */

void unpack_3d_permute1_n(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out,*instart,*begin,*end;
  register int iqty,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane,nqty;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+nqty*mid]);
      end = begin + nfast*nstride_plane;
      for (instart = begin; instart < end; instart += nstride_plane) {
	in = instart;
	for (iqty = 0; iqty < nqty; iqty++) *(in++) = *(out++);
      }
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, two axis permutation, 1 value/element */

void unpack_3d_permute2_1(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[slow+mid*nstride_plane]);
      end = begin + nfast*nstride_line;
      for (in = begin; in < end; in += nstride_line)
	*in = *(out++);
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, two axis permutation, 2 values/element */

void unpack_3d_permute2_2(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[2*slow+mid*nstride_plane]);
      end = begin + nfast*nstride_line;
      for (in = begin; in < end; in += nstride_line) {
	*in = *(out++);
	*(in+1) = *(out++);
      }
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, two axis permutation, nqty values/element */

void unpack_3d_permute2_n(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out,*instart,*begin,*end;
  register int iqty,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,nqty;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[nqty*slow+mid*nstride_plane]);
      end = begin + nfast*nstride_line;
      for (instart = begin; instart < end; instart += nstride_line) {
	in = instart;
	for (iqty = 0; iqty < nqty; iqty++) *(in++) = *(out++);
      }
    }
  }
}

#endif

/* ------------------------------------------------------------------- */
/* pack/unpack with pointers and memcpy function 
   no memcpy version of unpack_permute routines,
     just use PACK_POINTER versions */
/* ------------------------------------------------------------------- */

#ifdef PACK_MEMCPY

/* ------------------------------------------------------------------- */
/* pack from data -> buf */

void pack_3d(double *data, double *buf, struct pack_plan_3d *plan)

{
  register double *in,*out;
  register int mid,slow,size;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane,upto;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  size = nfast*sizeof(double);
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_plane;
    upto = slow*nmid*nfast;
    for (mid = 0; mid < nmid; mid++) {
      in = &(buf[upto+mid*nfast]);
      out = &(data[plane+mid*nstride_line]);
      memcpy(in,out,size);
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data */

void unpack_3d(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out;
  register int mid,slow,size;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane,upto;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  size = nfast*sizeof(double);
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_plane;
    upto = slow*nmid*nfast;
    for (mid = 0; mid < nmid; mid++) {
      in = &(data[plane+mid*nstride_line]);
      out = &(buf[upto+mid*nfast]);
      memcpy(in,out,size);
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, one axis permutation, 1 value/element */

void unpack_3d_permute1_1(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+mid]);
      end = begin + nfast*nstride_plane;
      for (in = begin; in < end; in += nstride_plane)
	*in = *(out++);
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, one axis permutation, 2 values/element */

void unpack_3d_permute1_2(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+2*mid]);
      end = begin + nfast*nstride_plane;
      for (in = begin; in < end; in += nstride_plane) {
	*in = *(out++);
	*(in+1) = *(out++);
      }
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, one axis permutation, nqty values/element */

void unpack_3d_permute1_n(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out,*instart,*begin,*end;
  register int iqty,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane,nqty;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+nqty*mid]);
      end = begin + nfast*nstride_plane;
      for (instart = begin; instart < end; instart += nstride_plane) {
	in = instart;
	for (iqty = 0; iqty < nqty; iqty++) *(in++) = *(out++);
      }
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, two axis permutation, 1 value/element */

void unpack_3d_permute2_1(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[slow+mid*nstride_plane]);
      end = begin + nfast*nstride_line;
      for (in = begin; in < end; in += nstride_line)
	*in = *(out++);
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, two axis permutation, 2 values/element */

void unpack_3d_permute2_2(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[2*slow+mid*nstride_plane]);
      end = begin + nfast*nstride_line;
      for (in = begin; in < end; in += nstride_line) {
	*in = *(out++);
	*(in+1) = *(out++);
      }
    }
  }
}

/* ------------------------------------------------------------------- */
/* unpack from buf -> data, two axis permutation, nqty values/element */

void unpack_3d_permute2_n(double *buf, double *data, struct pack_plan_3d *plan)

{
  register double *in,*out,*instart,*begin,*end;
  register int iqty,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,nqty;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[nqty*slow+mid*nstride_plane]);
      end = begin + nfast*nstride_line;
      for (instart = begin; instart < end; instart += nstride_line) {
	in = instart;
	for (iqty = 0; iqty < nqty; iqty++) *(in++) = *(out++);
      }
    }
  }
}

#endif
