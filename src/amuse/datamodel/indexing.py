"""
This module provides utility functions for handling
numpy indexing options.
"""

import numpy
import collections
try:
    from types import EllipsisType
except:
    EllipsisType = type(Ellipsis)
    
if hasattr(numpy, 'count_nonzero'):
    def count_nonzero(array):
        return numpy.count_nonzero(array)
else:
    def count_nonzero(array):
        return len(numpy.flatnonzero(array))

ceil=lambda x,y: (x//y+(x%y>0))
        
# unpack_slice: get start,stop step infsofar possible w/o length
def unpack_slice(s):
    step = 1 if s.step is None else s.step
    if step==0: raise ValueError("slice step can't be zero")
    if step>0:
        start = 0 if s.start is None else s.start 
        stop = s.stop
    else:
        start = -1 if s.start is None else s.start
        stop = s.stop
    return start, stop, step

# resolve_slice: determine consistent start,stop, step given length
# note the difference with slice().indices(length)
def resolve_slice(s,length):
    step = 1 if s.step is None else s.step
    if step==0: raise ValueError("slice step can't be zero")
    if step>0:
        start = 0 if s.start is None else s.start
        stop  = length if s.stop is None else s.stop
        if start<0: start=length+start
        if start<0: start=0
        if stop<0: stop=length+stop
        if stop>length: stop=length
        if stop<start: stop=start
    else:
        start = -1 if s.start is None else s.start
        stop  = -(length+1) if s.stop is None else s.stop
        if start>=0: start=-length+start
        if start>-1: start=-1  
        if stop>=0: stop=-length+stop
        if stop<-length: stop=-(length+1)
        if stop>start: stop=start
    return start,stop,step

# combine slices: s1, s2 must be resolved slices!! 
# note that it can return a 'unresolved slice'
def combine_slices(s1,s2):
    a1,b1,c1=unpack_slice(s1)
    a2,b2,c2=unpack_slice(s2)
  
    if b1 is None or b2 is None:
        raise Exception("combining slices not possible")
  
    c3=c1*c2
    imax= ceil( abs(b1-a1), abs(c1))
    if c2<0:
        a2=imax+a2
        b2=imax+b2
    
    a3=a1+a2*c1
      
    jmax=ceil( abs(b2-a2), abs(c2))
    b3=jmax*c3+a3
    if a3<0:
        if b3>-1: b3=None
    else:
        if b3<0: b3=None
    return a3,b3,c3
  
def combine_indices(index0, index1):
    if index1 is None or index0 is None:
        raise Exception("unhandled case, think about numpy")
    
    if isinstance(index0, tuple):
        if len(index0) == 1:
            index0 = index0[0]
        elif isinstance(index1, tuple):
            result=[]
            offset=0
            array_offset=None
            for i0 in index0:
                if isinstance(i0, (int,numpy.integer)):
                  result.append(i0)
                elif isinstance(i0, numpy.ndarray) and i0.dtype != "bool":
                  if array_offset is None:
                    result.append(combine_indices(i0,index1[offset]))
                    array_offset=offset
                    offset+=1
                  else:
                    result.append(combine_indices(i0,index1[array_offset]))                    
                else:
                  result.append(combine_indices(i0,index1[offset]))
                  offset+=1
            if offset<len(index1):
              result.extend(index1[offset:])
            return tuple(result)              
        else:
            index = 0
            for i, x in enumerate(index0):
                index = i
                if isinstance(x, (int,numpy.integer)):
                    continue
                elif isinstance(x, slice):
                    break
                elif isinstance(x, EllipsisType):
                    break
                else:
                    break
            result = []
            result.extend(index0[:index])
            continuation = combine_indices(index0[index], index1)
            if isinstance(continuation, collections.abc.Sequence):
                result.extend(continuation)
            else: 
                result.append(continuation)
            result.extend(index0[index+1:])
            return tuple(result)
        
    if isinstance(index0, (int,numpy.integer)):
        if isinstance(index1, tuple):
            return (index0,)+index1
        else:
            return (index0, index1)
    elif isinstance(index0, slice):
        if isinstance(index1, (int, numpy.integer)):
            start,stop,step = unpack_slice(index0)
            if index1>=0:
              return start + (index1 * step)
            else:
              imax= ceil( abs(stop-start), abs(step))
              stop=start+imax*step
              return stop + index1*step
        elif isinstance(index1, EllipsisType):
            return index0
        elif isinstance(index1, numpy.ndarray):
            start,stop,step = unpack_slice(index0)
            imax= ceil( abs(stop-start), abs(step))
            stop=start+imax*step
            return start+ (index1 *step)*(index1>=0)+(stop + index1*step)*(index1<0)
        else:
            if isinstance(index1, slice):
              start,stop,step = combine_slices(index0, index1)
              return numpy.s_[start:stop:step]
            else:
              return (combine_indices(index0, index1[0]),)+index1[1:]

    elif isinstance(index0, EllipsisType):
        if isinstance(index1, slice):
            return index1
        elif isinstance(index1, EllipsisType):
            return index0
        elif isinstance(index1, (int, numpy.integer)):
            return index1
        else:
            raise Exception("not handled yet")
    elif isinstance(index0, list) or isinstance(index0, numpy.ndarray):
        ndarray = numpy.asarray(index0)
        if ndarray.dtype == 'bool':
            ndarray1 = numpy.zeros_like(ndarray)
            ndarray2 = ndarray1[ndarray]
            ndarray2[index1] = True
            ndarray1[ndarray] = ndarray2
            return ndarray1
        else:
            return index0[index1]
                
    else:
        raise Exception("index must be a integer, slice or sequence")

def is_all_int(sequence):
    for x in sequence:
        if not (isinstance(x, (int, numpy.integer))):
            return False
    return True
    
def number_of_dimensions(array, index):
    return number_of_dimensions_after_index(array.ndim, index)

def number_of_dimensions_after_index(number_of_dimensions, index):
    if isinstance(index, EllipsisType):
        return number_of_dimensions
    elif isinstance(index, tuple):
        if is_all_int(index):
            return number_of_dimensions - len(index)
        else:
            result = number_of_dimensions
            arrays=[]
            for x in index:
                if isinstance(x, EllipsisType):
                    pass
                elif isinstance(x, slice):
                    pass 
                elif isinstance(x, numpy.ndarray):
                    arrays.append(x)
                    result-=1
                else:
                    result -= 1
            if arrays:
                b=numpy.broadcast(*arrays)
                result+=b.nd
            return result
                    
    elif isinstance(index, (int,numpy.integer)):
        return number_of_dimensions - 1
    elif isinstance(index, slice):
        return number_of_dimensions
    elif isinstance(index, list) or isinstance(index, numpy.ndarray):
        ndarray = numpy.asarray(index)
        if ndarray.dtype == 'bool':
            return number_of_dimensions - len(ndarray.shape) + 1
        else:
            if isinstance(index, list):
                raise Exception("indexing with lists is inconsistent with indexing numpy arrays, hence not permitted atm")
            return number_of_dimensions + len(ndarray.shape) - 1
    else:
        raise Exception("Not handled yet")
    
def normalize_slices(shape,index):
    """ returns index with slice expressions normalized, i.e.
    simplified using actual length, replace ellipsis, extend where necessary
    """
    if isinstance(index,EllipsisType):
        index=(Ellipsis,)
    if isinstance(index, tuple):
        if is_all_int(index):
            return tuple([ind if ind>=0 else shape[i]+ind for i,ind in enumerate(index)])
        else:
            result = []
            first_ellipsis = True
            arrays = [x for x in index if isinstance(x,numpy.ndarray)]
            for length,x in zip(shape,index):
                if isinstance(x, slice):
                    result.append( slice(*resolve_slice(x,length)) )
                elif isinstance(x,EllipsisType):
                    if first_ellipsis:
                      n=len(shape)-len(index)+1
                      result.extend([slice(0,shape[i+len(result)],1) for i in range(n)])
                      first_ellipsis=False
                    else:
                      result.append(slice(*resolve_slice(slice(None),length)))
                else:
                    result.append(x)
            n=len(shape)-len(result)
            result.extend([slice(0,shape[len(shape)-n+i],1) for i in range(n)])
            return tuple(result)
                    
    if isinstance(index, slice):
        if isinstance(shape,(int,numpy.integer)):
            return slice(*resolve_slice(index,shape))
        else:   
            return normalize_slices(shape,(index,))
    else:
        return index
    
    
def shape_after_index(shape, index):
    index=normalize_slices(shape,index)
    if isinstance(index, tuple):
        if is_all_int(index):
            return tuple(shape[len(index):])
        else:
            if len(index) != len(shape):
                raise Exception("should not be possible")

            shape_as_list = list(shape)
            result = []
            arrays = [x for x in index if isinstance(x,numpy.ndarray)]
            for i,x in enumerate(index):
                if isinstance(x, slice):
                    start,stop,step = resolve_slice(x, shape_as_list[i])
                    if step>0:
                      nitems = (stop - 1 - start) // step + 1
                    else:
                      nitems = (start- stop -1) // (-step) + 1
                    result.append(nitems)
                elif isinstance(x, numpy.ndarray):
                    if arrays:
                        b=numpy.broadcast(*arrays)
                        result+=b.shape
                        arrays=None
                else:
                    pass
            return tuple(result)
    elif isinstance(index, (int,numpy.integer)):
        return tuple(shape[1:])
    elif isinstance(index, slice):
        return shape_after_index(shape,(index,))
    elif isinstance(index, list) or isinstance(index, numpy.ndarray):
        ndarray = numpy.asarray(index)
        if ndarray.dtype == 'bool':
            if ndarray.shape == shape:
                return (count_nonzero(ndarray),)
            if len(ndarray.shape) < len(shape):
                if not ndarray.shape == shape[:len(ndarray.shape)]:
                    raise Exception("Shape is not compatible")
                    
                result = list(shape[len(ndarray.shape):])
                result.insert(0, count_nonzero(ndarray))
                return tuple(result)
            else:
                
                raise Exception("Not handled yet")

        else:
            if isinstance(index, list):
                raise Exception("indexing with lists is inconsistent with indexing numpy array, hence not permitted atm")

            return ndarray.shape+shape[1:]
            #~ return numpy.zeros(shape)[ndarray].shape # this is cheating a bit..

    else:
        raise Exception("Not handled yet")
        

def split_numpy_index_over_dimensions(index, dimension_values):
    """given a numpy index and a list of dimension values (the accepted values
    per dimension, return the selected values per dimension, values are always 
    arrays"""
    result = list(dimension_values)
    if isinstance(index, (int,numpy.integer)):
        result[0] = result[0][index]
        return result
    elif isinstance(index, slice):
        result[0] = result[0][index]
        return result
    elif isinstance(index, tuple):
        if is_all_int(index):
            for i, x in enumerate(index):
                result[i] = result[i][x]
            return result
        else:
            number_of_indices = len(index)
            i = 0
            for x in index:
                if isinstance(x, (int,numpy.integer)):
                    result[i] = result[i][x]
                elif x is Ellipsis:
                    result[i] = result[i]
                    for _ in range(len(dimension_values) - number_of_indices):
                        i += 1
                        result[i] = result[i]
                        number_of_indices += 1
                elif isinstance(x, slice):
                    result[i] = result[i][x]
                i += 1
            return result
    else:
        raise Exception("Not handled yet")
    
