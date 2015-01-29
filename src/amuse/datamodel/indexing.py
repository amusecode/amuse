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
        
def unpack_slice(slice):
    start = 0 if slice.start is None else slice.start
    stop = None if slice.stop is None else slice.stop
    step = 1 if slice.step is None else slice.step
    return start, stop, step

def combine_slices(slice0, slice1):
    start0,stop0,step0 = unpack_slice(slice0)
    start1,stop1,step1 = unpack_slice(slice1)

    if step0>0:
      newstart = start0 + (start1 * step0) if stop0 is None else min(start0 + (start1 * step0),stop0)
    else:
      newstart = max(start0 + (start1 * step0),stop0)
    if stop0 is None:
        if stop1 is None:
            newstop = None
        else:
            newstop = (stop1 * step0) + start0
    elif stop0<0:
        if stop1 is None:
          newstop = stop0
        elif stop1<0:
          newstop = stop0+stop1*step0
        else:
          newstop = stop1*step0+start0
    elif stop1 is None:
        newstop = stop0
    else:
        if step0>0:
          newstop = start0 + (stop1 * step0) if stop0 is None else min(start0 + (stop1 * step0),stop0)
        else:
          newstop = max((stop1 * step0) + start0, stop0)
    
    newstep = step0 * step1
    return newstart, newstop, newstep

def combine_indices(index0, index1):
    if index1 is None:
        return index0
    
    if isinstance(index0, tuple):
        if len(index0) == 1:
            index0 = index0[0]
        elif isinstance(index1, tuple):
            if len(index0) == len(index1):
                combined = [ combine_indices(p0, p1) for p0, p1 in zip(index0, index1)]
                return tuple(combined)
            else:
                raise Exception("unhandled case, two tuple one with different length")
        else:
            index = 0
            for i, x in enumerate(index0):
                index = i
                if isinstance(x, int) or isinstance(x, long):
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
            if isinstance(continuation, collections.Sequence):
                result.extend(continuation)
            else: 
                result.append(continuation)
            result.extend(index0[index+1:])
            return tuple(result)
        
    if isinstance(index0, int) or isinstance(index0, long):
        return (index0, index1)
    elif isinstance(index0, slice):
        if isinstance(index1, int) or isinstance(index1, long):
            start,stop,step = unpack_slice(index0)
            return start + (index1 * step)
        elif isinstance(index1, EllipsisType):
            return index0
        else:
            if isinstance(index1, slice):
              start,stop,step = combine_slices(index0, index1)
              return numpy.s_[start:stop:step]
            else:
              start,stop,step = combine_slices(index0, index1[0])
              return numpy.s_[start:stop:step]

    elif isinstance(index0, EllipsisType):
        if isinstance(index1, slice):
            return index1
        elif isinstance(index1, EllipsisType):
            return index0
        elif isinstance(index1, int) or isinstance(index1, long):
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
            raise Exception("index must be a integer, slice or sequence")
                
    else:
        raise Exception("index must be a integer, slice or sequence")

def is_all_int(sequence):
    for x in sequence:
        if not (isinstance(x, int) or isinstance(x, long)):
            return False
    return True
    
def number_of_dimensions(array, index):
    return number_of_dimensions_after_index(array.ndim, index)

def number_of_dimensions_after_index(number_of_dimensions, index):
    if isinstance(index, tuple):
        if is_all_int(index):
            return number_of_dimensions - len(index)
        else:
            result = number_of_dimensions
            for x in index:
                if isinstance(x, EllipsisType):
                    pass
                elif isinstance(x, slice):
                    pass 
                else:
                    result -= 1
            return result
                    
    elif isinstance(index, int) or isinstance(index, long):
        return number_of_dimensions - 1
    elif isinstance(index, slice):
        return number_of_dimensions
    elif isinstance(index, list) or isinstance(index, numpy.ndarray):
        ndarray = numpy.asarray(index)
        if ndarray.dtype == 'bool':
            return number_of_dimensions - len(ndarray.shape) + 1
        else:
            raise Exception("Not handled yet")
    else:
        raise Exception("Not handled yet")
    
def normalize_slices(shape,index):
    """ returns index with slice expressions normalized, i.e.
    simplified using actual length, replace ellipsis, extend where necessary
    """
    if isinstance(index, tuple):
        if is_all_int(index):
            return index
        else:
            result = []
            first_ellipsis = True
            for length,x in zip(shape,index):
                if isinstance(x, slice):
                    result.append( slice(*x.indices(length)))
                elif isinstance(x,EllipsisType):
                    if first_ellipsis:
                      n=len(shape)-len(index)+1
                      result.extend([slice(0,shape[i],1) for i in range(n)])
                      first_ellipsis=False
                    else:
                      result.append(slice(None))
                else:
                    result.append(x)
            n=len(shape)-len(result)
            result.extend([slice(0,shape[i],1) for i in range(n)])
            return tuple(result)
                    
    if isinstance(index, slice):
        if isinstance(shape,int) or isinstance(shape,long):
            return slice(*index.indices(shape))
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
            for i,x in enumerate(index):
                if isinstance(x, slice):
                    
                    start,stop,step = unpack_slice(x)
                    if start is None:
                        start = 0
                    if stop is None:
                        stop = shape_as_list[i]
                    if stop < 0:
                        stop = shape_as_list[i] + stop
                    if start < 0:
                        start = shape_as_list[i] + start
                    
                    if start < 0:
                        start = 0
                    if stop < 0:
                        stop = 0
                        
                    nmax =  min(stop, shape_as_list[i])
                    nitems = (nmax - start) // step
                    result.append(nitems)
                else:
                    pass
            return tuple(result)
    elif isinstance(index, int) or isinstance(index, long):
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
            raise Exception("Not handled yet")
    else:
        raise Exception("Not handled yet")
        

def split_numpy_index_over_dimensions(index, dimension_values):
    """given a numpy index and a list of dimension values (the accepted values
    per dimension, return the selected values per dimension, values are always 
    arrays"""
    result = list(dimension_values)
    if isinstance(index, int) or isinstance(index, long):
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
                if isinstance(x, int) or isinstance(x, long):
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
    
