"""
This module provides utility functions for handling
numpy indexing options.
"""

import numpy
import collections
from types import EllipsisType

def unpack_slice(slice):
    start = 0 if slice.start is None else slice.start
    stop = None if slice.stop is None else slice.stop
    step = 1 if slice.step is None else slice.step
    return start, stop, step

def combine_slices(slice0, slice1):
    start0,stop0,step0 = unpack_slice(slice0)
    start1,stop1,step1 = unpack_slice(slice1)
            
    newstart = start0 + (start1 * step0)
    
    if stop0 is None:
        if stop1 is None:
            newstop = None
        else:
            newstop = (stop1 * step0) + start0
    elif stop1 is None:
        newstop = stop0
    else:
        newstop = min((stop1 * step0) + start0, stop0)
    
    newstep = step0 * step1
    return newstart, newstop, newstep

def combine_indices(index0, index1):
    if index1 is None:
        return index0
        
    if isinstance(index0, tuple):
        if len(index0) == 1:
            index0 = index0[0]
        else:
            result = []
            result.extend(index0[:-1])
            continuation = combine_indices(index0[-1], index1)
            if isinstance(continuation, collections.Sequence):
                result.extend(continuation)
            else: 
                result.append(continuation)
            return tuple(result)
    
    if isinstance(index0, int) or isinstance(index0, long):
        return (index0, index1)
    elif isinstance(index0, slice):
        if isinstance(index1, int) or isinstance(index1, long):
            start,stop,step = unpack_slice(index0)
            return start + (index1 * step)
        else:
            start,stop,step = combine_slices(index0, index1)
            return numpy.s_[start:stop:step]
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
    
def shape_after_index(shape, index):
    if isinstance(index, tuple):
        if is_all_int(index):
            return tuple(shape[len(index):])
        else:
            if len(index) > len(shape):
                raise Exception("Not handled yet")
                
            shape_as_list = list(shape)
            result = []
            for i,x in enumerate(index):
                if isinstance(x, EllipsisType):
                    result.append(shape_as_list[i])
                else:
                    pass
            return tuple(result)
    elif isinstance(index, int) or isinstance(index, long):
        return tuple(shape[1:])
    elif isinstance(index, slice):
        result = list(shape)
        result[0] = len(index.indices(shape[0]))
        return tuple(result)
    elif isinstance(index, list) or isinstance(index, numpy.ndarray):
        ndarray = numpy.asarray(index)
        if ndarray.dtype == 'bool':
            raise Exception("Not handled yet")
        else:
            raise Exception("Not handled yet")
    else:
        raise Exception("Not handled yet")
    
    
    
