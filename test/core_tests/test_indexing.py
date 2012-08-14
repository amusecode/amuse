from amuse.test import amusetest

from amuse.support.interface import InCodeComponentImplementation

from amuse.datamodel.indexing import *
from amuse.datamodel import indexing

class TestIndexing(amusetest.TestCase):
    def test1(self):
        self.assertEquals(2, number_of_dimensions_after_index(3, 1))
        self.assertEquals(3, number_of_dimensions_after_index(3, numpy.s_[0:3]))
        self.assertEquals(1, number_of_dimensions_after_index(3, combine_indices(3,2)))
        self.assertEquals(0, number_of_dimensions_after_index(3, combine_indices(combine_indices(3,2),1)))
        self.assertEquals(3, indexing.number_of_dimensions_after_index(3, numpy.s_[1:2,...,...])) 
        self.assertEquals(3, indexing.number_of_dimensions_after_index(3, numpy.s_[1:2,:,:])) 
        
    def test2(self):
        a = numpy.arange(12).reshape(3,4)
        self.assertEquals(a[combine_indices(0,1)], a[0][1])
        self.assertEquals(a[combine_indices(1,0)], a[1][0])
        self.assertTrue(numpy.all(a[combine_indices(1,numpy.s_[0:2])] == a[1][0:2]))
        indirect = combine_indices(0,1)
        self.assertEquals(number_of_dimensions(a, indirect), 0)
        
        
    def test3(self):
        a = numpy.arange(12).reshape(3,4)
        print a, a[0:2][0], a[combine_indices(numpy.s_[0:2],0)]
        self.assertTrue(a[combine_indices(numpy.s_[0:2],0)].shape, a[0:2][0].shape)
        self.assertTrue(numpy.all(a[combine_indices(numpy.s_[0:2],0)] ==  a[0:2][0]))

    def test4(self):
        a = numpy.arange(12).reshape(3,4)
        direct =  a[1][:]
        indirect = a[combine_indices(1,numpy.s_[:])]
        print a, direct, indirect
        self.assertEquals(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))
        
    def test5(self):
        a = numpy.arange(12).reshape(3,4)
        direct =  a[0:2][:]
        indirect = a[combine_indices(numpy.s_[0:2],numpy.s_[:])]
        print a, direct, indirect
        self.assertEquals(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))

    def test6(self):
        a = numpy.arange(12).reshape(3,4)
        direct =  a[1:3][1:]
        indirect = a[combine_indices(numpy.s_[1:3],numpy.s_[1:])]
        print a, direct, indirect
        self.assertEquals(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))

    def test7(self):
        a = numpy.arange(30).reshape(6,5)
        direct =  a[1:5:2][1:]
        indirect = a[combine_indices(numpy.s_[1:5:2],numpy.s_[1:])]
        self.assertEquals(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))
        
    def test8(self):
        a = numpy.arange(30)
        direct =  a[2:14:3][1:5:2]
        indirect = a[combine_indices(numpy.s_[2:14:3],numpy.s_[1:5:2])]
        self.assertEquals(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))
        
    
    def test9(self):
        a = numpy.arange(100)
        for s in range(0,40):
            for e in range(40,101):
                for step in range(1,5):
                    direct =  a[s:e:step][1:5:2]
                    indirect = a[combine_indices(numpy.s_[s:e:step],numpy.s_[1:5:2])]
                    self.assertEquals(indirect.shape, direct.shape)
                    self.assertTrue(numpy.all(indirect ==  direct))
    
    
    def test10(self):
        a = numpy.arange(60).reshape(5,6,2)
        direct =  a[3][2][1]
        indirect = a[combine_indices(combine_indices(3,2),1)]
        print direct, indirect, combine_indices(combine_indices(3,2),1), a[(3,2,1)]
        self.assertEquals(indirect, direct)
        
    def test11(self):
        a = numpy.arange(60).reshape(5,6,2)
        direct =  a[3]
        indirect = a[combine_indices(3,None)]
        self.assertEquals(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))
    
    def test12(self):
        self.assertEquals((1,4,2), indexing.shape_after_index((5,4,2), numpy.s_[1:2,...,...])) 
        self.assertEquals((1,4,2), indexing.shape_after_index((5,4,2), numpy.s_[1:2,:,:])) 
        self.assertEquals((2,4,2), indexing.shape_after_index((5,4,2), numpy.s_[1:3,...,...])) 
        self.assertEquals((2,4,2), indexing.shape_after_index((5,4,2), numpy.s_[1:3,:,:])) 
        self.assertEquals((2,1,2), indexing.shape_after_index((5,4,2), numpy.s_[1:3,2:3,...])) 
        self.assertEquals((2,1,2), indexing.shape_after_index((5,4,2), numpy.s_[1:3,2:3,:])) 
        
    
    def test13(self):
        combined_indices = combine_indices(numpy.s_[1:3],numpy.s_[:])
        self.assertEquals(combined_indices, numpy.s_[1:3:1])
        combined_indices = combine_indices(numpy.s_[:], numpy.s_[1:3])
        self.assertEquals(combined_indices, numpy.s_[1:3:1])
        combined_indices = combine_indices((numpy.s_[:], numpy.s_[:]), (numpy.s_[1:3], numpy.s_[1:2]))
        self.assertEquals(combined_indices,(numpy.s_[1:3:1], numpy.s_[1:2:1]))
        
        combined_indices = combine_indices((numpy.s_[0:2], numpy.s_[:]), (numpy.s_[1:3], numpy.s_[1:2]))
        self.assertEquals(combined_indices,(numpy.s_[1:2:1], numpy.s_[1:2:1]))
        
        
    def test14(self):
        self.assertEquals((5,4,2), indexing.shape_after_index((5,4,2), numpy.s_[:10,...,...]))
        self.assertEquals((5,4,2), indexing.shape_after_index((5,4,2), numpy.s_[...,:10,...]))
        self.assertEquals((4,4,2), indexing.shape_after_index((5,4,2), numpy.s_[1:10,...,...]))
        self.assertEquals((1,4,2), indexing.shape_after_index((5,4,2), numpy.s_[-1:,...,...]))
        self.assertEquals((2,4,2), indexing.shape_after_index((5,4,2), numpy.s_[-2:,...,...]))
        self.assertEquals((1,4,2), indexing.shape_after_index((5,4,2), numpy.s_[-2:-1,...,...]))
        self.assertEquals((5,4,2), indexing.shape_after_index((5,4,2), numpy.s_[-10:,...,...]))
        

class TestSplitOverDimensions(amusetest.TestCase):
    
    def test1(self):
        dimension_values = [
            [3,4,5,6],
            ['a','b','c']
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions(0, dimension_values)
        self.assertEquals(len(split_dimension_values) ,2)
        self.assertEquals(split_dimension_values[0], 3)
        self.assertEquals(split_dimension_values[1], ['a','b','c'])
       
    def test2(self):
        dimension_values = [
            [3,4,5,6],
            ['a','b','c']
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions((1,2), dimension_values)
        self.assertEquals(len(split_dimension_values) ,2)
        self.assertEquals(split_dimension_values[0], 4)
        self.assertEquals(split_dimension_values[1], 'c')
    
    def test3(self):
        dimension_values = [
            [3,4,5,6],
            ['a','b','c']
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions(slice(0,2), dimension_values)
        self.assertEquals(len(split_dimension_values) ,2)
        self.assertEquals(split_dimension_values[0], [3,4])
        self.assertEquals(split_dimension_values[1], ['a','b','c'])
        
    def test4(self):
        dimension_values = [
            [0, 1, 2, 3, 4, 5, 6, 7, 8, ]
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions(slice(1,7,2), dimension_values)
        self.assertEquals(len(split_dimension_values) ,1)
        self.assertEquals(split_dimension_values[0], [1, 3, 5])
        
    
    def test5(self):
        dimension_values = [
            [0, 1, 2, 3, 4, 5, 6, 7, 8,9 ]
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions(slice(-2,10), dimension_values)
        self.assertEquals(split_dimension_values[0], [8, 9])
        split_dimension_values = indexing.split_numpy_index_over_dimensions(slice(-3,3,-1), dimension_values)
        self.assertEquals(split_dimension_values[0], [7, 6, 5, 4])
    
    
    def test6(self):
        dimension_values = [
            [0, 1, 2, 3, 4, 5, 6, 7, 8,9 ]
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions(slice(5,None), dimension_values)
        self.assertEquals(split_dimension_values[0], [5, 6, 7, 8, 9])
    
    
    def test7(self):
        dimension_values = [
            [0, 1],
            [0, 1, 2],
            [0]
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions(slice(1,2), dimension_values)
        self.assertEquals(split_dimension_values[0], [1])
        self.assertEquals(split_dimension_values[1], [0,1,2])
        self.assertEquals(split_dimension_values[2], [0])
    
    
    def test8(self):
        dimension_values = [
            [0, 1],
            [0, 1, 2],
            [0]
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions((Ellipsis,0), dimension_values)
        self.assertEquals(split_dimension_values[0], [0,1])
        self.assertEquals(split_dimension_values[1], [0,1,2])
        self.assertEquals(split_dimension_values[2], 0)
        
        split_dimension_values = indexing.split_numpy_index_over_dimensions((slice(None),slice(None),0), dimension_values)
        self.assertEquals(split_dimension_values[0], [0,1])
        self.assertEquals(split_dimension_values[1], [0,1,2])
        self.assertEquals(split_dimension_values[2], 0)
    
        split_dimension_values = indexing.split_numpy_index_over_dimensions((Ellipsis,0, Ellipsis), dimension_values)
        self.assertEquals(split_dimension_values[0], [0,1])
        self.assertEquals(split_dimension_values[1], 0)
        self.assertEquals(split_dimension_values[2], [0])