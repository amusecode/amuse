from amuse.test import amusetest

from amuse.support.interface import InCodeComponentImplementation

from amuse.datamodel.indexing import *
from amuse.datamodel import indexing

class TestIndexing(amusetest.TestCase):
    def test1(self):
        self.assertEqual(2, number_of_dimensions_after_index(3, 1))
        self.assertEqual(3, number_of_dimensions_after_index(3, numpy.s_[0:3]))
        self.assertEqual(1, number_of_dimensions_after_index(3, combine_indices(3,2)))
        self.assertEqual(0, number_of_dimensions_after_index(3, combine_indices(combine_indices(3,2),1)))
        self.assertEqual(3, indexing.number_of_dimensions_after_index(3, numpy.s_[1:2,...,...])) 
        self.assertEqual(3, indexing.number_of_dimensions_after_index(3, numpy.s_[1:2,:,:])) 
        
    def test2(self):
        a = numpy.arange(12).reshape(3,4)
        self.assertEqual(a[combine_indices(0,1)], a[0][1])
        self.assertEqual(a[combine_indices(1,0)], a[1][0])
        self.assertTrue(numpy.all(a[combine_indices(1,numpy.s_[0:2])] == a[1][0:2]))
        indirect = combine_indices(0,1)
        self.assertEqual(number_of_dimensions(a, indirect), 0)
        
        
    def test3(self):
        a = numpy.arange(12).reshape(3,4)
        self.assertTrue(a[combine_indices(numpy.s_[0:2],0)].shape, a[0:2][0].shape)
        self.assertTrue(numpy.all(a[combine_indices(numpy.s_[0:2],0)] ==  a[0:2][0]))

    def test4(self):
        a = numpy.arange(12).reshape(3,4)
        direct =  a[1][:]
        indirect = a[combine_indices(1, indexing.normalize_slices(a[1].shape,numpy.s_[:]))]
        self.assertEqual(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))
        
    def test5(self):
        a = numpy.arange(12).reshape(3,4)
        direct =  a[0:2][:]
        indirect = a[combine_indices(numpy.s_[0:2],indexing.normalize_slices(a[0:2].shape,numpy.s_[:]))]
        self.assertEqual(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))

    def test6(self):
        a = numpy.arange(12).reshape(3,4)
        direct =  a[1:3][1:]
        indirect = a[combine_indices(numpy.s_[1:3],indexing.normalize_slices(a[1:3].shape,numpy.s_[1:]))]
        self.assertEqual(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))

    def test7(self):
        a = numpy.arange(30).reshape(6,5)
        direct =  a[1:5:2][1:]
        indirect = a[combine_indices(numpy.s_[1:5:2],indexing.normalize_slices(a[1:5:2].shape,numpy.s_[1:]))]
        self.assertEqual(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))
        
    def test8(self):
        a = numpy.arange(30)
        direct =  a[2:14:3][1:5:2]
        indirect = a[combine_indices(numpy.s_[2:14:3],numpy.s_[1:5:2])]
        self.assertEqual(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))
        
    
    def test9(self):
        a = numpy.arange(100)
        for s in range(0,40):
            for e in range(40,101):
                for step in range(1,5):
                    direct =  a[s:e:step][1:5:2]
                    indirect = a[combine_indices(numpy.s_[s:e:step], 
                      indexing.normalize_slices(a[s:e:step].shape,numpy.s_[1:5:2]))]
                    self.assertEqual(indirect.shape, direct.shape)
                    self.assertTrue(numpy.all(indirect ==  direct))
    
    def test10(self):
        a = numpy.arange(60).reshape(5,6,2)
        direct =  a[3][2][1]
        indirect = a[combine_indices(combine_indices(3,2),1)]
        self.assertEqual(indirect, direct)
        
    def test11(self):
        a = numpy.arange(60).reshape(5,6,2)
        direct =  a[3]
        indirect = a[combine_indices(3,Ellipsis)]
        self.assertEqual(indirect.shape, direct.shape)
        self.assertTrue(numpy.all(indirect ==  direct))
    
    def test12(self):
        self.assertEqual((1,4,2), indexing.shape_after_index((5,4,2), numpy.s_[1:2,...,...])) 
        self.assertEqual((1,4,2), indexing.shape_after_index((5,4,2), numpy.s_[1:2,:,:])) 
        self.assertEqual((2,4,2), indexing.shape_after_index((5,4,2), numpy.s_[1:3,...,...])) 
        self.assertEqual((2,4,2), indexing.shape_after_index((5,4,2), numpy.s_[1:3,:,:])) 
        self.assertEqual((2,1,2), indexing.shape_after_index((5,4,2), numpy.s_[1:3,2:3,...])) 
        self.assertEqual((2,1,2), indexing.shape_after_index((5,4,2), numpy.s_[1:3,2:3,:])) 
        
    
    def xtest13(self):
        combined_indices = combine_indices(numpy.s_[1:3],numpy.s_[:])
        self.assertEqual(combined_indices, numpy.s_[1:3:1])
        combined_indices = combine_indices(numpy.s_[:], numpy.s_[1:3])
        self.assertEqual(combined_indices, numpy.s_[1:3:1])
        combined_indices = combine_indices((numpy.s_[:], numpy.s_[:]), (numpy.s_[1:3], numpy.s_[1:2]))
        self.assertEqual(combined_indices,(numpy.s_[1:3:1], numpy.s_[1:2:1]))
        
        combined_indices = combine_indices((numpy.s_[0:2], numpy.s_[:]), (numpy.s_[1:3], numpy.s_[1:2]))
        self.assertEqual(combined_indices,(numpy.s_[1:2:1], numpy.s_[1:2:1]))
        
        
    def test14(self):
        self.assertEqual((5,4,2), indexing.shape_after_index((5,4,2), numpy.s_[:10,...,...]))
        self.assertEqual((5,4,2), indexing.shape_after_index((5,4,2), numpy.s_[...,:10,...]))
        self.assertEqual((4,4,2), indexing.shape_after_index((5,4,2), numpy.s_[1:10,...,...]))
        self.assertEqual((1,4,2), indexing.shape_after_index((5,4,2), numpy.s_[-1:,...,...]))
        self.assertEqual((2,4,2), indexing.shape_after_index((5,4,2), numpy.s_[-2:,...,...]))
        self.assertEqual((1,4,2), indexing.shape_after_index((5,4,2), numpy.s_[-2:-1,...,...]))
        self.assertEqual((5,4,2), indexing.shape_after_index((5,4,2), numpy.s_[-10:,...,...]))
        

    
    def test15(self):
        a = numpy.arange(6).reshape(2,3)
        indices = numpy.asarray([[True, False, True],[True,False,True]])
        direct =  a[indices][list([1,3])]
        combined = combine_indices(indices,[1,3])
        indirect = a[combined]
        self.assertEqual(indirect, direct)
        
    
    def test16(self):
        self.assertEqual((4,), indexing.shape_after_index((2,3), [[True, False, True],[True,False,True]]))
        self.assertEqual((1,3), indexing.shape_after_index((2,3), [True, False]))
    
    def test17(self):
        a = numpy.arange(6).reshape(2,3)
        indices = numpy.asarray([True,False])
        direct =  a[indices, 1:][0,1:]
        combined = combine_indices(indexing.normalize_slices(a.shape,numpy.s_[indices,1:]),
            indexing.normalize_slices(a[indices,1:].shape,numpy.s_[0,1:]))
        indirect = a[combined]
        self.assertEqual(indirect, direct)
        
    def test18(self):
        a = numpy.arange(6).reshape(2,3)
        indices = numpy.asarray([True,False])
        direct =  a[indices, 1:][0]
        combined = combine_indices(numpy.s_[indices,1:],0)
        indirect = a[combined]
        self.assertEqual(indirect, direct)
        
    def test19(self):
        a = numpy.arange(6).reshape(2,3)
        indices = numpy.asarray([True,False])
        direct =  a[:1, 1:][0]
        combined = combine_indices(numpy.s_[:1,1:],0)
        indirect = a[combined]
        self.assertEqual(indirect, direct)

    def test20(self):
        combined=combine_indices( slice(1, 199, None), slice(3, 5, None) )
        self.assertEqual(combined,slice(4,6,1))
        
    def test21(self):
        shape=shape_after_index((200,), slice(4, 6, 1))
        self.assertEqual(shape,(2,))

    def xtest22(self):
        tiny=list(range(2))
        small=list(range(10))
        big=list(range(1000))
        
        # combining slicings w negative stops not possible! e.g. ((7,-1),(2,3),(9,10,1))
        # (without normalize)
        slicings=[ ((9,19),(5,9,2),(14,18,2)),
                   ((7,19,2),(5,9,1),(17,19,2)),
                   ((1,None),(1,10),(2,11,1)),
                   ((7,None),(1,10),(8,17,1)),
                   ((None,12),(3,5),(3,5,1)),
                   ((None,12),(3,15),(3,12,1)),
                   ((None,None),(3,5),(3,5,1)),
                   ((None,None),(3,15),(3,15,1)),
                   ((None,None),(None,15),(0,15,1)),
                   ((None,None),(None,None),(0,None,1)),
                   ((9,None),(None,None),(9,None,1)),
                   ((9,None),(6,None),(15,None,1)),
                   ((9,None),(None,40),(9,49,1)),
                   ((1,None),(None,40),(1,41,1)),
                   ((9,16),(None,40),(9,16,1)),
                   ((1,16),(None,40),(1,16,1)),
                   ((49,16),(None,40),(16,16,1)),
                   ((41,66),(None,40),(41,66,1)),
                 ]
        
        for t1,t2,t3 in slicings:
          s1=slice(*t1)
          s2=slice(*t2)
          s3=slice(*t3)
          self.assertEqual(combine_slices(s1,s2),t3)
          self.assertTrue(tiny[s1][s2]==tiny[s3])
          self.assertTrue(small[s1][s2]==small[s3])
          self.assertTrue(big[s1][s2]==big[s3])
        
    def xtest23(self):
        import random
        random.seed(123456)
        tiny=list(range(2))
        small=list(range(20))
        big=list(range(2000))
        
        Ntest=1000
        start0=[random.randint(0,20) for x in range(Ntest)]        
        stop0=[random.randint(15,50) for x in range(Ntest)]        
        step0=[random.randint(1,3) for x in range(Ntest)]        
        start1=[random.randint(0,10) for x in range(Ntest)]        
        stop1=[random.randint(5,25) for x in range(Ntest)]        
        step1=[random.randint(1,3) for x in range(Ntest)]        
        
        slicings=[]
        for x in zip(start0,stop0,step0,start1,stop1,step1):
          slicings.append(((x[0],x[1],x[2]),(x[3],x[4],x[5])))
        
        for t1,t2 in slicings:
          s1=slice(*t1)
          s2=slice(*t2)
          t3=combine_slices(s1,s2)
          s3=slice(*t3)
          self.assertTrue(tiny[s1][s2]==tiny[s3])
          self.assertTrue(small[s1][s2]==small[s3])
          self.assertTrue(big[s1][s2]==big[s3])

    def test24(self):
        import random
        random.seed(123456)
        tiny=list(range(2))
        small=list(range(20))
        big=list(range(2000))
        
        Ntest=1000
        stop0=[random.randint(0,20) for x in range(Ntest)]        
        start0=[random.randint(15,50) for x in range(Ntest)]        
        step0=[random.randint(-3,-1) for x in range(Ntest)]        
        start1=[random.randint(0,10) for x in range(Ntest)]        
        stop1=[random.randint(5,25) for x in range(Ntest)]        
        step1=[random.randint(1,3) for x in range(Ntest)]        
        
        slicings=[]
        for x in zip(start0,stop0,step0,start1,stop1,step1):
          slicings.append(((x[0],x[1],x[2]),(x[3],x[4],x[5])))
        
        for t1,t2 in slicings:
          s1=slice(*t1)
          s2=slice(*t2)

          t3=combine_slices(normalize_slices(len(tiny),s1),normalize_slices(len(tiny[s1]),s2))
          s3=slice(*t3)
          self.assertTrue(tiny[s1][s2]==tiny[s3])

          t3=combine_slices(normalize_slices(len(small),s1),normalize_slices(len(small[s1]),s2))
          s3=slice(*t3)
          self.assertTrue(small[s1][s2]==small[s3])
          t3=combine_slices(normalize_slices(len(big),s1),normalize_slices(len(big[s1]),s2))
          s3=slice(*t3)
          self.assertTrue(big[s1][s2]==big[s3])

    def test25(self):
        import random
        random.seed(123456)
        tiny=list(range(2))
        small=list(range(20))
        big=list(range(2000))
        
        Ntest=1000
        stop0=[random.randint(0,20) for x in range(Ntest)]        
        start0=[random.randint(15,50) for x in range(Ntest)]        
        step0=[random.randint(-3,-1) for x in range(Ntest)]        
        stop1=[random.randint(0,10) for x in range(Ntest)]        
        start1=[random.randint(5,25) for x in range(Ntest)]        
        step1=[random.randint(-3,-1) for x in range(Ntest)]        
        
        slicings=[]
        for x in zip(start0,stop0,step0,start1,stop1,step1):
            slicings.append(((x[0],x[1],x[2]),(x[3],x[4],x[5])))
        
        for t1,t2 in slicings:
            s1=slice(*t1)
            s2=slice(*t2)

            t3=combine_slices(normalize_slices(len(tiny),s1),normalize_slices(len(tiny[s1]),s2))
            s3=slice(*t3)
            #~ print(s1,s2,s3)
            self.assertTrue(tiny[s1][s2]==tiny[s3])

            t3=combine_slices(normalize_slices(len(small),s1),normalize_slices(len(small[s1]),s2))
            s3=slice(*t3)
            self.assertTrue(small[s1][s2]==small[s3])
            t3=combine_slices(normalize_slices(len(big),s1),normalize_slices(len(big[s1]),s2))
            s3=slice(*t3)
            self.assertTrue(big[s1][s2]==big[s3])

    def test26(self):
        oned=numpy.zeros(5)
        threed=numpy.zeros((4,5,6))
        for index in [0,[1],[1,2],[[1,2],[2,3]],[[2]],[[0,1]]]:
          i=numpy.array(index)
          self.assertEqual(len(oned[i].shape), number_of_dimensions_after_index(1, i ))
        for index in [0,[1],[1,2],[[1,2],[2,3]],[[2]],[[2,1]]]:
          i=numpy.array(index)
          self.assertEqual(len(threed[i].shape), number_of_dimensions_after_index(3, i ))

    def test27(self):
        oned=numpy.zeros(5)
        threed=numpy.zeros((4,5,6))
        for index in [0,[1],[1,2],[[1,2],[2,3]],[[2]],[[0,1]]]:
          i=numpy.array(index)
          self.assertEqual(oned[i].shape, shape_after_index(oned.shape, i ))
        for index in [0,[1],[1,2],[[1,2],[2,3]],[[2]],[[2,1]],[[[[0],[1],[1]]]]]:
          i=numpy.array(index)
          self.assertEqual(threed[i].shape, shape_after_index(threed.shape, i ))

    def test28(self):
        twod=numpy.zeros((5,6))
        threed=numpy.zeros((4,5,6))
        for _i,_j in [([0],[1]),([0,2],[1,3]),([0,2],[1,3]),([[0,1],[1,2]],[[2,3],[3,4]])]:
          i=numpy.array(_i)
          j=numpy.array(_j)
          self.assertEqual(len(twod[i,j].shape), number_of_dimensions_after_index(2, (i,j) ))
        for _i,_j in [([0],[1]),([0,2],[1,3]),([0,2],[1,3]),([[0,1],[1,2]],[[2,3],[3,4]])]:
          i=numpy.array(_i)
          j=numpy.array(_j)
          self.assertEqual(len(threed[i,j].shape), number_of_dimensions_after_index(3, (i,j) ))

    def test29(self):
        twod=numpy.zeros((5,6))
        for _i,_j in [([0],[1]),([0,2],[1,3]),([0,2],[1,3]),([[0,1],[1,2]],[[2,3],[3,4]])]:
          i=numpy.array(_i)
          j=numpy.array(_j)
          self.assertEqual(twod[i,j].shape, shape_after_index(twod.shape, (i,j) ))

        threed=numpy.zeros((4,5,6))
        for _i,_j in [([0],[1]),([0,2],[1,3]),([0,2],[1,3]),([[0,1],[1,2]],[[2,3],[3,4]])]:
          i=numpy.array(_i)
          j=numpy.array(_j)
          self.assertEqual(threed[i,j].shape, shape_after_index(threed.shape, (i,j) ))

        fourd=numpy.zeros((4,5,6,7))
        for _i,_j in [([0],[1]),([0,2],[1,3])]:
          i=numpy.array(_i)
          j=numpy.array(_j)
          self.assertEqual(fourd[i,j].shape, shape_after_index(fourd.shape, (i,j) ))
        for _i,_j in [([0,2],[1,3])]:
          i=numpy.array(_i)
          j=numpy.array(_j)
          self.assertEqual(fourd[i,Ellipsis,j].shape, shape_after_index(fourd.shape, (i,Ellipsis,j) ))
        for _i,_j in [([0,2],[1,3])]:
          i=numpy.array(_i)
          j=numpy.array(_j)
          self.assertEqual(fourd[i,slice(None),j].shape, shape_after_index(fourd.shape, (i,slice(None),j) ))


        
class TestSplitOverDimensions(amusetest.TestCase):
    
    def test1(self):
        dimension_values = [
            [3,4,5,6],
            ['a','b','c']
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions(0, dimension_values)
        self.assertEqual(len(split_dimension_values) ,2)
        self.assertEqual(split_dimension_values[0], 3)
        self.assertEqual(split_dimension_values[1], ['a','b','c'])
       
    def test2(self):
        dimension_values = [
            [3,4,5,6],
            ['a','b','c']
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions((1,2), dimension_values)
        self.assertEqual(len(split_dimension_values) ,2)
        self.assertEqual(split_dimension_values[0], 4)
        self.assertEqual(split_dimension_values[1], 'c')
    
    def test3(self):
        dimension_values = [
            [3,4,5,6],
            ['a','b','c']
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions(slice(0,2), dimension_values)
        self.assertEqual(len(split_dimension_values) ,2)
        self.assertEqual(split_dimension_values[0], [3,4])
        self.assertEqual(split_dimension_values[1], ['a','b','c'])
        
    def test4(self):
        dimension_values = [
            [0, 1, 2, 3, 4, 5, 6, 7, 8, ]
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions(slice(1,7,2), dimension_values)
        self.assertEqual(len(split_dimension_values) ,1)
        self.assertEqual(split_dimension_values[0], [1, 3, 5])
        
    
    def test5(self):
        dimension_values = [
            [0, 1, 2, 3, 4, 5, 6, 7, 8,9 ]
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions(slice(-2,10), dimension_values)
        self.assertEqual(split_dimension_values[0], [8, 9])
        split_dimension_values = indexing.split_numpy_index_over_dimensions(slice(-3,3,-1), dimension_values)
        self.assertEqual(split_dimension_values[0], [7, 6, 5, 4])
    
    
    def test6(self):
        dimension_values = [
            [0, 1, 2, 3, 4, 5, 6, 7, 8,9 ]
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions(slice(5,None), dimension_values)
        self.assertEqual(split_dimension_values[0], [5, 6, 7, 8, 9])
    
    
    def test7(self):
        dimension_values = [
            [0, 1],
            [0, 1, 2],
            [0]
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions(slice(1,2), dimension_values)
        self.assertEqual(split_dimension_values[0], [1])
        self.assertEqual(split_dimension_values[1], [0,1,2])
        self.assertEqual(split_dimension_values[2], [0])
    
    
    def test8(self):
        dimension_values = [
            [0, 1],
            [0, 1, 2],
            [0]
        ]
        split_dimension_values = indexing.split_numpy_index_over_dimensions((Ellipsis,0), dimension_values)
        self.assertEqual(split_dimension_values[0], [0,1])
        self.assertEqual(split_dimension_values[1], [0,1,2])
        self.assertEqual(split_dimension_values[2], 0)
        
        split_dimension_values = indexing.split_numpy_index_over_dimensions((slice(None),slice(None),0), dimension_values)
        self.assertEqual(split_dimension_values[0], [0,1])
        self.assertEqual(split_dimension_values[1], [0,1,2])
        self.assertEqual(split_dimension_values[2], 0)
    
        split_dimension_values = indexing.split_numpy_index_over_dimensions((Ellipsis,0, Ellipsis), dimension_values)
        self.assertEqual(split_dimension_values[0], [0,1])
        self.assertEqual(split_dimension_values[1], 0)
        self.assertEqual(split_dimension_values[2], [0])
