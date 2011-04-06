import unittest
from amuse.support.codes import lit
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.community.bhtree.interface import BHTree
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.sse.interface import SSE
from amuse.community.bse.interface import BSE
from amuse.community.mesa.interface import MESA
from amuse.community.evtwin.interface import EVtwin

class TestLiteratureRefs(unittest.TestCase):
    def setUp(self):
        lit.LiteratureReferencesMixIn.literature_list = []

    def test1(self):
        class ClassLitrefs(lit.LiteratureReferencesMixIn):
            """ some comments with added lit refs, i.e.  [#]_ and [3]_ etc...

                .. [#] Gratia, Exempli, *Journal of I.M.*, **1**, 1--100 (2009).
                .. [3] Dude, John, *The Intern. Foo Journal of Bars*, **51**, 1647--1751 (2009)  
                .. [4] Hat, John, *The Ex. Foo Journal of Bars*, **51**, 1647--1751 (2009)  
            """
            def __init__(self):
                lit.LiteratureReferencesMixIn.__init__(self)

        classnames = lit.LiteratureReferencesMixIn.names_of_classes_with_references()
        self.assertFalse("ClassLitrefs" in classnames)

        instance = ClassLitrefs()
        
        classnames = lit.LiteratureReferencesMixIn.names_of_classes_with_references()
        self.assertTrue("ClassLitrefs" in classnames)

    def test2(self):
        class ClassLitrefs(lit.LiteratureReferencesMixIn):
            """ some comments with added lit refs, i.e.  [#]_ and [3]_ etc...

                .. [#] Gratia, Exempli, *Journal of I.M.*, **1**, 1--100 (2009).
                .. [3] Dude, John, *The Intern. Foo Journal of Bars*, **51**, 1647--1751 (2009)  
                .. [4] Hat, John, *The Gal. Foo Journal of Bars*, **3**, 16--51 (2009)  
            """
            def __init__(self):
                lit.LiteratureReferencesMixIn.__init__(self)
        
        lit.TrackLiteratureReferences.default().registered_classes = set([])
        string = lit.LiteratureReferencesMixIn.all_literature_references_string()
        self.assertTrue("AMUSE" in string)
        self.assertTrue("multiphysics and multiscale software environment" in string)

        instance = ClassLitrefs()
        
        string = lit.LiteratureReferencesMixIn.all_literature_references_string()
        lit.LiteratureReferencesMixIn.print_literature_references()
        self.assertFalse(".. [#] Gratia, Exem" in string)
        self.assertFalse(".. [3] Dude" in string)
        self.assertFalse(".. [4] Hat" in string)

        self.assertTrue("Gratia, Exem" in string)
        self.assertTrue("Dude" in string)
        self.assertTrue("Hat" in string)

        self.assertTrue("ClassLitrefs" in string)
        
    def test3(self):
        class ClassLitrefs(lit.LiteratureReferencesMixIn):
            """ some comments with added lit refs, i.e.  [#]_ and [3]_ etc...

                .. [#] Gratia, Exempli, *Journal of I.M.*, **1**, 1--100 (2009).
                .. [3] Dude, John, *The Intern. Foo Journal of Bars*, **51**, 1647--1751 (2009)  
                .. [4] Hat, John, *The Gal. Foo Journal of Bars*, **3**, 16--51 (2009)  
            """
            def __init__(self):
                lit.LiteratureReferencesMixIn.__init__(self)
        
        string = lit.LiteratureReferencesMixIn.export2html()
        print string
        #from docutils import core
        #print core.publish_string(source = string)
        
    def test4(self):
        print "This test shows how the references to currently used legacy codes can be obtained."
       
        gravity = BHTree()
        gravity2 = Hermite()
        gravity3 = PhiGRAPE()
        stellar_evolution1 = SSE()
        stellar_evolution2 = BSE()
        stellar_evolution3 = EVtwin()
        #stellar_evolution4 = MESA()
        print "Each legacy code instance collects all references."
        print "They can be retrieved with the print_refs method:"
        print ">>> gravity.print_refs()"
        gravity.print_literature_references()
        all_refs_as_returned_by_code1 = gravity.all_literature_references_string()
        all_refs_as_returned_by_code2 = stellar_evolution2.all_literature_references_string()
        print "Checking whether each of them returns the same set of references... ",
        self.assertEquals(all_refs_as_returned_by_code1, all_refs_as_returned_by_code2)
        print "ok."
        self.assertTrue('Barnes, J., Hut, P., A Hierarchical O(N log N) force-calculation algorithm, *Nature*, **4**, 324 (1986)' 
            in all_refs_as_returned_by_code1)
        self.assertTrue('Hurley J.R., Pols O.R., Tout C.A., 2000, MNRAS, 315, 543' 
            in all_refs_as_returned_by_code1)
        self.assertTrue('Eggleton, P.P. 1971, MNRAS, 151, 351: "The evolution of low mass stars"' 
            in all_refs_as_returned_by_code1)
        gravity.stop()
        stellar_evolution1.stop()
        stellar_evolution2.stop()
        stellar_evolution3.stop()
        #stellar_evolution4.stop()
