import unittest
import docutils.core

from amuse.support import literature

class TestLiteratureRefs(unittest.TestCase):
    def setUp(self):
        literature.LiteratureReferencesMixIn.literature_list = []

    def test1(self):
        class ClassLitrefs(literature.LiteratureReferencesMixIn):
            """ some comments with added lit refs, i.e.  [#]_ and [3]_ etc...

                .. [#] Gratia, Exempli, *Journal of I.M.*, **1**, 1--100 (2009).
                .. [3] Dude, John, *The Intern. Foo Journal of Bars*, **51**, 1647--1751 (2009)  
                .. [4] Hat, John, *The Ex. Foo Journal of Bars*, **51**, 1647--1751 (2009)  
            """
            def __init__(self):
                literature.LiteratureReferencesMixIn.__init__(self)

        classnames = literature.LiteratureReferencesMixIn.names_of_classes_with_references()
        self.assertFalse("ClassLitrefs" in classnames)

        instance = ClassLitrefs()
        
        classnames = literature.LiteratureReferencesMixIn.names_of_classes_with_references()
        self.assertTrue("ClassLitrefs" in classnames)

    def test2(self):
        class ClassLitrefs(literature.LiteratureReferencesMixIn):
            """ some comments with added lit refs, i.e.  [#]_ and [3]_ etc...

                .. [#] Gratia, Exempli, *Journal of I.M.*, **1**, 1--100 (2009).
                .. [3] Dude, John, *The Intern. Foo Journal of Bars*, **51**, 1647--1751 (2009)  
                .. [4] Hat, John, *The Gal. Foo Journal of Bars*, **3**, 16--51 (2009)  
            """
            def __init__(self):
                literature.LiteratureReferencesMixIn.__init__(self)
        
        literature.TrackLiteratureReferences.default().registered_classes = set([])
        string = literature.LiteratureReferencesMixIn.all_literature_references_string()
        self.assertTrue("AMUSE" in string)
        self.assertTrue("multiphysics and multiscale software environment" in string)

        instance = ClassLitrefs()
        
        string = literature.LiteratureReferencesMixIn.all_literature_references_string()
        literature.LiteratureReferencesMixIn.print_literature_references()
        self.assertFalse(".. [#] Gratia, Exem" in string)
        self.assertFalse(".. [3] Dude" in string)
        self.assertFalse(".. [4] Hat" in string)

        self.assertTrue("Gratia, Exem" in string)
        self.assertTrue("Dude" in string)
        self.assertTrue("Hat" in string)

        self.assertTrue("ClassLitrefs" in string)
