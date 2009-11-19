import unittest
from amuse.legacy.support import lit

class TestLiteratureRefs(unittest.TestCase):
    def setUp(self):
        lit.LiteratureRefs.literature_list = []

    def test1(self):
        class ClassLitrefs(lit.LiteratureRefs):
            """ some comments with added lit refs, i.e.  [#]_ and [3]_ etc...

                .. [#] Gratia, Exempli, *Journal of I.M.*, **1**, 1--100 (2009).
                .. [3] Dude, John, *The Intern. Foo Journal of Bars*, **51**, 1647--1751 (2009)  
                .. [4] Hat, John, *The Ex. Foo Journal of Bars*, **51**, 1647--1751 (2009)  
            """
            def __init__(self):
                lit.LiteratureRefs.__init__(self)

        classnames = lit.LiteratureRefs.names_of_classes_with_references()
        self.assertFalse("ClassLitrefs" in classnames)

        instance = ClassLitrefs()
        
        classnames = lit.LiteratureRefs.names_of_classes_with_references()
        self.assertTrue("ClassLitrefs" in classnames)

    def test2(self):
        class ClassLitrefs(lit.LiteratureRefs):
            """ some comments with added lit refs, i.e.  [#]_ and [3]_ etc...

                .. [#] Gratia, Exempli, *Journal of I.M.*, **1**, 1--100 (2009).
                .. [3] Dude, John, *The Intern. Foo Journal of Bars*, **51**, 1647--1751 (2009)  
                .. [4] Hat, John, *The Gal. Foo Journal of Bars*, **3**, 16--51 (2009)  
            """
            def __init__(self):
                lit.LiteratureRefs.__init__(self)

        string = lit.LiteratureRefs.all_literature_references_string()
        self.assertEquals("", string)

        instance = ClassLitrefs()
        
        string = lit.LiteratureRefs.all_literature_references_string()
        lit.LiteratureRefs.print_refs()
        self.assertFalse(".. [#] Gratia, Exem" in string)
        self.assertFalse(".. [3] Dude" in string)
        self.assertFalse(".. [4] Hat" in string)

        self.assertTrue("Gratia, Exem" in string)
        self.assertTrue("Dude" in string)
        self.assertTrue("Hat" in string)

        self.assertTrue("ClassLitrefs" in string)
        
    def test3(self):
        class ClassLitrefs(lit.LiteratureRefs):
            """ some comments with added lit refs, i.e.  [#]_ and [3]_ etc...

                .. [#] Gratia, Exempli, *Journal of I.M.*, **1**, 1--100 (2009).
                .. [3] Dude, John, *The Intern. Foo Journal of Bars*, **51**, 1647--1751 (2009)  
                .. [4] Hat, John, *The Gal. Foo Journal of Bars*, **3**, 16--51 (2009)  
            """
            def __init__(self):
                lit.LiteratureRefs.__init__(self)
        
        string = lit.LiteratureRefs.export2html()
        print string
        #from docutils import core
        #print core.publish_string(source = string)
        
