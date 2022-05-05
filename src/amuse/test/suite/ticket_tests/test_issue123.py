from amuse.test import amusetest

from amuse.community.sse.interface import SSE
from amuse.community.bhtree.interface import BHTree

class TestsForIssue123(amusetest.TestCase):
  def test1(self): # doesn't trigger recursion error
    self.assertRaises(Exception, BHTree, name_of_the_worker="bogus",expected_message= 
      "__init__() got multiple values for keyword argument 'name_of_the_worker'")
  def test2(self): # does
    self.assertRaises(Exception, SSE, name_of_the_worker="bogus",expected_message= 
      "__init__() got multiple values for keyword argument 'name_of_the_worker'")
