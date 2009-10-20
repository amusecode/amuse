import unittest

from amuse.legacy.support.core import is_mpd_running

class TestWithMPI(unittest.TestCase):

    def setUp(self):
        self.assertTrue(is_mpd_running(), "MPICH2 mpd deamon process not running, cannot run this test as it requires MPI")
            
            
            