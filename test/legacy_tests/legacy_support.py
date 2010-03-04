from amuse.test import amusetest
from amuse.legacy.support.core import is_mpd_running, stop_interfaces

class TestWithMPI(amusetest.TestCase):

    def setUp(self):
        self.assertTrue(is_mpd_running(), "MPICH2 mpd deamon process not running, cannot run this test as it requires MPI")
            
    def tearDown(self):
        stop_interfaces()
            
