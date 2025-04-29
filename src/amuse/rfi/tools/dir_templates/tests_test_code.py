from amuse.support.testing.amusetest import TestWithMPI

from amuse.community.{code}.interface import {interface_class}, {user_class}

class {interface_class}Tests(TestWithMPI):

    def test_echo_int(self):
        instance = {interface_class}()
        result,error = instance.echo_int(12)
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()
