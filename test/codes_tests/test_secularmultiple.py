from amuse.test.amusetest import TestWithMPI
import os
import sys
import numpy

from amuse.community.secularmultiple.interface import SecularMultipleInterface, SecularMultiple

class TestSecularMultipleInterface(TestWithMPI):
    def test1(self):
        instance = SecularMultipleInterface()

class TestSecularMultiple(TestWithMPI):
    def test1(self):
        instance = SecularMultiple()
        
        print instance.parameters
