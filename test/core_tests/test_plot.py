import numpy

from amuse.test import amusetest
from amuse.units import units, quantities

from amuse.support.console import set_printing_strategy

import amuse.plot as aplot

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

class TestPlot(amusetest.TestCase):

    def test1(self):
        """ Test a basic plot with units and labels"""
        if not HAS_MATPLOTLIB:
            return self.skip()
        pyplot.clf()

        x = numpy.linspace(0, 100, 100) | units.yr
        y = numpy.linspace(0, 200, 100) | units.RSun

        aplot.plot(x, y)
        self.assertEquals("[yr]", self.xaxis().get_label_text())
        self.assertEquals("[RSun]", self.yaxis().get_label_text())

        aplot.xlabel("time")
        aplot.ylabel("radius")
        self.assertEquals("time [yr]", self.xaxis().get_label_text())
        self.assertEquals("radius [RSun]", self.yaxis().get_label_text())

    def test2(self):
        """ Test a basic plot with and without units and labels"""
        if not HAS_MATPLOTLIB:
            return self.skip()
        pyplot.clf()

        x = numpy.linspace(0, 100, 100) | units.yr
        y = numpy.linspace(0, 200, 100)

        aplot.plot(x, y)
        self.assertEquals("[yr]", self.xaxis().get_label_text())
        self.assertEquals("", self.yaxis().get_label_text())

        aplot.xlabel("time")
        aplot.ylabel("radius")
        self.assertEquals("time [yr]", self.xaxis().get_label_text())
        self.assertEquals("radius ", self.yaxis().get_label_text())

    def test3(self):
        """ Test a plot with preferred units """
        if not HAS_MATPLOTLIB:
            return self.skip()
        pyplot.clf()

        x = numpy.linspace(0, 100, 100) | units.yr
        y = numpy.linspace(0, 200, 100) | units.RSun

        set_printing_strategy('custom', preferred_units=[units.Myr, units.AU])
        aplot.plot(x, y)

        self.assertEquals("[Myr]", self.xaxis().get_label_text())
        self.assertEquals("[AU]", self.yaxis().get_label_text())
        self.assertEquals((0., 0.0001), pyplot.xlim())

    def test4(self):
        """ Test text in a plot """
        if not HAS_MATPLOTLIB:
            return self.skip()
        pyplot.clf()
        set_printing_strategy('default')

        x = numpy.linspace(0, 100, 100) | units.yr
        y = numpy.linspace(0, 200, 100) | units.RSun

        aplot.plot(x, y)

        text = aplot.text(50|units.yr, 0.5|units.AU, "test text")

        self.assertEquals(50., text.get_position()[0])
        self.assertAlmostEquals(107.546995464, text.get_position()[1])

    def test5(self):
        """ Test errorbar plot """
        if not HAS_MATPLOTLIB:
            return self.skip()
        pyplot.clf()
        set_printing_strategy('default')

        x = numpy.linspace(0, 100, 100) | units.yr
        y = numpy.linspace(0, 200, 100) | units.RSun
        yerr = [2e5]*len(y) | units.km

        line = aplot.errorbar(x, y, yerr=yerr)

        points, caps, bars = line
        bottoms, tops = caps
        error_height = tops.get_ydata()[0] - bottoms.get_ydata()[0]

        self.assertAlmostEquals(0.575125808, error_height)

    def test6(self):
        """ Test setting the x limits on a plot """
        if not HAS_MATPLOTLIB:
            return self.skip()
        pyplot.clf()
        set_printing_strategy('default')

        x = numpy.linspace(0, 100, 100) | units.yr
        y = numpy.linspace(0, 200, 100) | units.RSun

        line = aplot.plot(x, y)

        aplot.xlim(0|units.yr, 2e9|units.s)

        self.assertAlmostEquals(0, pyplot.xlim()[0])
        self.assertAlmostEquals(63.37752924, pyplot.xlim()[1])

    def test7(self):
        """ Test setting the x and y limits in various ways"""
        if not HAS_MATPLOTLIB:
            return self.skip()
        pyplot.clf()
        set_printing_strategy('default')

        x = numpy.linspace(0, 100, 100)
        y = numpy.linspace(0, 200, 100) | units.RSun

        line = aplot.plot(x, y)

        aplot.xlim(-10, 80)
        self.assertEquals(-10, pyplot.xlim()[0])
        self.assertEquals(80, pyplot.xlim()[1])

        aplot.xlim(xmax=90)
        self.assertEquals(-10, pyplot.xlim()[0])
        self.assertEquals(90, pyplot.xlim()[1])

        aplot.ylim([-12, 110]|units.RSun)
        self.assertEquals(-12, pyplot.ylim()[0])
        self.assertEquals(110, pyplot.ylim()[1])

        aplot.ylim(ymin=1e6|units.km)
        self.assertAlmostEquals(1.43781452, pyplot.ylim()[0])
        self.assertEquals(110, pyplot.ylim()[1])

    def xaxis(self):
        return pyplot.gca().get_xaxis()

    def yaxis(self):
        return pyplot.gca().get_yaxis()

    def skip(self):
        print "Matplotlib not installed. Skipping test."
