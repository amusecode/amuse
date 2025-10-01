import pytest
from amuse.support.testing import amusetest

from amuse.units import units
from amuse.units.trigo import arctan2


class TestsForIssue1180(amusetest.TestCase):

    def test_arctan2_without_units(self):
        "Test when input has no units"
        x = 1.0
        y = 2.0
        result = arctan2(x, y).value_in(units.rad)
        assert result == pytest.approx(0.4636476, 1e-7)

    def test_arctan2_with_units(self):
        "Test when input has units"
        x = 1.0 | units.m
        y = 2.0 | units.m
        result = arctan2(x, y).value_in(units.rad)
        assert result == pytest.approx(0.4636476, 1e-7)

    def test_arctan2_with_units_and_no_units(self):
        "Test when input has units and no units (should fail)"
        x = 1.0 | units.m
        y = 2.0
        with pytest.raises(AttributeError):
            result = arctan2(x, y).value_in(units.rad)

    def test_arctan2_with_different_units(self):
        "Test when input has different units"
        x = 1.0 | units.m
        y = 2.0e-3 | units.km
        result = arctan2(x, y).value_in(units.rad)
        assert result == pytest.approx(0.4636476, 1e-7)
