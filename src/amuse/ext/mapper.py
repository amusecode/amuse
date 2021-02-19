#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make 2D maps from SPH particles
"""

import logging
from amuse.units import units, constants, nbody_system
try:
    from amuse.community.fi.interface import FiMap
except ImportError:
    FiMap = False


def gas_mean_molecular_weight(h2ratio=1):
    "Return mean molecular weight of hydrogen gas"
    gmmw = (
        (
            2.0 * h2ratio
            + (1. - 2. * h2ratio)
            + 0.4
        ) /
        (
            0.1 + h2ratio + (1. - 2. * h2ratio)
        )
    ) | units.amu
    return gmmw


def u_to_temperature(
        internal_energy,
        gmmw=gas_mean_molecular_weight(),
):
    """
    Convert internal energy to temperature, by default assumes all gas is
    molecular hydrogen.
    """
    temperature = (
        2/3*internal_energy/(constants.kB/gmmw)
    )
    return temperature


class MapHydro():
    class Maps():
        "Store all calculated maps"

        def __init__(self):
            self.counts = None
            self.mass = None
            self.temperature = None
            self.vx = None
            self.vy = None
            self.vz = None

    def __init__(self, gas=None, stars=None, sinks=None):
        """
        Sets up MapHydro instance. Maps are created as-needed.
        Usage:
        mapper = MapHydro(gas=gas, stars=None, sinks=None)

        mapper.width = 10 | units.pc
        mapper.bins = 800
        mapper.offset = (0, 0, 0) | units.pc

        counts_map = mapper.counts
        density_map = mapper.density
        vx_map = mapper.vx
        vy_map = mapper.vy
        temperature_map = mapper.temperature
        mapper.stop()
        """
        self.__bins_x = 800
        self.__bins_y = 800
        self.__width = 1 | units.pc
        self.__axis_x = 'x'
        self.__axis_y = 'y'
        self.__axis_z = 'z'
        self.__offset_x = 0 | units.pc
        self.__offset_y = 0 | units.pc
        self.__offset_z = 0 | units.pc
        self.__maps = self.Maps()
        self.__unit_mass = units.MSun
        self.__unit_length = units.pc
        self.__unit_time = units.Myr
        self.__unit_speed = units.kms
        self.__unit_temperature = units.K
        self.__mapper = None

        self.__gas = gas
        self.__stars = stars
        self.__sinks = sinks
        self.__state = "INIT"

        self.__constant_mass = False
        if not gas.is_empty():
            if self.__gas.mass.std() == 0 | self.__unit_mass:
                self.__constant_mass = self.__gas[0].mass

    def stop(self):
        "Clean up"
        if self.__mapper is not None:
            self.__mapper.stop()
            self.__mapper = None
        del self.__gas
        del self.__maps
        self.__maps = self.Maps()
        self.__state = "STOP"

    def __new_gas_mapper(self):
        "Start a new mapper instance"
        converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.pc)
        gas = self.__gas
        x_axis = self.__axis_x
        y_axis = self.__axis_y

        mapper = FiMap(converter, mode="openmp", redirection="none")
        if not hasattr(gas, "radius"):
            gas.radius = gas.h_smooth
        mapper.particles.add_particles(gas)

        # positive x = up
        if y_axis == 'x':
            mapper.parameters.upvector = [1, 0, 0]
            if x_axis == 'y':
                # negative z = top layer
                mapper.parameters.projection_direction = [0, 0, 1]
            elif x_axis == 'z':
                # positive y = top layer
                mapper.parameters.projection_direction = [0, -1, 0]
            else:
                raise ValueError(
                    'Incorrect value for x_axis or y_axis'
                )

        # positive y = up
        if y_axis == 'y':
            mapper.parameters.upvector = [0, 1, 0]
            if x_axis == 'x':
                # positive z = top layer
                mapper.parameters.projection_direction = [0, 0, -1]
            elif x_axis == 'z':
                # negative x = top layer
                mapper.parameters.projection_direction = [1, 0, 0]
            else:
                raise ValueError(
                    'Incorrect value for x_axis or y_axis'
                )

        # positive z = up
        if y_axis == 'z':
            mapper.parameters.upvector = [0, 0, 1]
            if x_axis == 'x':
                # negative y = top layer
                mapper.parameters.projection_direction = [0, 1, 0]
            elif x_axis == 'y':
                # positive x = top layer
                mapper.parameters.projection_direction = [-1, 0, 0]
            else:
                raise ValueError(
                    'Incorrect value for x_axis or y_axis'
                )

        mapper.parameters.target_x = self.__offset_x
        mapper.parameters.target_y = self.__offset_y
        mapper.parameters.target_z = self.__offset_z
        mapper.parameters.image_width = self.__width
        mapper.parameters.image_size = [self.__bins_x, self.__bins_y]
        self.__mapper = mapper
        self.__state = "RUN"

    @property
    def axes(self):
        "Get axes"
        return (self.__axis_x, self.__axis_y, self.__axis_z)

    @axes.setter
    def axes(self, axes):
        "Set axes"
        self.__axis_x = axes[0]
        self.__axis_y = axes[1]
        self.__axis_z = axes[2]

    @property
    def xmin(self):
        "Return smallest value on x axis"
        return self.__offset_x - (self.width/2)

    @property
    def xmax(self):
        "Return largest value on x axis"
        return self.__offset_x + (self.width/2)

    @property
    def ymin(self):
        "Return smallest value on x axis"
        return self.__offset_y - (self.height/2)

    @property
    def ymax(self):
        "Return largest value on x axis"
        return self.__offset_y + (self.height/2)

    @property
    def extent(self):
        "Return image extent"
        return (
            self.xmin.value_in(self.__unit_length),
            self.xmax.value_in(self.__unit_length),
            self.ymin.value_in(self.__unit_length),
            self.ymax.value_in(self.__unit_length),
        )

    @property
    def offset(self):
        "Get origin offset"
        return self.__offset_x, self.__offset_y, self.__offset_z

    @offset.setter
    def offset(self, offset=(0, 0, 0) | units.pc):
        "Set origin offset"
        self.__offset_x = offset[0]
        self.__offset_y = offset[1]
        self.__offset_z = offset[2]
        self.__state = "EDIT"

    @property
    def height(self):
        "Get height of map"
        return self.__width * (self.__bins_y / self.__bins_x)

    @property
    def width(self):
        "Get width of map"
        return self.__width

    @width.setter
    def width(self, width=1 | units.pc):
        "Set width of map"
        self.__width = width
        self.__state = "EDIT"

    @property
    def bins(self):
        "Get number of bins"
        return self.__bins_x, self.__bins_y

    @bins.setter
    def bins(self, bins):
        "Set number of bins"
        if isinstance(bins, (int)):
            self.__bins_x = bins
            self.__bins_y = bins
        elif isinstance(bins, (list, tuple)):
            self.__bins_x = bins[0]
            self.__bins_y = bins[1]
        else:
            raise TypeError("bins needs to be 'int', 'tuple' or 'list'")
        self.__state = "EDIT"

    @property
    def counts(self):
        "Return a counts map"
        if (
            (self.__maps.counts is not None)
            and (self.__state != "EDIT")
        ):
            return self.__maps.counts

        if self.__state != "RUN":
            self.__new_gas_mapper()

        self.__mapper.particles.weight = 1
        self.__maps.counts = self.__mapper.image.pixel_value.transpose()
        return self.__maps.counts

    @property
    def mass(self):
        "Return a mass map"
        if self.__constant_mass:
            return self.__maps.counts * self.__constant_mass

        if (
            (self.__maps.mass is not None)
            and (self.__state != "EDIT")
        ):
            return self.__maps.mass

        if self.__state != "RUN":
            self.__new_gas_mapper()

        self.__mapper.particles.weight = self.__gas.mass.value_in(
            self.__unit_mass)
        self.__maps.mass = self.__mapper.image.pixel_value.transpose(
        ) | self.__unit_mass
        return self.__maps.mass

    @property
    def column_density(self):
        "Return a column density map"
        length_x = self.__width / self.__bins_x
        length_y = length_x * (self.__bins_y / self.__bins_x)
        pixel_size = length_x * length_y
        map_mass = self.mass

        return map_mass / pixel_size

    @property
    def temperature(self):
        "Return a temperature map"
        if (
            (self.__maps.temperature is not None)
            and (self.__state != "EDIT")
        ):
            return self.__maps.temperature

        if self.__state != "RUN":
            self.__new_gas_mapper()

        gas = self.__gas
        if hasattr(gas, "h2ratio"):
            gmmw = gas_mean_molecular_weight(gas.h2ratio)
        else:
            gmmw = gas_mean_molecular_weight()
        temperature = u_to_temperature(gas.u, gmmw=gmmw)
        self.__mapper.particles.weight = temperature.value_in(
            self.__unit_temperature)
        counts = self.counts
        self.__maps.temperature = (
            self.__mapper.image.pixel_value.transpose(
            ) | self.__unit_temperature
        ) / counts
        return self.__maps.temperature

    @property
    def vx(self):
        "Return a vx map"
        if (
            (self.__maps.vx is not None)
            and (self.__state != "EDIT")
        ):
            return self.__maps.vx

        if self.__state != "RUN":
            self.__new_gas_mapper()

        gas = self.__gas

        self.__mapper.particles.weight = gas.vx.value_in(self.__unit_speed)
        self.__maps.vx = self.__mapper.image.pixel_value.transpose(
        ) | self.__unit_speed
        return self.__maps.vx

    @property
    def vy(self):
        "Return a vy map"
        if (
            (self.__maps.vy is not None)
            and (self.__state != "EDIT")
        ):
            return self.__maps.vy

        if self.__state != "RUN":
            self.__new_gas_mapper()

        gas = self.__gas

        self.__mapper.particles.weight = gas.vy.value_in(self.__unit_speed)
        self.__maps.vy = self.__mapper.image.pixel_value.transpose(
        ) | self.__unit_speed
        return self.__maps.vy

    @property
    def vz(self):
        "Return a vz map"
        if (
            (self.__maps.vz is not None)
            and (self.__state != "EDIT")
        ):
            return self.__maps.vz

        if self.__state != "RUN":
            self.__new_gas_mapper()

        gas = self.__gas

        self.__mapper.particles.weight = gas.vz.value_in(self.__unit_speed)
        self.__maps.vz = self.__mapper.image.pixel_value.transpose(
        ) | self.__unit_speed
        return self.__maps.vz

    @property
    def gas(self):
        "Return gas"
        return self.__gas

    @property
    def stars(self):
        "Return stars"
        return self.__stars

    @property
    def sinks(self):
        "Return sinks"
        return self.__sinks