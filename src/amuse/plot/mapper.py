#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make 2D maps from SPH particles using FiMap
"""

import logging
import numpy
from amuse.units import units, constants, nbody_system
from amuse.datamodel.rotation import new_rotation_matrix, rotate
from amuse.community.fi import FiMap


def gas_mean_molecular_weight(
    h2ratio=0.5,
):
    "Return mean molecular weight of hydrogen gas"
    meanmwt = (
        (
            2.0 * h2ratio
            + (1. - 2. * h2ratio)
            + 0.4
        ) /
        (
            0.1 + h2ratio + (1. - 2. * h2ratio)
        )
    ) | units.amu
    return meanmwt


def u_to_temperature(
    internal_energy,
    meanmwt=gas_mean_molecular_weight(),
    gamma=5/3,
):
    """
    Convert internal energy to temperature, by default assumes all gas is
    molecular hydrogen and gamma = 5/3.
    """
    temperature = (
        (gamma - 1) * meanmwt * internal_energy
        / constants.kB
    )
    return temperature


def temperature_to_u(
    temperature,
    meanmwt=gas_mean_molecular_weight(),
    gamma=5/3,
):
    """
    Convert temperature to internal energy, by default assumes all gas is
    molecular hydrogen and gamma = 5/3.
    """
    internal_energy = (
        constants.kB * temperature
        / ((gamma - 1) * meanmwt)
    )
    return internal_energy


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
        mapper.origin = (0, 0, 0) | units.pc

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
        self.__origin = [0., 0., 0.] | units.pc
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

        self.__phi = 0 | units.deg
        self.__theta = 0 | units.deg
        self.__psi = 0 | units.deg

        self.__constant_mass = False
        if not gas.is_empty():
            if self.__gas.mass.std() == 0 | self.__unit_mass:
                self.__constant_mass = self.__gas[0].mass

    def set_unit_length(self, unit=units.pc):
        # need to assert this is indeed a length
        self.__unit_length = unit

    def rotate(
            self,
            phi=0 | units.deg,
            theta=0 | units.deg,
            psi=0 | units.deg,
    ):
        dphi = phi - self.__phi
        dtheta = theta - self.__theta
        dpsi = psi - self.__psi
        self.__phi = phi
        self.__theta = theta
        self.__psi = psi
        if self.__stars is not None:
            if not self.__stars.is_empty():
                rotate(self.__stars, dphi, dtheta, dpsi)
        if self.__sinks is not None:
            if not self.__sinks.is_empty():
                rotate(self.__sinks, dphi, dtheta, dpsi)

        if self.__state == "RUN":
            rotation_matrix = new_rotation_matrix(
                self.__phi, self.__theta, self.__psi,
            ).transpose()
            self.__mapper.parameters.upvector = numpy.inner(
                rotation_matrix,
                [0, 1, 0],
            )
            self.__mapper.parameters.projection_direction = numpy.inner(
                rotation_matrix,
                [0, 0, -1],
            )
            del self.__maps
            self.__maps = self.Maps()

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
        converter = nbody_system.nbody_to_si(100 | units.MSun, 1000 | units.pc)
        gas = self.__gas
        # x_axis = self.__axis_x
        # y_axis = self.__axis_y

        mapper = FiMap(converter, mode="openmp", redirection="none")
        # if not hasattr(gas, "radius"):
        # gas.radius = gas.h_smooth
        mapper.particles.add_particles(gas)

        rotation_matrix = new_rotation_matrix(
            self.__phi, self.__theta, self.__psi,
        ).transpose()
        mapper.parameters.upvector = numpy.inner(
            rotation_matrix,
            [0, 1, 0],
        )
        mapper.parameters.projection_direction = numpy.inner(
            rotation_matrix,
            [0, 0, -1],
        )
        mapper.parameters.image_target = self.__origin
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
    def phi(self):
        return self.__phi

    @property
    def theta(self):
        return self.__theta

    @property
    def psi(self):
        return self.__psi

    @property
    def xmin(self):
        "Return smallest value on x axis"
        return self.origin[0] - (self.width/2)

    @property
    def xmax(self):
        "Return largest value on x axis"
        return self.origin[0] + (self.width/2)

    @property
    def ymin(self):
        "Return smallest value on x axis"
        return self.origin[1] - (self.height/2)

    @property
    def ymax(self):
        "Return largest value on x axis"
        return self.origin[1] + (self.height/2)

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
    def origin(self):
        "Get origin (rotated)"
        rotation_matrix = new_rotation_matrix(
            self.__phi, self.__theta, self.__psi,
        ).transpose()
        return self.__origin.dot(rotation_matrix)

    @origin.setter
    def origin(self, origin=[0., 0., 0.] | units.pc):
        "Set origin offset"
        delta_origin = -self.__origin
        for i, o_i in enumerate(origin):
            delta_origin[i] += o_i
        self.__origin = origin
        if self.__state in ["RUN"]:
            self.__mapper.parameters.image_target = self.__origin

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
            return self.counts * self.__constant_mass

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
        if hasattr(gas, "mu"):
            print(f"h2ratio: {gas.mu.mean()}")
            meanmwt = gas.mu
        elif hasattr(gas, "h2ratio"):
            print(f"h2ratio: {gas.h2ratio.mean()}")
            meanmwt = gas_mean_molecular_weight(gas.h2ratio)
        else:
            meanmwt = gas_mean_molecular_weight()
        temperature = u_to_temperature(gas.u, meanmwt=meanmwt)
        print(f"temperature range: {min(temperature)} - {max(temperature)}")
        self.__mapper.particles.weight = temperature.value_in(
            self.__unit_temperature)
        # counts = self.counts
        self.__maps.temperature = numpy.nan_to_num(
            self.__mapper.image.pixel_value.transpose(),
            # / counts,
            nan=0,
        ) | self.__unit_temperature
        print(
            f"min/max: {self.__maps.temperature.min()} "
            f"- {self.__maps.temperature.max()}"
        )
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
