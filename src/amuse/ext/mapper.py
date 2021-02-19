#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make 2D maps from SPH particles
"""

from amuse.units import units, constants, nbody_system
from amuse.community.fi.interface import FiMap


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
            self.temp = None
            self.vx = None
            self.vy = None
            self.vz = None

    def __init__(self, gas=None, stars=None, sinks=None):
        self.set_bins()
        self.set_width()
        self.set_axes()
        self.set_offset()
        self.__maps = self.Maps()
        self.__unit_mass = units.MSun
        self.__unit_length = units.pc
        self.__unit_time = units.Myr
        self.__unit_speed = units.kms
        self.__unit_temp = units.K
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
                print('Wrong input for x_axis or y_axis: please check!')
                return None

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
                print('Wrong input for x_axis or y_axis: please check!')
                return None

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
                print('Wrong input for x_axis or y_axis: please check!')
                return None

        mapper.parameters.target_x = self.__offset_x
        mapper.parameters.target_y = self.__offset_y
        mapper.parameters.target_z = self.__offset_z
        mapper.parameters.image_width = self.__width
        mapper.parameters.image_size = [self.__bins_x, self.__bins_y]
        self.__mapper = mapper
        self.__state = "RUN"

    def set_axes(self, x='x', y='y', z='z'):
        "Set axes"
        self.__axis_x = x
        self.__axis_y = y
        self.__axis_z = z
        self.__state = "EDIT"

    def set_offset(self, x=0 | units.pc, y=0 | units.pc, z=0 | units.pc):
        "Set origin offset"
        self.__offset_x = x
        self.__offset_y = y
        self.__offset_z = z
        self.__state = "EDIT"

    def get_offset(self):
        "Get origin offset"
        return self.__offset_x, self.__offset_y, self.__offset_z

    def set_width(self, width=1 | units.pc):
        "Set width of map"
        self.__width = width
        self.__state = "EDIT"

    def get_width(self):
        "Get width of map"
        return self.__width

    def set_bins(self, x=800, y=None):
        "Set number of bins"
        if y is None:
            y = x
        self.__bins_x = x
        self.__bins_y = y
        self.__state = "EDIT"

    def get_bins(self):
        "Get number of bins"
        return self.__bins_x, self.__bins_y

    def get_counts(self):
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

    def get_mass(self):
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

    def get_column_density(self):
        "Return a column density map"
        length_x = self.__width / self.__bins_x
        length_y = length_x * (self.__bins_y / self.__bins_x)
        pixel_size = length_x * length_y
        map_mass = self.get_mass()

        return map_mass / pixel_size

    def get_temperature(self):
        "Return a temperature map"
        if (
            (self.__maps.temp is not None)
            and (self.__state != "EDIT")
        ):
            return self.__maps.temp

        if self.__state != "RUN":
            self.__new_gas_mapper()

        gas = self.__gas
        if hasattr(gas, "h2ratio"):
            gmmw = gas_mean_molecular_weight(gas.h2ratio)
        else:
            gmmw = gas_mean_molecular_weight()
        temperature = u_to_temperature(gas.u, gmmw=gmmw)
        self.__mapper.particles.weight = temperature.value_in(
            self.__unit_temp)
        counts = self.get_counts()
        self.__maps.temp = (
            self.__mapper.image.pixel_value.transpose(
            ) | self.__unit_temp
        ) / counts
        return self.__maps.temp

    def get_vx(self):
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

    def get_vy(self):
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

    def get_vz(self):
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
