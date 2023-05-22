#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot hydro and/or stars
"""
import os
import sys
import logging
import numpy
import copy
import argparse

import matplotlib
import matplotlib.pyplot as plt
# import matplotlib.cm as cm

from amuse.datamodel import Particles
from amuse.io import read_set_from_file
from amuse.units import units, nbody_system
from amuse.units.quantities import as_vector_quantity
from amuse.datamodel.rotation import new_rotation_matrix

from amuse.plot.mapper import MapHydro


logger = logging.getLogger(__name__)


def velocity_divergence(vx_field, vy_field, dx, velocity_unit=units.kms):
    # Return divergence of velocity fields

    div = (numpy.ufunc.reduce(
        numpy.add,
        [
            -1 * numpy.gradient(vx_field.value_in(velocity_unit), axis=1),
            -1 * numpy.gradient(vy_field.value_in(velocity_unit), axis=0),
        ]
    ) | velocity_unit) / dx

    return div


def plot_column_density(
    plot_axes,
    maps,
    unit_col_density=units.MSun * units.pc**-2,
    vmin=-1,
    vmax=4,
    cmap='viridis',
):
    cmap = copy.copy(matplotlib.colormaps[cmap])
    cmap.set_bad('k', alpha=1.0)

    column_density = maps.column_density
    logscale = numpy.log10(column_density.value_in(unit_col_density))
    return plot_axes.imshow(
        logscale,
        extent=maps.extent,
        vmin=vmin,
        vmax=vmax,
        origin='lower',
        cmap=cmap,
    )


def plot_temperature(
    plot_axes,
    maps,
    vmin=0,
    vmax=5,
    cmap="inferno",
):
    cmap = copy.copy(matplotlib.colormaps[cmap])
    cmap.set_bad('k', alpha=1.0)
    temperature_map = maps.temperature
    logscale_temperature_map = numpy.log10(
        temperature_map.value_in(units.K)
    )

    return plot_axes.imshow(
        logscale_temperature_map,
        extent=maps.extent,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        origin="lower",
    )


def plot_divergence(
    plot_axes,
    maps,
    x_axis="vx",
    y_axis="vy",
    contours=False,
    normalise_by_counts=True,
    vmin=None,
    vmax=None,
    div_unit=units.Myr**-1,
):
    dx = maps.width / maps.bins[0]
    dy = maps.width / maps.bins[1]
    assert (dx == dy), "dx/dy are spaced differently!"
    x_range = numpy.arange(
        maps.xmin.value_in(units.pc), maps.xmax.value_in(units.pc),
        (maps.width / maps.bins[0]).value_in(units.pc),
    )
    y_range = numpy.arange(
        maps.ymin.value_in(units.pc), maps.ymax.value_in(units.pc),
        (maps.height / maps.bins[1]).value_in(units.pc),
    )
    x_grid, y_grid = numpy.meshgrid(x_range, y_range)
    div = velocity_divergence(
        getattr(maps, x_axis) / maps.counts,
        getattr(maps, y_axis) / maps.counts,
        dx,
    )
    if not normalise_by_counts:
        div = div * maps.mass

    if contours:
        levels = 1 / ([-4, -10, -100, 100, 10, 4] | units.Myr)
        return plot_axes.contour(
            x_grid, y_grid,
            div.value_in(div_unit),
            levels=levels.value_in(div_unit),
            # colors=["blue", "red"],
            # cmap="coolwarm",
            cmap='bwr',
        )
    minmax = max(-div.min(), div.max())
    if not vmin:
        vmin = -minmax
    if not vmax:
        vmax = minmax
    # minmax = 500
    # print(minmax)
    # print(vmin, vmax)
    return plot_axes.imshow(
        div.value_in(div_unit),
        extent=maps.extent,
        vmin=vmin.value_in(div_unit),
        vmax=vmax.value_in(div_unit),
        cmap='bwr',
        origin='lower',
    )


def new_figure(
    aspect=1,
    colorbar_width=0 | units.cm,
    title_margin=0 | units.cm,
    xlabel_margin=0.9 | units.cm,
    ylabel_margin=1.25 | units.cm,
    right_padding=0.18 | units.cm,
    top_padding=0.15 | units.cm,
    # plot_style="default",
):
    "Return a new figure and a new axes object"
    # /plt.style.use(plot_style)
    fig = plt.figure()

    # Get the width of the figure - this is the fixed value
    # height will depend on aspect and other attributes
    figure_width = (fig.get_figwidth() | units.inch)

    # Calculate in 'figure units'
    left_frac = (ylabel_margin) / figure_width
    right_frac = (colorbar_width + right_padding) / figure_width
    if colorbar_width > 0 | units.cm:
        right_frac += (ylabel_margin) / figure_width

    left = left_frac
    width = 1 - left - right_frac

    # now we set the figure height to get the aspect ratio right
    width_cm = aspect * width * figure_width

    figure_height = width_cm + top_padding + title_margin + xlabel_margin
    fig.set_figheight(figure_height.value_in(units.inch))

    bottom_frac = (xlabel_margin) / figure_height
    top_frac = (title_margin + top_padding) / figure_height

    bottom = bottom_frac
    height = 1 - bottom - top_frac

    axes = fig.add_axes(
        [
            left,
            bottom,
            width,
            height,
        ]
    )
    if colorbar_width > 0 | units.cm:
        cax_width = (colorbar_width) / figure_width
        colorbar_axes = fig.add_axes(
            [
                left+width,
                bottom,
                cax_width,
                height,
            ]
        )
        return fig, axes, colorbar_axes

    return fig, axes


def plot_hydro_and_stars(
    maps,
    fig=None,
    ax=None,
    title="",
    plot="column density",
    length_unit=units.parsec,
    colorbar=True,
    **kwargs
):
    "Plot gas and stars"
    maps.set_unit_length(length_unit)

    if fig:
        if ax is None:
            ax = fig.gca()
        ax.cla()
        # for artist in ax.lines + ax.collections:
        #     artist.remove()
        cax = False
    else:
        if colorbar:
            fig, ax, cax = new_figure(
                colorbar_width=0.2 | units.cm
            )
        else:
            fig, ax = new_figure()
            cax = False
    if "stars" in plot:
        ax.set_facecolor('black')

    gasplot = False
    if plot == "column density":
        gasplot = plot_column_density(ax, maps)
        gasplot_unit = units.MSun * units.pc**-2
    elif plot == "temperature":
        gasplot = plot_temperature(ax, maps)
        gasplot_unit = units.K
    elif plot == "divergence":
        gasplot = plot_divergence(
            ax,
            maps,
            normalise_by_counts=True,
            div_unit=units.Myr**-1,
        )
    elif plot in ["hubble", "hubbledust"]:
        from amuse.plot.fresco.fresco import make_image as fresco_image
        if plot == "hubbledust":
            extinction = True
        else:
            extinction = False

        converter = nbody_system.nbody_to_si(1 | units.pc, 1 | units.MSun)
        image = fresco_image(
            stars=maps.stars,
            gas=maps.gas if extinction is True else None,
            image_width=maps.width,
            image_size=maps.bins,
            return_vmax=False,
            converter=converter,
            extinction=extinction,
            # **kwargs
        )
        plt.imshow(
            image,
            extent=maps.extent,
            origin="lower",
        )

    cmap = copy.copy(matplotlib.colormaps["coolwarm"])

    if (
            (maps.stars is not None)
            and not maps.stars.is_empty()
    ):
        s = 1 * (
            (maps.stars.mass / (7 | units.MSun))**(3.5 / 2)
        )
        # s = 0.5
        x = getattr(maps.stars, 'x').value_in(length_unit)
        y = getattr(maps.stars, 'y').value_in(length_unit)
        z = getattr(maps.stars, 'z').value_in(length_unit)
        c = (
            "cyan" if plot in [
                "temperature",
            ]
            else "white"
        )
        ax.scatter(x, y, s=s, c=c, lw=0)

    ax.set_xlabel("%s (%s)" % (maps.axes[0], length_unit))
    ax.set_ylabel("%s (%s)" % (maps.axes[1], length_unit))
    ax.set_xlim(maps.extent[0], maps.extent[1])
    ax.set_ylim(maps.extent[2], maps.extent[3])

    if cax:
        cbar = plt.colorbar(gasplot, cax=cax)
        if gasplot_unit == units.MSun * units.pc**-2:
            cbar.set_label("log10(M$_{\odot}$ / parsec$^2$)")
        elif gasplot_unit == units.Myr**-1:
            cbar.set_label("Myr$^{-1}$")
        elif gasplot_unit == units.K:
            cbar.set_label("log10(temp [K])")
        else:
            cbar.set_label("%s" % gasplot_unit)

    return gasplot


def new_argument_parser():
    "Parse command line arguments"
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        '-s',
        dest='starsfilename',
        default='',
        help='file containing stars (optional)',
    )
    parser.add_argument(
        '-g',
        dest='gasfilename',
        default='',
        help='file containing gas (optional)',
    )
    parser.add_argument(
        '-o',
        dest='imagefilename',
        default='test',
        help='write image to this file',
    )
    parser.add_argument(
        '-n',
        dest='bins',
        default=800,
        type=int,
        help='number of bins',
    )
    parser.add_argument(
        '-x',
        dest='x',
        default=0 | units.pc,
        type=units.pc,
        help='Central X coordinate',
    )
    parser.add_argument(
        '-y',
        dest='y',
        default=0 | units.pc,
        type=units.pc,
        help='Central Y coordinate',
    )
    parser.add_argument(
        '-z',
        dest='z',
        default=0 | units.pc,
        type=units.pc,
        help='Central Z coordinate',
    )
    parser.add_argument(
        '-w',
        dest='w',
        default=10 | units.pc,
        type=units.pc,
        help='Width',
    )
    parser.add_argument(
        '--com',
        dest='use_com',
        action='store_true',
        default=False,
        help='Center on center of mass',
    )
    parser.add_argument(
        '-X',
        dest='x_axis',
        default='x',
        help='Horizontal axis',
    )
    parser.add_argument(
        '-Y',
        dest='y_axis',
        default='y',
        help='Vertical axis',
    )
    parser.add_argument(
        '-Z',
        dest='z_axis',
        default='z',
        help='Line-of-sight axis',
    )
    parser.add_argument(
        '-f',
        dest='followfilename',
        default=None,
        help=(
            'file containing star keys to center on (optional)\n'
            '  implies --com'
        ),
    )
    # parser.add_argument(
    #     '--length',
    #     dest='length_unit',
    #     default='parsec',
    #     help='Length unit (default: parsec)',
    # )
    return parser.parse_args()


def main():
    args = new_argument_parser()
    # length_unit = getattr(units, args.length_unit)
    length_unit = units.pc
    gasfilename = args.gasfilename
    starsfilename = args.starsfilename
    imagefilename = args.imagefilename
    followfilename = args.followfilename
    bins = args.bins
    finish = as_vector_quantity([args.x, args.y, args.z])
    width_finish = args.w
    x_axis = args.x_axis
    y_axis = args.y_axis
    z_axis = args.z_axis
    bins = args.bins

    plots = []
    if followfilename is not None:
        use_com = True
    else:
        use_com = args.use_com

    if os.path.isfile(starsfilename):
        stars = read_set_from_file(
            starsfilename,
            "amuse",
        ) if starsfilename != "" else Particles()
    else:
        stars = Particles()
    if gasfilename:
        gas = read_set_from_file(
            gasfilename,
            "amuse",
        )
        if hasattr(gas, "itype"):
            gas = gas[gas.itype == 1]
        # gas.h_smooth = gas.h
    else:
        gas = Particles()

    if not gas.is_empty():
        plots += ["column density"]
        if hasattr(gas, 'u'):
            plots.append("temperature")

    if not gas.is_empty():
        mtot = gas.total_mass()
        com = mtot * gas.center_of_mass()
        time = gas.get_timestamp()
    else:
        mtot = 0 | units.MSun
        com = [0, 0, 0] | units.pc * units.MSun
        time = False
    if not stars.is_empty():
        mtot += stars.total_mass()
        com += stars.total_mass() * stars.center_of_mass()
        if not time:
            time = stars.get_timestamp()
    com = com / mtot
    if use_com:
        if followfilename is not None:
            followstars = read_set_from_file(
                followfilename, close_file=True,
            )
            center_on_these_stars = followstars.get_intersecting_subset_in(
                stars,
            )
            if center_on_these_stars.is_empty():
                print("Stars not found")
                sys.exit()
            com = center_on_these_stars.center_of_mass()
        finish = com

    if time is None:
        print('Unable to get timestamp, set time to 0 Myr.')
        time = 0.0 | units.Myr

    maps = MapHydro(gas, stars)
    maps.axes = (x_axis, y_axis, z_axis)
    maps.bins = bins

    def xyz_to_xzy():
        maps.rotate(90 | units.deg, 0 | units.deg, 0 | units.deg)

    def xyz_to_zyx():
        maps.rotate(0 | units.deg, 90 | units.deg, 0 | units.deg)

    def xyz_to_yzx():
        maps.rotate(0 | units.deg, 90 | units.deg, 90 | units.deg)

    maps.origin = finish
    maps.width = width_finish

    if x_axis == "x" and y_axis == "z":
        xyz_to_xzy()
    elif x_axis == "y" and y_axis == "z":
        xyz_to_yzx()

    for plot in plots:
        gasplot = plot_hydro_and_stars(
            maps,
            title="time = %06.2f %s" % (
                time.value_in(units.Myr),
                units.Myr,
            ),
            plot=plot,
            length_unit=length_unit,
        )
        plotname = (
            plot.replace(" ", "_")
            + "-" + imagefilename + ".png"
        )
        plt.savefig(
            plotname,
        )


if __name__ == "__main__":
    main()
