"""
Plot a Hertzsprung-Russell diagram of a star cluster, for a number of
snapshots, and make a movie.

Uses data from
https://astronomy.stackexchange.com/questions/39994/what-is-the-rgb-curve-for-blackbodies
to convert temperature to RGB.

Shows a progress bar.
"""

import sys
import os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from amuse.io import read_set_from_file
from amuse.units import units, constants

# package for progress bar
from tqdm import tqdm


def lumrad_to_temp(luminosity, radius):
    temperature = ((
        luminosity
        / (constants.four_pi_stefan_boltzmann * radius**2)
    )**0.25).in_(units.K)
    return temperature


def lumtemp_to_rad(luminosity, temperature):
    radius = ((luminosity / (constants.four_pi_stefan_boltzmann * temperature**4))**0.5).in_(units.RSun)
    return radius


def templum_to_xyz(temperature, luminosity):
    log_temperature = np.nan_to_num(np.log10(temperature.value_in(units.K)))
    log_luminosity = np.nan_to_num(np.log10(luminosity.value_in(units.LSun)))
    color = temp_to_rgb(temperature)
    return log_temperature, log_luminosity, color


def temp_to_rgb(temperature):
    temp = temperature.value_in(units.K)
    logT = np.log(temp)
    logT1000 = np.log(temp - 1000.0)
    rgb = np.zeros((len(temp), 3))
    rgb[:, 0] = 1.0
    rgb[:, 1] = 0.390081972 * logT - 2.427925631
    rgb[:, 2] = 0.543206396 * logT1000 - 3.698136688

    t6600 = temp > 6600
    rgb[t6600, 0] = 2.4054 * (temp[t6600]-6000)**(-0.1332047592)
    rgb[t6600, 1] = 1.6 * (temp[t6600]-6000)**(-0.0755148492)
    rgb[t6600, 2] = 1.0
   
    rgb = np.clip(rgb, 0, 1)
    return rgb



class StarHRPlotter:
    def __init__(self, name, extension="amuse"):
        self.fig = plt.figure(figsize=(10, 10))
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel("log(Teff)")
        self.ax.set_ylabel("log(L)")
        self.temperature_range = [5.5, 3]
        self.luminosity_range = [-5, 9]
        self.ax.set_xlim(self.temperature_range)
        self.ax.set_ylim(self.luminosity_range)
        self.ax.set_facecolor("k")
        self.scatter = None
        self.scatter2 = None
        self.name = name
        self.extension = extension
        self.ndigit = 6


    def make_movie(self, start, end):
        "Find all snapshots, and make a movie"
        def update(frame):
            """
            update frame and update progress bar. The progress bar doesn't work
            yet so printing dots too.
            """
            print(".", end="", flush=True)
            i = start + frame
            
            filename = f"{self.name}{i:0{self.ndigit}d}.{self.extension}"
            stars = read_set_from_file(filename)
            size = 4 * stars.radius.value_in(units.RSun)**0.5
            x, y, color = templum_to_xyz(stars.temperature, stars.luminosity)
            try:
                age = set(stars.age.value_in(units.Myr))
                if len(age) == 1:
                    age = age.pop()
                self.ax.set_title(f"Snapshot {i} ({age:.2f} Myr)")
            except:
                self.ax.set_title(f"Snapshot {i}")
            self.scatter.set_offsets(np.array([x, y,]).T)
            self.scatter.set_sizes(size)
            self.scatter.set_facecolors(color)


        i = start
        filename = f"{self.name}{i:0{self.ndigit}d}.{self.extension}"
        stars = read_set_from_file(filename)
        x, y, color = templum_to_xyz(stars.temperature, stars.luminosity)
        # radius_is_zero = stars.radius == 0 | units.RSun
        # print(radius_is_zero)
        # stars[radius_is_zero] = lumtemp_to_rad(
        #     stars[radius_is_zero].luminosity,
        #     stars[radius_is_zero].temperature,
        # )
        size = 4 * stars.radius.value_in(units.RSun)**0.5
        self.scatter = self.ax.scatter(x, y, s=size, c=color, edgecolor="none")
        anim = animation.FuncAnimation(
            self.fig,
            update,
            frames=tqdm(range(end - start), position=0, file=sys.stdout),
            interval=30,
            repeat=False,
        )
        anim.save(f"{self.name}.mp4", dpi=150, writer=animation.FFMpegWriter(fps=25))


    def templum_to_xy(self):
        log_temperature = np.nan_to_num(np.log10(stars.temperature.value_in(units.K)))
        log_luminosity = np.nan_to_num(np.log10(stars.luminosity.value_in(units.LSun)))


    def plot_hr(self, stars, stars2=None):
        log_temperature = np.nan_to_num(np.log10(stars.temperature.value_in(units.K)))
        log_luminosity = np.nan_to_num(np.log10(stars.luminosity.value_in(units.LSun)))
        col = temp_to_rgb(stars.temperature)
        if not self.scatter:
            self.scatter = self.ax.scatter(
                log_temperature,
                log_luminosity,
                s=42,
                c=temp_to_rgb(stars.temperature),
                edgecolor="none",
            )
        else:
            self.scatter.set_offsets(
                np.array(
                    [
                        log_temperature,
                        log_luminosity,
                    ]
                ).T
            )
            self.scatter.set_facecolors(col)
        if stars2:

            log_temperature2 = np.nan_to_num(np.log10(stars2.temperature.value_in(units.K)))
            log_luminosity2 = np.nan_to_num(np.log10(stars2.luminosity.value_in(units.LSun)))
            if not self.scatter2:
                self.scatter2 = self.ax.scatter(
                    log_temperature2,
                    log_luminosity2,
                    s=2,
                )
            else:
                self.scatter2.set_offsets(
                    np.array(
                        [
                            log_temperature2,
                            log_luminosity2,
                        ]
                    ).T
                )


    def savefig(self, filename):
        self.fig.savefig(filename)



def new_argument_parser():
    "Parse command line arguments, show defaults"
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--infile",
        type=str,
        default="",
        help="The first snapshot to plot",
    )
    parser.add_argument(
        "-I", 
        "--infile2",
        type=str,
        default="",
        help="Second set of snapshots to plot",
    )
    parser.add_argument(
        "-n",
        "--number",
        type=int,
        default=1,
        help="The number of snapshots to plot",
    )
    return parser    


def main():
    args = new_argument_parser().parse_args()
    filename_template = args.infile
    # check extension of the file
    filename_template = filename_template.split(".")
    extension = filename_template[-1]
    name = filename_template[0]
    if args.infile2:
        filename_template2 = args.infile2
        filename_template2 = filename_template2.split(".")
        extension2 = filename_template2[-1]
        name2 = filename_template2[0]
    # check if the name ends with a number, if so, check how many characters
    # are there and store it and strip it
    i = 0
    while name[-i - 1].isdigit():
        i += 1 
    if i > 1:
        name = name[:-i]
        if args.infile2:
            name2 = name2[:-i]
        ndigit = i
    else:
        raise ValueError("The name of the first snapshot should end with a number")

    # set up plotter
    plotter = StarHRPlotter(name)
    plotter.make_movie(0, args.number)
    sys.exit()

    # read the snapshots
    for i in range(args.number):
        filename = f"{name}{i:0{ndigit}}.{extension}"
        snapshot = read_set_from_file(filename)
        if args.infile2:
            filename2 = f"{name2}{i:0{ndigit}}.{extension}"
            snapshot2 = read_set_from_file(filename2)
            print(filename, filename2)
            plotter.plot_hr(snapshot, snapshot2)
            plotter.savefig(f"{name}{name2}{i:0{ndigit}}.png")
        else:
            print(filename)
            plotter.plot_hr(snapshot)
            plotter.savefig(f"rgb{name}{i:0{ndigit}}.png")


if __name__ == "__main__":
    main()
