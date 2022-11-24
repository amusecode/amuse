import numpy
import time

from amuse.datamodel import Particle
from amuse.units import units
from amuse.community.genec import Genec
from amuse.io import write_set_to_file, read_set_from_file
from amuse.support.console import set_printing_strategy

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from plot_models import StellarModelPlot

from amuse.community.genec.interface import SPECIES_NAMES

class NeedToSave:
    def __init__(self):
        self.model = 0
        self.current = {
            'luminosity': 0 | units.LSun,
            'temperature': 0 | units.K,
            'central_temperature': 0 | units.K,
            'central_density': 0 | units.g * units.cm**-3,
            'surface_velocity': 0 | units.kms,
            'mass': 0 | units.MSun,
            'time_step': 0 | units.s,
        }
        for species in SPECIES_NAMES.keys():
            self.current[species] = []
        self.previous = self.current.copy()
        self.derivs = {}
        for key in self.current.keys():
            self.derivs[key] = 1
        return

    def get_values(self):
        return self.current

    @property
    def save_next(self):
        return True

    def update(self, star):
        self.previous = self.current.copy()
        for key in self.current.keys():
            self.current[key] = getattr(star, key)


def read_saved_star_timeline(star_key):
    star = read_set_from_file(f'star-{star_key}.amuse')[0]
    age, radius = star.get_timeline_of_attribute_as_vector('radius')
    print(age.in_(units.yr))
    print(radius.in_(units.RSun))


def write_backup(
    step,
    star,
    abundances,
    append=True,
):
    # GENEC reads the following on a restore/resume:
    #   gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,ab,dm_lost,m,(q(i),p(i),t(i),r(i),s(i),x(i),y(i),&
    #           xc12(i),vp(i),vt(i),vr(i),vs(i),xo16(i),vx(i),vy(i),vxc12(i),vxo16(i),i=1,m),drl,drte,dk,drp,drt,drr,rlp,&
    #           rlt,rlc,rrp,rrt,rrc,rtp,rtt,rtc,tdiff,vsuminenv,(CorrOmega(i),i=1,npondcouche),xLtotbeg,dlelexprev,&
    #           zams_radius
    #
    #   read (51) (y3(i),xc13(i),xn14(i),xn15(i),xo17(i),xo18(i),vy3(i),vxc13(i),vxn14(i),vxn15(i),vxo17(i),vxo18(i),xne20(i), &
    #           xne22(i),xmg24(i),xmg25(i),xmg26(i),vxne20(i),vxne22(i),vxmg24(i),vxmg25(i),vxmg26(i),omegi(i),vomegi(i),i=1,m)
    #
    #   read(51) (xf19(i),xne21(i),xna23(i),xal26(i),xal27(i),xsi28(i),vxf19(i),vxne21(i),vxna23(i),vxal26g(i),vxal27(i), &
    #           vxsi28(i),xneut(i),xprot(i),xc14(i),xf18(i),xbid(i),xbid1(i),vxneut(i),vxprot(i),vxc14(i),vxf18(i),vxbid(i), &
    #           vxbid1(i),i=1,m)

    #   do ii=1,nbelx
    #    read(51) (abelx(ii,i),vabelx(ii,i),i=1,m)
    #   enddo
    #
    #   if (isugi >= 1) then
    #     read(51) nsugi
    #   endif
    #
    #   if (bintide) then
    #     read(51) period,r_core,vna,vnr
    #   endif

    if append:
        filename = f'star-{star.key}.amuse'
    else:
        filename = f'star-{star.key}-{step}.amuse'
    write_set_to_file(
        star.as_set(),
        filename,
        timestamp=star.age if append else None,
        append_to_file=append,
        compression=True,
    )

    # For now, abundances aren't part of the single star particle
    # numpy.savez_compressed(
    #     f'star-abundances-{star.key}-{step}.npz',
    #     abundances=abundances,
    # )
    return

ROTATING_STAR = True
MASS_UNIT = units.MSun
LENGTH_UNIT = units.RSun
SPEED_UNIT = units.kms
TIME_UNIT = units.Myr
MASSLOSS_UNIT = units.MSun / units.yr
TEMPERATURE_UNIT = units.K
LUMINOSITY_UNIT = units.LSun
SPEEDUP_UNIT = units.Myr / units.minute
set_printing_strategy(
    "custom",
    preferred_units=[
        MASS_UNIT, LENGTH_UNIT, TIME_UNIT, SPEED_UNIT, MASSLOSS_UNIT,
        TEMPERATURE_UNIT, LUMINOSITY_UNIT, SPEEDUP_UNIT
    ],
    precision=6,
    prefix="",
    separator=" ",
    # separator=" [",
    suffix="",
    # suffix="]",
)

star_difficult = Particle(
    mass=60 | units.MSun,
    metallicity=0.014,
)
star = Particle(mass=60.0 | units.MSun, metallicity=0.014)
evo = Genec(redirection="none")
# evo = Genec()
if ROTATING_STAR:
    evo.parameters.vwant = 0.7
    evo.parameters.ianiso = 1
    # evo.parameters.irot = 1
    # evo.parameters.starname = "AmuseDifficultStar"
else:
    evo.parameters.vwant = 0.0
evo.parameters.ipoly = 0
# evo.parameters.idebug = 1
star_in_evo = evo.fullparticles.add_particle(star)

# for p in params.items():
#     setattr(evo.parameters, p[0], p[1])
evo.commit_particles()
print(evo.parameters)
# print(star_in_evo)
font = {
    'size': 8,
}
plt.rc('font', **font)
# iplt.switch_backend('macosx')
plt.ion()

save_every = 10
store_every = 1
plot_time = 10 | units.s
plot_models = 20
step = 0

model_of_last_save = 0
model_of_last_plot = 0
time_start = time.time() | units.s
time_of_last_plot = 0 | units.s
age_of_last_plot = star_in_evo.age

plotting = None
plotting = StellarModelPlot(star_in_evo)

# evo.parameters.nzmod = 100
evo.parameters.n_snap = 0

print("age   mass   radius   temp   lum   phase   vequat   h0   vwant  xcn")
while True:
    time_elapsed = (time.time() | units.s) - time_start
    star = star_in_evo.copy()
    # number_of_zones = star_in_evo.get_number_of_zones()
    # density_profile = star_in_evo.get_density_profile()
    # radius_profile = star_in_evo.get_radius_profile()
    # temperature_profile = star_in_evo.get_temperature_profile()
    # luminosity_profile = star_in_evo.get_luminosity_profile()
    # pressure_profile = star_in_evo.get_pressure_profile()
    chemical_abundance_profile = star_in_evo.get_chemical_abundance_profiles()

    # print(evo.fullparticles[0])
    # print(evo.fullparticles[0].get_number_of_species())
    # print(evo.fullparticles[0].get_names_of_species())
    # print(evo.fullparticles[0].get_mass_profile())
    # exit()
    print(
        star.age.in_(units.Myr),
        star.mass.in_(units.MSun),
        star.radius.in_(units.RSun),
        star.temperature.in_(units.K),
        star.luminosity.in_(units.LSun),
        evo.parameters.phase,
        star.surface_velocity,
        star.abundance_h[0],
        evo.parameters.vwant,
        evo.parameters.xcn,
    )
    print(f"step: {step} time: {star.age} timestep: {star.time_step} xcn: {evo.parameters.xcn}")    
    if (step % store_every == 0) and plotting is not None:
        plotting.update(star_in_evo, phase=evo.parameters.phase)
        if (
            (time_elapsed - time_of_last_plot) > plot_time
            or step - model_of_last_plot > plot_models
        ):
            speed = (
                (star.age - age_of_last_plot).value_in(units.Myr)
                / (time_elapsed - time_of_last_plot).value_in(units.minute)
            ) | units.Myr / units.minute
            plotting.plot_all(speed=speed, step=step)
            model_of_last_plot = step
            time_of_last_plot = time_elapsed
            age_of_last_plot = star.age

    if step % save_every == 0:
        write_backup(
            step,
            star,
            # density_profile,
            # radius_profile,
            # temperature_profile,
            # luminosity_profile,
            # pressure_profile,
            chemical_abundance_profile,
        )

    age_previous = star_in_evo.age
    star_in_evo.evolve_one_step()
    # print(f"Condition: {evo.parameters.stopping_condition}")
    # if evo.parameters.stopping_condition != "none":
    #     # star_in_evo.age == age_previous:
    #     if step > 1:
    #         print("Stopping - not evolving!")
    #         print(f"Condition: {evo.parameters.stopping_condition}")
    #         break
    step += 1

runtime = (time.time() | units.s) - time_start
print(
    f"Running {step} models took {runtime.value_in(units.minute)} minutes"
)
