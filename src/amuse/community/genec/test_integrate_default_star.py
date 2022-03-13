import numpy

from amuse.datamodel import Particle
from amuse.units import units
from amuse.community.genec import Genec
from amuse.io import write_set_to_file


def write_backup(
    step,
    star,
    # density,
    # radius,
    # temperature,
    # luminosity,
    # pressure,
    abundances,
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

    write_set_to_file(
        star.as_set(),
        f'star-{star.key}-{step}.amuse'
    )
    # units = dict()
    # units['density'] = f'{density.unit}'
    # units['radius'] = f'{radius.unit}'
    # units['temperature'] = f'{temperature.unit}'
    # units['luminosity'] = f'{luminosity.unit}'
    # units['pressure'] = f'{pressure.unit}'
    
    numpy.savez_compressed(
        f'star-abundances-{star.key}-{step}.npz',
        # units=units,
        # density=density.value_in(density.unit),
        # radius=radius.value_in(radius.unit),
        # temperature=temperature.value_in(temperature.unit),
        # luminosity=luminosity.value_in(luminosity.unit),
        # pressure=pressure.value_in(pressure.unit),
        abundances=abundances,
    )
    return

star = Particle(mass=7 | units.MSun, metallicity=0.014)
evo = Genec()
star_in_evo = evo.fullparticles.add_particle(star)

plot_every = -1
save_every = 100
step = 0

while True:
    star = star_in_evo.copy()
    # number_of_zones = star_in_evo.get_number_of_zones()
    # density_profile = star_in_evo.get_density_profile()
    # radius_profile = star_in_evo.get_radius_profile()
    # temperature_profile = star_in_evo.get_temperature_profile()
    # luminosity_profile = star_in_evo.get_luminosity_profile()
    # pressure_profile = star_in_evo.get_pressure_profile()
    chemical_abundance_profile = star_in_evo.get_chemical_abundance_profiles()

    print(evo.particles[0])
    print(star.age.in_(units.Myr), star.mass.in_(units.MSun))
    if step%save_every == 0:
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

    star_in_evo.evolve_one_step()
    step += 1
