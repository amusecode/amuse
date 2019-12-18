import numpy
from amuse.lab import *
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.ext.spherical_model \
    import new_uniform_spherical_particle_distribution

from matplotlib import pyplot
from prepare_figure import single_frame, figure_frame
from distinct_colours import get_distinct

def binned_mean_data(r, x):
    R = numpy.arange(0, r[-1], 0.1)
    X = numpy.zeros(len(R))
    N = numpy.zeros(len(R))
    for i in range(len(R)-1):
        for j in range(len(r)):
            if r[j]>=R[i] and r[j]<=R[i+1]:
                X[i] += x[j]
                N[i] += 1.
    for i in range(len(X)):
        if X[i]>0 and N[i]>0:
            X[i] = X[i]/float(N[i])

    return R, X

def get_ionization_fraction1d(ism):
    r = [] | units.parsec
    x = []
    for p in ism:
        r.append(p.position.length())
        x.append(p.xion)
    r, x = list(zip(*sorted(zip(r.value_in(units.parsec), x))))
    R, X = binned_mean_data(r, x)

    return r, x, R, X

def get_ionization_fraction2d(ism):
    x = []
    y = []
    xi = []
    h = []
    l = 0
    for p in ism:
        xp = p.x.value_in(units.parsec)
        yp = p.y.value_in(units.parsec)
        zp = p.z.value_in(units.parsec)
        if abs(zp) <= 0.5:
            x.append(xp)
            if abs(xp) > l: l = abs(xp)
            y.append(yp)
            if abs(yp) > l: l = abs(yp)
            xi.append(p.xion)
            h.append(p.h_smooth.value_in(units.parsec))
    return x, y, xi, h, l

def plot_ionization(ism, N, which, rmax, t_end, rS, ip):

    r, x, R, X = get_ionization_fraction1d(ism)
    xx, yy, xi, hh, ll = get_ionization_fraction2d(ism)

    w = 14
    if ip > 1:
        pyplot.close('all')
        pyplot.clf()
    pyplot.figure(figsize=(w,6))
    pyplot.rcParams.update({'font.size': 22})

    colors = get_distinct(3)
#    rscolor = 'green'
    
    ax1 = pyplot.subplot(1,2,1)
    #ax1.scatter(r, x, c=colors[0], lw=0, s=10)
    ax1.scatter(r, x, c=x, lw=0, s=20, alpha=0.5, cmap="jet")
    ax1.plot(R[1:], X[1:], c=colors[1], lw=3)
    ax1.plot([rS, rS], [-1,2], color=colors[2], linestyle='dashed', lw=3)
    ax1.text(rS+0.08, 1.06, r'$R_s$', color=colors[2], fontsize=20)
    ax1.set_xlim(0, rmax)
    ax1.set_ylim(-0.04, 1.19)
    ax1.set_xlabel("r [pc]", fontsize=20)
    ax1.set_ylabel(r'$\xi_{\rm ion}$', fontsize=20)

    ax2 = pyplot.subplot(1,2,2)
    h = numpy.array(hh)
    h *= 72*w/(6.*ll)		# approximate scaling
    ax2.set_aspect(1)
    sc2 = ax2.scatter(xx, yy, c=xi, s=numpy.pi*h**2,
                      alpha=0.05, edgecolor=None, cmap="jet")
    if ll < rmax: ll = rmax
    ax2.set_xlim(-ll,ll)
    ax2.set_ylim(-ll,ll)
    ax2.set_xlabel("x [pc]", fontsize=20)
    ax2.set_ylabel("y [pc]", fontsize=20)
    #for axi in ax2.flat:
    ax2.xaxis.set_major_locator(pyplot.MaxNLocator(6))
    ax2.yaxis.set_major_locator(pyplot.MaxNLocator(6))
    #pyplot.colorbar(sc2, ax=ax2, fraction=0.046, pad=0.04)	# magic numbers!
    circle = pyplot.Circle((0, 0), rS, color=colors[2],
                           fill=False, linestyle='dashed', lw=3)
    ax2.add_artist(circle)
    
    if which == 0:
        id_string = 'Plummer_model'
    else:
        id_string = 'Homogeneous_sphere'

    param_string = "_N=%d_t=%.3f_Myr" % (N, t_end)

    #pyplot.suptitle(id_string+param_string, y=0.96, fontsize=16)

    savefile = 'fig_ionization_of_GMC_'+id_string+param_string+'.png'
    pyplot.savefig(savefile, dpi=300)
    print('Figure saved in file', savefile)
    
    if ip == 0: pyplot.show()

def plot_density(r, rho, rho0, a):
    pyplot.figure(figsize=(10,8))
    pyplot.scatter(r, rho, lw=0, s=10)
    rhop = rho0*(1+(r/a)**2)**(-2.5)
    pyplot.scatter(r, rhop, lw=0, s=10)
    pyplot.xlim(0, 6)
    pyplot.ylim(0.005, 10.)
    pyplot.semilogy()
    pyplot.show()

###BOOKLISTSTART1###
def generate_ism_initial_conditions(N, M=10|units.MSun, R=3|units.parsec,
                                    boxsize=10|units.parsec, which=1):
    
    converter = nbody_system.nbody_to_si(M, R)

    # Option 0: Plummer model with mass M and virial radius R.
    # Option 1: Homogeneous sphere with radius R and same central
    #		density as the Plummer model.

    a = R/1.7				# Plummer parameter
    rho0 = 3*M/(4*numpy.pi*a**3)	# Plummer model central density
    rhp = 1.3*a				# Plummer model half-mass radius

    print('a =', a.in_(units.parsec))
    print('rhp =', rhp.in_(units.parsec))
    print('rho0 =', rho0.in_(units.MSun/units.parsec**3))

    if which == 0:
        rmax = boxsize
        ism = new_plummer_gas_model(N, converter)
    else:
        Ru = R
        rmax = Ru
        ism = new_uniform_spherical_particle_distribution(N, Ru,
                                                      4*numpy.pi*rho0*Ru**3/3,
                                                          type="random")
        rhu = Ru/3**0.5
        print('rhu =', rhu.in_(units.parsec))
        ism.vx = 0.|units.kms
        ism.vy = 0.|units.kms
        ism.vz = 0.|units.kms
        ism.u = (0.075|units.kms)**2
        
        print('M =', ism.mass.sum().in_(units.MSun))

    #print 'mean u =', ism.u.sum().in_(units.kms**2)/N
    #print 'max u =', ism.u.max().in_(units.kms**2)

    rr = ((ism.position)**2).sum(axis=1).sqrt().value_in(units.parsec)
    ii = numpy.argsort(rr)
    print('rh =', rr[ii[N/2]], 'parsec')

    ism.flux = 0. | units.s**-1
    ism.xion = 0.0

    hydro = Fi(converter)
    hydro.gas_particles.add_particles(ism)
    hydro.evolve_model(1|units.yr)
    hydro.gas_particles.new_channel_to(ism).copy()
    hydro.stop()

    '''
    plot_density(rr, ism.rho.value_in(units.MSun/units.parsec**3),
                 rho0.value_in(units.MSun/units.parsec**3),
                 a.value_in(units.parsec))
    '''

    #print 'max density =', ism.rho.max().in_(units.MSun/units.parsec**3)
    rc = a/3.
    cent = ism.select(lambda r: r.length() < rc, ["position"])
    print('approximate central density =', \
        (3*cent.mass.sum()/(4*numpy.pi*rc**3)).in_(units.MSun/units.parsec**3))

    ism = ism.select(lambda r: r.length() < 0.5*boxsize, ["position"])
    #print "max density in box =", \
    #      ism.rho.max().in_(units.MSun/units.parsec**3), \
    #      ism.rho.max().in_(units.amu/units.cm**3)
    if rmax > boxsize/2: rmax = boxsize/2
    
    return ism,rho0,a,M,rmax
###BOOKLISTSTOP1###

def radius_containing(mass, ism):
    rr = ((ism.position)**2).sum(axis=1).sqrt().value_in(units.parsec)
    isort = numpy.argsort(rr)
    msum = zero
    rp = []
    mp = []
    i = 0
    #print msum, mass.in_(units.MSun)
    while msum < mass:
        msum += ism.mass[isort[i]]
        rp.append(rr[isort[i]])
        mp.append(msum.value_in(units.MSun))
        if i >= len(rr)-1: break
        i += 1

    #pyplot.plot(rp, mp)
    #pyplot.show()
                  
    return rp[-1]|units.parsec
    
###BOOKLISTSTART###
def main(N, Lstar, boxsize, t_end, which, np):
    ism,rho0,a,M,rmax \
        = generate_ism_initial_conditions(N, boxsize=boxsize, which=which)
    print('particles in box =', len(ism))

    source = Particle()
    source.position = (0, 0, 0) |units.parsec
    source.flux = Lstar/(20. | units.eV)
    source.rho = ism.rho.max()
    source.xion = ism.xion.max()
    source.u = (9.|units.kms)**2	# = kT/m

    S = source.flux.in_(units.s**-1)
    print('ionizing flux =', S)
    print('central density =', rho0.in_(units.MSun/units.parsec**3))
    n = rho0/(1|units.amu)
    print('number density =', n.in_(units.cm**-3))
    T = source.u*(1|units.amu)/constants.kB
    print('temperature =', T.in_(units.K))
    alpha = (2.e-16*T.value_in(units.K)**-0.75)|(units.m**3/units.s)
    print('alpha =', alpha.in_(units.m**3/units.s))
    trec = 1/(n*alpha)
    print('recombination time =', trec.in_(units.Myr))
    MS = (S*trec)*1|units.amu
    if which == 0:
        rS = a/((M/MS)**(2./3) - 1.)**0.5		# Plummer
    else:
        rS = (3*S/(4*numpy.pi*n**2*alpha))**(1./3)	# Homogeneous

    print('Stromgren radius (analytic) =', rS.in_(units.parsec))
    rS1 = radius_containing(MS, ism)
    print('Stromgren radius (measured) =', rS1.in_(units.parsec))

    radiative = SimpleX(redirection='none')
    radiative.parameters.box_size = boxsize
    radiative.parameters.timestep = 0.125*trec
    #radiative.parameters.timestep = 0.05|units.Myr
    radiative.particles.add_particle(source)
    radiative.particles.add_particles(ism)

    te = t_end.value_in(units.Myr)
    tlist = numpy.linspace(te/np, te, np)|units.Myr

    ip = 0
    for t in tlist:
        radiative.evolve_model(t)
        print('t =', t, \
              'model_time =', radiative.model_time.in_(units.Myr))

        radiative.particles.new_channel_to(ism).copy()

        print('')
        print("min ionization:", ism.xion.min())
        print("average ionization:", ism.xion.mean())
        print("max ionization:", ism.xion.max())

        if np > 1: ip += 1
        plot_ionization(ism, N, which,
                        rmax.value_in(units.parsec), t.value_in(units.Myr),
                        rS1.value_in(units.parsec), ip)

    radiative.stop()
###BOOKLISTSTOP###

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", type="int",
                      dest="np", default = 1,
                      help="number of plots [%default]")
    result.add_option("-N", type="int",
                      dest="N", default = 10000,
                      help="number of gas particles [%default]")
    result.add_option("-t", unit=units.Myr, type="float",
                      dest="t_end", default = 1.0|units.Myr,
                      help="radiation time [%default]")
    result.add_option("-L", unit=units.LSun, type="float",
                      dest="Lstar", default = 1.e2|units.LSun,
                      help="luminosity of ionizing source [%default]")
    result.add_option("-d", unit=units.parsec,
                      dest="boxsize", default = 10|units.parsec,
                      help="size of the density box [%default]")
    result.add_option("-w", type='int',
                      dest="which", default = 1,
                      help="which model [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    numpy.random.seed(12345)
    main(**o.__dict__)

