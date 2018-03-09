import numpy
from matplotlib import pyplot
from amuse.lab import *
from prepare_figure import single_frame

def make_e_map(sph,N=100,L=1):

    x,y=numpy.indices( ( N+1,N+1 ))

    x=L*(x.flatten()-N/2.)/N
    y=L*(y.flatten()-N/2.)/N
    z=x*0.
    vx=0.*x
    vy=0.*x
    vz=0.*x

    x=units.parsec(x)
    y=units.parsec(y)
    z=units.parsec(z)
    vx=units.kms(vx)
    vy=units.kms(vy)
    vz=units.kms(vz)

    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
    rho=rhoe.reshape((N+1,N+1))

    return rho

def setup_grid(N, L):
    x,y=numpy.indices( ( N+1,N+1 ))
    x=L*(x.flatten()-N/2.)/N
    y=L*(y.flatten()-N/2.)/N
    z=0.*x 
    vx=0.*x
    vy=0.*x
    vz=0.*x
    x=units.parsec(x)
    y=units.parsec(y)
    z=units.parsec(z)
    vx=units.kms(vx)
    vy=units.kms(vy)
    vz=units.kms(vz)
    return x, y, z, vx, vy, vz

def make_hydromap_and_show_picture(sph_particles, N=100, L=10):
    x_label = "x [R$_\odot$]"
    y_label = "y [R$_\odot$]"
    fig = single_frame(x_label, y_label, logx=False, logy=False, xsize=12, ysize=12)

    hydro = Gadget2(converter)
    hydro.gas_particles.add_particles(sph_particles)
    x, y, z, vx, vy, vz = setup_grid(N, L)
    rho,rhovx,rhovy,rhovz,rhoe=hydro.get_hydro_state_at_point(x,y,z,vx,vy,vz)
    rho=rhoe.reshape((N+1,N+1))
    rho_e=make_map(hydro,N=50,L=L)
    hydro.stop()
    print("extrema:", rho_e.value_in(units.erg/units.RSun**3).min(), rho_e.value_in(units.erg/units.RSun**3).max())
    cax = pyplot.imshow(numpy.log10(rho_e.value_in(units.erg/units.RSun**3)), extent=[-L/2,L/2,-L/2,L/2],vmin=4,vmax=11)
    cbar = fig.colorbar(cax, ticks=[4, 7.5, 11], orientation='vertical', fraction=0.045)
    cbar.ax.set_yticklabels(['Low', ' ', 'High'])  # horizontal colorbar
    cbar.set_label('mid-plane energy-density', rotation=270)
    
    t = int(0.5+gas.get_timestamp().value_in(units.s))
    filename = "supernova_sph_T"+str(t)+".pdf"
    pyplot.savefig(filename)
#    pyplot.show()

def make_map(sph,N=100,L=1):

    x,y=numpy.indices( ( N+1,N+1 ))

    x=L*(x.flatten()-N/2.)/N
    y=L*(y.flatten()-N/2.)/N
    z=x*0.
    vx=0.*x
    vy=0.*x
    vz=0.*x

    x=units.RSun(x)
    y=units.RSun(y)
    z=units.RSun(z)
    vx=units.kms(vx)
    vy=units.kms(vy)
    vz=units.kms(vz)

    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
    rho=rho.reshape((N+1,N+1))

    return rho

def plot_e_sph(sph, time):
#    pyplot.rcParams.update({'font.size': 30})
#    figure = pyplot.figure(figsize=(12, 12))

    L = 10
    max_dens = sph.gas_particles.rho.value_in(units.g/units.cm**3).max()
    print("Density extrema:", max_dens)

    x_label = "X [R$_\odot$]"
    y_label = "Y [R$_\odot$]"
    figure = single_frame(x_label, y_label, logx=False, logy=False, xsize=12, ysize=12)


    rho_e=make_e_map(sph,N=20,L=L)
    cax = pyplot.imshow(rho_e.value_in(units.erg/units.MSun**3), extent=[-L/2,L/2,-L/2,L/2], interpolation='bicubic', origin = 'lower', cmap="hot")
    cbar = figure.colorbar(cax, ticks=[1.e-8, 0.5*max_dens, max_dens], orientation='vertical', fraction=0.045)
    rmin = 0.0 
    rmid = "%.1f" % (0.5*max_dens)
    rmax = "%.1f" % (max_dens)
    cbar.ax.set_yticklabels([rmin, ' ', rmax])  # horizontal colorbar
    cbar.set_label('mid-plane density [$erg/M_\odot^3$]', rotation=270)
    pyplot.xlabel("x [R$_\odot$]")
    pyplot.ylabel("y [R$_\odot$]")

    t = int(0.5+time.value_in(units.s))
    filename = "supernova_sph_T"+str(t)+".pdf"
    pyplot.savefig(filename)
    pyplot.show()

def plot_sph(sph, time):
    L = 10
    unit = units.g/units.cm**3
    max_dens = sph.gas_particles.rho.value_in(unit).max()
    print("Density extrema:", max_dens)

    x_label = "X [R$_\odot$]"
    y_label = "Y [R$_\odot$]"
    figure = single_frame(x_label, y_label, logx=False, logy=False, xsize=12, ysize=12)

    rho_e=make_map(sph,N=50,L=L)
    max_dens = rho_e.value_in(unit).max()
    print("extrema:", rho_e.value_in(unit).min(), max_dens)
    cax = pyplot.imshow(rho_e.value_in(unit), extent=[-L/2,L/2,-L/2,L/2], interpolation='bicubic', origin='lower', cmap="hot", vmin=0.0, vmax=max_dens)    
#    cbar = figure.colorbar(cax, ticks=[-0.4*max_dens, 0.0*max_dens, 0.5*max_dens], orientation='vertical', fraction=0.045)
    cbar = figure.colorbar(cax, ticks=[0.0, 0.5*max_dens, 0.99*max_dens], orientation='vertical', fraction=0.045)

    rmin = "0.0"
    rmid = "%.1f" % (0.5*max_dens)
    rmax = "%.1f" % (max_dens)
    print("min/max=", rmin, rmid, rmax)
    cbar.ax.set_yticklabels([rmin, ' ', rmax])  # horizontal colorbar
#    cbar.ax.set_yticklabels(['a', 'b', 'c'])  # horizontal colorbar
    cbar.set_label('mid-plane density [$g/cm^3$]', rotation=270)
    pyplot.xlabel("x [R$_\odot$]")
    pyplot.ylabel("y [R$_\odot$]")

    t = int(0.5+time.value_in(units.s))
    filename = "supernova_sph_T"+str(t)+".pdf"
    pyplot.savefig(filename)
    pyplot.show()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-t", unit=units.s,
                      dest="tplot", type="float", default = 300|units.s,
                      help="plotting time [%default]")
    result.add_option("-f", 
                      dest="filename", default = "supernova_sph_gadget.amuse",
                      help="input filename [%default]")
    return result
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    converter = nbody_system.nbody_to_si(1|units.MSun, 1|units.RSun)
    bodies = read_set_from_file(o.filename, "amuse")
    t_end = o.tplot
    for gas in bodies.history:
        time = gas.get_timestamp()
        print("time=", time.in_(units.s), "N=", len(gas))
        if int(0.5+time.value_in(units.s)) >= t_end.value_in(units.s):
            gas.move_to_center()
            hydro = Gadget2(converter)
            hydro.gas_particles.add_particles(gas)
            plot_sph(hydro, time)
            hydro.stop()
            break

