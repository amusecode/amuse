"""
equilibrium_tides_averaged.py

Written by Dr. Alessandro Alberto Trani; adapted for testing Tsunami in AMUSE.

Secular tidal evolution integrator for a two-body system using the equilibrium tide model and an RK45 solver.
Defines the EqTides class, which encapsulates tidal parameters and orbital/spin vectors, computes derivatives, and evolves the system.
"""

import numpy as np
from scipy.integrate import RK45
import tsunami
import matplotlib.pyplot as plt


class EqTides:
    """
    Equilibrium tide secular evolution for a two-body system.

    Attributes:
        m1 (float): Mass of body 1 [mass units].
        m2 (float): Mass of body 2 [mass units].
        R1 (float): Radius of body 1 [length units].
        R2 (float): Radius of body 2 [length units].
        k1 (float): Love number of body 1 [dimensionless].
        k2 (float): Love number of body 2 [dimensionless].
        tau1 (float): Tidal time lag for body 1 [time units].
        tau2 (float): Tidal time lag for body 2 [time units].
        rg1 (float): Gyration radius factor for body 1 [dimensionless].
        rg2 (float): Gyration radius factor for body 2 [dimensionless].
        rotdist (bool): Include rotational distortion (default True).
        usespin (bool): Evolve spin state (default True).
        KU (KeplerUtils): Utility for unit conversions.
        I1, I2 (float): Moments of inertia for bodies 1 and 2 [mass*length^2].
        Q1, Q2 (float): Eggleton parameter related to Love numbers [dimensionless].
        A1, A2, Tv1, Tv2 (float): Derived tidal coefficients [various units].
        noT1, noT2 (bool): Flags to disable tides if parameters are zero.
        integrator (RK45): Scipy integrator instance.
    """

    def __init__(self, m1, m2, R1, R2, k1, k2, tau1, tau2, rg1, rg2,
                 rotdist=True, usespin=True):
        """
        Initialize tidal and rotational parameters.

        Args:
            m1 (float): Mass of body 1 [mass units].
            m2 (float): Mass of body 2 [mass units].
            R1 (float): Radius of body 1 [length units].
            R2 (float): Radius of body 2 [length units].
            k1 (float): Love number of body 1 [dimensionless].
            k2 (float): Love number of body 2 [dimensionless].
            tau1 (float): Tidal time lag for body 1 [time units].
            tau2 (float): Tidal time lag for body 2 [time units].
            rg1 (float): Gyration radius factor for body 1 [dimensionless].
            rg2 (float): Gyration radius factor for body 2 [dimensionless].
            rotdist (bool, optional): Include rotational distortion. Defaults to True.
            usespin (bool, optional): Evolve spin state. Defaults to True.
        """
        self.KU = tsunami.KeplerUtils()
        self.k1, self.k2 = k1, k2
        self.tau1, self.tau2 = tau1, tau2
        self.m1, self.m2 = m1, m2
        self.mtot = m1 + m2
        self.mred = m1 * m2 / self.mtot
        self.R1, self.R2 = R1, R2
        self.I1 = m1 * (R1 * rg1) ** 2
        self.I2 = m2 * (R2 * rg2) ** 2

        # Eggleton parameters
        self.Q1 = 2 * k1 / (1 + 2 * k1)
        self.Q2 = 2 * k2 / (1 + 2 * k2)
        self.A1 = R1 ** 5 * self.Q1 / (1 - self.Q1)
        self.A2 = R2 ** 5 * self.Q2 / (1 - self.Q2)
        self.Tv1 = 3 * R1 ** 3 / m1 / (tau1 * k1) / (1 - self.Q1) ** 2
        self.Tv2 = 3 * R2 ** 3 / m2 / (tau2 * k2) / (1 - self.Q2) ** 2

        self.noT1 = tau1 == 0 or R1 == 0 or k1 == 0
        self.noT2 = tau2 == 0 or R2 == 0 or k2 == 0
        self.rotdist = rotdist
        self.usespin = usespin

        print("k1, k2 = {:e}, {:e}".format(self.k1, self.k2))
        print("tau1, tau2 = {:e}, {:e}".format(self.tau1, self.tau2))

        print("A1, A2 = {:e}, {:e}".format(self.A1, self.A2))
        print("Tv1, Tv2 = {:e}, {:e}".format(self.Tv1, self.Tv2))

        sigma1 = 2 / 3 * R1 ** 5 / self.A1 ** 2 * k1 * tau1
        sigma2 = 2 / 3 * R2 ** 5 / self.A2 ** 2 * k2 * tau2

        print("sigma1, sigma2 = {:e}, {:e}".format(sigma1, sigma2))

        print("Rotational distortion:", self.rotdist)
        if self.noT1: print("No tides on particle 1")
        if self.noT2: print("No tides on particle 2")

    def setup_orbital_vectors(self, a0, e0, i0=0.0, ome0=0.0, Ome0=0.0):
        """
         Define initial orbital vectors and compute derived parameters.

         Args:
             a0 (float): Semi-major axis [length units].
             e0 (float): Orbital eccentricity [dimensionless].
             i0 (float, optional): Inclination [radians]. Defaults to 0.0.
             ome0 (float, optional): Argument of periapse [radians]. Defaults to 0.0.
             Ome0 (float, optional): Longitude of ascending node [radians]. Defaults to 0.0.
         """
        self.P0 = 2 * np.pi * (a0 ** 3 / self.mtot) ** 0.5
        self.n0 = 2 * np.pi / self.P0
        rp = a0 * (1 - e0)
        print("P0 = {:g} yr".format(self.P0 * self.KU.Tscale))
        print("n0 = {:g} 1/yr".format(self.n0 / self.KU.Tscale))
        print("Rp, R1/Rp, R2/Rp = {:g}, {:e}, {:e}".format(rp, self.R1 / rp, self.R2 / rp))

        self.h = (self.mtot * (1 - e0 * e0) * a0) ** 0.5
        self.hvec = self.h * np.array([np.sin(i0) * np.sin(Ome0),
                                       -np.sin(i0) * np.cos(Ome0),
                                       np.cos(i0)])
        self.e = e0
        self.evec = self.e * np.array([np.cos(ome0) * np.cos(Ome0) - np.sin(ome0) * np.cos(i0) * np.sin(Ome0),
                                       np.cos(ome0) * np.sin(Ome0) + np.sin(ome0) * np.cos(i0) * np.cos(Ome0),
                                       np.sin(ome0) * np.sin(i0)])

        self.qvec = np.cross(self.hvec, self.evec)
        self.q = (self.qvec * self.qvec).sum() ** 0.5
        self.qhat = self.qvec / self.q

    def setup_spin_vectors(self, spin1, spin2, obl1, obl2, azim1=0.0, azim2=0.0):
        """
        Define initial spin rates and orientation vectors.

        Args:
            spin1 (float): Spin rate of body 1 [rad/time units].
            spin2 (float): Spin rate of body 2 [rad/time units].
            obl1 (float): Obliquity of body 1 [radians].
            obl2 (float): Obliquity of body 2 [radians].
            azim1 (float, optional): Azimuth angle of spin1 [radians]. Defaults to 0.0.
            azim2 (float, optional): Azimuth angle of spin2 [radians]. Defaults to 0.0.
        """
        print("Pspin1 = {:g} yr".format(2 * np.pi / spin1 * self.KU.Tscale))
        print("spin1 = {:g} 2pi/yr".format(spin1 / self.KU.Tscale))
        print("spin1/n0 = {:e}".format(spin1 / self.n0))
        print("Pspin2 = {:g} yr".format(2 * np.pi / spin2 * self.KU.Tscale))
        print("spin2 = {:g} 2pi/yr".format(spin2 / self.KU.Tscale))
        print("spin2/n0 = {:e}".format(spin2 / self.n0))

        self.spin1 = spin1
        self.spin1_vec = self.spin1 * np.array([np.sin(obl1) * np.sin(azim1),
                                                -np.sin(obl1) * np.cos(azim1),
                                                np.cos(obl1)])

        self.spin2 = spin2
        self.spin2_vec = self.spin2 * np.array([np.sin(obl2) * np.sin(azim2),
                                                -np.sin(obl2) * np.cos(azim2),
                                                np.cos(obl2)])

    def setup_vectors(self, y0=None):
        """
        Initialize integrator state vector.

        Args:
            y0 (ndarray, optional): Initial state array of length 12 [evec,hvec,spin1_vec,spin2_vec]. Defaults to None.
        """
        self.y = np.concatenate((self.evec, self.hvec, self.spin1_vec, self.spin2_vec)) if y0 is None else y0

    def calculate_quantities(self):
        """
        Update scalar orbital and spin quantities from state vector.

        No arguments.
        """
        self.evec, self.hvec, self.spin1_vec, self.spin2_vec = np.split(self.y, 4)

        e2 = (self.evec * self.evec).sum()
        self.e = e2 ** 0.5
        h2 = (self.hvec * self.hvec).sum()
        self.h = h2 ** 0.5
        omecc = (1 - e2)
        self.a = h2 / (self.mtot * omecc)
        self.spin1 = (self.spin1_vec * self.spin1_vec).sum() ** 0.5
        self.spin2 = (self.spin2_vec * self.spin2_vec).sum() ** 0.5
        self.obl1 = np.arccos((self.hvec * self.spin1_vec).sum() / (self.h * self.spin1))
        self.obl2 = np.arccos((self.hvec * self.spin2_vec).sum() / (self.h * self.spin2))

    def derivative_dydt(self, t, y):
        """
        Compute time derivatives for the integrator.

        Args:
            t (float): Current time [time units].
            y (ndarray): State vector [12] [evec,hvec,spin1_vec,spin2_vec].

        Returns:
            ndarray: Time derivative of state [12].
        """
        self.evec, self.hvec, self.spin1_vec, self.spin2_vec = np.split(y, 4)

        e2 = (self.evec * self.evec).sum()
        self.e = e2 ** 0.5
        self.ehat = self.evec / self.e

        h2 = (self.hvec * self.hvec).sum()
        self.h = h2 ** 0.5
        self.hhat = self.hvec / self.h

        self.qvec = np.cross(self.hvec, self.evec)
        self.q = (self.qvec * self.qvec).sum() ** 0.5
        self.qhat = self.qvec / self.q

        spin1_e = (self.spin1_vec * self.ehat).sum()
        spin1_h = (self.spin1_vec * self.hhat).sum()
        spin1_q = (self.spin1_vec * self.qhat).sum()

        spin2_e = (self.spin2_vec * self.ehat).sum()
        spin2_h = (self.spin2_vec * self.hhat).sum()
        spin2_q = (self.spin2_vec * self.qhat).sum()

        omecc = (1 - e2)
        self.a = h2 / (self.mtot * omecc)
        a3 = self.a * self.a * self.a
        n = (self.mtot / a3) ** 0.5
        a5 = a3 * self.a * self.a

        omecc2 = omecc * omecc
        omecc3 = omecc2 * omecc
        omecc5 = omecc2 * omecc2 * omecc

        fecc4 = (1.5 + e2 * 0.125) * e2 + 1
        fecc6 = (4.5 + e2 * 0.625) * e2 + 1
        y_diss_fac = 0.5 * fecc4 / (n * omecc5)
        x_diss_fac = 0.5 * fecc6 / (n * omecc5)
        fecc2 = ((e2 * 0.3125 + 5.625) * e2 + 7.5) * e2 + 1
        fecc3 = ((e2 * 7.8125e-2 + 1.875) * e2 + 3.75) * e2 + 1
        omecc3h = omecc3 ** 0.5
        v_diss_fac = fecc3 / omecc3h
        w_diss_fac = fecc2 / omecc3h
        fecc5 = (e2 * 0.375 + 3) * e2 + 1

        if not self.usespin:
            spin2_h = spin1_h = EqTides.pseudosync_spin(self.e, n)
            spin1_e = spin2_e = spin1_q = spin2_q = 0.0

        y1 = y2 = x1 = x2 = z1 = z2 = zgr = 0.0
        if not self.noT1:
            xyz_prefac1 = self.m2 * self.A1 / (2 * self.mred * n * a5 * omecc2)
            if self.rotdist:
                rotdistx1 = - xyz_prefac1 * spin1_h * spin1_e
                rotdisty1 = - xyz_prefac1 * spin1_h * spin1_q
                rotdistz1 = 0.5 * (2 * spin1_h * spin1_h - spin1_e * spin1_e - spin1_q * spin1_q)
            else:
                rotdistx1 = rotdisty1 = rotdistz1 = 0.0

            Tf1 = self.Tv1 / 9 * (self.a / self.R1) ** 8 * (self.m1 * (1 - self.Q1)) ** 2 / (self.mtot * self.m2)
            x1 = rotdistx1 - spin1_q / Tf1 * x_diss_fac
            y1 = rotdisty1 + spin1_e / Tf1 * y_diss_fac
            z1 = xyz_prefac1 * (rotdistz1
                                + 15 * self.m2 / a3 * fecc4 / omecc3)
            v1 = 9 / (Tf1 * omecc5) * (v_diss_fac - 11 * spin1_h / (18 * n) * fecc4)
            w1 = 1 / (Tf1 * omecc5) * (w_diss_fac - spin1_h / n * fecc5)

        if not self.noT2:
            xyz_prefac2 = self.m1 * self.A2 / (2 * self.mred * n * a5 * omecc2)
            if self.rotdist:
                rotdistx2 = - xyz_prefac2 * spin2_h * spin2_e
                rotdisty2 = - xyz_prefac2 * spin2_h * spin2_q
                rotdistz2 = 0.5 * (2 * spin2_h * spin2_h - spin2_e * spin2_e - spin2_q * spin2_q)
            else:
                rotdistx2 = rotdisty2 = rotdistz2 = 0.0
            Tf2 = self.Tv2 / 9 * (self.a / self.R2) ** 8 * (self.m2 * (1 - self.Q2)) ** 2 / (self.mtot * self.m1)
            x2 = rotdistx2 - spin2_q / Tf2 * x_diss_fac
            y2 = rotdisty2 + spin2_e / Tf2 * y_diss_fac
            z2 = xyz_prefac2 * (rotdistz2
                                + 15 * self.m1 / a3 * fecc4 / omecc3)
            v2 = 9 / (Tf2 * omecc5) * (v_diss_fac - 11 * spin2_h / (18 * n) * fecc4)
            w2 = 1 / (Tf2 * omecc5) * (w_diss_fac - spin2_h / n * fecc5)

        dydt = np.zeros_like(y)
        devec, dhvec, dspin1, dspin2 = np.split(dydt, 4)

        devec[...] = (z1 + z2 + zgr) * self.qhat - (y1 + y2) * self.hhat - (v1 + v2) * self.ehat
        devec[...] = self.e * devec

        dhvec[...] = (y1 + y2) * self.ehat - (x1 + x2) * self.qhat - (w1 + w2) * self.hhat
        dhvec[...] = self.h * dhvec

        dspin1[...] = self.mred * self.h / self.I1 * (-y1 * self.ehat + x1 * self.qhat + w1 * self.hhat)
        dspin2[...] = self.mred * self.h / self.I2 * (-y2 * self.ehat + x2 * self.qhat + w2 * self.hhat)

        return dydt

    def initialize_integrator(self, t0=None):
        """
        Initialize RK45 integrator with current state.

        Args:
            t0 (float, optional): Start time [time units]. Defaults to current time (0).
        """
        self.rtol, self.atol = 1e-8, 1e-7
        self.t = 0.0 if t0 is None else t0
        # BDF LSODA RK45 Radau
        self.integrator = RK45(self.derivative_dydt, self.t, self.y,
                               t_bound=1.3e10,
                               rtol=self.rtol,
                               atol=self.atol)

    def evolve(self, tfin, dt_out):
        """
        Evolve system to final time, outputting at fixed intervals.

        Args:
            tfin (float): Final integration time [time units].
            dt_out (float): Output timestep [time units].

        Returns:
            ndarray: Array of shape (N,7) with columns [t,a,e,spin1,spin2,obl1,obl2].
        """
        print("t0 = {:g} ".format(self.t))
        print("tfin = {:g}  ".format(tfin))
        print("dt_out = {:g}  ".format(dt_out))

        dt_out_next = dt_out
        dataout = []
        self.calculate_quantities()
        if self.t == 0.0: dataout.append((self.t, self.a, self.e, self.spin1, self.spin2, self.obl1, self.obl2))

        while self.t < tfin:
            self.integrator.step()

            if self.integrator.status == "failed":
                print("Integrator failed at t={:e}".format(self.t))
                break

            self.t = self.integrator.t
            self.y = self.integrator.y

            self.calculate_quantities()

            if self.a * (1 - self.e) < (self.R1 + self.R2):
                print("System has collided at t={:e}".format(self.t))
                break

            if self.t >= dt_out_next:
                dataout.append((self.t, self.a, self.e, self.spin1, self.spin2, self.obl1, self.obl2))
                dt_out_next = self.t + dt_out
            print("Progress: {:2.2%} ".format(self.t / tfin), end="\033[K\r")

        return np.array(dataout)

    def pseudosync_spin(e, n):
        """
        Compute pseudosynchronous spin frequency.

        Args:
            e (float): Eccentricity [dimensionless].
            n (float): Mean motion [rad/time units].

        Returns:
            float: Pseudosynchronous spin rate [rad/time units].
        """
        e2 = e * e
        fecc2 = ((e2 * 0.3125 + 5.625) * e2 + 7.5) * e2 + 1
        fecc5 = (e2 * 0.375 + 3) * e2 + 1
        omecc3h = (1 - e2) ** 1.5

        pseudos = fecc2 / (fecc5 * omecc3h) * n
        return pseudos

    def eccentricity_spin(e, n):
        """
        Compute eccentricity-driven equilibrium spin frequency.

        Args:
            e (float): Eccentricity [dimensionless].
            n (float): Mean motion [rad/time units].

        Returns:
            float: Equilibrium spin rate [rad/time units].
        """
        e2 = e * e
        fecc3 = ((e2 * 7.8125e-2 + 1.875) * e2 + 3.75) * e2 + 1
        fecc4 = (1.5 + e2 * 0.125) * e2 + 1
        omecc3h = (1 - e2) ** 1.5

        eccs = 18 / 11 * fecc3 / (fecc4 * omecc3h) * n
        return eccs

    def plot(self, outdata):
        """
        Plot evolution results.

        Args:
            outdata (ndarray): Output from evolve(), shape (N,7) [t,a,e,spin1,spin2,obl1,obl2].
        """
        f, ax = plt.subplots(4, 1, figsize=(7, 7), sharex=True, gridspec_kw=dict(hspace=0), tight_layout=True)
        plt.rc('lines', linewidth=2)
        plt.rc('text', usetex=True)
        t = outdata[:, 0] / self.P0
        a = outdata[:, 1]
        e = outdata[:, 2]
        n = (self.mtot / a ** 3) ** 0.5 / self.KU.Tscale

        spin1 = outdata[:, 3] / self.KU.Tscale
        spin2 = outdata[:, 4] / self.KU.Tscale
        obl1 = np.degrees(outdata[:, 5])
        obl2 = np.degrees(outdata[:, 6])
        nsync = EqTides.pseudosync_spin(e, n)
        necc = EqTides.eccentricity_spin(e, n)
        ax[0].plot(t, a)
        ax[0].set_ylabel("a [au]")

        ax[1].plot(t, e)
        ax[1].set_ylabel("ecc")

        ax[2].plot(t, spin1, alpha=0.66)
        ax[2].plot(t, spin2, alpha=0.66)
        ax[2].plot(t, n, color="black", ls="--", label="$n$")
        ax[2].plot(t, nsync, color="grey", ls=":", label=r"$n_{\\rm p\hbox{-}s}$")
        ax[2].plot(t, necc, color="olivedrab", ls=":", label="$n_{\\rm ecc}$")
        ax[2].set_ylabel(r"spin, n [2$\pi$/yr]")
        ax[2].set_yscale("log")
        ax[2].legend()

        ax[3].plot(t, obl1, alpha=0.66)
        ax[3].plot(t, obl2, alpha=0.66)
        ax[3].set_ylabel("obl [deg]")

        ax[3].set_xlabel("time / $P_0$")
        for axx in ax:
            axx.margins(x=0.01)
            axx.grid()

        plt.show()


def make_tsunami_ic(m1, m2, R1, R2, a0, e0, spin1, spin2):
    KU = tsunami.KeplerUtils()
    m = np.array((m1, m2), dtype=np.float64)
    R = np.array((R1, R2), dtype=np.float64)

    pos_vel1 = np.zeros(6, dtype=np.float64)
    pos_vel2 = KU.kepl_to_cart(pos_vel1, m1, m2, a0, e0, 0.0, 0.0, 0.0, 0.0)

    pos_vel_com = m1 * pos_vel1 + m2 * pos_vel2
    pos_vel_com /= m1 + m2

    pos_vel1 -= pos_vel_com
    pos_vel2 -= pos_vel_com

    p = np.array((pos_vel1[:3], pos_vel2[:3]), dtype=np.float64)
    v = np.array((pos_vel1[3:], pos_vel2[3:]), dtype=np.float64)
    pt = np.zeros_like(m, dtype=np.int64)
    spin = np.array((spin1, spin2))

    return m, R, p, v, spin, pt


if __name__ == "__main__":
    KU = tsunami.KeplerUtils()

    m1 = 1  # MSun
    m2 = 1  # MSun
    R1 = 1 * KU.RSun2au
    R2 = 1 * KU.RSun2au

    kap = 0.05
    k1 = kap
    k2 = kap

    tausec = 1e3
    tau1 = tausec / KU.yr2sec / KU.Tscale
    tau2 = tausec / KU.yr2sec / KU.Tscale

    rg = 0.28
    rg1 = rg
    rg2 = rg

    a0, e0 = 0.09, 0.3
    Ps1 = 100 / 365.25 / KU.Tscale
    spin1 = 2 * np.pi / Ps1

    Ps2 = 0.0015 / KU.Tscale
    spin2 = 2 * np.pi / Ps2

    ET = EqTides(m1, m2, R1, R2, k1, k2, tau1, tau2, rg1, rg2, rotdist=False, usespin=True)

    ET.setup_orbital_vectors(a0, e0)

    obl1 = np.radians(30)
    obl2 = 0.0

    ET.setup_spin_vectors(spin1, spin2, obl1, obl2)
    ET.setup_vectors()
    ET.initialize_integrator()

    save_txt = False

    tfin = 5e5 * ET.P0
    dt = 5e1 * ET.P0

    m, R, p, v, spin, pt = make_tsunami_ic(m1, m2, R1, R2, a0, e0, ET.spin1_vec, ET.spin2_vec)
    if save_txt:
        dat = np.hstack((p, v, m.reshape(-1, 1), R.reshape(-1, 1), pt.reshape(-1, 1), spin))
        np.savetxt("testspintides.dat", dat,
                   fmt="%1.16e %1.16e %1.16e  %1.16e %1.16e %1.16e  %1.16e  %1.16e  %4d  %1.16e %1.16e %1.16e")

    print("rg:", rg, "tausec", tausec, "k", kap, "ft", tfin)

    outdata = ET.evolve(tfin, dt)
    ET.plot(outdata)
    if save_txt:
        np.savetxt("testspintides_secular.txt", outdata)
