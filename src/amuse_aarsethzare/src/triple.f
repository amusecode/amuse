*	triple.f
*
*********************************************************************************
*
*	Calculates evolution of 3-body system.
*	--------------------------------------
*	Inner orbit at arbitrary inclination to outer orbit.
*	---------------------------------------------------
*	Self-consistent treatment: total angular momentum conserved.
*	-----------------------------------------------------------
*
*********************************************************************************

      implicit real*8 (a-h,m,o-z)
      real*8 mass(3),pos(3,3),vel(3,3), Etot, DEtot
c      open(98,file='exampletriple.in')
      real*8 dmdt, r(3), v(3), elements(6), dm, ttend, tend, tstart
      integer NR, NRtot
      integer AarsethZare
      real*8 DegenareteCoreMass, mc
      open(98,file='triple.in')

c     index i over particle number 1 to 3
c     pos(3,i): z-coordinate of particle "i"
      read (98,*) tstart
      do i=1,3
         read(98,*) id, mass(i), pos(1,i), pos(2,i), pos(3,i), 
     &                           vel(1,i), vel(2,i), vel(3,i)
      enddo
      read (98,*) dmdt, ttend

c      write (*,*) "Input number of orbits to integrate:"
c      read (*,*) norb

      write (*,*) "Initial Configurations"
c      call print_binary_combinations(time, mass, pos, vel)
      call print_inner_binary_parameters(tstart, mass,pos,vel,
     &                                   Etot,ai,ao)

c
c      time_acceleration_factor = 1.e+5
c      time_acceleration_factor = 10
c     mass loss rate in Msun/year
c     dmdt = dm/(tend-tstart)

      m123 = mass(1)+mass(2)+mass(3)
      call get_outer_elements(mass, pos, vel, elements)
      Porb = orbital_period(elements(1), m123)
      Tunit = orbital_time_unit(elements(1), m123)
c      time_acceleration_factor = max(1., 0.1*Tunit)
c      time_acceleration_factor = min(50.0, 5000./Porb)
c      time_acceleration_factor = 10.0
c      time_acceleration_factor = 50.0
      time_acceleration_factor = 1.0
      year_per_step = time_acceleration_factor * Porb/Tunit
      write (*,*) "Time step constant ", Porb, Tunit, Porb/Tunit, 
     1                                   time_acceleration_factor

c      mcore = min(0.50, 0.17 + 0.0657 * dlog10(365.25*Porb))
      mcore = DegenerateCoreMass(mass, pos, vel)

      tend = ttend / year_per_step
      write (*,*) "Porbit", Porb, Tunit, year_per_step, tstart, 
     &                      tend, dmdt, mcore
      ttme = tstart
c     Output time step in years
      tstep = 1.e+4
      dtev = 0
      time = 0
      ttime = 0
      dt = 0
      cont_run = 0
      mass_transferred = 0
      NRtot = 0
      DEtot = 0
      do while (ttime<tend)
         time = 0
         cont_run = AarsethZare(time,mass,pos,vel,Etot,NR)
         DEtot = DEtot + Etot
         NRtot = NRtot + NR
c         call print_inner_binary_parameters(ttme, mass, pos, vel)
         ttime = ttime + time
         ttme = ttme + time * year_per_step
c         write (*,*) "Time=", time, ttime, ttme, NRtot
         dt = dt + time * year_per_step
         if (ttme.gt.dtev) then
            dtev  = ttme + tstep
            write (*,*) "Time = ", ttime, ttime/Tunit, ttme, 
     1                             dt, dtev, "[year], NR=", NRtot, DEtot
            call print_inner_binary_parameters(ttme, mass,pos,vel,
     1                                         Etot, ai, ao)
            call ApplyMassTransfer(dt, dmdt, mass, pos, vel, mcore)
            if (ai.lt.0.and.ao.lt.0) then
               call terminate_run(ttme, mass, pos,vel,Etot)
            endif
            mc = DegenerateCoreMass(mass, pos, vel)
            if (mc.gt.0) then
               mcore = mc
            endif
            call write_restart_file(ttme, mass, pos, vel, dmdt, ttend)
            dt = 0
            if (mass(1).lt.mcore) then
               dmdt = 0
               tend = min(tend, (ttme+1.e+7)/ year_per_step)
               write (*,*) "Reached end of mass transfer:", dmdt, tend
               call print_binary_combinations(ttme, mass, pos,vel,Etot)
               write (*,*) "Continue until termination time is reached."
c               STOP
            endif
         endif
         if (cont_run.eq.0)  then
            write (*,*) "Triple resolved"
            call print_binary_combinations(ttme, mass, pos, vel, Etot)
            STOP
         endif

      enddo
c      write (*,*) "Final Configurations"
c      call print_binary_combinations(time, mass, pos, vel)
      end

      subroutine terminate_run(ttme, mass, pos,vel,Etot)
      implicit real*8 (a-h,m,o-z)
      real*8 ttme, mass(3), pos(3,3),vel(3,3),Etot
      write (*,*) "Triple unrsolved: terminate run"
      call print_binary_combinations(ttme, mass, pos,vel,Etot)
      STOP
      end

      subroutine write_restart_file(time, mass, pos, vel, dm, tend)
      implicit real*8 (a-h,m,o-z)
      real*8 r(3), v(3), mass(3), pos(3,3), vel(3,3), elements(6)

      open(99,file='triple.rst')
      write(99,*) time
      do i=1,3
         write(99,*) i, mass(i), pos(1,i), pos(2,i), pos(3,i), 
     &                           vel(1,i), vel(2,i), vel(3,i)
      enddo
      write (99,*) dm, tend
      close(99)
      end

      subroutine transfer_mass(ain, aout, md, ma, dmd, dma) 
      implicit real*8 (a-h,m,o-z)
      real*8 ain, aout, md, ma, dmd, dma
      aout = ain * (md*ma/((md-dm)*(ma+dm)))**2

      mt = md + ma
      mdf = md - dmd
      maf = ma + dma
      mtf = mdf + maf
      alpha = (mt-mtf)/dmd
      aout = ain * (((mdf/md) * (maf/ma))**(1/(1-alpha)) )**(-2) 
     &     * (mt/mtf)
      md = mdf
      ma = maf
      end

      subroutine get_outer_elements(mass, pos, vel, elements)
      implicit real*8 (a-h,m,o-z)
      real*8 r(3), v(3), mass(3), pos(3,3), vel(3,3), elements(6)
      m3 = mass(3)
      m12 = mass(1)+mass(2) 
      m123 = m12 + m3
      q123 = m123/(1.0*m12)
      do i=1,3
         r(i)= q123*pos(i,3)
         v(i)= q123*vel(i,3)
      enddo
      call transform1new(m3,m12,r,v,elements)
      end

      subroutine get_inner_elements(mass, pos, vel, elements)
      implicit real*8 (a-h,m,o-z)
      real*8 r(3), v(3), mass(3), pos(3,3), vel(3,3), elements(6)
      do i=1,3
         r(i)=pos(i,2)-pos(i,1)
         v(i)=vel(i,2)-vel(i,1)
      enddo
      call transform1new(mass(1),mass(2),r,v,elements)
      end
     
      subroutine print_inner_binary_parameters(time,mass,pos,vel,
     1                                         Etot,ai,ao)
      implicit real*8 (a-h,m,o-z)
      real*8 time, r(3), v(3), mass(3), pos(3,3), vel(3,3), elements(6)
      real*8 ai, ao
      aout = ain * (md*ma/((md-dm)*(ma+dm)))**2

      write(*,*) ' TIME= ', time/1.e+6, ' DE/E=', Etot

      do i=1,3
         r(i)=pos(i,2)-pos(i,1)
         v(i)=vel(i,2)-vel(i,1)
      enddo
      call transform1new(mass(1),mass(2),r,v,elements)
      write(*,*) 'Inner (1,2) time=',1.e-6 * time,
     $        ' [Myear] a=',elements(1),
     $        'e=',elements(2),	mass(1), mass(2),
     $        elements(3),elements(4),elements(5),elements(6)

      ai = elements(1)
      m3 = mass(3)
      m12 = mass(1)+mass(2) 
      m123 = m12 + m3
      q123 = m123/(1.0*m12)
      do i=1,3
         r(i)= q123*pos(i,3)
         v(i)= q123*vel(i,3)
      enddo
      call transform1new(m3, m12,r,v,elements)
      write(*,*) 'Outer (12,3) time=', 1.e-6*time,
     $     ' [Myear] a=',elements(1),' e=',elements(2), m12, m3,
     $     elements(3),elements(4),elements(5),elements(6)
      ao = elements(1)
      end

c     orbital period in years      
      function orbital_period(a, mtot)
      implicit real*8 (a-h,m,o-z)
      real*8 a, mtot
      orbital_period = 365.25*sqrt((0.004652*a)**3/mtot)
      end

      function orbital_time_unit(a, mtot)
      implicit real*8 (a-h,m,o-z)
      real*8 a, mtot, pi
      pi=4.0d0*atan(1.0d0)
      orbital_time_unit = 2*pi*sqrt(a**3/mtot)
      end

      subroutine center_of_mass(mass, pos, vel, rcom, vcom)
      implicit real*8 (a-h,m,o-z)
      real*8 time, dmdt, mass(3),pos(3,3),vel(3,3)
      real*8 rcom(3), vcom(3)

      do i=1,3
         rcom(i) = mass(1)*pos(i,1) + mass(2)*pos(i,2)
         vcom(i) = mass(1)*vel(i,1) + mass(2)*vel(i,2)
      enddo
      mtot = mass(1)+mass(2)
      do i=1,3
         rcom(i) = rcom(i)/mtot
         vcom(i) = vcom(i)/mtot
      enddo
      end

      function RocheRadius(ain, md, ma)
      implicit real*8 (a-h,m,o-z)
      real*8 ain, md, ma, q, qcrt, qcrt2, RR
      q = md/ma;
      qcrt = q**(1./3.)
      qcrt2 = qcrt*qcrt
      RR = ain*0.49*qcrt2/(0.6*qcrt2 + log(1 + qcrt))
      RocheRadius = RR
      end

      subroutine circularize(a, e)
      implicit real*8 (a-h,m,o-z)
      real*8 a, e
      a = a*(1-e*e)
      e = 0
      end

      function FitToDegenerateCoreMass(Porb)
      implicit real*8 (a-h,m,o-z)
      real*8 Porb, a, b, c, mmax
      mmax = 0.4
      a = 5.00
      b = 1.0e+5
      c = 0.110
      FitToDegenerateCoreMass = min(mmax, (Porb/b)**(1./a) + c)
      end

      function DegenerateCoreMass(mass, pos, vel)
      implicit real*8 (a-h,m,o-z)
      real*8 mass(2), pos(3,3), vel(3), elements(6)
      real*8 a, mtot, Porb 
      call get_inner_elements(mass, pos, vel, elements)
      a = elements(1)
      if (a.gt.0) then
         mtot = mass(1)+mass(2)
         Porb = orbital_period(a, mtot)
         DegenerateCoreMass = FitToDegenerateCoreMass(Porb)
      else 
         DegenerateCoreMass = -1
      endif
      write (*,*) "Degenerate Core Mass=", 
     &            Porb, elements(1), mass(1), DegenerateCoreMass

      end

      subroutine ApplyMassTransfer(dt, dmdt, mass, pos, vel, mcore)
      implicit real*8 (a-h,m,o-z)
      real*8 time, dt, dmdt, mass(3),pos(3,3),vel(3,3), mcore
      real*8 rcom(3), vcom(3), dmd, dma
      real*8 md, ma, r(3), v(3), elements(6)
      real*8 AccretionEfficiency
c      write (*,*) "apply mass transfer"
c      write(*,*) dt, dmdt, mass, pos, vel, dmdt, dm
c      write(*,*) dt, dmdt, mass, dmdt, dm

      call center_of_mass(mass, pos, vel, rcom, vcom)
c      write(*,*) "CoM:", rcom, vcom
      
c      AccretionEfficiency = 1./2.0
      AccretionEfficiency = 1.
      dmd =  dmdt * dt
      dMEddington = dt*1.5e-8
      dma = min(AccretionEfficiency*dmd, dMEddington)
      md = mass(1)
      ma = mass(2)
c      mcore = 0.01
      if (md-dmd.lt.mcore) then
         fm = max(0., (md-mcore)/dmd)
         dmd = fm * dmd
         dma = fm * dma
         write (*,*) 'Minimum mass reached for mass transfer', 
     &                fm, dmd, dma, mass
      endif
      do i=1,3
         r(i)=pos(i,2)-pos(i,1)
         v(i)=vel(i,2)-vel(i,1)
      enddo
c      write(*,*) "Pre Transformation: r, v=", r, v
      call transform1new(md,ma,r,v,elements)
c      write(*,*) "Pre Elements: ", elements
c     Apply mass transfer
      ain = elements(1)
      RRoche = RocheRadius(ain, md, ma)
      if (RRoche.lt.150) then
c         write(*,*) "PreMassTransfer: ", ain, md, ma, dmd, dma
c         write(*,*) "Pre mass transfer: r, v=", r, v
         call transfer_mass(ain, aout, md, ma, dmd, dma)
c         write(*,*) "PostMassTransfer: ", ain, aout, md, ma
         elements(1) = aout
         mass(1) = md
         mass(2) = ma
         mtot = mass(1)+mass(2)
c         write (*,*) "Pre circularize:", elements(1), elements(2)
c         call circularize(elements(1), elements(2))
c         write (*,*) "Post circularize:", elements(1), elements(2)
         call transform2(mtot,elements,r,v)
c        write(*,*) "Post Transformation: r, v=", r, v
c        write(*,*) "pre, Massposvel:", mass, pos, vel
         do i=1,3
            pos(i, 1) = rcom(i) - r(i)* mass(2)/mtot
            vel(i, 1) = vcom(i) - v(i)* mass(2)/mtot

            pos(i, 2) = rcom(i) + r(i)* mass(1)/mtot
            vel(i, 2) = vcom(i) + v(i)* mass(1)/mtot
         enddo
c        write(*,*) "post, Massposvel:", mass, pos, vel
      endif
      end

      subroutine print_binary_combinations(time, mass, pos, vel, Etot)
      implicit real*8 (a-h,m,o-z)
      real*8 time, mass(3),pos(3,3),vel(3,3), Etot
      real*8 m1, m2, r(3), v(3), elements(6)

      write(*,*) ' TIME= ', time/1.e+6, ' DE/E=', Etot
      
c      write (*,*) "Print binary parameters"
c     Inserted from RM (by SPZ)
c     try (1,2)
      m12 = mass(1)+mass(2)
      do i=1,3
         r(i)=pos(i,2)-pos(i,1)
         v(i)=vel(i,2)-vel(i,1)
      enddo
      call transform1new(mass(1),mass(2),r,v,elements)
      if (elements(1) .gt.0.and.elements(2)<1) then
         write(*,*) 'Inner (1,2) time=', 1.e-6*time,
     $        ' [Myear]  a=',elements(1),
     $        'e=',elements(2),	mass(1), mass(2),
     $        elements(3),elements(4),elements(5),elements(6)
      endif
c     try (12,3)
      m3 = mass(3)
      m123 = m12 + m3
      q123 = m123/(1.0*m12)
      do i=1,3
         r(i)= q123*pos(i,3)
         v(i)= q123*vel(i,3)
      enddo
      call transform1new(m3, m12,r,v,elements)
      if (elements(1) .gt.0.and.elements(2)<1) then
         write(*,*) 'Outer (12,3) time=', 1.e-6*time,
     $        ' [Myear]  a=',elements(1),'e=',elements(2), m12, m3,
     $        elements(3),elements(4),elements(5),elements(6)
         write (*,*), m3, m12, r, v
      endif

c     try (1,3)
      do i=1,3
         r(i)=pos(i,3)-pos(i,1)
         v(i)=vel(i,3)-vel(i,1)
      enddo
      call transform1new(mass(1), mass(3), r, v, elements)
      if (elements(1).gt.0.and.elements(2)<1) then
         write(*,*) 'Binary (1,3) time=',time*1.e-6,' [Myear]  a=',
     $        elements(1), 'e=',elements(2), mass(1), mass(3), 
     $        elements(3),elements(4),elements(5),elements(6)
      endif
c     try (2,3)
      do i=1,3
         r(i)=pos(i,2)-pos(i,3)
         v(i)=vel(i,2)-vel(i,3)
      enddo
      call transform1new(mass(2), mass(3), r, v, elements)
      if (elements(1) .gt.0.and.elements(2)<1) then
         write(*,*) 'Binary (2,3) time=',time*1.e-6,' [Myear]  a=',
     $        elements(1),'e=',elements(2),	mass(2), mass(3), 
     $        elements(3),elements(4),elements(5),elements(6)
      endif

      do is=1,3
         write (*,*) "PosVel: of star(", is, ")", time/1.e+6, 
     &        mass(is), pos(1, is), pos(2, is), pos(3, is), 
     &        vel(1, is), vel(2, is), vel(3, is)
      enddo

      end
c      include 'AarsethZare.f'
      include 'Aarseth.f'
      include 'transform.f'

