
      program SPH

c------------------------------------------------------------------
c     This is a three dimensional SPH code. the followings are the 
c     basic parameters needed in this codeor calculated by this code

c     ntotal-- total particle number ues                             [in]
      
c     dt--- Time step used in the time integration                   [in]
c     time -- time of the snapshot
c     itimestep -- number of the time step.
c     mspace(1,i) = id -- label of each particle                   [out]      
c     mspace(2 to 4,i) = x-- coordinates of particles              [in/out]
c     mspace(5 to 7,i) = vx-- velocities of particles              [in/out]
c     mspace(8,i) = mass-- mass of particles                       [in]
c     mspace(9,i) = rho-- dnesities of particles                   [in/out]
c     mspace(10,i) = p-- pressure of particles                     [in/out]
c     mspace(11,i) = u-- internal energy of particles              [in/out]
c     mspace(12,i) = itype-- types of particles                    [in]
c     mspace(13,i) = hsml-- smoothing lengths of particles         [in/out]
c     mspace(14,i) = c-- sound velocity of particle                [out]
c     mspace(15,i) = s-- entropy of particles                      [out]
c     mspace(16,i) = e-- total energy of particles                 [out]      
c     mspace(17 to 19,i) = dx  : dx = vx = dx/dt                   [out]
c     mspace(20 to 22,i) = dvx = dvx/dt, force per unit mass       [out]
c     mspace(23,i) = du        : du = du/dt                        [out]
c     mspace(24,i) = ds        : ds = ds/dt                        [out]
c     mspace(25,i) = drho      : drho = drh,o/dt                   [out]

      
      implicit none
      include 'param.inc'
      
      integer ntotal,nfluid,nvirt,ntype(2) !, itype(nmax),id(nmax)
c      double precision x(dim, nmax), vx(dim, nmax), mass(nmax)
c     +     ,rho(nmax), p(nmax), u(nmax), c(nmax), s(nmax),
c     +     e(nmax), hsml(nmax), dt
      double precision dt,mspace(25,nmax)

      real :: start_time, end_time, total_time

! Obtener el tiempo inicial
      call CPU_TIME(start_time)

      dt = dt0
      
      call input(mspace,ntotal,nfluid,nvirt,0)
c     ntotal = nfluid + nvirt
c     nfluid = ntotal - nvirt

      if(ntotal.ne.nmax)then
         write(*,*)'Different ntotal and nmax, 
     +change nmax in param file'
         stop
      endif
      
      ntype(1) = nmax!nfluid
      ntype(2) = nvirt
      write(*,*)' **************************************************'
      write(*,*)'        The  maximal time steps = ', maxtimestep
      write(*,*)'        Time steps integration  = ', dt
      write(*,*)' **************************************************' 

      call time_integration(mspace,ntotal,ntype,dt)


! Obtener el tiempo final
      call CPU_TIME(end_time)
! Calcular el tiempo de CPU consumido
      total_time = end_time - start_time
      open(83,file='time_CPU.dat')
      write(83,*)'START TIME=', start_time
      write(83,*)'END TIME=', end_time
      write(83,*)'cpu_time = end_time - start_time=', total_time
      close(83)
      
      end      
