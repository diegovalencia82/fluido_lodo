
      subroutine viscous_bingham(mspace,ntype,epsilon2,mud)

c     mspace(10,i) = p-- pressure of particles                     [in]

      implicit none
      include 'param.inc'

      integer i,ntype(2)
      double precision mspace(25,nmax),epsilon2(ntype(1)),gammap
      double precision c,phi,nc,mud(ntype(1)),nume,dem,K,tauB,mu_inf,
     +     mu00

      double precision concentracion, p
      
      
      c = 10.
      phi = 30 * pi / 180.
      mu_inf = mu
      mu00 = 50.*mu_inf
      
      do i=1,ntype(1)
         mud(i) = mu * 50
      enddo

      ! Modelo de Binghman
c      do i = 1,ntype(1)!-ntype(2)
c         gammap = sqrt(0.5*epsilon2(i))
c         tauB = ( c + mspace(10,i)*tan(phi) ) * 0.5
c         K = mu00 / tauB
c         nume = K*mu_inf*gammap + mu00
c         dem = K*gammap + 1
c         mud(i) = nume / dem
c         write(*,*)i,mspace(12,i),mspace(10,i),K,mu_inf,gammap,mu00,
c     +        mud(i),tauB
c      enddo

c      nc = mu / 1000.
c      concentracion = 90 !%
c     p = (-1./100.)*concentracion + 1
      p = 2.!1./100000.
      nc = mu * p
      
      do i = 1,ntype(1)!-ntype(2)
         gammap = sqrt(0.5*epsilon2(i))
         mud(i) = nc * gammap**2
c         write(*,*)i,mspace(12,i),mspace(10,i),K,mu_inf,gammap,mu00,
c     +        mud(i),tauB
      enddo
      
      
      end


