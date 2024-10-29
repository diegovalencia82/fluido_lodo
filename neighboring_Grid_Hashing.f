

      subroutine neighboring_Grid_Hashing(rdomain,mspace,ntype,npairs,
     +     pairs,nfilas,mrij,mxij,mvij,mvxij)

c----------------------------------------------------------------------
c     Subroutine to determine the neighboring of each particle

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
c     mspace(14,i) = c-- sound velocity of particle                [in]
c     mspace(15,i) = s-- entropy of particles                      [in]
c     mspace(16,i) = e-- total energy of particles                 [in]
c     mspace(17 to 19,i) = dx  : dx = vx = dx/dt                   [in]
c     mspace(20 to 22,i) = dvx = dvx/dt, force per unit mass       [in]
c     mspace(23,i) = du        : du = du/dt                        [in]
c     mspace(24,i) = ds        : ds = ds/dt                        [in]
c     mspace(25,i) = drho      : drho = drh,o/dt                   [in]
c     mspace(26,i) = eta_c     : Coeficiente de Viscosidad No Lineal [in]      

c     t         : Temperature                                      [in]
c     tdsdt     : Production of viscous entropy t*ds/dt            [in]
c     av        : Monaghan average velocity                        [in]
c     rdomain   : smoothing length                                 [in]      

c     npairs    : maximum number of pairs interaction             [out]
c     nfilas    : number of interaction for each fluid particle   [out]
c     mrij      : matrix of ri-rj for all fluid particles for each interaction [out]
c     mxiij      : matrix of xi-xj for all fluid particles for each interaction [out]            
      
      implicit none
      include 'param.inc'

      integer ntype(2),sel,npairs
      integer pairs(npairs,ntype(1)),nfilas(ntype(1)),i,j,d
      double precision mspace(26,nmax)
      double precision rij,xij,yij,zij
      double precision mrij(npairs,ntype(1)),mxij(3,npairs,ntype(1))
      double precision vij,vxij,vyij,vzij,rdomain
      double precision mvij(npairs,ntype(1)),mvxij(3,npairs,ntype(1))


      double precision :: cell_size
      double precision dxmin, dxmax, dymin, dymax, dzmin, dzmax
      integer :: ncelda_x, ncelda_z, ndomain


c      mrij = 0.
c      mxij = 0.
c      mvij = 0.
c      mvxij = 0.
c      nfilas = 0

      dxmin = -0.25
      dxmax = 2.
c      ymin = xmin
c      ymax = 0.3
      dzmin = dxmin
      dzmax = 1.

      ndomain = 4
      cell_size = rdomain / (ndomain + 0.0)
      
      ncelda_x = int((dxmax - dxmin) / cell_size)
c      ncelda_y = int((dymax - dymin) / cell_size)
      ncelda_z = int((dzmax - dzmin) / cell_size)

c      write(*,*)'ncelda_x, ncelda_y, ncelda_z',ncelda_x, ncelda_y,
c     +     ncelda_z


      call search(rdomain,mspace,ntype,npairs,pairs,nfilas,mrij,mxij,
     +     mvij,mvxij,ncelda_x,ncelda_z,cell_size,ndomain+1)


      
      end

c     ---------------------------------------------------------------------

      subroutine search(rdomain,mspace,ntype,npairs,pairs,nfilas,mrij,
     +   mxij,mvij,mvxij,ncelda_x,ncelda_z,cell_size,ndomain)

      implicit none
      include 'param.inc'

      integer ntype(2),sel,npairs
      integer pairs(npairs,ntype(1)),nfilas(ntype(1)),i,j,d
      double precision mspace(26,nmax)
      double precision rij,xij,yij,zij
      double precision mrij(npairs,ntype(1)),mxij(3,npairs,ntype(1))
      double precision vij,vxij,vyij,vzij,rdomain
      double precision mvij(npairs,ntype(1)),mvxij(3,npairs,ntype(1))

      
!     Crear un arreglo para almacenar las partículas en cada celda

      double precision dxmin, dxmax, dymin, dymax, dzmin, dzmax
      integer :: ncelda_x, ncelda_z, nceldas, k
      parameter (nceldas = 20)
      integer, dimension(ncelda_x, ncelda_z, nceldas) :: celdas
      integer, dimension(ncelda_x, ncelda_z) :: sceldas
      double precision :: cell_size, xmin_inv, ymin_inv, zmin_inv
      integer :: cell_x, cell_y, cell_z, ndomain

      integer :: dx, dy, dz, neighbor_x, neighbor_y, neighbor_z

      celdas = -1
      sceldas = 0

      dxmin = -0.25
      dxmax = 2.
c      ymin = xmin
c      ymax = 0.3
      dzmin = dxmin
      dzmax = 1.

      
c      write(*,*)'0000'

      xmin_inv = 1.0 / cell_size
      ymin_inv = 1.0 / cell_size
      zmin_inv = 1.0 / cell_size
      
      do i = 1, ntype(1)
         ! Calcular las celdas correspondientes para la partícula i
         cell_x = int((mspace(2,i) - dxmin) * xmin_inv) + 1
         cell_z = int((mspace(4,i) - dzmin) * zmin_inv) + 1
         if (cell_x > 0 .and. cell_x <= ncelda_x .and. 
     +        cell_z > 0 .and. cell_z <= ncelda_z) then
              ! Añadir la partícula a la celda correspondiente
            sceldas(cell_x, cell_z) = sceldas(cell_x, cell_z) + 1
            celdas(cell_x, cell_z, sceldas(cell_x, cell_z)) =
     +           int(mspace(1,i))
            if(sceldas(cell_x, cell_z).gt.20)then
               write(*,*)'There are cells with more than 20 particle'
               write(*,*)'Number ',sceldas(cell_x, cell_z)
               write(*,*)'Stop, neighboring_Grid_Hashing.f'
               stop
            endif
         endif
      enddo

c      write(*,*)'1111'
      
      mrij = 0.
      mxij = 0.
      mvij = 0.
      mvxij = 0.
      nfilas = 0




      do i = 1, ntype(1) - ntype(2)
         ! Calcular la celda de la partícula
         cell_x = int((mspace(2,i) - dxmin) / cell_size) + 1
         cell_z = int((mspace(4,i) - dzmin) / cell_size) + 1
         ! Buscar en las celdas vecinas
         do dx = -ndomain, ndomain
            do dz = -ndomain, ndomain
               neighbor_x = cell_x + dx
               neighbor_z = cell_z + dz
               
                  ! Verificar que la celda está dentro de los límites
               if (neighbor_x >= 1 .and. neighbor_x <= ncelda_x.and.
     +              neighbor_z >= 1 .and. neighbor_z <= ncelda_z)then
                     
!     Buscar partículas en la celda vecina
                  do k = 1,sceldas(neighbor_x,neighbor_z)
                     if (celdas(neighbor_x, neighbor_z, k)
     +                    .ne. -1) then
!     Calcular la distancia entre las partículas i y j
                        j = celdas(neighbor_x, neighbor_z, k)
                        if(i.ne.j)then
                           call radioij(mspace(2,i), mspace(3,i),
     +                          mspace(4,i),mspace(2,j), mspace(3,j),
     +                          mspace(4,j), rij, xij, yij, zij)
                           call velij(mspace(5,i),mspace(6,i),
     +                          mspace(7,i),mspace(5,j),mspace(6,j),
     +                          mspace(7,j),vij,vxij,vyij,vzij)
                           if (rij.lt.rdomain) then
!     Es vecino, almacenar el resultado como antes
                              nfilas(i) = nfilas(i) + 1
                              if(nfilas(i).gt.npairs)then
                                 write(*,*)
     +                                'neighboring_search.f STOP, i=',i
     +                                ,'j=',j,'nfilas = ',nfilas(i)
                                 stop
                              endif
                              pairs(nfilas(i),i) = int(mspace(1,j))
                              mrij(nfilas(i),i) = rij
                              mxij(1,nfilas(i),i)  = xij
                              mxij(2,nfilas(i),i)  = yij
                              mxij(3,nfilas(i),i)  = zij
                              mvij(nfilas(i),i) = vij
                              mvxij(1,nfilas(i),i)  = vxij
                              mvxij(2,nfilas(i),i)  = vyij
                              mvxij(3,nfilas(i),i)  = vzij
                           endif
                        endif
                     endif
                  enddo
               endif
            enddo
         enddo
c      write(*,*)i,nfilas(i)
      enddo
      
      end
      
c     ----------------------------------------------------------------------

      subroutine radioij(xi,yi,zi,xj,yj,zj,rij,xij,yij,zij)

      implicit none
      include 'param.inc'

      double precision xi,yi,zi,xj,yj,zj,xij,yij,zij,rij

      
      if(dim.eq.2)then
         xij = xi-xj
         yij = 0.0d0
         zij = zi-zj
         rij = sqrt(xij*xij + zij*zij)
      endif

      if(dim.eq.3)then
         xij = xi-xj
         yij = yi-yj
         zij = zi-zj
         rij = sqrt(xij*xij + yij*yij + zij*zij)
      endif
      
      end 

c     ----------------------------------------------------------------------

      subroutine velij(vxi,vyi,vzi,vxj,vyj,vzj,vij,vxij,vyij,vzij)

      implicit none
      include 'param.inc'

      double precision vxi,vyi,vzi,vxj,vyj,vzj,vxij,vyij,vzij,vij

      
      if(dim.eq.2)then
         vxij = vxi-vxj
         vyij = 0.0
         vzij = vzi-vzj
         vij = sqrt(vxij*vxij + vzij*vzij)
      endif

      if(dim.eq.3)then
         vxij = vxi-vxj
         vyij = vyi-vyj
         vzij = vzi-vzj
         vij = sqrt(vxij*vxij + vyij*vyij + vzij*vzij)
      endif
      
      end       
      
