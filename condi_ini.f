
      program condini
      implicit none

      integer maxn,dim

c     itype = 1: Fluid particles
c     itype = -1: Walls particles or virtual particles 
c     nt1 = 1 : Fluid particles
c     nt2 = -1: Walls particles or virtual particles 
      
      parameter(maxn=100000,dim=2)
      integer itype(maxn),itypev(maxn), ntotal
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn),
     &     p(maxn), u(maxn), hsml(maxn), rho(maxn) !, PI
      double precision xv(dim, maxn), vxv(dim, maxn), massv(maxn),
     &     pv(maxn), uv(maxn), hsmlv(maxn), rhov(maxn) !, PI
      integer i, d, im

      real x1(maxn),y1(maxn),xmin,xmax,ymin,ymax,r,dr,xt,yt
      double precision xa(dim, maxn),volum,PI,rhoa,masa,dmasa,drr,dv

      integer j,k,m,n,nvirt,mp,np,nt,nt1,nt2
      integer itimestep
      double precision xx1,yy1,dx,dy,r0,bc(dim,4),vxi(dim),xybox(5)
      double precision ht,gamma,rho0,b,lengsx,lengsy,delt,time,mu,beta
      double precision C,vmax,Dt,g 

      open(10,file="condiciones_iniciales.txt")

! Box position
      xybox(1) = 0.0d0          ! xybox(1): x Initial Position of the box
      xybox(2) = 0.0d0          ! xybox(2): y Initial Position of the box
      xybox(3) = 0.61d0          ! xybox(3): High of the box
      xybox(4) = 1.61d0          ! xybox(4): lower side of the box
      xybox(5) =0.0          ! xybox(5): angle position of the box
!-------------------end box position
      
      PI=4.*atan(1.0)
      g = 9.82
      rho0=1000.
      mu = 1.002e-3
      beta = 100.

      delt = 0.008
      
      lengsx = 600.
      lengsy = 300.
      xx1 = lengsx*1.0e-3
      yy1 = lengsy*1.0e-3

c      m=nint(yy1/delt) + 1
      
c      delt = yy1/(2.*m)

      bc(1,1)=0.0 + delt ! Initial position "x" of the box fluid (bottom left)
      bc(2,1)=0.0 + delt ! Initial position "y" of the box fluid (bottom left)
      bc(1,2)=xx1 + delt ! Initial position "x" of the box fluid (bottom rigth)
      bc(2,2)=0.0 + delt ! Initial position "y" of the box fluid (bottom rigth)
      bc(1,3)=0.0 + delt ! Initial position "x" of the box fluid (top left)
      bc(2,3)=yy1 + delt ! Initial position "y" of the box fluid (top left)
      bc(1,4)=xx1 + delt ! Initial position "x" of the box fluid (top rigth)
      bc(2,4)=yy1 + delt ! Initial position "y" of the box fluid (top rigth)
      
      ht=yy1
      vxi(1)=0.0
      vxi(2)=0.0
      nt=0
      write(*,*)'1 box, nt=',nt,'box=',bc
      
      call boxparticles01(dim,m,maxn,rho0,bc,ht,vxi,dx,nt,x,vx,mass,
     +     rho,p,u,itype,hsml,delt,xybox,mu,beta)
      
      nt1 = nt
      do i=1,nt
         x1(i)=x(1,i)
         y1(i)=x(2,i)
      enddo

      write(*,*)' *******************************************'
      write(10,*)' *******************************************'
      write(*,*)' Initial particle configuration generated   '
      write(10,*)' Initial particle configuration generated   '
      write(*,*)' Total number of particles   ', nt
      write(10,*)' Total number of particles   ', nt
      write(*,*)'                      hsml = ',hsml(1),delt
      write(10,*)'                      hsml = ',hsml(1)
      r0 = y1(200)-y1(199)
      write(*,*)'The cutoff distance for external forces r0 = ',r0
      write(10,*)'The cutoff distance for external forces r0 = ',r0
      write(*,*)'r0 = ',(y1(2)-y1(1))/2.
      write(10,*)'r0 = ',(y1(2)-y1(1))/2.
      write(*,*)' *******************************************'
      write(10,*)' *******************************************'      
      
      call pgbeg(0,'pgp_condini.ps/cps',1,1)
      call pgscf(2)
      call pgsch(1.0)
      xmin = -0.1!-10*dx
      xmax = 1.62!7*xx1+10*dx
      ymin = -0.10*dx
      ymax = 0.65!y1+10*dx
      call pgenv(xmin,xmax,ymin,ymax,1,2)
      call pglab('x','y','')
c      call pgenv(-0.1,3.5,-0.3,1.5,1,0)
      call pgsci(2)
      call pgpt(nt,x1,y1,1)
c      call pgline(nt,x1,y1)
      call pgsci(1)
      call pgtext(0.4,0.4,'T=0.00')
      
c      call virtpartfallwater(nt,maxn,dim,dx,hsml(1),mass(1),
c     +     itypev,hsmlv,massv,xv,vxv,rhov,uv,pv,nvirt,ht,rho0)

      call virtbox(beta,nt,maxn,dim,dx,hsml(1),mass(1),
     +     itypev,hsmlv,massv,xv,vxv,rhov,uv,pv,nvirt,ht,rho0,xybox,0)
      
      nt2 = nvirt
     
      call pgend

      
      write(*,*)'=======================Open files ini'
      write(10,*)'=======================Open files ini'
      write(*,*)'nfluid =',nt1
      write(10,*)'nfluid =',nt1
      write(*,*)'nvirt =',nt2
      write(10,*)'nvirt =',nt2
      write(*,*)'ntotal =',nt1+nt2
      write(10,*)'ntotal =',nt1+nt2

      open(1,file="snapshot_000")

      itimestep = 0
      time = 0.0
      write(1,*)itimestep,time,nt1+nt2, nt1, nt2
      do i = 1, nt1
         write(1,*) i, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim),
     +        mass(i), rho(i), p(i), u(i), itype(i), hsml(i)
      enddo
      do i = nt1+1, nt1+nt2
         write(1,*) i, (xv(d, i),d = 1, dim), (vxv(d, i),d = 1, dim),
     +        massv(i), rhov(i), pv(i), uv(i), itypev(i), hsmlv(i)
      enddo
      
      close(1)

! Calculos of Cournt condition
      C = 0.2
      vmax = sqrt(2*g*yy1)
      Dt = C * yy1 / vmax
      write(*,*)'Time step of Couran conditio = ',Dt
      write(10,*)'Time step of Couran conditio = ',Dt

      
      close(10)
      end

c==============================

      subroutine virtpartfallwater(ntotal,maxn,dim,dx,hsml1,dmasa,
     +     itype,hsml,mass,x,vx,rho,u,p,nvirt,ht,rho0)
      implicit none

      integer itimestep, ntotal, nvirt, maxn, dim
c      parameter (maxn=1000,dim=3)
      integer itype(maxn)
      double precision hsml(maxn),mass(maxn),x(dim,maxn),vx(dim, maxn),
     &     rho(maxn), u(maxn), p(maxn)
      integer i, j, d, im, mp
      double precision leng, dx,v_inf,hsml1,dl,dmasa,xy(2),vxy(2),ang,PI

      double precision gamma,beta,c,b,ht,rho0
      real x1(maxn),y1(maxn),xmin,xmax,ymin,ymax,lengs

      PI=4.*atan(1.0)

      gamma=7.
      beta = 5.0!0.003
      c = beta*sqrt(9.8*ht)
      b = rho0*c*c/gamma

      v_inf = 1.e-3
      nvirt = 0

      lengs = 600.

      ! Left side
      leng = 1000*1.0e-3 
      ang=90* (PI/180.)
      xy(1)=0.0
      xy(2)=0.0
      vxy(1)=0.0
      vxy(2)=0.0!v_inf
      call walls(dx,leng,xy,vxy,ang,maxn,ntotal,dim,x,vx,nvirt)
      write(*,*)'nvirt=',nvirt
      
      !lower side
      leng = 1500*1.0e-3-dx/2.
      ang=0.0
      xy(1)=dx/2.
      xy(2)=0.0
      vxy(1)=0.0!-v_inf
      vxy(2)=0.0
      call walls(dx,leng,xy,vxy,ang,maxn,ntotal,dim,x,vx,nvirt)
      write(*,*)'nvirt=',nvirt

      !Rigth side
      leng = 1000*1.0e-3-dx/2.
      ang=90* (PI/180.)
      xy(1)=1500*1.0e-3
      xy(2)=dx/2.
      vxy(1)=0.0
      vxy(2)=0.0!-v_inf
      call walls(dx,leng,xy,vxy,ang,maxn,ntotal,dim,x,vx,nvirt)
      write(*,*)'nvirt=',nvirt

      !Upper side
c      leng = 100*1.0e-3-dx/2.
c      ang=0.0
c      xy(1)=1000*1.0e-3
c      xy(2)=dx/2.
c      vxy(1)=0.0
c      vxy(2)=0.0
c      call walls(dx,leng,xy,vxy,45.d0,maxn,ntotal,dim,x,vx,nvirt)
c      write(*,*)'nvirt=',nvirt

      !Rigth side
      leng = 400*1.0e-3-dx/2.
      ang=90* (PI/180.)
      xy(1)=200*1.0e-3
      xy(2)=200*1.0e-3+dx/2.
      vxy(1)=0.0
      vxy(2)=0.0!-v_inf
c      call walls(dx,leng,xy,vxy,ang,maxn,ntotal,dim,x,vx,nvirt)
      write(*,*)'nvirt=',nvirt

            !Upper side
      leng = 200*1.0e-3-2*dx/2.
      ang=0.0
      xy(1)=dx/2.
      xy(2)=600*1.0e-3
      vxy(1)=0.0!v_inf
      vxy(2)=0.0
c      call walls(dx,leng,xy,vxy,ang,maxn,ntotal,dim,x,vx,nvirt)
      write(*,*)'nvirt=',nvirt
      
      do i = 1, nvirt
         x1(i)=x(1, ntotal + i)
         y1(i)=x(2, ntotal + i)
      enddo
      call pgsci(4)
      call pgpt(nvirt,x1,y1,16)
      
      
      mp=nvirt
      
      write(*,*)'r_vmin=',x1(2)-x1(1)
      write(10,*)'r_vmin=',x1(2)-x1(1)
      
c      nvirt = nvirt-1
      write(*,*)'nvirt=',nvirt,'ntotal=',ntotal
      write(10,*)'nvirt=',nvirt,'ntotal=',ntotal
      
      do i = 1, nvirt
         rho(ntotal + i) = 1000.
         mass(ntotal + i) = rho(ntotal + i) * dx * dx
         p(ntotal + i) = b*((rho(ntotal+i)/rho0)**gamma - 1.0)
         u(ntotal + i) = 357.1
         itype(ntotal + i) = -1
         hsml(ntotal + i) = hsml1 !dx
c         write(*,*)ntotal + i,'nvirt=',nvirt,x(1, ntotal + i)
      enddo
      write(*,*)'hsml1, dx',hsml1,dx
      write(10,*)'hsml1, dx',hsml1,dx

c      call pgbeg(0,'pgp_virt_part.ps/cps',1,0)
c      call pgscf(2)
c      call pgsch(1.5)
      
c      call pgsci(1)
c      call pgpt(nvirt,x1,y1,-1)
      
      return
      end

ccccccccccccccccccccccccccccccccccccccccc

      subroutine walls(dx,leng,xy,vxy,ang,maxn,ntotal,dim,x,vx,nvirt)
      implicit none

      integer i,maxn,nvirt,mp,dim,ntotal
      double precision dl,dx,leng,xy(dim),vxy(dim),ang,x(dim,maxn),
     +     vx(dim, maxn),ang1,PI

      PI=4.*atan(1.0)
      
      dl = 2.
      mp=nint(dl*leng/dx)
      ang1 = ang*PI/180.
      write(*,*)'Write mp=',mp      
      
      do i=1,mp+1
         nvirt = nvirt + 1
         x(1, ntotal + nvirt) = xy(1) + (i-1)*dx/dl*cos(ang)
         x(2, ntotal + nvirt) = xy(2) + (i-1)*dx/dl*sin(ang)
         vx(1, ntotal + nvirt) = vxy(1)
         vx(2, ntotal + nvirt) = vxy(2)
         
      enddo      
      return
      end
ccccccccccccccccccccccccccccccccccccccccc

      subroutine boxparticles(dim,m,maxn,rho0,bc,ht,vxi,dx,nt,x,vx,
     +     mass,rho,p,u,itype,hsml,xybox)
      implicit none

      integer dim,nt,ni,m,maxn,mp,n,np,ntotal,i,j,k
      integer itype(maxn)
      double precision rho0,bc(dim,4),dx,dy,ht,xx1,yy1
      double precision b,gamma,rhoa
      double precision xa(dim, maxn),vxi(2)
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn),
     &     p(maxn), u(maxn), hsml(maxn), rho(maxn),dxy
      double precision c,beta,g, xybox(5),PI,ang

      PI = 4.*atan(1.0)
      ang = xybox(5)*PI/180.
      g = 9.8
      
      xx1 = bc(1,2)-bc(1,1) 
      yy1 = bc(2,3)-bc(2,1)
      
      np = m-1
      dy = yy1/np
      dx = dy
      m = nint(xx1/dx) + 1
      mp = m-1
      ntotal =  mp * np
      ni = nt
      dxy = dx
      
c     nt=0
      gamma=7.
      beta = 5.0!0.003
      c = beta*sqrt(9.8*ht)
c      c=1480.
c      b = c*c*rho0/gamma
c      b = beta*g*ht*rho0/(gamma)

c      c = 1400
      b = rho0*c*c/gamma
      
      do i = 1, mp
         do j = 1 , np
            nt = nt + 1
            x(1, nt) = ((i-1)*dx + dx/2. + bc(1,1) )
            x(2, nt) = ((j-1)*dy + dy/2. + bc(2,1) )
            vx(1,nt)=vxi(1)
            vx(2,nt)=vxi(2)
            rho (nt) = rho0*( 1+ rho0*9.8*(ht-x(2,j))/b  )**(1./gamma) !1000.
c            write(*,*)i,j,nt,rho(nt)
            mass(nt) = dx*dy*rho(i)
c            b = 200*9.8*ht/(rho(nt)*gamma)
            p(nt)= b*((rho(nt)/rho0)**gamma - 1.0)
            u(nt)=357.1
            itype(nt) = 1
            hsml(nt) = dxy
         enddo
      enddo      

      write(*,*)'m=',m,'n=',n,'dx=',dx,'dy=',dy
     +     ,'ntotal=',ntotal,'nt=',nt,'dxy=',dxy
            
      return
      end


c=================================================================

! this subrutain make a box
      subroutine virtbox(beta,ntotal,maxn,dim,dx,hsml1,dmasa,
     +     itype,hsml,mass,x,vx,rho,u,p,nvirt,ht,rho0,xybox,sel)
      implicit none

! xybox(1): x Initial Position of the box
! xybox(2): y Initial Position of the box
! xybox(3): High of the box
! xybox(4): lower side of the box
! xybox(5): angle position of the box
      
      
      integer ntotal, nvirt, maxn, dim, sel
      integer itype(maxn)
      double precision hsml(maxn),mass(maxn),x(dim,maxn),vx(dim, maxn),
     &     rho(maxn), u(maxn), p(maxn)
      integer i, j, mp
      double precision leng, dx,v_inf,hsml1,dl,dmasa,xy(2),vxy(2),ang,PI

      double precision gamma,beta,c,b,ht,rho0,xybox(5)
      real x1(maxn),y1(maxn),xmin,xmax,ymin,ymax

      PI=4.*atan(1.0)

      gamma=7.
c      beta = 5.0!0.003
      c = beta*sqrt(9.8*ht)
      b = rho0*c*c/gamma

      nvirt = 0

      ! Left side
      leng   = xybox(3)
      ang    = (xybox(5) + 90.0d0) * (PI/180.0d0)! 90 * (PI/180.)
      xy(1)  = xybox(1)
      xy(2)  = xybox(2)
      vxy(1) = 0.0
      vxy(2) = 0.0
      call walls(dx,leng,xy,vxy,ang,maxn,ntotal,dim,x,vx,nvirt)
      write(*,*)'nvirt=',nvirt
      
      !lower side
      leng   = xybox(4) - dx
      ang    = xybox(5) * (PI/180.)
      xy(1)  = dx/2. + xybox(1)
      xy(2)  = xybox(2)
      vxy(1) = 0.0
      vxy(2) = 0.0
      call walls(dx,leng,xy,vxy,ang,maxn,ntotal,dim,x,vx,nvirt)
      write(*,*)'nvirt=',nvirt

      !Rigth side
      leng   = xybox(3)
      ang    = (xybox(5) + 90) * (PI/180.)
      xy(1)  = xybox(1) + xybox(4)*dcos(xybox(5)*PI/180.)
      xy(2)  = xybox(2) + xybox(4)*dsin(xybox(5)*PI/180.)
      vxy(1) = 0.0
      vxy(2) = 0.0
      call walls(dx,leng,xy,vxy,ang,maxn,ntotal,dim,x,vx,nvirt)
      write(*,*)'nvirt=',nvirt

      if(sel.eq.1)then
      !Upper side
         leng   = xybox(4) - dx/2.
         ang    = xybox(5) * (PI/180.)
         xy(1)  = xybox(1) - xybox(3)*dsin(xybox(5)*PI/180.)
         xy(2)  = xybox(2) + dx/2. + xybox(3)*dcos(xybox(5)*PI/180.)
         vxy(1)=0.0
         vxy(2)=0.0
         call walls(dx,leng,xy,vxy,ang,maxn,ntotal,dim,x,vx,nvirt)
         write(*,*)'nvirt=',nvirt
      endif
         
      
      do i = 1, nvirt
         x1(i)=x(1, ntotal + i)
         y1(i)=x(2, ntotal + i)
      enddo
      call pgsci(4)
      call pgpt(nvirt,x1,y1,1)
      
      
      mp=nvirt
      
      write(*,*)'r_vmin=',x1(2)-x1(1)
      write(10,*)'r_vmin=',x1(2)-x1(1)
      
      write(*,*)'nvirt=',nvirt,'ntotal=',ntotal
      write(10,*)'nvirt=',nvirt,'ntotal=',ntotal
      
      do i = 1, nvirt
         rho(ntotal + i) = 1000.
         mass(ntotal + i) = rho(ntotal + i) * dx * dx
         p(ntotal + i) = 0.0!b*((rho(ntotal+i)/rho0)**gamma - 1.0)
         u(ntotal + i) = 357.1
         itype(ntotal + i) = -1
         hsml(ntotal + i) = hsml1 !dx
c         write(*,*)ntotal + i,'nvirt=',nvirt,x(1, ntotal + i)
      enddo
      write(*,*)'hsml1, dx',hsml1,dx
      write(10,*)'hsml1, dx',hsml1,dx
     
      return
      end

c========================================================================      
      subroutine boxparticles01(dim,m,maxn,rho0,bc,ht,vxi,dx,nt,x,vx,
     +     mass,rho,p,u,itype,hsml,delt,xybox,mu,beta)
      implicit none

      integer dim,nt,ni,m,maxn,mp,n,np,ntotal,i,j,k
      integer itype(maxn)
      double precision rho0,bc(dim,4),dx,dy,ht,htt,xx1,yy1,delt
      double precision b,gamma,rhoa
      double precision xa(dim, maxn),vxi(2)
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn),
     &     p(maxn), u(maxn), hsml(maxn), rho(maxn),dxy
      double precision c,beta,g,mu,xybox(5),PI,ang,vxy(2),xy(2)
      double precision delttiempo,alpha,vis

      PI = 4.*atan(1.0)
      ang = xybox(5)*PI/180.
      vxy(1) = 0.0
      vxy(2) = 0.0
      xx1 = bc(1,2)-bc(1,1) 
      yy1 = bc(2,3)-bc(2,1)
      htt = ht*cos(ang)!bc(1,3)*cos(ang) - bc(2,3)*sin(ang) + xybox(2)
      write(*,*)'PARMA ARCHIV, ht = ',htt
      
      g = 9.82
     
c      np = m-1
c      dy = yy1/np
c      dx = dy
c      m = nint(xx1/dx) + 1
c      mp = m-1
c      ntotal =  mp * np
c      ni = nt
c      dxy = dx

      dy = delt
      np = nint(yy1/dy)
      dx = dy
      m = nint(xx1/dx) + 1
      mp = m-1
      ntotal =  mp * np
      ni = nt
      dxy = dx
      
c     nt=0
      gamma=7.
c      beta = 5.0!0.003
      c = sqrt(beta*g*htt)
c      c=1480.
c      b = c*c*rho0/gamma
c      b = beta*g*ht*rho0/(gamma)

c      c = 1400
      b = rho0*c*c/gamma
      
      do i = 1, mp
         do j = 1 ,np
            nt = nt + 1
            xy(1) = (i-1)*dx + dx/2. + bc(1,1) 
            xy(2) = (j-1)*dy + dy/2. + bc(2,1) 
            x(1, nt) = xy(1)*cos(ang) - xy(2)*sin(ang) + xybox(1)
            x(2, nt) = xy(1)*sin(ang) + xy(2)*cos(ang) + xybox(2)
            vx(1,nt)=vxi(1)
            vx(2,nt)=vxi(2)
c     rho (nt) = rho0*( 1+ rho0*9.8*(ht-x(2,j))/b  )**(1./gamma)
c            rho (nt) = rho0*( 1+ rho0*9.8*(ht-
c     +           x(2,nt) )/b  )**(1./gamma)
c            rho (nt) = rho0*( 1+ rho0*9.8*((ht+xybox(2)) -
c     +           abs(x(2,nt)) )/b  )**(1./gamma)
            rho (nt) = rho0*( 1 + rho0 * g * ( htt -
     +           abs(x(2,nt) ) )/b  )**(1./gamma)
            mass(nt) = dx*dy*rho(i)
            p(nt)= b*((rho(nt)/rho0)**gamma - 1.0)
            u(nt)=357.1
            itype(nt) = 1
            hsml(nt) = dxy
c            write(*,*)nt,ht,x(2,nt),rho(nt),p(nt)
         enddo
      enddo      

      write(*,*)'m=',m,'n=',n,'dx=',dx,'dy=',dy
     +     ,'ntotal=',ntotal,'nt=',nt,'dxy=',dxy


!     Calculos del paso de tiempo

!CFL condition

      delttiempo = 0.25 * htt / ( c * sqrt(2*htt*g))
      write(*,*)'Time step CFL condition',delttiempo

!     viscous condtion
      alpha = 0.1
      vis = 1./(2.*(dim+2.)) * alpha*delt*c
      delttiempo = 0.125 * htt*htt / vis
      write(*,*)'Time step viscous condtion',delttiempo

!     physical kinematic viscosity
     
c      delttiempo = 0.125 * htt*htt / v
c      write(*,*)'Time step physical kinematic viscosity',delttiempo
      
      return
      end
