
      program plot_snaps
      implicit none

      integer i,j,k,np,rangg(2),npp,ip,nt,ntype(2),nf,nfluid
      parameter(nt=401,np=20000)
      integer id(np),itype(np)
      real t(np),x(np),y(np),z(np),vx(np),vy(np),vz(np),mass(np),p(np)
     +     ,rho(np),u(np),hsml(np),xf(np),zf(np)
      real xmin,xmax,ymin,ymax
      character infilebas*80,outfile(10000)*80,tch*5

c --  Make an array of outfiles using a base outfile
      infilebas = 'snapshot'
      rangg(1)=1
      rangg(2)=10000
      call array_infilebase(infilebas,rangg,outfile,1,nf)

      call pgbeg(0,'pgp_snpas.ps/cps',1,1)
      call pgscf(2)
      call pgsch(1.0)
      xmin = -0.05              
      xmax =  1.7
      ymin =  -0.05
      ymax =  0.7
    
      do i=1,730
         write(*,*)'loading file = ',outfile(i)
         open(1,file=outfile(i))
         ip = 1
         read(1,*)npp,t(i),ntype(1),nfluid,ntype(2)
 3       read(1,*,end=4)id(ip),x(ip),z(ip),vx(ip),vz(ip),mass(ip),p(ip)
     +        ,rho(ip),u(ip),itype(ip),hsml(ip)

         ip = ip + 1
         goto 3
 4       ip = ip -1
         write(*,*)ip,npp,t(i),ntype
         
         nf = nfluid+1!ntype(1)-ntype(2)+1
         k=0
         do j=nf,ntype(1)
            k=k+1
            xf(k) = x(j)
            zf(k) = z(j)            
         enddo

         write(tch,'(F5.3)')t(i)

         call pgsci(1)
         call pgenv(xmin,xmax,ymin,ymax,1,0)
         call pgsci(4)         
         call pgpt(nf-1,x,z,1)
         call pgsci(2)
         call pgpt(k-1,xf,zf,17)
         call pgsci(1)
         call pgtext(0.1, -0.05, 't = '//tch)
      enddo

      call pgend
      
      end

      
c======================================================================
c THIS SUBROUTINE MAKE AN ARRAY OF INFILES USING A BASE OF INFILE.

      SUBROUTINE array_infilebase(infilebas,rangg,infile,step,nf)
      implicit none

      integer n,rangg(2),i,l,j,step,nf,k
      character infilebas*80,infile(10000)*80,ch1,ch2*2,ch3*3,ch4*4
      character ch5*5,ch6*6

      l = len_trim(infilebas)
      k = 0
      
      j=rangg(1)-1
      DO 10 i=1,10000
         j=j+1
c         if(mod(i,step).eq.1)then
            if(j.le.rangg(2))then
               k=k+1
               if(j.lt.10)write(ch1,'(I1)')j
               if(j.lt.10)infile(k)=infilebas(1:l)//'_000'//ch1
               
               if(j.ge.10.and.j.lt.100)write(ch2,'(I2)')j
               if(j.ge.10.and.j.lt.100)
     +              infile(k)=infilebas(1:l)//'_00'//ch2
               
               if(j.ge.100.and.j.lt.1000)write(ch3,'(I3)')j
               if(j.ge.100.and.j.lt.1000)
     +              infile(k)=infilebas(1:l)//'_0'//ch3
               
               if(j.ge.1000.and.j.lt.10000)write(ch4,'(I4)')j
               if(j.ge.1000.and.j.lt.10000)
     +              infile(k)=infilebas(1:l)//'_'//ch4
               
               if(j.ge.10000.and.j.lt.100000)write(ch5,'(I5)')j
               if(j.ge.10000.and.j.lt.100000)
     +              infile(k)=infilebas(1:l)//'_'//ch5
               
               if(j.ge.100000.and.j.lt.1000000)write(ch6,'(I6)')j
               if(j.ge.100000.and.j.lt.1000000)
     +              infile(k)=infilebas(1:l)//'_'//ch6
               
c     write(*,*)i,infile(i),j
            endif
c         endif
 10   CONTINUE
      
      nf = k
      
      RETURN
      END      
