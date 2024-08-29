      program ordenar_particulas
      implicit none
      
      integer n
      parameter (n = 50)
      double precision x(n), y(n), x_temp(n), y_temp(n)
      integer indx_x(n), indx_y(n)
      integer i
      
! Inicializar posiciones (aquí se usa un ejemplo de posiciones distribuidas uniformemente)
      call random_number(x)
      call random_number(y)
      x = x * 0.6               ! Escalar a 0.6 m en x
      y = y * 0.3               ! Escalar a 0.3 m en y
      write(*,*)x,y
      
! Llamar a la subrutina indexx para ordenar por x
      call indexx(n, x, indx_x)
      
! Reordenar las posiciones usando los índices ordenados por x
      do i = 1, n
         x_temp(i) = x(indx_x(i))
         y_temp(i) = y(indx_x(i))
      end do
      
! Copiar los valores reordenados de vuelta a los arreglos originales
      x = x_temp
      y = y_temp
      
! Llamar a la subrutina indexx para ordenar por y
      call indexx(n, y, indx_y)
      
! Reordenar las posiciones usando los índices ordenados por y
      do i = 1, n
         x_temp(i) = x(indx_y(i))
         y_temp(i) = y(indx_y(i))
      end do
      
! Copiar los valores reordenados de vuelta a los arreglos originales
      x = x_temp
      y = y_temp
      
! Imprimir los primeros 10 valores para verificar
      print *, "Primeros 10 valores ordenados:"
      do i = 1, 10
         print *, x(i), y(i)
      end do
      end
      
      
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      double precision arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)stop! 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *!^3#!0Y..
