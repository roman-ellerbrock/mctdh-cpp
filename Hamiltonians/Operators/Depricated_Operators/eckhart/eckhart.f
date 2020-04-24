c######################################################################
c MCTDH-Modul Operatoren
c
c Contains subroutines "h init" and "h".
c "h init" initializes the arrays "koeff" and "diag" which
c specify the coefficients employed in the Hamiltonian series.
c "h" calculates the action of a single-particle operator (on the
c primitive grid). "mode" and "teil" define the number of the
c coordinate and of the term in the series, respectively.
c
c----------------------------------------------------------------------
c 3D-kinetic energy operator in Jacobi-coordinates:
c 1st coordinate: FFT
c 2nd coordinate: FFT
c 3rd coordinate: Legendre-DVR
c parameters: 
c             m(1), m(2) (dissociative and diatomic reduced masses) 
c             r0   (dividing surface)
c             rabs (beginning of absorber)
c             rend (end of absorber)
c
c V2.0, 25.3.1996
c######################################################################

      subroutine h initeckhart (koeff,diag,dim)

      implicit none

      complex*16  koeff(*)
      integer     dim, k zahl
      logical*4   diag(3)

      k zahl=4
      koeff(1)=1
      koeff(2)=1
      koeff(3)=1
      koeff(4)=1
      koeff(5)=1

      koeff(k zahl+1)=0

c Bestimmung der diagonalen Komponenten des Hamilton-Operators

      call diagoneckhart(diag,dim,k zahl)
      return
      end

c--------------------------------------------------------------

      subroutine diagoneckhart(diag,dim,k zahl)

      integer     dim, k zahl, i, k
      logical*4   diag(dim,k zahl)

      do k=1,k zahl
         do i=1,dim
            diag(i,k)=.false.
         enddo
      enddo
      return
      end

c#####################################################################

      subroutine heckhart (mode,teil,h psi,psi,dim,matrix,trafo,ort)

      implicit none

      integer       mode, teil, dim, i, j
      complex*16    h psi(dim), psi(dim), h psi2(dim)
      complex*16    work(1024), work2(1024), p psi(1024)
      real*8        matrix(dim,dim), trafo(dim,dim), ort(dim)
      real*8        x

      if (teil.eq.-1) goto (11,1) mode
      if (teil.eq.-2) goto (12,1) mode

      goto (100,200) mode
 100  goto (2,1,4,1,6) teil
 200  goto (1,3,1,5,1) teil
      return

c H = 1

 1    do i=1,dim
         h psi(i)=psi(i)
      enddo
      return

c H = -(d/dx)**2

 2    do i=1,dim
         h psi(i)=psi(i)
      enddo
      call fft(h psi,work,trafo(1,3),trafo(1,5),dim,1,1,.true.)
      do i=1,dim
         h psi(i)=matrix(i,2)*h psi(i)
      enddo
      call fft(h psi,work,trafo(1,1),trafo(1,5),dim,1,1,.false.)
      return

c     H = -0.5*(d/dx)**2

 3    do i=1,dim
         h psi(i)=0
         do j=1,dim
            h psi(i)=h psi(i)+matrix(i,j)*psi(j)
         enddo
      enddo
      return

c 

 4    do i=1,dim
         x=ort(i)/dsqrt(1060.d0)
         h psi(i)=psi(i)*(0.425/27.21)*(2/(exp(x)+exp(-x)))**2
      enddo
      return

c 

 5    do i=1,dim
         h psi(i)=psi(i)*0.5*(0.2/27.21)**2*ort(i)**2
      enddo
      return

 6    do i=1,dim
	   if (abs(ort(i)).gt.500) then
	      h psi(i)= 5/27.2114*psi(i)
         endif
      enddo

c H = 0

 9    do i=1,dim
         h psi(i)=0
      enddo
      return

c h=theta(r0-r)

 11   do i=1,dim
         if (ort(i).lt.-400) then
            h psi(i)=0
         else
            h psi(i)=psi(i)
         endif
      enddo
      return

c i[H,h(r0-r)]

 12   do i=1,dim
         h psi(i)=psi(i)
      enddo
      call fft(h psi,work,trafo(1,1),trafo(1,5),dim,1,1,.false.)
      call rvecvec(h psi,h psi,matrix(1,2),dim,2,1)
      call fft(h psi,work,trafo(1,3),trafo(1,5),dim,1,1,.true.)
      do i=1,dim
         if (ort(i).lt.-300) h psi(i)=0
      enddo

      do i=1,dim
         if (ort(i).lt.-300) then
            work2(i)=0
         else
            work2(i)=psi(i)
         endif
      enddo
      call fft(work2,work,trafo(1,3),trafo(1,5),dim,1,1,.true.)
      call rvecvec(work2,work2,matrix(1,2),dim,2,1)
      call fft(work2,work,trafo(1,1),trafo(1,5),dim,1,1,.false.)
      do i=1,dim
         h psi(i)=(work2(i)-h psi(i))!*(0,1)
      enddo
      return

      end


