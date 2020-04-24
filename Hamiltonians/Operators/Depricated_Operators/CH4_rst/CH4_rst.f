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
c Operator for CH_4 by Gerd Schiffel
c
c######################################################################

      subroutine h initch4rst(koeff,diag,dim)
      use iso_c_binding
      implicit none

      integer     dim, k zahl, i
      integer     diag(3)
      complex*16  koeff(*),c(200)
      real*8      dph2,dth2,dph,dchi

      dph2 = 0.0d0
      dth2 = 0.0d0
      dph  = 0.0d0
      dchi = 0.0d0

      dph2 = 0.0030d0
      dth2 = 0.0045d0
      dph  = 0.0100d0
      dchi = 0.0300d0

c      call potinich4
c      call potini
c     call koordread

c [1,1],[2,2],[3,3]
      do i=1,4
         c(i)=1.0d0
      end do
c [4,4]
      do i=5,6
         c(i)=1.0d0+ 0.25d0*
     .   (8.0d0*dph2+16.0d0/3.0d0*dth2+4.0d0/3.0d0*dph+4.0d0/9.0d0*dchi)
      end do
c [5,5]
      c(7)=2.0d0+ 0.25d0*
     .    (16.0d0*dph2+32.0d0/3.0d0*dth2+4.0d0*dph+20.0d0/9.0d0*dchi)
      c(8)=-1.0d0+ 0.25d0*
     .     (-11.0d0*dph2-8.0d0/3.0d0*dth2-4.0d0*dph-20.0d0/9.0d0*dchi)
c [6,6]
      c(9)=2.0d0+ 0.25d0*
     . (16.0d0*dph2+32.0d0/3.0d0*dth2+20.0d0/3.0d0*dph+4.0d0/3.0d0*dchi)
      c(10)=-1.0d0+ 0.25d0*
     . (-5.0d0*dph2-8.0d0*dth2-20.0d0/3.0d0*dph-4.0d0/3.0d0*dchi)
c [4,5]
      do i=11,12
       c(i)=1.0d-99+ 0.25d0*
     .   (-9.0d0*dph2+8.0d0*dth2+7.0d0*dph-7.0d0/3.0d0*dchi)/sqrt(3.0d0)
      end do
c [7,7]
      do i=13,22
        c(i)=2.0d0+ 0.25d0*
     .         ( 10.0d0*dph2+16.0d0*dth2+8.0d0*dph+16.0d0/3.0d0*dchi)
      end do
c [8,8]
      do i=23,32
        c(i)=2.0d0+ 0.25d0*
     . (22.0d0*dph2+16.0d0/3.0d0*dth2+16.0d0/3.0d0*dph-8.0d0/9.0d0*dchi)
      end do
c [9,9]
      do i=33,36
         c(i)=5.0d0+ 0.25d0*
     .       (52.0d0*dph2+16.0d0*dth2+8.0d0*dph+16.0d0*dchi)
      end do
      do i=37,40
         c(i)=-2.0d0+0.25d0*
     .       (-10.0d0*dph2-16.0d0*dth2-8.0d0*dph-16.0d0*dchi)
      end do
c [6,9]
      do i=41,42
         c(i)=2.0d0+ 0.25d0*
     .      (19.0d0*dph2+8.0d0*dth2+5.0d0*dph+17.0d0/9.0d0*dchi)
      end do
      do i=43,44
         c(i)=-1.0d0+ 0.25d0*
     .      (-5.0d0*dph2-8.0d0*dth2-5.0d0*dph-17.0d0/9.0d0*dchi)
      end do
c [4,8]
      do i=45,47
        c(i)=1.0d-99+ 0.25d0*
     .  (6.0d0*dph2-16.0d0/3.0d0*dth2-10.0d0/3.0d0*dph+2.0d0/9.0d0*dchi)
      end do
c [6,7]
      do i=48,50
        c(i)=-3.0d0+ 0.25d0*
     .   (-15.0d0*dph2-24.0d0*dth2-15.0d0*dph-17.0d0/3.0d0*dchi)
      end do
c [5,8]
      do i=51,53
       c(i)=-sqrt(3.0d0)*0.25d0
     .     *(4.0d0+11.0d0*dph2+8.0d0/3.0d0*dth2+3.0d0*dph+dchi/9.0d0)
      end do
c [7,9]
      do i=54,60
        c(i)=2.0d0+ 0.25d0*
     .    (10.0d0*dph2+16.0d0*dth2+8.0d0*dph+16.0d0/3.0d0*dchi)
      end do
c Ha, Hb
      do i=61,69
        c(i)=1.0d0
      end do

      koeff( 1)=-2.000d0
      koeff( 2)=-2.000d0
      koeff( 3)=-2.000d0
      koeff( 4)= 0.250d0
      koeff( 5)= 0.750d0
      koeff( 6)=-2.000d0
      koeff( 7)= 1.500d0
      koeff( 8)= 1.500d0
      koeff( 9)= 4.500d0
      koeff(10)= 4.500d0
      koeff(11)= 1.000d0
      koeff(12)= 1.000d0
      koeff(13)= 0.250d0
      koeff(14)= 0.250d0
      koeff(15)= 0.250d0
      koeff(16)=-2.000d0
      koeff(17)=-0.250d0
      koeff(18)= 0.250d0
      koeff(19)= 0.250d0
      koeff(20)=-0.500d0
      koeff(21)=-1.000d0
      koeff(22)=-0.500d0
      koeff(23)= 0.250d0
      koeff(24)= 0.250d0
      koeff(25)= 0.250d0
      koeff(26)=-0.500d0
      koeff(27)=-1.000d0
      koeff(28)=-0.500d0
      koeff(29)= 0.250d0
      koeff(30)= 0.250d0
      koeff(31)=-0.250d0
      koeff(32)=-2.000d0
      koeff(33)= 0.500d0
      koeff(34)=-2.000d0
      koeff(35)=-2.000d0
      koeff(36)=-0.500d0
      koeff(37)= 0.500d0
      koeff(38)=-2.000d0
      koeff(39)=-2.000d0
      koeff(40)=-0.500d0
      koeff(41)= 6.000d0
      koeff(42)=-6.000d0
      koeff(43)= 6.000d0
      koeff(44)=-6.000d0
      koeff(45)=-1.000d0
      koeff(46)=-1.000d0
      koeff(47)=-1.000d0
      koeff(48)= 1.000d0
      koeff(49)= 1.000d0
      koeff(50)= 1.000d0
      koeff(51)=-1.000d0
      koeff(52)=-1.000d0
      koeff(53)=-1.000d0
      koeff(54)= 0.500d0
      koeff(55)= 2.000d0
      koeff(56)=-2.000d0
      koeff(57)= 2.000d0
      koeff(58)=-1.500d0
      koeff(59)= 0.500d0
      koeff(60)= 0.500d0
      koeff(61)=-2.000d0
      koeff(62)=-0.500d0
      koeff(63)=-0.500d0
      koeff(64)=-1.000d0
      koeff(65)=-1.000d0
      koeff(66)=-0.500d0
      koeff(67)=-0.500d0
      koeff(68)= 0.500d0
      koeff(69)= 0.500d0

      kzahl=69

      koeff(kzahl+1)=0
      do i=1,kzahl
         koeff(i)=-0.5d0*koeff(i)*c(i)
      end do


c Bestimmung der diagonalen Komponenten des Hamilton-Operators
      koeff(kzahl+1)=0
      call diagonch4rst(diag,dim,kzahl)

      return
      end

c--------------------------------------------------------------

      subroutine diagonch4rst(diag,dim,k zahl)
      use iso_c_binding

      integer     dim, k zahl, i, k
      integer     diag(dim,k zahl)


      real*8      ort,matrix(2)
      integer     trafo(12)
      complex*16  hpsi,psi

      do k=1,k zahl
         do i=1,dim
c            diag(i,k)=.false.
            diag(i,k)=0
         enddo
      enddo

      hpsi     =(1.0d0,0.0d0)
      psi      =(1.0d0,0.0d0)
      matrix(1)= 1.0d0
      matrix(2)= 1.0d0
      do i=1,6
         trafo(i) = 1
      end do

      do k=1,kzahl
         do i=1,dim
            ort=998.0d0
            call h ch4 rst(i,k,hpsi,psi,1,matrix,trafo,ort)
            if (abs(ort-999.0d0).lt.1.0d-2) then
c               diag(i,k)=.true.
               diag(i,k)=1
            endif
         enddo
      enddo

      return
      end

c#####################################################################

      subroutine h ch4rst (mode,teil,h psi,psi,dim,matrix,trafo,ort)

      implicit      none

      integer       mode, teil, dim, i, j
      complex*16    h psi(dim), psi(dim)
      complex*16    work(1024),work2(1024),work3(1024)
      real*8        matrix(dim,dim), trafo(dim,dim), ort(dim), r0

      real*8        xphich4,xchich4

      real*8        pi,pi3
      parameter     (pi=3.14159265358979323844d0)
      parameter     (pi3=3.14159265358979323844d0/3.0d0)

      goto (101,102,103,104,105,106,107,108,109) mode
      print*,'unknown mode'
      stop

101   goto (
     &30, 6, 6, 6, 6,   6, 6, 6, 6, 6,    6, 6, 6, 6, 6,  6, 6, 6, 6, 6,
     & 6, 6,
     &    6, 6, 6, 6,   6, 6, 6, 6, 6,    6,
     &                                    6, 6, 6, 6,
     &       6, 6, 6,   6,                            6,  6,       6, 6,
     &
     &       6, 6, 6,            6, 6,    6,          6,  6, 6,
     & 6, 6, 6, 6, 6,   6, 6,
     &                     1, 1, 1, 1,    1, 1, 1, 1, 1,
     &             1) teil
      if (teil.eq.-1) goto 500
      if (teil.eq.-2) goto 500
      print *, 'unknown teil in 1'
      stop

102   goto (
     & 1,30, 8, 9, 1,   1, 1, 1, 1, 1,    1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
     & 1, 1,
     &    1, 1, 1, 1,   1, 1, 1, 1, 1,    1,
     &                                    1, 1, 1, 1,
     &       1, 1, 1,   1,                            1,  1,       1, 1,
     &
     &       1, 1, 1,            1, 1,    1,          1,  1, 1,
     & 1, 1, 1, 1, 1,   1, 1,
     &                     1, 1, 1, 1,    1, 1, 1, 1, 1,
     &             1) teil
      if (teil.eq.-1) goto 500
      if (teil.eq.-2) goto 500
      print *, 'unknown teil in 2'
      stop

103   goto (
     & 1, 1,30, 1, 1,   1, 1, 1, 1, 1,    1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
     & 1, 1,
     &    1, 1, 1, 1,   1, 1, 1, 1, 1,    1,
     &                                    1, 1, 1, 1,
     &       1, 1, 1,   1,                            1,  1,       1, 1,
     &
     &       1, 1, 1,            1, 1,    1,          1,  1, 1,
     & 1, 1, 1, 1, 1,   1, 1,
     &                     1, 1, 1, 1,    1, 1, 1, 1, 1,
     &             1) teil
      if (teil.eq.-1) goto 500
      if (teil.eq.-2) goto 500
      print *, 'unknown teil in 3'
      stop

104   goto (
     & 1, 1, 1, 1,10,  30, 8, 1, 8, 1,   24,11, 1, 1, 1,  1, 1, 1, 1, 1,
     & 1, 1,
     &    1, 1, 1, 1,   1, 1, 1, 1, 1,    1,
     &                                    8, 8, 8, 8,
     &       1, 1, 1,   1,                            8,  8,       1, 1,
     &
     &      20,20,20,            7, 7,    7,          7,  7, 7,
     & 7, 7, 7, 7, 7,   7, 7,
     &                     1, 1, 1, 1,    1, 1, 1, 1, 1,
     &             1) teil
      if (teil.eq.-1) goto 500
      if (teil.eq.-2) goto 500
      print *, 'unknown teil in 4'
      stop

105   goto (
     & 1, 1, 1, 1, 1,   1,25,25, 1, 1,   20,27, 1, 1, 1,  1, 1, 1, 1, 1,
     & 1, 1,
     &    1, 1, 1, 1,   1, 1, 1, 1, 1,    1,
     &                                    1, 1, 1, 1,
     &       1, 1, 1,   1,                            1,  1,       1, 1,
     &
     &       1, 1, 1,            1, 1,    1,         20, 20,20,
     & 1, 1, 1, 1, 1,   1, 1,
     &                     1, 1, 1, 1,    1, 1, 1, 1, 1,
     &             1) teil
      if (teil.eq.-1) goto 500
      if (teil.eq.-2) goto 500
      print *, 'unknown teil in 5'
      stop

106   goto (
     & 1, 1, 1, 1, 1,   1, 1, 1,26,26,    1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
     & 1, 1,
     &    1, 1, 1, 1,   1, 1, 1, 1, 1,    1,
     &                                    1, 1, 1, 1,
     &       1, 1, 1,   1,                           20, 20,      20,20,
     &
     &       1, 1, 1,           20,20,   20,          1,  1, 1,
     & 1, 1, 1, 1, 1,   1, 1,
     &                     1, 1, 1, 1,    1, 1, 1, 1, 1,
     &             1) teil
      if (teil.eq.-1) goto 500
      if (teil.eq.-2) goto 500
      print *, 'unknown teil in 6'
      stop

107   goto (
     & 1, 1, 1, 1, 1,   1, 1, 1, 1, 1,    1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
     & 1, 1,
     &    1, 1, 1, 1,   1, 1, 1, 1, 1,    1,
     &                                    1, 1, 1, 1,
     &       1, 1, 1,   1,                            1,  1,       1, 1,
     &
     &       1, 1, 1,            1, 1,    1,          1,  1, 1,
     & 1, 1, 1, 1, 1,   1, 1,
     &                    30, 6, 6, 6,    6, 6, 6, 6, 6,
     &             1) teil
      if (teil.eq.-1) goto 201
      if (teil.eq.-2) goto 202
      print *, 'unknown teil in 7'
      stop

108   goto (
     & 1, 1, 1, 1, 1,   1, 1, 1, 1, 1,    1, 1, 1, 3, 1, 31,23,22,22, 5,
     & 4, 1,
     &    1, 3, 1,32,  31,30,23,22,22,    3,
     &                                    1,30, 3,22,
     &       1,30, 3,  22,                            2, 20,       2,20,
     &
     &       2,20,21,           22, 4,    1,          2, 20,21,
     & 2,12,33, 2,21,  20,20,
     &                     1,30,14,31,   13,32, 1,13, 1,
     &             1) teil
      if (teil.eq.-1) goto 500
      if (teil.eq.-2) goto 500
      print *, 'unknown teil in 8'
      stop

109   goto (
     & 1, 1, 1, 1, 1,   1, 1, 1, 1, 1,    1, 1, 1, 1, 3,  3,22,22,23,30,
     &31,32,
     &    1, 1, 3, 1,   4, 5,22,22,23,   31,
     &                                    1, 3,30,22,
     &       1, 3,30,  22,                           20,  2,      20, 2,
     &
     &      22, 4, 1,            2,20,   21,         22,  4, 1,
     & 1,30, 3,31,22,  22,23,
     &                     1,14,30,13,   31, 1,32, 1,13,
     &             1) teil
      if (teil.eq.-1) goto 500
      if (teil.eq.-2) goto 500
      print *, 'unknown teil in 9'
      stop

      print *, 'gar nicht gut'
      stop


c --------------------------------------------------------------------
c operators:
c
c  1:  1
c  2:  1/x**2
c --------------------------------------------------------------------

c H = 1

500   do i=1,dim
         h psi(i)=psi(i)
      enddo
      return

c H = 1, set diagon

 1    if (abs(ort(1)-998.0d0).gt.1.0d-2) then
         write(*,*) "Error in subroutine h."
         write(*,*) "Do not use label 1 for identity. Use label 500."
         stop
      endif
      ort(1)=999.0d0
      return

c H = x

 2    do i=1,dim
         h psi(i)=psi(i)*ort(i)
      enddo
      return

c H = x^2

 3    do i=1,dim
         h psi(i)=psi(i)*ort(i)*ort(i)
      enddo
      return

c H = (1-x^2)

 4    do i=1,dim
         h psi(i)=psi(i)*(-ort(i)*ort(i)+1.0d0)
      enddo
      return

c H = (1-x^2)^2

 5    do i=1,dim
c        h psi(i)=psi(i)*((ort(i)*ort(i)-2.0d0)*ort(i)*ort(i)+1.0d0)
         h psi(i)=psi(i)*(ort(i)**2-1.0d0)**2
      enddo
      return

c H = 1/x**2

 6    do i=1,dim
         h psi(i)=psi(i)/ort(i)**2
      end do
      return

c H = cot(x)

 7    do i=1,dim
         h psi(i)=psi(i)/tan(ort(i))
      end do
      return

c H = 1/sin(x)**2

 8    do i=1,dim
         h psi(i)=psi(i)/sin(ort(i))**2
      end do
      return


c H = 1 + 1/sin(x)**2

 9    do i=1,dim
         h psi(i)=psi(i)*(1.0d0+1.0d0/sin(ort(i))**2)
      end do
      return

c H = 3 - 1/sin(x)**2

10    do i=1,dim
         h psi(i)=psi(i)*(3.0d0-1.0d0/sin(ort(i))**2)
      end do
      return

c H = 1/sin(x)**2 - 3/2

11    do i=1,dim
         h psi(i)=psi(i)*(1.0d0/sin(ort(i))**2-1.5d0)
      end do
      return

c H = x*(1-x^2)

12    do i=1,dim
         h psi(i)=-psi(i)*(ort(i)*ort(i)-1.0d0)*ort(i)
      enddo
      return

c H = 1+x^2

13    do i=1,dim
         h psi(i)=psi(i)*(ort(i)*ort(i)+1.0d0)
      enddo
      return

c H = (1+x^2)^2

14    do i=1,dim
         h psi(i)=psi(i)*(ort(i)**2+1.0d0)**2
      enddo
      return

c H = d/dx

20    call ddx (mode,h psi,psi,dim,matrix,trafo)
      return

c H = x d/dx x

21    do i=1,dim
         work(i)=ort(i)*psi(i)
      end do
      call ddx (mode,h psi,work,dim,matrix,trafo)
      do i=1,dim
         hpsi(i)=hpsi(i)*ort(i)
      end do
      return

c H =  d/dx x  +  x d/dx

22    do i=1,dim
         work(i)=psi(i)*ort(i)
      end do
      call ddx (mode,h psi,work,dim,matrix,trafo)
      call ddx (mode,work,psi,dim,matrix,trafo)
      do i=1,dim
         hpsi(i)=hpsi(i)+work(i)*ort(i)
      end do
      return

c H =  d/dx x  +  x d/dx

23    do i=1,dim
         work(i)=psi(i)*ort(i)**3
      end do
      call ddx (mode,h psi,work,dim,matrix,trafo)
      call ddx (mode,work,psi,dim,matrix,trafo)
      do i=1,dim
         hpsi(i)=hpsi(i)+work(i)*ort(i)**3
      end do
      return

c H =  d/dx cot(x)  +  cot(x) d/dx

24    do i=1,dim
         work(i)=psi(i)*cos(ort(i))/sin(ort(i))
      end do
      call ddx (mode,h psi,work,dim,matrix,trafo)
      call ddx (mode,work,psi,dim,matrix,trafo)
      do i=1,dim
         hpsi(i)=hpsi(i)+work(i)*cos(ort(i))/sin(ort(i))
      end do
      return

c H = 1/sqrt( y(x) ) * d/dx * y(x) * d/dx * 1/sqrt( y(x) ) , y=x(phi)

25    do i=1,dim
         work(i)=psi(i)/sqrt(xphich4(ort(i)))
      end do
      call ddx (mode,h psi,work,dim,matrix,trafo)
      do i=1,dim
         work(i)=hpsi(i)*xphich4(ort(i))
      end do
      call ddx (mode,h psi,work,dim,matrix,trafo)
      do i=1,dim
         hpsi(i)=hpsi(i)/sqrt(xphich4(ort(i)))
      end do
      return

c H = 1/sqrt( y(x) ) * d/dx * y(x) * d/dx * 1/sqrt( y(x) ) , y=x(chi)

26    do i=1,dim
         work(i)=psi(i)/sqrt(xchich4(ort(i)))
      end do
      call ddx (mode,h psi,work,dim,matrix,trafo)
      do i=1,dim
         work(i)=hpsi(i)*xchich4(ort(i))
      end do
      call ddx (mode,h psi,work,dim,matrix,trafo)
      do i=1,dim
         hpsi(i)=hpsi(i)/sqrt(xchich4(ort(i)))
      end do
      return

c sqrt( x(phi) * d/dx * 1/sqrt( x(phi) )  -  1/sqrt( x(phi) * d/dx * sqrt( x(phi) )

27    do i=1,dim
         work(i)=psi(i)/sqrt(xphich4(ort(i)))
      end do
      call ddx (mode,h psi,work,dim,matrix,trafo)
      do i=1,dim
         hpsi(i)=hpsi(i)*sqrt(xphich4(ort(i)))
      end do
      do i=1,dim
         work(i)=sqrt(xphich4(ort(i)))*psi(i)
      end do
      call ddx (mode,work2,work,dim,matrix,trafo)
      do i=1,dim
         hpsi(i)=hpsi(i)-work2(i)/sqrt((xphich4(ort(i))))
      enddo
      return

c H = -1/2 * (d/dx)^2

30    call kin(mode,h psi,psi,dim,matrix,trafo)
      return

c H =  -1/2 x (d/dx)^2 x

31    do i=1,dim
         work(i)=ort(i)*psi(i)
      end do
      call kin(mode,h psi,work,dim,matrix,trafo)
      do i=1,dim
         hpsi(i)=hpsi(i)*ort(i)
      end do
      return


c H =  -1/2 x^2 (d/dx)^2 x^2

32    do i=1,dim
         work(i)=ort(i)**2*psi(i)
      end do
      call kin(mode,h psi,work,dim,matrix,trafo)
      do i=1,dim
         hpsi(i)=hpsi(i)*ort(i)**2
      end do
      return


c H =  -1/2 (d/dx)^2 x  -1/2 x (d/dx)^2

33    do i=1,dim
         work(i)=ort(i)*psi(i)
      end do
      call kin(mode,h psi,work,dim,matrix,trafo)
      call kin(mode,work,psi,dim,matrix,trafo)
      do i=1,dim
         hpsi(i)=hpsi(i)+work(i)*ort(i)
      end do
      return


c -----------------------------------------------------------------

c h=theta(r-r0)

201   r0=120.0d0
      do i=1,dim
         if (ort(i).lt.r0) then
            h psi(i)=0.0d0
         else
            h psi(i)=psi(i)
         endif
      enddo
      return

c [H,h(r-r0)]

202   r0=120.0d0
      call kin(mode,h psi,psi,dim,matrix,trafo)
      do i=1,dim
         if (ort(i).lt.r0) h psi(i)=0.0d0
      enddo

      do i=1,dim
         if (ort(i).lt.r0) then
            work(i)=0.0d0
         else
            work(i)=psi(i)
         endif
      enddo
      call kin(mode,work2,work,dim,matrix,trafo)
      do i=1,dim
         h psi(i)=work2(i)-h psi(i)
      enddo
      return

      end

c -------------------------------------------------------------------
      real*8 function xphich4(phi)
      implicit none

      real*8     phi
      real*8     pi,pi3
      parameter (pi=3.14159265358979323844d0,pi3=pi/3.0d0)

      xphich4=(1.0d0/( 1.0d0 +                   (phi-pi3)**2
     .                    - sqrt(3.0d0)/9.0d0*(phi-pi3)**3
     .                    +       3.0d0/4.0d0*(phi-pi3)**4  ))

      end
c -------------------------------------------------------------------
      real*8 function xchich4(chi)
      implicit none

      real*8     chi
      real*8     pi
      parameter (pi=3.14159265358979323844d0)

      xchich4=( 1.0d0/( 1.0d0 + ((chi-pi)**2)/3.0d0
     .                     + ((chi-pi)**4)/12.0d0 ))

      end
c -------------------------------------------------------------------

      subroutine IntToCartch4rst(q,x)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine transforms the internal coordinates       c
c to cartesian ones.                                        c
c Subroutine koordread must be called in main program       c
c GS / 9 Feb 2006                                           c
c rst, RST                                                  c
c GS / 25 July 2008                                         c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer    na,na3,dimq,dimch3
      parameter  (na=5,na3=3*na,dimq=na3-6,dimch3=12)

      real*8  q(dimq),x(na3)
      real*8  r(3),xcm(3)
      real*8  radau(dimch3)
      real*8  mch3,mch4
      real*8     hmass
      parameter (hmass=1836.109d0)

      real*8   rmch3,rmch4,srmch31,srmch41
      real*8   dist1,dist2,nenn1,nenn2

      integer i,j,k

      real*8 rctrans(4,4),m(na),sm1(na)
      common /rctrans/ rctrans,m,sm1

c mass of CH3-group
      mch3=0.0d0
      do i=1,4
         mch3=mch3+m(i)
      end do

c mass of CH4-group
      mch4=0.0d0
      do i=1,5
         mch4=mch4+m(i)
      end do

c reduced mass between CH3 and 4th H

      rmch3=mch3*m(5)/(mch3+m(5))
      srmch31=1.0d0/sqrt(rmch3)


c reduced mass between CH4 and 5th H

c      rmch4=mch4*m(6)/(mch4+m(6))
c      srmch41=1.0d0/sqrt(rmch4)


c order of q's:
c  1 rho          4 theta     7 r           10 R
c  2 theta_rho    5 phi       8 s           11 S
c  3 phi_rho      6 chi       9 t           12 T

c order of radau's
c  1 center of mass    4 x H_1     7 x H_2   10 x H_3
c  2 "                 5 y H_1     8 y H_2   11 y H_3
c  3 "                 6 z H_1     9 z H_2   12 z H_3

c order of cartesians
c 1-3 C (x,y,z)
c 4-6, 7-9, 10-12 methyl-H
c 13-15, 16-18 H2-H

c First step: calculate r1,r2,r3 from rho_2, phi_2, theta_2

      r(1) = q(1)*cos(q(2))
      r(2) = q(1)*sin(q(2))*cos(q(3))
      r(3) = q(1)*sin(q(2))*sin(q(3))

c Second step: calculate mass weighted Radau's
c              from r1,r2,r3,theta,phi,chi

c     print*,"r's in Ruecktransformation:"
c     print'(F10.5)',(r(i),i=1,3)

      radau( 1)=0.0d0
      radau( 2)=0.0d0
      radau( 3)=0.0d0

      radau( 4)=r(1)*sin(q(4))
      radau( 5)=0.0d0
      radau( 6)=r(1)*cos(q(4))

      radau( 7)=r(2)*sin(q(4))*cos(q(5)+q(6))
      radau( 8)=r(2)*sin(q(4))*sin(q(5)+q(6))
      radau( 9)=r(2)*cos(q(4))

      radau(10)=r(3)*sin(q(4))*cos(-q(5)+q(6))
      radau(11)=r(3)*sin(q(4))*sin(-q(5)+q(6))
      radau(12)=r(3)*cos(q(4))

c     print*,"Radau in Ruecktransformation"
c     do i=1,4
c        print'(3F12.6)',(radau(j+(i-1)*3),j=1,3)
c     end do

c Third step: calculate cartesians for methyl group
c             from Radau's
      do i=1,dimch3,3
         x(i  )=0.0d0
         x(i+1)=0.0d0
         x(i+2)=0.0d0
         do j=1,dimch3,3
           x(i  )=x(i  )+rctrans(i/3+1,j/3+1)*radau(j  )
           x(i+1)=x(i+1)+rctrans(i/3+1,j/3+1)*radau(j+1)
           x(i+2)=x(i+2)+rctrans(i/3+1,j/3+1)*radau(j+2)
         end do
      end do

c  unmassweight cartesians, correct sign
      do i=1,dimch3,3
         x(i  )=-x(i  )*sm1(i/3+1)
         x(i+1)=-x(i+1)*sm1(i/3+1)
         x(i+2)=-x(i+2)*sm1(i/3+1)
      end do

c Fourth step: add 4th hydrogen atom

c  unmassweight r
      dist1=q(7)*srmch31
      nenn1=1.0d0/(1.0d0+q(8)**2+q(9)**2)

c  calculate cartesians
      x(13)=2*dist1*q(8)*nenn1
      x(14)=2*dist1*q(9)*nenn1
      x(15)=dist1*(1.0d0-q(8)**2-q(9)**2)*nenn1

      end
c ------------------------------------------------------------------------
      subroutine koordreadch4rst
      implicit none

      integer    na,na3,dimch3
      parameter (na=5,na3=3*na,dimch3=12)

      real*8     mch3,hmass
      parameter (hmass=1836.109d0)

      integer i,j

      real*8 work(4,4),m(na),sm1(na),alpha(4)
      common /rctrans/ work,m,sm1

c masses
      m(1)=11.907d0
      m(2)=1.9984d0
      m(3)=1.9984d0
      m(4)=1.9984d0
      m(5)=1.0d0
!      m(1)=11.907d0
!      m(2)=1.0d0
!      m(3)=1.0d0
!      m(4)=1.0d0
!      m(5)=1.0d0
      do i=1,na
         m(i)=m(i)*hmass
         sm1(i)=1.0d0/sqrt(m(i))
      end do

      mch3=0.0d0
      do i=1,4
         mch3=mch3+m(i)
      end do

c calculate matrix for transformation between mass weighted
c cartesians and Radaus

      do i=1,4
         alpha(i)=sqrt(m(i)/mch3)
      end do

      do i=1,4
         work(i,1)=alpha(i)
         work(1,i)=alpha(i)
      end do

      do i=2,4
         do j=2,4
            work(j,i)=alpha(j)*alpha(i)/(alpha(1)+1.0d0)
         end do
      end do

      do i=2,4
         work(i,i)=work(i,i)-1.0d0
      end do

      end
c ------------------------------------------------------------------------
