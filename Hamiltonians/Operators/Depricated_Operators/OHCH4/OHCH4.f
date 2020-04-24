c######################################################################
c MCTDH-Modul Operatoren
c
c Contains subroutines "h init" and "h".
c "h init" initializes the arrays "coeff" and "diag" which
c specify the coefficients employed in the Hamiltonian series.
c "h" calculates the action of a single-particle operator (on the
c primitive grid). "mode" and "teil" define the number of the
c coordinate and of the term in the series, respectively.
c
c generic KEO
c DS May 2018
c######################################################################
        
        include 'set_coefficients.f'
        include 'set_diag.f'
        include 'get_index.f'


      subroutine h init(coeff,diag,dim)

      implicit none

      integer                           dim, number_of_coefficients
      integer                           N, k, maxN
      parameter                         (maxN=20)
      integer                           i, coeff_per_part(4:maxN)
      character*20                      coord(5:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      logical*4                         diag(3)
      complex*16                        coeff(*)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart
      SAVE      /coordinate_choice/
      integer   njx(5:maxN), njy(5:maxN), njz(5:maxN)
      common    /nj/                    njx, njy, njz
      SAVE      /nj/

c     TODO: adapt this number of atoms to your system!
      
        N = 7     !number of atoms in your system
        
        if (N.gt.maxN) then
        print*, 'ERROR: maximum number of atoms is 20'
        print*, 'yours is: ', N
        stop
        endif


c     TODO: adapt to your choice of coordinates for the atoms 5 to N
c       KEY             COORDINATES
c       stereo          stereographic coordinates (r,s,t)
c       spherical       spherical coordinates (r, theta, phi)


        coord(5) = 'stereo'
        coord(6) = 'stereo'
        coord(7) = 'spherical'

c      initalize potential
c       TODO: change the intialization routine according to your
c       potential
c       koordread is called later so don't call it here!
       call pes_init
        
c       generic part starts here
c       don't change anything of the rest of the file unless you want to
c       change the definition of the kinetic energy operator
        call koordread

        do i=5, N
                if (coord(i).eq.'stereo') then
                        print*, 'streographic coordinates for atom ', i
                else if (coord(i).eq.'spherical') then
                        print*, 'spherical coordinates for atom', i
                else
                        print*, 'ERROR: not supported: ', coord(i)
                        print*, 'coordinates for atom:', i
                        stop
                endif
        enddo
        

        number_of_coefficients = 0

c      Methyl-Potential
       call set_coefficients_methyl(coeff, number_of_coefficients)
       
       coeff_per_part(4) = number_of_coefficients

c      potential for atom 5 (no mixed terms between atoms not belonging
c      to ch3)

       call set_coefficients_atom(coeff, number_of_coefficients, 5)

       coeff_per_part(5) = number_of_coefficients - coeff_per_part(4)
c      potential for atom i (i=6 to N) (mixed terms to atoms 5 to i-1
c      included)
        
       do i = 6, N
                call set_coefficients_atom_mix(coeff,
     .                  number_of_coefficients, i)
                coeff_per_part(i) = number_of_coefficients
               do k = 4, i-1       
               coeff_per_part(i) = coeff_per_part(i) - coeff_per_part(k)
               enddo
       enddo

        coeff(number_of_coefficients+1)=0.d0

c Bestimmung der diagonalen Komponenten des Hamilton-Operators
      call diagon(diag,dim,number_of_coefficients)


      end

c----------------------------------------------------------------
      subroutine diagon(diag,dim,number_of_coefficients)

      implicit none

      integer     dim, number_of_coefficients, i, k, N, part
      logical*4   diag(dim,number_of_coefficients)


        do i = 1, number_of_coefficients
           do k= 1, dim
                diag(k,i)=.true.
           enddo
        enddo
        
      N = (dim - 6)/3 + 4


      call set_diag_methyl(diag, dim, i)
      call set_diag_atom(diag, dim, i, 5)

      do part = 6, N
        call set_diag_atom_mix(diag, dim, i, part)
      enddo

      end


c#####################################################################

      subroutine h (mode,teilin,h psi,psi,dim,matrix,trafo,ort)

      implicit none

      integer       mode, teilin, dim, i
      complex*16    h psi(dim), psi(dim)
c      complex*16    work(1024),work2(1024),work3(1024)
      complex*16    work(dim), work2(dim), work3(dim)
      real*8        matrix(dim,dim), trafo(dim,dim), ort(dim)

      real*8        xphi,xchi
      REAL*8        R0

      integer       term, teil, part, index

      integer                           N, maxN
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      character*20                      coord(5:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart

      integer         sys(4)
      real*8          para(4,15)
c     complex*16      statpsi(maxdim), faktor
      complex*16      statpsi(dim), faktor
      common /sys/    sys, para

      complex*16       tphase 
      real*8           genau, CDVRgenau, SDgenau
      integer          ireal
      logical*4        save, fixed, matcalc, bottomLayer_, TDDVRinit
      logical*4        h2alloc, h2calc, SPFinit, ldummy
      common /steuer/  tphase, genau, save, fixed, matcalc, bottomLayer_,
     .                 CDVRgenau, SDgenau, TDDVRinit, ireal,
     .                 h2alloc, h2calc, SPFinit, ldummy


! h=theta(r-r0)  dividing surface
      if ((teilin.eq.-1)) then
        r0=700.0 
        if (mode.eq.10) then
          do i=1,dim
            if (ort(i).lt.r0) then
              h psi(i)=0.0d0
            else
              h psi(i)=psi(i)
            endif
          enddo
        else
          do i=1,dim
            h psi(i)=psi(i)
          enddo
        endif
        return
      endif

! [H,h(r-r0)] flux operator
      if ((teilin.eq.-2)) then
        if (mode.eq.10) then
          r0=700.0 
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
        else if (mode.eq.7) then
          r0= 120.0
          do i=1,dim
            if (ort(i).gt.r0) then
              h psi(i)=0.0d0
            else
              h psi(i)=psi(i)
            endif
          enddo
        else
          do i=1,dim
            hpsi(i)=psi(i)
          enddo
        endif
        return
      endif

! statistical sampling
      if ((teilin.eq.-19)) then
         if (mode.eq.15) then
           call genstatist (statpsi, ireal, dim, trafo,
     .                      2, para(1,mode), 1d99)
           faktor=(0.0d0,0.0d0)
           do i=1,dim
              faktor=faktor+dconjg(statpsi(i))*psi(i)
           enddo
c           if (abs(faktor/real(faktor)-1).gt.1d-4) then
c              write (6,*) 'wavefunction not real'
c              stop
c           endif
           do i=1,dim
              h psi(i)=real(faktor)*statpsi(i)
           enddo
        else
          do i=1,dim
            hpsi(i)=psi(i)
          enddo
        endif
        return
      endif

! Here continues as Daniela's code
        i = 4
        teil = teilin

        do while (teil.gt.0)
                teil = teil-coeff_per_part(i)
                i = i +1
        enddo

        part = i-1
        term = teil + coeff_per_part(part)

        if (part.eq.4) then
                call get_index_methyl(index, term, mode)
        else if (part.eq.5) then
                call get_index_atom(index, term, mode, part)
        else if (part.gt.5) then
                call get_index_atom_mix(index, term, mode, part)
        else
                stop
        endif
        
c        print*, 'Teil: ', teilin, ' Mode: ', mode, ' Index: ', index
        

        SELECT CASE(index)
                CASE( 2) ! H = 1/x**2
                        do i=1,dim
                                h psi(i)=psi(i)/ort(i)**2
                        end do
                        return
                CASE( 3) ! H = sin(x)
                        do i=1,dim
                                h psi(i)=psi(i)*sin(ort(i))
                        end do
                        return
                CASE( 4) ! H = cos(x)
                        do i=1,dim
                                h psi(i)=psi(i)*cos(ort(i))
                        end do
                        return
                CASE( 5) ! H = cot(x)
                        do i=1,dim
                                h psi(i)=psi(i)*cos(ort(i))/sin(ort(i))
                        end do
                        return
                CASE( 6) ! H = sin(x)**2
                        do i=1,dim
                                h psi(i)=psi(i)*sin(ort(i))**2
                        end do
                        return
                CASE( 7) ! H = cos(x)**2
                        do i=1,dim
                                h psi(i)=psi(i)*cos(ort(i))**2
                        end do
                        return
                CASE( 8) ! H = cot(x)**2
                        do i=1,dim
                                h psi(i)=psi(i)/tan(ort(i))**2
                        end do
                        return
                CASE( 9) ! H = cos(x)**2 - sin(x)**2
                        do i=1,dim
                        h psi(i)=psi(i)*(cos(ort(i))**2-sin(ort(i))**2)
                        end do
                        return
                CASE(10) ! H = 1/sin(x)**2
                        do i=1,dim
                                h psi(i)=psi(i)/sin(ort(i))**2
                        end do
                        return
                CASE(11) ! H = 5 - 3/sin(x)**2
                        do i=1,dim
                        h psi(i)=psi(i)*(1.0d0+1.0d0/sin(ort(i))**2)
                        end do
                        return
                CASE(12) ! H = 3 - 1/sin(x)**2
                        do i=1,dim
                        h psi(i)=psi(i)*(3.0d0-1.0d0/sin(ort(i))**2)
                        end do
                        return
                CASE(13) ! H = 1/sin(x)**2 - 3/2
                        do i=1,dim
                        h psi(i)=psi(i)*(1.0d0/sin(ort(i))**2-1.5d0)
                        end do
                        return
                CASE(14) ! H = d/dx 
                        call ddx (mode,h psi,psi,dim,matrix,trafo)
                        return
                CASE(15) ! H = -1/2 * (d/dx)**2
                        call kin(mode,h psi,psi,dim,matrix,trafo)
                        return
                CASE(16) ! H = cot(x) d/dx + 1/sin(x) d/dx cos(x)
                         ! Legendre
                        do i=1,dim
                                work(i)=psi(i)*cos(ort(i))/sin(ort(i))
                        end do
                        call ddx (mode,work2,work,dim,matrix,trafo)
                        do i=1,dim
                                work2(i)=work2(i)/sin(ort(i))
                        end do
                        do i=1,dim
                                work(i)=psi(i)/sin(ort(i))
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        do i=1,dim
                        hpsi(i)=work2(i)+hpsi(i)*cos(ort(i))/sin(ort(i))
                        end do
                        return
                CASE(17) ! H = 1/sin(x) d/dx sin(x)  +  d/dx
                         !  Legendre 
                        do i=1,dim
                                work(i)=psi(i)/sin(ort(i))
                        end do
                        call ddx (mode,h psi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)/sin(ort(i))
                        end do
                        return
                CASE(18) ! H = d/dx sin(x)  +  sin(x) d/dx
                        do i=1,dim
                                work(i)=sin(ort(i))*psi(i)
                        end do
                        call ddx (mode,h psi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)*sin(ort(i))
                        end do
                        return
                CASE(19) ! H = d/dx cos(x)  +  cos(x) d/dx
                        do i=1,dim
                                work(i)=cos(ort(i))*psi(i)
                        end do
                        call ddx (mode,h psi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)*cos(ort(i))
                        end do
                        return
                CASE(20) ! H =  d/dx cot(x)  +  cot(x) d/dx
                        do i=1,dim
                                work(i)=psi(i)*cos(ort(i))/sin(ort(i))
                        end do
                        call ddx (mode,h psi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                        hpsi(i)=hpsi(i)+work(i)*cos(ort(i))/sin(ort(i))
                        end do
                        return
                CASE(21) ! H = d/dx sin(x) cos(x) + sin(x) cos(x) d/dx
                        do i=1,dim
                                work(i)=cos(ort(i))*sin(ort(i))*psi(i)
                        end do
                        call ddx (mode,h psi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                        hpsi(i)=hpsi(i)+work(i)*cos(ort(i))*sin(ort(i))
                        end do
                        return
                CASE(22) ! H = 1/sqrt( y(x) ) * d/dx * y(x) * d/dx * 1/sqrt( y(x) ) , y=x(phi)
                        do i=1,dim
                                work(i)=psi(i)/sqrt(xphi(ort(i)))
                        end do
                        call ddx (mode,h psi,work,dim,matrix,trafo)
                        do i=1,dim
                                work(i)=hpsi(i)*xphi(ort(i))
                        end do
                        call ddx (mode,h psi,work,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)/sqrt(xphi(ort(i)))
                        end do
                        return 
                CASE(23) ! H = 1/sqrt( y(x) ) * d/dx * y(x) * d/dx * 1/sqrt( y(x) ) , y=x(chi)
                        do i=1,dim
                                work(i)=psi(i)/sqrt(xchi(ort(i)))
                        end do
                        call ddx (mode,h psi,work,dim,matrix,trafo)
                        do i=1,dim
                                work(i)=hpsi(i)*xchi(ort(i))
                        end do
                        call ddx (mode,h psi,work,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)/sqrt(xchi(ort(i)))
                        end do
                        return 
                CASE(24) ! H =  1/sin(x) d/dx sin(x) d/dx  + 2 cot(x)**2  - 1/sin(x) (d/dx)**2 sin(x)  - (d/dx)**2
                         ! Legendre
                    call kin(mode,h psi,psi,dim,matrix,trafo)
                    do i=1,dim
                    hpsi(i)=-2.0d0*hpsi(i)+2.0d0*psi(i)/(tan(ort(i))**2)
                    end do
                    call ddx(mode,work,psi,dim,matrix,trafo)
                    do i=1,dim
                        work(i)=work(i)/sin(ort(i))
                    end do
                    call ddx(mode,work2,work,dim,matrix,trafo)
                    do i=1,dim
                        hpsi(i)=hpsi(i)-work2(i)/sin(ort(i))
                    end do
                    do i=1,dim
                        work(i)=psi(i)/sin(ort(i))
                    end do
                    call ddx(mode,work2,work,dim,matrix,trafo)
                    do i=1,dim
                        work(i)=work2(i)/sin(ort(i))
                    end do
                    call ddx(mode,work2,work,dim,matrix,trafo)
                    do i=1,dim
                        hpsi(i)=hpsi(i)-work2(i)
                    end do
                    return    
                CASE(25) ! H =  -1/2 (d/dx)**2  +  1/2 (d/dx)**2 sin(x)**2  +  1/2 sin(x)**2 (d/dx)**2
                        do i=1,dim
                                work(i)=-sin(ort(i))**2*psi(i)
                        end do
                        call kin(mode,h psi,work,dim,matrix,trafo)
                        call kin(mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)-work(i)*sin(ort(i))**2
                        end do
                        call kin(mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)
                        end do
                        return
                CASE(26) ! H =  -1/2 (d/dx)**2 cos(x)  -1/4 cos(x)  -1/2 cos(x) (d/dx)**2
                        do i=1,dim
                                work(i)=cos(ort(i))*psi(i)
                        end do
                        call kin(mode,h psi,work,dim,matrix,trafo)
                        call kin(mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)*cos(ort(i))
     .                                  -0.25d0*cos(ort(i))*psi(i)
                        end do
                        return
               CASE(27) ! H = sqrt( x(phi) * d/dx * 1/sqrt( x(phi) )  -  1/sqrt( x(phi) * d/dx * sqrt( x(phi) ) 
                        do i=1,dim
                                work(i)=psi(i)/sqrt(xphi(ort(i)))
                        end do
                        call ddx (mode,h psi,work,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)*sqrt(xphi(ort(i)))
                        end do
                        do i=1,dim
                                work(i)=sqrt(xphi(ort(i)))*psi(i)
                        end do
                        call ddx (mode,work2,work,dim,matrix,trafo)
                        do i=1,dim
                           hpsi(i)=hpsi(i)-work2(i)/sqrt(xphi(ort(i)))
                        enddo
                        return
                CASE(28) ! H = x
                        do i=1,dim
                                h psi(i)=psi(i)*ort(i)
                        enddo
                        return
                CASE(29) ! H = x^2
                        do i=1,dim
                                h psi(i)=psi(i)*ort(i)*ort(i)
                        enddo
                        return
                CASE(30) ! H = (1-x^2)
                        do i=1,dim
                                h psi(i)=psi(i)*(-ort(i)*ort(i)+1.0d0)
                        enddo
                        return
                CASE(31) ! H = (1-x^2)^2
                        do i=1,dim
                                h psi(i)=psi(i)*(ort(i)**2-1.0d0)**2
                        enddo
                        return
                CASE(32) ! H = x*(1-x^2)
                        do i=1,dim
                           h psi(i)=-psi(i)*(ort(i)*ort(i)-1.0d0)*ort(i)
                        enddo
                        return
                CASE(33) ! H = 1+x^2
                        do i=1,dim
                                h psi(i)=psi(i)*(ort(i)*ort(i)+1.0d0)
                        enddo
                        return

                CASE(35) ! H = (1+x^2)^2
                        do i=1,dim
                                h psi(i)=psi(i)*(ort(i)**2+1.0d0)**2
                        enddo
                        return
                CASE(36) ! H = x d/dx x
                        do i=1,dim
                                work(i)=ort(i)*psi(i)
                        end do
                        call ddx (mode,h psi,work,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)*ort(i)
                        end do
                        return
                CASE(37) ! H =  d/dx x  +  x d/dx
                        do i=1,dim
                                work(i)=psi(i)*ort(i)
                        end do
                        call ddx (mode,h psi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)*ort(i)
                        end do
                        return
                CASE(38) ! H =  d/dx x^3  +  x^3 d/dx
                        do i=1,dim
                                work(i)=psi(i)*ort(i)**3
                        end do
                        call ddx (mode,h psi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)*ort(i)**3
                        end do
                        return
                CASE(39) ! H =  -1/2 x (d/dx)^2 x
                        do i=1,dim
                                work(i)=ort(i)*psi(i)
                        end do
                        call kin(mode,h psi,work,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)*ort(i)
                        end do
                        return
                CASE(40) ! H =  -1/2 x^2 (d/dx)^2 x^2
                        do i=1,dim
                                work(i)=ort(i)**2*psi(i)
                        end do
                        call kin(mode,h psi,work,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)*ort(i)**2
                        end do
                        return
                CASE(41) ! H =  -1/2 (d/dx)^2 x  -1/2 x (d/dx)^2
                        do i=1,dim
                                work(i)=ort(i)*psi(i)
                        end do
                        call kin(mode,h psi,work,dim,matrix,trafo)
                        call kin(mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)*ort(i)
                        end do
                        return
                CASE DEFAULT
                        print*, 'ERROR: unknown index in h'
                        print*, 'teilin: ', teilin
                        print*, 'part: ', part
                        print*, 'term: ', term
                        print*, 'mode: ', mode
                        print*, 'index: ', index
                        stop
        END SELECT
                        


      end subroutine
c -------------------------------------------------------------------
      real*8 function xphi(phi)
      implicit none
 
      real*8     phi
      real*8     pi,pi3
      parameter (pi=3.14159265358979323844d0,pi3=pi/3.0d0)

      xphi=(1.0d0/( 1.0d0 +                   (phi-pi3)**2
     .                    - sqrt(3.0d0)/9.0d0*(phi-pi3)**3
     .                    +       3.0d0/4.0d0*(phi-pi3)**4  ))

      end 
c -------------------------------------------------------------------
      real*8 function xchi(chi)
      implicit none

      real*8     chi
      real*8     pi
      parameter (pi=3.14159265358979323844d0)
 
      xchi=( 1.0d0/( 1.0d0 + ((chi-pi)**2)/3.0d0
     .                     + ((chi-pi)**4)/12.0d0 ))

      end
c-------------------------------------------------------
        subroutine print_error(name, error, term, mode)

        implicit none
        integer         term, mode
        character*30    name,error

        print*, 'ERROR in KEO: ', name
        print*, error
        print*, 'term: ', term
        print*, 'mode: ', mode

        stop

        end subroutine
