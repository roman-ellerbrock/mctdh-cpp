c *********************************************************************
c generic kinetic energy opearator for N atoms (e.g. for CH(N-1) )
c
c built up_ from three-atomic system (r1,r2,theta) defining molecular
c frame (theta -> Hermite-DVR, FFT)
c
c subsequent atoms described in stereographic projection coordinates
c (ri,si,ti) with north or south being the pole of projection
c
c last coordinate refers to rotational part (J.ne.0)
c *********************************************************************

      subroutine hinitch4rad(koeff,diag)
      use iso_c_binding
      implicit none
      
      integer dim, k zahl, N
      integer i,j,k,l
      parameter (N=5)
      integer sys(4)
      integer vect(N-3)
      real*8 para(4,3*N-6)
      integer diag(3)
      real*8 mu(N-1),par_,si,sj
      complex*16 hilf,hilf1
      complex*16 koeff(*)
      
      common /sys/ sys,para

      dim=3*N-6
      
c      call potini
c      call koordread

c decide whether to project the N-3 vectors from north (.ne.1) or south (.eq.1)

      vect(1)=0
      vect(2)=0
	if (N.gt.2) then
	      vect(3)=0
	endif
      
      k zahl=4+34*(N-3)+38*(N-3)*(N-4)+15*(N-3)+1+4
      
      do i=1,N-1
         mu(i)=1.0d0
      enddo
      
      koeff(1)=-0.5d0/mu(1)
      koeff(2)=-0.5d0/mu(2)
      koeff(3)=-0.5d0/mu(1)
      koeff(4)=-0.5d0/mu(2)
      
c     (d/dr)**2      
      
      l=4
      do i=1,N-3
         koeff(l+i)=-0.5d0/mu(i+2)
      enddo
      
      l=l+(N-3)
      hilf=-0.25d0*0.5d0

c     L**2
      
      do i=1,N-3
         do j=1,6
             koeff(l+6*(i-1)+j)=hilf/mu(i+2)
         enddo
      enddo

      l=l+(N-3)*6
      hilf=-0.25d0*0.5d0/mu(1)
      
c     Lx**2+Ly**2      
      
      do i=1,N-3
         do j=1,6
            koeff(l+10*(i-1)+j)=2.0d0*hilf
         enddo
         do j=7,10
            koeff(l+10*(i-1)+j)=hilf
         enddo
      enddo

      l=l+(N-3)*10
      hilf=-0.5d0*0.5d0/mu(1)

c     cot(theta)**2*Lz**2

      do i=1,N-3
         koeff(l+4*(i-1)+1)=hilf
         koeff(l+4*(i-1)+2)=2.0d0*hilf
         koeff(l+4*(i-1)+3)=2.0d0*hilf
         koeff(l+4*(i-1)+4)=-hilf
      enddo

      l=l+(N-3)*4
      hilf=-0.5d0*0.5d0/mu(2)

c     csc(theta)**2*Lz**2

      do i=1,N-3
         koeff(l+4*(i-1)+1)=hilf
         koeff(l+4*(i-1)+2)=2.0d0*hilf
         koeff(l+4*(i-1)+3)=2.0d0*hilf
         koeff(l+4*(i-1)+4)=-hilf
      enddo
      
      l=l+(N-3)*4
      hilf=0.5d0*0.5d0/mu(1)

c modified for variable projection pole, term 1

c     Lx*Lz+Lz*Lx

      do i=1,N-3
         si=1.0d0
       if (vect(i).eq.1) si=-1.0d0
       koeff(l+6*(i-1)+1)=hilf*si
         koeff(l+6*(i-1)+2)=-2.0d0*hilf*si
         koeff(l+6*(i-1)+3)=2.0d0*hilf*si
         koeff(l+6*(i-1)+4)=-2.0d0*hilf*si
         do j=5,6
            koeff(l+6*(i-1)+j)=hilf*si
         enddo
      enddo

      
      l=l+(N-3)*6
      hilf=0.5d0/mu(1)

c     Ly      
      
      do i=1,N-3
         si=1.0d0
       if (vect(i).eq.1) si=-1.0d0
       do j=1,3
            koeff(l+3*(i-1)+j)=hilf*si
         enddo
      enddo
      
      l=l+(N-3)*3
      hilf=-0.125d0/mu(1)

c     off-diagonal Elemente

c modified for variable projection pole, term 2

c     Li,x*Lj,x+Li,y*Lj,y   

      do i=1,N-3
         do j=1,N-3
          if (i.ne.j) then
             si=1.0d0
             sj=1.0d0
             if (vect(i).eq.1) si=-1.0d0
             if (vect(j).eq.1) sj=-1.0d0
             do k=1,18
                koeff(l+k)=hilf*si*sj
             enddo
             l=l+18
          endif
       enddo
      enddo
      
      hilf=-0.5d0/mu(1)

c     cot(theta)**2*Li,z*Lj,z

      do i=1,(N-3)*(N-4)
         koeff(l+4*(i-1)+1)=hilf
         koeff(l+4*(i-1)+2)=-hilf
         koeff(l+4*(i-1)+3)=-hilf
         koeff(l+4*(i-1)+4)=hilf
      enddo

      l=l+(N-3)*(N-4)*4
      hilf=-0.5d0/mu(2)

c     csc(theta)**2*Li,z*Lj,z      
      
      do i=1,(N-3)*(N-4)
         koeff(l+4*(i-1)+1)=hilf
         koeff(l+4*(i-1)+2)=-hilf
         koeff(l+4*(i-1)+3)=-hilf
         koeff(l+4*(i-1)+4)=hilf
      enddo

      l=l+(N-3)*(N-4)*4 
      hilf=-0.5d0*0.5d0/mu(1)

c modified for variable projection pole, term 3

c     Li,x*Lj,z+Li,z*Lj,x

      do i=1,N-3
         do j=1,N-3
            if (i.ne.j) then
             par_=1.d0
             sj=1.d0
             if (vect(i).eq.1) sj=-1.d0
             do k=1,6
                  koeff(l+k)=par_*hilf*sj
                par_=-par_
             enddo
               si=1.d0
             if (vect(j).eq.1) si=-1.d0
             do k=7,12
                koeff(l+k)=par_*hilf*si
              par_=-par_
             enddo
             l=l+12         
          endif
       enddo
      enddo
      
      hilf=(0.0d0, 0.5d0)/mu(1)*(0.5d0, 0.d0)

c ------- Koeffizienten fuer J-linear ------- 

c modified for variable projection pole, term 4

c Jx*Li,x

      do i=1,N-3
         si=1.d0
       if (vect(i).eq.1) si=-1.d0
       do j=1,3
            koeff(l+3*(i-1)+j)=hilf*si
         enddo
      enddo
      
      l=l+(N-3)*3
      hilf=(0.0d0, -0.5d0)/mu(1)*(0.0d0, -0.5d0)
c J+ = Jx-iJy, J- = Jx+iJy in body-fixed frame !      
      hilf=-hilf

c modified for variable projection pole, term 5

c Jy*Li,y 

      do i=1,N-3
         si=1.d0
       if (vect(i).eq.1) si=-1.d0
       do j=1,3
            koeff(l+3*(i-1)+j)=hilf*si
         enddo
      enddo

      l=l+(N-3)*3
      hilf=(0.0d0, 1.0d0)/mu(1)

c cot(theta)**2*Jz*Li,z

      do i=1,N-3
         koeff(l+2*(i-1)+1)=hilf
         koeff(l+2*(i-1)+2)=-1.0d0*hilf
      enddo

      l=l+(N-3)*2
      hilf=(0.0d0, 1.0d0)/mu(2)

c csc(theta)**2*Jz*Li,z

      do i=1,N-3
         koeff(l+2*(i-1)+1)=hilf
         koeff(l+2*(i-1)+2)=-1.0d0*hilf
      enddo
      
      l=l+(N-3)*2
      hilf=(0.0d0, 1.0d0)/mu(1)*(0.5d0, 0.0d0)

c Jx*Li,z

      do i=1,N-3
         koeff(l+2*(i-1)+1)=hilf
         koeff(l+2*(i-1)+2)=-1.0d0*hilf
      enddo

      l=l+(N-3)*2
      hilf=(0.0d0, 0.5d0)/mu(1)

c Jz*Li,x

c modified for variable projection pole, term 6

      do i=1,N-3
         si=1.d0
       if (vect(i).eq.1) si=-1.d0
       do j=1,3
            koeff(l+3*(i-1)+j)=hilf*si
         enddo
      enddo

      l=l+(N-3)*3
      
c Jy      

      koeff(l+1)=(0.0d0, 1.0d0)/mu(1)*(0.0d0, -0.5d0)
c J+ = Jx-iJy, J- = Jx+iJy in body-fixed frame !     
      koeff(l+1)=-koeff(l+1)

      l=l+1
      hilf=0.5d0/mu(1)

c ------- Koeffizienten für J-quadratisch -------
      
      koeff(l+1)=hilf
      koeff(l+2)=hilf
      koeff(l+3)=0.5d0/mu(2)
      koeff(l+4)=hilf

c ----------- must be used for J=0 -------------      

      k zahl=k zahl-15*(N-3)-4-1

c ---------------------------------------------
      
      koeff (k zahl+1)=0
      
c ------ diagonale Komponenten des Hamilton-Operators -------

      call diagon ch4rad (diag, dim, k zahl)
      
      return
      end
      
c ------------------------------------------------

      subroutine diagon ch4rad (diag, dim, k zahl)
      
      integer dim, k zahl, i, k
      integer diag(dim, k zahl)

      do k=1,k zahl
         do i=1,dim
!            diag(i,k)=.false.
            diag(i,k)=0
         enddo
      enddo
      return
      end
c -------------------------------------------------

      subroutine hch4rad(mode, teil, h psi, psi, dim, matrix, trafo,
     .  ort)
      use iso_c_binding
      
      implicit none
      
      integer   maxdim
      parameter (maxdim=1024)


      integer    mode, teil, dim, N
      integer    jmin, jmax, kmin, kmax
      integer    i,j,k,l,u
      integer    teildiag, teiloffdiag, teilmax, modemax
      integer    teiljlin
      integer    term, vec
      integer    blan, blanoff, blanjlin, jquad
      integer    rz, sz, tz
      parameter  (blan=7, blanoff=4, blanjlin=7, jquad=4)
      integer    block(blan), blockoff(blanoff), blockjlin(blanjlin)
      logical*4  rmode, smode, tmode
      parameter  (N=5)
      
      complex*16 h psi(dim), psi(dim), work(maxdim)
      real*8     matrix(dim,dim,*), trafo(dim,dim), ort(dim)
      real*8     para(4,3*N-6),wurz
      integer    sys(4)

      common /sys/ sys, para

c ---------------------------------------      
c      if (mode.eq.3) then
c        do i=1, dim
c          ort(i)=acos(ort(i))
c        enddo
c      endif
          
      block(1)=1
      block(2)=6
      block(3)=10
      block(4)=4
      block(5)=4
      block(6)=6
      block(7)=3
      
      blockoff(1)=18
      blockoff(2)=4
      blockoff(3)=4
      blockoff(4)=12
      
      blockjlin(1)=3
      blockjlin(2)=3
      blockjlin(3)=2
      blockjlin(4)=2
      blockjlin(5)=2
      blockjlin(6)=3
      blockjlin(7)=1
      
      teildiag=0
   
      do i=1,blan
         teildiag=teildiag+(N-3)*block(i)
      enddo
      
      teiloffdiag=0

      do i=1,blanoff
         teiloffdiag=teiloffdiag+(N-3)*(N-4)*blockoff(i)
      enddo
      
      teiljlin=0

      do i=1,blanjlin-1
         teiljlin=teiljlin+(N-3)*blockjlin(i)
      enddo
            
      teilmax=4+teildiag+teiloffdiag+teiljlin+blockjlin(7)+jquad
      modemax=3*N-5

c    ----------- must be used for J=0 -------

      teilmax=teilmax-teiljlin-blockjlin(7)-jquad
      modemax=modemax-1
      
c    ----------------------------------------

      if ((mode.lt.0).or.(mode.gt.modemax)) then
         print*, 'mode out of range'
         read(5, *)
         stop
      endif

      if ((teil.lt.0).or.(teil.gt.teilmax)) then
         print*, 'teil out of range'
         read(5, *)
         stop
      endif
      
      if (N.lt.3) then
         print*, 'N out of range'
         read(5, *)
         stop
      endif

c     Fehler korrigiert (23.03.09): .le. statt .lt. !!

      rmode=(mod(mode-4,3).eq.0).and.(mode.le.modemax).and.(mode.gt.3)
      smode=(mod(mode-5,3).eq.0).and.(mode.le.modemax).and.(mode.gt.3)
      tmode=(mod(mode-6,3).eq.0).and.(mode.le.modemax).and.(mode.gt.3)
      
      rz=int((mode-4)/3)+3
      sz=int((mode-5)/3)+3
      tz=int((mode-6)/3)+3

c -------- Teile 1-4 im KEO --------------

      if (teil.lt.5) then
         if (mode.lt.4) then
            goto (200,201,202) mode
         read(5, *)
            stop

 200        goto (10, 1, 4, 1) teil
         read(5, *)
            stop

 201        goto (1, 10, 1, 4) teil
         read(5, *)
            stop

c     ------- 12 --> HO-DVR, 35 --> Legendre-DVR ----------

c202         goto (1, 1, 35, 35) teil
202        goto (1, 1, 12, 12) teil
         read(5, *)
            stop
         else
            goto 1
         read(5, *)
            stop
         endif
      endif

c --------- Diagonalterme ----------------

         if (teil.lt.(teildiag+5)) then
            u=teil-4   
            i=0
            do while (u.gt.0)
               i=i+1
               u=u-(N-3)*block(i)
            enddo           
            
            u=u+(N-3)*block(i)
            term=mod(u,block(i))
            vec=int(u/block(i))+3
            
            if (term.eq.0) then
               term=term+block(i)
               vec=vec-1
            endif
        
            goto (210, 211, 212, 213, 214, 215, 216) i
         read(5, *)
            stop 'bei i'

c (d/dri)**2

 210        if ((rmode).and.(vec.eq.rz)) then
               goto 10
         read(5, *)
               stop 'bei (d/dri)**2'
            else
               goto 1
         read(5, *)
               stop '2'
            endif

c (Li)**2

 211        if ((rmode).and.(vec.eq.rz)) then
               goto 4
         read(5, *)
               stop '3'
            else if((smode).and.(vec.eq.sz)) then
               goto (10, 14, 15, 13, 16, 1) term
         read(5, *)
               stop '4'
            else if ((tmode).and.(vec.eq.tz)) then
               goto (14, 10, 13, 15, 1, 16) term
         read(5, *)
               stop '5'
            else
               goto 1
         read(5, *)
               stop '6'
            endif

c (Li,x)**2+(Li,y)**2

 212        if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '7'
            else if ((smode).and.(vec.eq.sz)) then
               goto (1, 3, 1, 17, 13, 11, 19, 10, 16, 1) term
         read(5, *)
               stop '8'
            else if ((tmode).and.(vec.eq.tz)) then
               goto (1, 1, 3, 13, 17, 11, 10, 19, 1, 16) term
         read(5, *)
               stop '9'
            else
               goto 1
         read(5, *)
               stop '10'
            endif

c (cot th )**2 (Li,z)**2

 213        if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '11'
            else if (mode.eq.3) then
               goto 20
         read(5, *)
               stop '12'
            else if ((smode).and.(vec.eq.sz)) then
               goto (1, 10, 3, 11) term
         read(5, *)
               stop '13'
            else if ((tmode).and.(vec.eq.tz)) then
               goto (1, 3, 10, 11) term
         read(5, *)
               stop '14'
            else
               goto 1
         read(5, *)
               stop '15'
            endif

c (csc th)**2 (Li,z)**2

 214        if (mode.eq.2) then 
               goto 4
         read(5, *)
               stop '16'
            else if (mode.eq.3) then
               goto 21
         read(5, *)
               stop '17'
            else if ((smode).and.(vec.eq.sz)) then
               goto (1, 10, 3, 11) term
         read(5, *)
               stop '18'
            else if ((tmode).and.(vec.eq.tz)) then
               goto (1, 3, 10, 11) term
         read(5, *)
               stop '19'
            else
               goto 1
         read(5, *)
               stop '20'
            endif

c Li,x * Li,z + Li,z * Li,x

 215        if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '21'
            else if (mode.eq.3) then
               goto 7
         read(5, *)
               stop '22'
            else if ((smode).and.(vec.eq.sz)) then
               goto (2, 22, 23, 2, 25, 9) term
         read(5, *)
               stop '23'
            else if ((tmode).and.(vec.eq.tz)) then
               goto (1, 10, 3, 17, 11, 24) term
         read(5, *)
               stop '24'
            else
               goto 1
         read(5, *)
               stop '25'
            endif

c Li,y

 216        if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '26'
            else if (mode.eq.3) then
c ------------ 9 --> HO-DVR, 34 --> Legendre-DVR ------------ 
c             goto 34
            goto 9
         read(5, *)
               stop '27'
            else if ((smode).and.(vec.eq.sz)) then
               goto (2, 9, 26) term
         read(5, *)
               stop '28'
            else if ((tmode).and.(vec.eq.tz)) then 
               goto (11, 18, 1) term
         read(5, *)
               stop '29'
            else
               goto 1
         read(5, *)
               stop '30'
            endif

         endif

 
c -------------- off-Diagonalterme ---------------------

         if (teil.lt.(teildiag+teiloffdiag+5)) then
            u=teil-teildiag-4
            i=0
            do while(u.gt.0) 
               i=i+1
               u=u-(N-3)*(N-4)*blockoff(i)
            enddo
            u=u+(N-3)*(N-4)*blockoff(i)
            term=mod(u,blockoff(i))
            
            l=i

            i=int(u/(blockoff(l)*(N-4)))+3
            if (mod(u,blockoff(l)*(N-4)).eq.0) i=i-1
            
            vec=int(u/blockoff(l))-(i-3)*(N-4)+1
           
            if (term.eq.0) then
               term=term+blockoff(l)
               vec=vec-1
            endif

            j=vec+2
            if (j.ge.i) j=j+1
            
            goto (221, 222, 223, 224) l
         read(5, *)
            stop 'bei l'


c Li,x Lj,x + Li,y Lj,y

 221        if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '31'
            else if (smode.and.((i.eq.sz).or.(j.eq.sz))) then
               if (i.eq.sz) then
                  goto (11, 11, 11, 18, 18, 18,  1,  1,  1,
     .                   2,  2,  2,  9,  9,  9, 26, 26, 26) term
         read(5, *)
                  stop '32'
               endif
               if (j.eq.sz) then
                  goto (11, 18,  1, 11, 18,  1, 11, 18,  1,
     .                   2,  9, 26,  2,  9, 26,  2,  9, 26) term
         read(5, *)
                  stop '33'
               endif
            else if (tmode.and.((i.eq.tz).or.(j.eq.tz))) then
               if (i.eq.tz) then
                  goto ( 2,  2,  2,  9,  9,  9, 26, 26, 26,
     .                  11, 11, 11, 18, 18, 18,  1,  1,  1) term 
         read(5, *)
                  stop '34'
               endif
               if (j.eq.tz) then
                  goto ( 2,  9, 26,  2,  9, 26,  2,  9, 26,
     .                  11, 18,  1, 11, 18,  1, 11, 18,  1) term
         read(5, *)
                  stop '35'
               endif
            else
               goto 1
         read(5, *)
               stop '36'
            endif

c (cot th)**2 Li,z Lj,z

 222       if (mode.eq.1) then
              goto 4
         read(5, *)
              stop '37'
           else if (mode.eq.3) then
              goto 20
         read(5, *)
              stop '38'
           else if (smode.and.((i.eq.sz).or.(j.eq.sz))) then 
              if (i.eq.sz) then
                 goto (2, 2, 9, 9) term
                 print*, 'Term', term
         read(5, *)
                 stop '39'
              endif
              if (j.eq.sz) then
                 goto (2, 9, 2, 9) term
         read(5, *)
                 stop '40'
              endif
           else if (tmode.and.((i.eq.tz).or.(j.eq.tz))) then
              if (i.eq.tz) then
                 goto (9, 9, 2, 2) term
         read(5, *)
                 stop '41'
              endif
              if (j.eq.tz) then
                 goto (9, 2, 9, 2) term
         read(5, *)
                 stop '42'
              endif
           else
              goto 1
         read(5, *)
              stop '43'
           endif

c (csc th)**2 Li,z Lj,z

 223       if (mode.eq.2) then
              goto 4
         read(5, *)
              stop '44'
           else if (mode.eq.3) then
              goto 21
         read(5, *)
              stop '45'
           else if (smode.and.((i.eq.sz).or.(j.eq.sz))) then
               if (i.eq.sz) then
                 goto (2, 2, 9, 9) term
         read(5, *)
                 stop '46'
              endif
              if (j.eq.sz) then
                 goto (2, 9, 2, 9) term
         read(5, *)
                 stop '47'
              endif
           else if (tmode.and.((i.eq.tz).or.(j.eq.tz))) then
              if (i.eq.tz) then
                 goto (9, 9, 2, 2) term
         read(5, *)
                 stop '48'
              endif
              if (j.eq.tz) then
                 goto (9, 2, 9, 2) term
         read(5, *)
                 stop '49'
              endif
           else
              goto 1
         read(5, *)
              stop '50'
           endif

c Li,x Lj,z + Li,z Lj,x

 224       if (mode.eq.1) then
              goto 4
         read(5, *)
              stop '51'
           else if (mode.eq.3) then
              goto 7
         read(5, *)
              stop '52'
           else if (smode.and.((i.eq.sz).or.(j.eq.sz))) then
              if (i.eq.sz) then
                 goto (11, 11, 18, 18, 1, 1, 2, 9, 2, 9, 2, 9) term
         read(5, *)
                 stop '53'
              endif
              if (j.eq.sz) then
                 goto (2, 9, 2, 9, 2, 9, 11, 11, 18, 18, 1, 1) term
         read(5, *)
                 stop '54'
              endif
           else if (tmode.and.((i.eq.tz).or.(j.eq.tz))) then
              if (i.eq.tz) then
                 goto (2, 2, 9, 9, 26, 26, 9, 2, 9, 2, 9, 2) term
         read(5, *)
                 stop '55'
              endif
              if (j.eq.tz) then
                 goto (9, 2, 9, 2, 9, 2, 2, 2, 9, 9, 26, 26) term
         read(5, *)
                 stop '56'
              endif
           else
              goto 1
         read(5, *)
              stop '57'
           endif

         endif


c -------------- Terme mit Groß-J linear ---------------

         if (teil.lt.(teildiag+teiloffdiag+teiljlin+5)) then
            u=teil-teildiag-teiloffdiag-4
            i=0
            do while (u.gt.0)
               i=i+1
               u=u-(N-3)*blockjlin(i)
            enddo
                      
            u=u+(N-3)*blockjlin(i)
            term=mod(u,blockjlin(i))
            vec=int(u/blockjlin(i))+3
                        
            if(term.eq.0) then
               term=term+blockjlin(i)
               vec=vec-1
            endif
                        
            goto (231, 232, 233, 234, 235, 236) i
         read(5, *)
            stop 'bei i, J linear'

c Jx sum Li,x

 231        if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '58'
            else if (mode.eq.modemax) then
               goto 27
         read(5, *)
               stop '59'
            else if ((smode).and.(vec.eq.sz))then
               goto (11, 18, 1) term
         read(5, *)
               stop '60'
            else if ((tmode).and.(vec.eq.tz)) then
               goto (2, 9, 26) term
         read(5, *)
               stop '61'
            else
               goto 1
         read(5, *)
               stop '62'
            endif

c Jy sum Li,y
 
232        if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '63'
            else if (mode.eq.modemax) then
               goto 28
         read(5, *)
               stop '64'
            else if ((smode).and.(vec.eq.sz)) then
               goto (2, 9, 26) term
         read(5, *)
               stop '65'
            else if ((tmode).and.(vec.eq.tz)) then
               goto (11, 18, 1) term
         read(5, *)
               stop '66'
            else
               goto 1
         read(5, *)
               stop '67'
            endif

c (cot th)**2 Jz sum Li,z

 233        if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '68'
            else if (mode.eq.3) then
               goto 20
         read(5, *)
               stop '69'
            else if (mode.eq.modemax) then
               goto 29
         read(5, *)
               stop '70'
            else if ((smode).and.(vec.eq.sz)) then
               goto (2, 9) term
         read(5, *)
               stop '71'
            else if ((tmode).and.(vec.eq.tz)) then
               goto (9, 2) term
         read(5, *)
               stop '72'
            else
               goto 1
         read(5, *)
               stop '73'
            endif

c (csc th)**2 Jz sum Li,z

 234        if (mode.eq.2) then
               goto 4
         read(5, *)
               stop '74'
            else if (mode.eq.3) then
               goto 21
         read(5, *)
               stop '75'
            else if (mode.eq.modemax) then
               goto 29
         read(5, *)
               stop '76'
            else if ((smode).and.(vec.eq.sz)) then
               goto (2, 9) term
         read(5, *)
               stop '77'
            else if ((tmode).and.(vec.eq.tz)) then
               goto (9, 2) term
         read(5, *)
               stop '78'
            else
               goto 1
         read(5, *)
               stop '79'
            endif

c Jx sum Li,z

 235        if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '80'
            else if (mode.eq.3) then
               goto 7
         read(5, *)
               stop '81'
            else if (mode.eq.modemax) then
               goto 27
         read(5, *)
               stop '82'
            else if ((smode).and.(vec.eq.sz)) then
               goto (2, 9) term
         read(5, *)
               stop '83'
            else if ((tmode).and.(vec.eq.tz)) then
               goto (9, 2) term
         read(5, *)
               stop '84'
            else
               goto 1
         read(5, *)
               stop '85'
            endif

c Jz sum Li,x

 236        if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '86'
            else if (mode.eq.3) then
               goto 7
         read(5, *)
               stop '87'
            else if (mode.eq.modemax) then
               goto 29
         read(5, *)
               stop '88'
            else if ((smode).and.(vec.eq.sz)) then
               goto (11, 18, 1) term
         read(5, *)
               stop '89'
            else if ((tmode).and.(vec.eq.tz)) then
               goto (2, 9, 26) term
         read(5, *)
               stop '90'
            else 
               goto 1
         read(5, *)
               stop '91'
            endif
       
         endif

         if (teil.eq.(teildiag+teiloffdiag+teiljlin+5)) then

c Jy
           if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '92'
            else if (mode.eq.3) then
               
c ------------ 9 --> HO-DVR, 34 --> Legendre-DVR ------------ 

c               goto 34
               goto 9
         read(5, *)
               stop '93'
            else if (mode.eq.modemax) then
               goto 28
         read(5, *)
               stop '94'
            else 
               goto 1
         read(5, *)
               stop '95'
            endif
         
         endif
       

c ------------ Terme mit Groß-J quadratisch -----------------

         if (teil.le.teilmax) then
            u=teil-teildiag-teiloffdiag-teiljlin-blockjlin(7)-4
            
            goto (241, 243, 244, 245) u
         read(5, *)
            stop 'bei jquad, u'

c Jx**2 + Jy**2

 241        if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '96'
            else if (mode.eq.modemax) then
               goto 30
         read(5, *)
               stop '97'
            else
               goto 1
         read(5, *)
               stop '98'
            endif

c (cot th)**2 Jz**2

 243        if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '102'
            else if (mode.eq.3) then
               goto 20
         read(5, *)
               stop '103'
            else if (mode.eq.modemax) then
               goto 32
         read(5, *)
               stop '104'
            else
               goto 1
         read(5, *)
               stop '105'
            endif

c (csc th)**2 Jz**2

 244        if (mode.eq.2) then
               goto 4
         read(5, *)
               stop '106'
            else if (mode.eq.3) then
               goto 21
         read(5, *)
               stop '107'
            else if (mode.eq.modemax) then
               goto 32
         read(5, *)
               stop '108'
            else
               goto 1
         read(5, *)
               stop '109'
            endif

c JxJz+JzJx

 245        if (mode.eq.1) then
               goto 4
         read(5, *)
               stop '110'
            else if (mode.eq.3) then
               goto 7
         read(5, *)
               stop '111'
            else if (mode.eq.modemax) then
               goto 33
         read(5, *)
               stop '112'
            else 
               goto 1
         read(5, *)
               stop '113'
            endif
            
         endif
     








c ------- Operatoren -------------

c H = 1

c 1    write(*,999) '1'

 1    do i=1,dim
         h psi(i)=psi(i)
      enddo

      return


c H = x

c 2    write(*,999) 'x'

 2    do i=1,dim
         h psi(i)=psi(i)*ort(i)
      enddo

      return

c H = x**2

c 3    write(*,999) 'x**2'

 3    do i=1,dim
         h psi(i)=psi(i)*ort(i)**2
      enddo

      return

c H = 1/x**2

c 4    write(*,999) '1/x**2'

 4    do i=1,dim
         h psi(i)=psi(i)/ort(i)**2
      enddo

      return

c H = sin x

c 5    write(*,999) 'sin x'
      
 5    do i=1,dim
         h psi(i)=psi(i)*sin(ort(i))
      enddo
      
      return

c H = cos x

c 6    write(*,999) 'cos x'
      
 6    do i=1,dim
         h psi(i)=psi(i)*cos(ort(i))
      enddo

      return

c H = cot x

c 7    write(*,999) 'cot x'
      
 7    do i=1,dim
         h psi(i)=psi(i)*cos(ort(i))/sin(ort(i))
      enddo

      return

c H = csc x

c 8    write(*,999) 'csc x'
     
 8    do i=1,dim
         h psi(i)=psi(i)/sin(ort(i))
      enddo

      return

c H = d/dx

c 9    write(*,999) 'd/dx'
      
 9    call ddx (mode, h psi, psi, dim, matrix, trafo)

      return

c H = (d/dx)**2

c 10   write(*,999) '(d/dx)**2'
      
 10   call kin (mode, h psi, psi, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=-2.0d0*h psi(i)
      enddo

      return

c H = x d/dx + d/dx x

c 11   write(*,999) 'x d/dx + d/dx x'
      
 11   do i=1,dim
         work(i)=psi(i)*ort(i)
      enddo
      call ddx(mode, h psi, work, dim, matrix, trafo)
      call ddx(mode, work, psi, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=h psi(i)+work(i)*ort(i)
      enddo
      
      return

c H = 1/sqrt(sin x) d/dx sin x d/dx 1/sqrt(sin x)

c 12   write(*,999) '1/sqrt(sin x) d/dx sin x d/dx 1/sqrt(sin x)'
      
 12   do i=1,dim
         work(i)=psi(i)/sqrt(sin(ort(i)))
      enddo
      call ddx (mode, h psi, work, dim, matrix, trafo)
      do i=1,dim
         work(i)=h psi(i)*sin(ort(i))
      enddo
      call ddx (mode, h psi, work, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=h psi(i)/sqrt(sin(ort(i)))
      enddo

      return

c H=1+x**2

c 13   write(*,999) '1+x**2'
      
 13   do i=1,dim
         h psi(i)=psi(i)*(1.0d0+ort(i)**2)
      enddo

      return

c H=(1+x**2)**2

c 14   write(*,999) '(1+x**2)**2'
      
 14   do i=1,dim
         h psi(i)=psi(i)*(1.0d0+ort(i)**2)**2
      enddo

c       do i=1,dim
c          h psi(i)=psi(i)
c       enddo
      return
      

c H=2*(x (d/dx)**2 x + 1)

c 15   write(*,999) '2*(x (d/dx)**2 x + 1)'
      
 15   do i=1,dim
         work(i)=psi(i)*ort(i)
      enddo
      call kin (mode, h psi, work, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=-4.0d0*(h psi(i)*ort(i)-0.5d0*psi(i))
      enddo

      return

c H=x**2 (d/dx)**2 x**2

c 16   write(*,999) 'x**2 (d/dx)**2 x**2'
      
 16   do i=1,dim
         work(i)=psi(i)*ort(i)**2
      enddo
      call kin (mode, h psi, work, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=-2.0d0*h psi(i)*ort(i)**2
      enddo

      return
      
c H=x (d/dx)**2 x

c 17   write(*,999) 'x (d/dx)**2 x'
      
 17   do i=1,dim
         work(i)=psi(i)*ort(i)
      enddo
      call kin (mode, h psi, work, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=-2.0d0*h psi(i)*ort(i)
      enddo

      return

c H=1-x**2

c 18   write(*,999) '1-x**2'
      
 18   do i=1,dim
         h psi(i)=psi(i)*(1.0d0-ort(i)**2)
      enddo

      return

c H=(1-x**2)**2

c 19   write(*,999) '(1-x**2)**2'
      
 19   do i=1,dim
         h psi(i)=psi(i)*(1.0d0-ort(i)**2)**2
      enddo

      return
 
c H=(cot x)**2

c 20   write(*,999) '(cot x)**2'
      
 20   do i=1,dim
         h psi(i)=psi(i)*(cos(ort(i))/sin(ort(i)))**2
      enddo

      return

c H=(csc x)**2

c 21   write(*,999) '(csc x)**2'
      
 21   do i=1,dim
         h psi(i)=psi(i)/sin(ort(i))**2
      enddo

      return

c H=x * (1-x**2)

c 22   write(*,999) 'x*(1-x**2)'
      
 22   do i=1,dim
         h psi(i)=psi(i)*ort(i)*(1.0d0-ort(i)**2)
      enddo

      return

c H=x (d/dx)**2 + (d/dx)**2 x

c 23   write(*,999) 'x (d/dx)**2 + (d/dx)**2 x'
      
 23   do i=1,dim
         work(i)=psi(i)*ort(i)
      enddo
      call kin(mode, h psi, work, dim, matrix, trafo)
      call kin(mode, work, psi, dim, matrix, trafo)
      do i=1, dim
         h psi(i)=-2.0d0*(h psi(i)+work(i)*ort(i))
      enddo

      return

c H=x**3 d/dx + d/dx x**3

c 24   write(*,999) 'x**3 d/dx + d/dx x**3'
      
 24   do i=1,dim
         work(i)=psi(i)*ort(i)**3
      enddo
      call ddx(mode, h psi, work, dim, matrix, trafo)
      call ddx(mode, work, psi, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=h psi(i)+work(i)*ort(i)**3
      enddo

      return

c H=d/dx - 3*x d/dx x

c 25   write(*,999) 'd/dx - 3*x d/dx x'
      
 25   do i=1,dim
         work(i)=psi(i)*ort(i)
      enddo
      call ddx(mode, h psi, work, dim, matrix, trafo)
      call ddx(mode, work, psi, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=-3.0d0*ort(i)*h psi(i)+work(i)
      enddo

      return

c H=x d/dx x

c 26   write(*,999) 'x d/dx x'
      
 26   do i=1,dim
         work(i)=psi(i)*ort(i)
      enddo
      call ddx(mode, h psi, work, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=ort(i)*h psi(i)
      enddo

      return

c H=2*Jx = J+ + J-

c 27   write(*,999) 'Jx'
      
 27   jmin=int(para(1,mode))
      jmax=int(para(2,mode))
      kmin=int(para(3,mode))
      kmax=int(para(4,mode))
      i=1
      do j=jmin,jmax
         h psi(i)=0.0d0
         do k=max(-j,kmin),min(j,kmax)-1
            wurz=dsqrt(dfloat(j*(j+1)-k*(k+1)))
            h psi(i)=h psi(i)+wurz*psi(i+1)
            h psi(i+1)=wurz*psi(i)
            i=i+1
         enddo
         i=i+1
      enddo

      return

c H=2i*Jy = J+ - J-
c H=-2i*Jy in body-fixed frame !

c 28   write(*,999) 'Jy'
      
 28   jmin=int(para(1,mode))
      jmax=int(para(2,mode))
      kmin=int(para(3,mode))
      kmax=int(para(4,mode))
      i=1
      do j=jmin,jmax
         h psi(i)=0.0d0
         do k=max(-j,kmin),min(j,kmax)-1
            wurz=dsqrt(dfloat(j*(j+1)-k*(k+1)))
            h psi(i)=h psi(i)-wurz*psi(i+1)
            h psi(i+1)=wurz*psi(i)
            i=i+1
         enddo
         i=i+1
      enddo

      return

c H=Jz

c 29   write(*,999) 'Jz'
      
 29   jmin=int(para(1,mode))
      jmax=int(para(2,mode))
      kmin=int(para(3,mode))
      kmax=int(para(4,mode))
      i=1
      do j=jmin,jmax
         do k=max(-j,kmin),min(j,kmax)
            h psi(i)=psi(i)*k
            i=i+1
         enddo
      enddo

      return

c H=Jx**2 + Jy**2 = J**2 - Jz**2

c 30   write(*,999) 'Jx**2 + Jy**2'
      
 30   jmin=int(para(1,mode))
      jmax=int(para(2,mode))
      kmin=int(para(3,mode))
      kmax=int(para(4,mode))
      i=1
      do j=jmin,jmax
         do k=max(-j,kmin),min(j,kmax)
            h psi(i)=psi(i)*(j*(j+1)-k**2)
            i=i+1
         enddo
      enddo

      return

c H=Jz**2

c 32   write(*,999) 'Jz**2'
      
 32   jmin=int(para(1,mode))
      jmax=int(para(2,mode))
      kmin=int(para(3,mode))
      kmax=int(para(4,mode))
      i=1
      do j=jmin,jmax
         do k=max(-j,kmin),min(j,kmax)
            h psi(i)=psi(i)*k**2
            i=i+1
         enddo
      enddo

      return

c H=JxJz+JzJx

c 33   write(*,999) 'JxJz+JzJx'
      
 33   jmin=int(para(1,mode))
      jmax=int(para(2,mode))
      kmin=int(para(3,mode))
      kmax=int(para(4,mode))
      i=1
      do j=jmin,jmax
         h psi(i)=0.0d0
         do k=max(-j,kmin),min(j,kmax)-1
            wurz=dsqrt(dfloat(j*(j+1)-k*(k+1)))
            h psi(i)=h psi(i)+wurz*psi(i+1)*(k+0.5d0)
            h psi(i+1)=wurz*psi(i)*(k+0.5d0)
            i=i+1
         enddo
         i=i+1
      enddo
            
      return

c H=d/dx + 0.5 cot x
c Legendre-DVR
c 34   write(*,999) 'd/dx + 0.5 cot x, Legendre'

 34   do i=1,dim
         work(i)=psi(i)/sin(ort(i))
      enddo
      call ddx (mode, h psi, work, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=h psi(i)+0.5d0*cos(ort(i))/sin(ort(i))*psi(i)
      enddo
      return

c H=1/sin x d/dx sin x d/dx
c Legendre-DVR
c 35   write(*,999) '1/sin x d/dx sin x d/dx, Legendre'

 35   call kin (mode, h psi, psi, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=-2.0d0*h psi(i)
      enddo 
      return

c 999  format (A, '   ', $)

      end
      

      subroutine IntToCartCH4Radau(q, xout)
      use iso_c_binding
	implicit none
      
      integer N
      parameter (N=5)
      real*8    q(3*N-6), v
      real*8 cart(N,3)
	real*8 xout(3*N)
      integer i, j, k
      logical*4 vm
      real*8 mh
      real*8 ar
      
c     Massenparameter, Handbook of Atomic, Molecular,
c     and Optical Physics, Gordon W. F. Drake
      parameter (ar=1822.88848d0)

      real*8 r1, r2, theta
      real*8 x(0:N-1,3), rad(0:N-1,3)
      real*8 Cadj (0:N-1,0:N-1)
      real*8 m(0:N-1), alpha(0:N-1), mg
      real*8 dist(N-1,N-1)
      
      r1=q(1)
      r2=q(2)
c     Adjustment for c++-Version
c      q(3)=acos(q(3))
      theta=q(3)

c     12C-Isotopenmasse      
      
      m(0)=12.0d0*ar

      mg=m(0)
      
c     Isotopenmasse, Handbook of Atomic, Molecular,
c     and Optical Physics, Gordon W. F. Drake
      mh=1.007825047d0*ar      
       
      do i=1,N-1
         m(i)=mh
         mg=mg+m(i)
      enddo

      do i=0,N-1
         alpha(i)=dsqrt(m(i)/mg)
      enddo

c     Umrechnung interne -> kart. Radaus

c     Massenschwerpunkt -> (0,0,0)T
      
      rad(0,1)=0.0d0
      rad(0,2)=0.0d0
      rad(0,3)=0.0d0

      rad(1,1)=0.0d0
      rad(1,2)=0.0d0
      rad(1,3)=r1

      rad(2,1)=r2*dsin(theta)
      rad(2,2)=0.0d0
      rad(2,3)=r2*dcos(theta)
      
      do i=3,N-1
         rad(i,1)=2.0d0*q(3*(i-2)+1)*q(3*(i-2)+2)
     .            /(1+q(3*(i-2)+2)**2+q(3*(i-2)+3)**2)

         rad(i,2)=2.0d0*q(3*(i-2)+1)*q(3*(i-2)+3)
     .            /(1+q(3*(i-2)+2)**2+q(3*(i-2)+3)**2)

c     negatives Vorzeichen z-Komponente bei Projektion vom Nordpol       
       
       rad(i,3)=-1.0d0*q(3*(i-2)+1)*(1-q(3*(i-2)+2)**2-q(3*(i-2)+3)**2)
     .            /(1+q(3*(i-2)+2)**2+q(3*(i-2)+3)**2)
      enddo

c     Transformationsmatrix aufbauen
            
      do i=0,N-1
         Cadj(0,i)=alpha(i)
         Cadj(i,0)=alpha(i)
      enddo

      do i=1,N-1
         do j=1, N-1
           Cadj(i,j)=alpha(i)*alpha(j)/(alpha(0)+1.0d0)
         enddo
      enddo

c     Diagonalelemente korrigieren

      do i=1,N-1
         Cadj(i,i)=Cadj(i,i)-1.0d0              
      enddo
      
c     Transformation der Radau-Vektoren

      do i=0,N-1
         x(i,1)=0.0d0
         x(i,2)=0.0d0
         x(i,3)=0.0d0
      enddo

      do i=0,N-1
         do j=0,N-1
            x(i,1)=x(i,1)+Cadj(i,j)*rad(j,1)
            x(i,2)=x(i,2)+Cadj(i,j)*rad(j,2)
            x(i,3)=x(i,3)+Cadj(i,j)*rad(j,3)
         enddo
      enddo

c     Massengewichtung rueckgaengig, Vorzeichenwechsel

      do i=0,N-1
         x(i,1)=-1.0d0*x(i,1)/dsqrt(m(i))
         x(i,2)=-1.0d0*x(i,2)/dsqrt(m(i))
         x(i,3)=-1.0d0*x(i,3)/dsqrt(m(i))
      enddo

c     Koordinaten Potentialroutine      
      
      do i=0,N-1
         cart(i+1,1)=x(i,1)
         cart(i+1,2)=x(i,2)
         cart(i+1,3)=x(i,3)
      enddo
	k=1
      do i=1,N
         do j=1, 3
		xout(k)=cart(i, j)
	   	k=k+1
         enddo
      enddo

c     test for unphysical interproton distances
c      call checkdist (cart,N,vm)
     
	end subroutine

      subroutine checkdist (x,N,vm)

      implicit none
      
      real*8 hh2,hh3,hh4,hh5
      parameter (hh2=1.12d0,hh3=1.7d0,hh4=1.9d0,hh5=3.d0)
      integer nh
      parameter (nh=10)

      integer N,i,j,l
      logical*4 vm
      real*8 x(N,3),dist(nh),h

c calculate interproton distances
      l=0
      do i=2,N
         do j=i+1,N
            l=l+1
            dist(l)=(x(j,1)-x(i,1))**2
     .             +(x(j,2)-x(i,2))**2
     .             +(x(j,3)-x(i,3))**2        
            dist(l)=dsqrt(dist(l))
         enddo
      enddo

c     sort interproton distances
      do i=1,nh
         do j=i+1,nh
            if (dist(j).lt.dist(i)) then
               h=dist(i)
               dist(i)=dist(j)
               dist(j)=h
            endif
         enddo
      enddo

c     check interproton distances
      vm=.true.
c      print*, 'x='
c      do i=1,6
c      print*, x(i,:)
c      enddo
c      print*, 'dist='
c      print*, dist(:)
c      read(5,*)
      if (dist(2).le.hh2) return 
      if (dist(3).le.hh3) return 
      if (dist(4).le.hh4) return 
      if (dist(5).le.hh5) return 
      vm=.false.

      return
      end
