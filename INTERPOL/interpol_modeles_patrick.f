      program interpol_model

 
****************************************************************************
* interpolate a model for babsma between max 8 models given in input at the*
* point defined by the requested parameters (Teff,logg,z)                  *
* TM 07/2003                                                               *
****************************************************************************
! to compile with Fortran 90 or 95 otherwise risk of incompatibility ! 
!
! Feb05: Adapted to New MARCS models - P.deLaverny
!        Not tested for all kind of models
!        Compilation with lf95



      implicit none
      integer :: file,k,ndp,ndepth_ref,out,intryc
      parameter (ndp=100)
      logical :: verif,test,extrapol
      real :: temp_ref,logg_ref,z_ref,x,y,z,xinf,
     &        xsup,yinf,ysup,zinf,zsup,xlr,alpha_ref,scale
      character*256, dimension (10) :: FILE_IN
      real, dimension (10,ndp):: taus,T,Pe,Pg,xit
      integer, dimension (10) :: ndepth
      real, dimension (10) :: teff,logg,metal,alpha
      character*1, dimension (10) :: geometry
      character*1 :: geo_ref
      external :: extract,blend_103
      real, external :: inf,sup

      verif=.true.

      write(*,*) '*****************************'
      write(*,*) '* begining of interpolation *'
      write(*,*) '*****************************'  
      write(*,*)
      write(*,*)'Modeles doivent avoir nombres de couches identiques,'
      write(*,*)'        meme [alpha/Fe] et meme geometrie (pp/sph)'
      write(*,*)

******  read 8 models, put in tables ******
      out=9
      write(*,*) 'Interpolation between :'
      do file=1,9
         read(*,*) FILE_IN(file)
      end do 
       read(*,*) temp_ref
       read(*,*) logg_ref
       read(*,*) z_ref
       read(*,*) test
      do file=1,8
      call extract(FILE_IN(file),teff(file),logg(file),metal(file),
     & ndepth(file),taus(file,:),T(file,:),Pe(file,:),
     & Pg(file,:),xit(file,:),geometry(file),alpha(file))
      ndepth_ref=ndepth(1)
      geo_ref = geometry(1)
      alpha_ref = alpha(1)

      write(*,*)'tester verif avec alpha et geometry'

      verif=verif.and.(ndepth(file).eq.ndepth_ref).and.
     &      (geometry(file).eq.geo_ref).and.(alpha(file).eq.alpha_ref)
      write(*,78) file,teff(file),logg(file),metal(file)      
c      write(*,*)'    Ndepth ',ndepth(file)
 78   format('model',i2,'  Teff=',f8.0,'  logg=',f5.2,'  z=',f6.2)
      end do
      write(*,*)
      write(*,79) temp_ref,logg_ref,z_ref
 79   format('Interpolation point : Teff=',f8.0,'  logg=',f5.2,
     &                            '  z=',f6.2)

************** check if files are length and depth ref compatible *******
      if (.not.(verif)) then
         write(*,*) 'All the models do not have the same' 
          write(*,*) 'number of Ross depths or geometry or [alpha/Fe]'
          write(*,*) 'no interpolation done'
      
********* calculation of the interpolation point(x,y,z) for routine blend_103 ******************

       else    
         xinf=inf(teff)
         yinf=inf(logg)
         zinf=inf(metal)
         xsup=sup(teff)
         ysup=sup(logg)
         zsup=sup(metal)
         if (xsup.eq.xinf) then
            x=0
         else
          x=(temp_ref-xinf)/(xsup-xinf)
         end if 
         if (ysup.eq.yinf) then
            y=0
         else      
          y=(logg_ref-yinf)/(ysup-yinf)
         end if
         if (zsup.eq.zinf) then
            z=0
         else   
          z=(z_ref-zinf)/(zsup-zinf)
           if (z.lt.0) then 
                     z=-((abs(z))**((temp_ref/4000)**3))   ! empiric calibration by myself
                     else 
                        z=z**((temp_ref/4000)**3)
           end if   
         end if
         extrapol=((x.lt.0).or.(x.gt.1).or.(y.lt.0).or.(y.gt.1)
     &         .or. (z.lt.0).or.(z.gt.1))
         if (extrapol) then
         write(*,*) '!!!  WARNING : extrapolation  !!!'
         end if
         
****** interpolation for each column (taus,teff,Pe,Pg,microt) and at each Rosseland depth *****************
      
         do k=1,ndepth_ref
          call blend_103(x,y,z,taus(1,k),taus(2,k),
     &     taus(3,k),taus(4,k),taus(5,k),taus(6,k),taus(7,k),taus(8,k)
     &     ,taus(out,k))
          
          call blend_103(x,y,z,T(1,k),T(2,k),T(3,k),T(4,k)
     &     ,T(5,k),T(6,k),T(7,k),T(8,k),T(out,k))
          
          call blend_103(x,y,z,Pe(1,k),Pe(2,k),Pe(3,k),Pe(4,k)
     &     ,Pe(5,k),Pe(6,k),Pe(7,k),Pe(8,k),Pe(out,k))
          
          call blend_103(x,y,z,Pg(1,k),Pg(2,k),Pg(3,k),Pg(4,k)
     &     ,Pg(5,k),Pg(6,k),Pg(7,k),Pg(8,k),Pg(out,k))
          
          call blend_103(x,y,z,xit(1,k),xit(2,k),
     &     xit(3,k),xit(4,k),xit(5,k),xit(6,k),xit(7,k),xit(8,k)
     &     ,xit(out,k))        
       end do

       ndepth(out)=ndepth_ref
       xlr= 5000.0

******** write interpolated model in file nber out ***********

       open(unit=21,file=FILE_IN(out))
c       write(21,1966) ndepth_ref,xlr,temp_ref,logg_ref,z_ref
        intryc = 0        ! lu par babsma mais utilite ????
        scale = 0.d0      ! babsma: t(k)=(1.000+scale)*t(k)
        write(21,1966) ndepth_ref,xlr,logg_ref,intryc,scale
1966     format('''INTERPOL''',x,i3,f8.0,2x,f8.2,2x,i2,2x,f8.2)
         do k=1,ndepth_ref
           write(21,1965) taus(out,k),T(out,k),Pe(out,k),
     &                    Pg(out,k),xit(out,k)
1965       format(f8.4,x,f8.2,3(x,f8.4))
        enddo
      close(21)
********** write files for check plotting (see shell script interpol.com)************
      open (unit=34,file='modele.sm')
      do file=1,8
       do k=1,ndepth_ref
        write(34,*) taus(file,k),T(file,k),Pe(file,k),Pg(file,k)
       end do
      end do

!case of a 10th check file
      read(*,*) FILE_IN(10)
      if (test) then 
         file=10
      call extract(FILE_IN(file),teff(file),logg(file),metal(file),
     & ndepth(file),taus(file,:),T(file,:),Pe(file,:),
     & Pg(file,:),xit(file,:),geometry(file),alpha(file))
      do k=1,ndepth_ref
        write(34,*) taus(10,k),T(10,k),Pe(10,k),Pg(10,k)
      end do
      end if


      write (*,*) 
      if (extrapol) then
          write (*,*) 'Extrapolation done'
          else
          write (*,*) 'Interpolation done'
      end if
      end if

      end

c--------------------------------------------------------------------------------

      subroutine extract(FILE,Teff,grav,metal,ndepth,tau,T,Pe,Pg,
     &                   xit,geometry,alpha)

!     Version  grandement modifiee en fevrier 2005 pour modeles NMARCS
!              lecture modele comme dans babsma
!              Patrick

!extracted from osplot.f 07/2003, to get tau,T,Pe,Pg,microturb from a model

       CHARACTER*256 FILE
       CHARACTER*50 MOCODE,blabla
       CHARACTER*1 geometry
       real :: Teff,grav,metal,xic,alpha
       integer :: ndepth
       real, dimension (100):: tau,T,Pe,Pg,xit

          imod = 10  
          blabla=''
          OPEN(UNIT=imod,FILE=FILE,STATUS='OLD')
          read(imod,'(a)') mocode
          print*
c          print*,mocode
          if (mocode(1:1).eq.'s') then
c            print*,' this model is SPHERICALLY SYMMETRIC'
          else if (mocode(1:1).eq.'p') then
c            print*,' this model is PLANE PARALLEL'
          else
            print*,' This model may not be a NewMARCS model!'
          endif
          geometry=mocode(1:1)
          read(imod,*)Teff
          read(imod,*)
          read(imod,*)grav
          grav=log10(grav)
          read(imod,*)xic
          read(imod,*)
          read(imod,*)metal,alpha
          do while (blabla.ne.'Model structure')
            read(imod,'(a)') blabla
          enddo
          backspace(imod)
          backspace(imod)
          read(imod,*)ndepth
          read(imod,*)
          read(imod,*)          
          do k=1,ndepth
* the first dum is tauross, the second is Depth, the third is Prad, 
*     the fourth is Pturb
            read(imod,*) idum,dum,tau(k),dum,T(k),
     &                   Pe(k),Pg(k),dum,dum
c           tau(k)=10.**tau(k)
            Pe(k)  = log10(Pe(k))
            Pg(k)  = log10(Pg(k))
            xit(k) = xic
          enddo
         close(imod)
         RETURN
         END

c-------------------------------------------------------------------------

      real function inf(tab)
      implicit none
      integer :: n
      real,dimension(9) :: tab
      inf=tab(1)
      do n=2,8
         if (tab(n).lt.inf) then
            inf=tab(n)
         end if
      end do
      end function

      real function sup(tab)
      implicit none
      integer :: n
      real,dimension(9) :: tab
      sup=tab(1)
      do n=2,8
         if (tab(n).gt.sup) then
            sup=tab(n)
         end if
      end do
      end function



c---------------------------------------------------------------------------------
      subroutine blend_103 (r,s,t,x000,x001,x010,x011,x100,x101,x110, 
     & x111, x )
!
!*******************************************************************************
!from http://www.math.iastate.edu/burkardt/f_src/
!
!! BLEND_103 extends scalar point data into a cube.
!
!
!  Diagram:
!
!    011--------------111
!      |               |
!      |               |
!      |               |
!      |               |
!      |               |
!    001--------------101
!
!
!      *---------------*
!      |               |
!      |               |
!      |      rst      |
!      |               |
!      |               |
!      *---------------*
!
!
!    010--------------110
!      |               |
!      |               |
!      |               |
!      |               |
!      |               |
!    000--------------100
!
!
!  Formula:
!
!    Written as a polynomial in R, S and T, the interpolation map has the
!    form:
!
!      X(R,S,T) =
!        1         * ( + x000 )
!      + r         * ( - x000 + x100 )
!      +     s     * ( - x000        + x010 )
!      +         t * ( - x000               + x001 )
!      + r * s     * ( + x000 - x100 - x010                       + x110 )
!      + r     * t * ( + x000 - x100        - x001        + x101 )
!      +     s * t * ( + x000        - x010 - x001 + x011 )
!      + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, T, the coordinates where an interpolated value
!    is desired.
!
!    Input, real X000, X001, X010, X011, X100, X101, X110, X111, the
!    data values at the corners.
!
!    Output, real X, the interpolated data value at (R,S,T).
!
      implicit none
!
      real r
      real s
      real t
      real x
      real x000
      real x001
      real x010
      real x011
      real x100
      real x101
      real x110
      real x111
!
!  Interpolate the interior point.
!
      x = 
     & 1.0E+00     * ( + x000 ) 
     & + r         * ( - x000 + x100 ) 
     & +     s     * ( - x000        + x010 ) 
     & +         t * ( - x000               + x001 ) 
     & + r * s     * ( + x000 - x100 - x010                      + x110) 
     & + r     * t * ( + x000 - x100        - x001        + x101 ) 
     & +     s * t * ( + x000        - x010 - x001 + x011 ) 
     & + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110+
     &                                                            x111 )

      return
      end
