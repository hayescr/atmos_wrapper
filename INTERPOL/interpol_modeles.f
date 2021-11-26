      program interpol_modeles


****************************************************************************
* interpolation of model atmosphere
* parameter space of interpolation {Teff,logg,z}
* 8 MARCS binary models as input required:
* ! the order of input models matters !
* ! only standard composition has been tested, no guaranty for peculiar composition !
* turbospectrum/babsma format compatible output

* Interpolation scheme of model atmosphere(@)
* in {Teff,logg,z} stellar parameter space
*         ^ z
*         |_ _ _ _ _
*         |        /|
*        /|       /
*       _ |_ _ _ /  |
*       | |      |
*         |         |
*       | |    @ |
*         |--------------> logg
*       |/       | /
*       /_ _ _ _ ./
*      /
*     /
*    Teff
*
* Each structure component of each model (T,Pe,Pg,xit,optical depth, kappaross)
* is first resampled on a common optical depth basis tau(5000 or ross) (see resample routine).
*
* Then each component of the atmosphere(@) (T,Pe,Pg,xit,kappa,geom depth +tau5000 and tauross))
* is interpolated individually along each dimension between each input model value(*)
* (physically it is impossible to ensure a simple relationship between models in more than one dimension at the same time,
* that's why the input models MUST form a "cube" in the stellar parameter space {Teff,logg,z})
* The interpolation is successively done at each optical depth (tau)
* + interpolation point weighted by an empirical value (see TM thesis + manual)
*
*
*                 ^ T,Pe,Pg,xit or radius
*                 |
*                 |            *
*                 |           /|
*                 |          @
*                 |         /| |
*                 |   /    /
*                 |  /    /  | |
*         ^       | /    *
*         |       |/     |   | |
*         |       -----------------------> Teff, logg and z
*         |      /      low ref up
*         |     /*      /   / /
*         |    / |\
*         |   /    \  /   / /
*         |  /   |  \
*         | /       /@  / /
*         |/     |   |\
*         / _  _ /_ _ /*_ _ _ _ _
*        /      low ref up
*       /
*   tau{5000 or Ross}
*
***************************************************************************
* TM 07/2003

c  07/2004 resampling of each model on a common optical depth basis
c  06/2006 works for spherical geometry models
c  09/2007 new calibration of free parameter alpha
c           + modified to read Uppsala ascii models
c           + kappa interpolation
c           + rhox calculation
c           + 2 outputs (babsma and ATLAS/MOOG)
c  10/2007 error estimates
c  10/2011 two non crucial bugs fixed
c             -unformatted->formatted read for ascii models because there is a couple of trouble makers in the grid
c             -dimension of taubas was not matching tauR (emo@astro.indiana.edu)
c  12/2012 bug fixed whiled reading MARCS interpolated models
c  12/2012 include Phoenix reading  (+ Kurucz model to be finished)
c        (->modtype become a character)
c  01/2013 choose output format Turbo MOOG Atlas + (+MARCS web format to be finished)
c  11/2018 interpolate rhox; comment the rhox computation that does not seem to work
****************************************************************************
!  compile with Fortran 90 or 95


      implicit none
      integer :: file,nfile,k,ndp,ndepth_ref,out,nlinemod,extension
      integer :: signteffpoint,signloggpoint,signmetpoint
      parameter (ndp=200)
      parameter (nfile=10)
      logical :: verif,check,test,extrapol,binary,optimize
      real :: lambda_ref,temp_ref,logg_ref,z_ref,x,y,z,xinf,mass,
     &        xsup,yinf,ysup,zinf,zsup,teffpoint,loggpoint,metpoint
      character*256, dimension (nfile) :: FILE_IN
      real, dimension (:,:), allocatable:: taus,tauR,T,Pe,Pg,xit,rr,
     &xkapref,rhox,taus_aux,tauR_aux,T_aux,Pe_aux,Pg_aux,xit_aux,rr_aux,
     &   xkapref_aux,Prad_aux,Prad,Pturb_aux,Pturb,rhox_aux,ro,ro_aux
      integer, dimension (nfile) :: ndepth
      real, dimension (nfile) :: xlr,teff,logg,metal
      logical, dimension (nfile) :: sph
      external :: blend_103
      real, external :: inf,sup
      real, dimension(12,3) :: lin_dif,power
      real, dimension(92,nfile) :: abu
      character*20 :: modtype,outtype
      INTERFACE reec
        subroutine resample(taus,tauR,T,Pe,Pg,Prad,Pturb,xit,rr,xkapref
     &                              ,rhox,ro)
        real, dimension (:,:) :: taus,tauR,T,Pe,Pg,Prad,Pturb,xit,rr,
     &                            xkapref,rhox,ro
        end
      END INTERFACE reec

      write(*,*) '*****************************'
      write(*,*) '* begining of interpolation *'
      write(*,*) '*****************************'

******* you can choose here to switch of the "optimization" and prefer simple linear interpolation

      optimize = .true.

******  read 8 models, put in tables,
****** check number of layer, reference optical depth, and geometry compatibility ******

      out=9
      write(*,*) 'Interpolation between :'
      do file=1,9
         read(*,*) FILE_IN(file)
      end do
       read(*,*) modtype
       read(*,*) temp_ref
       read(*,*) logg_ref
       read(*,*) z_ref
       read(*,*) outtype
       read(*,*) test

      verif=.true.
      check=.true.
      nlinemod=ndp

      allocate(taus_aux(nlinemod,nfile),tauR_aux(nlinemod,nfile),
     & T_aux(nlinemod,nfile),Pe_aux(nlinemod,nfile),
     & Pg_aux(nlinemod,nfile),xit_aux(nlinemod,nfile),
     & rr_aux(nlinemod,nfile),xkapref_aux(nlinemod,nfile),
     & Prad_aux(nlinemod,nfile),Pturb_aux(nlinemod,file),
     & rhox_aux(nlinemod,file),ro_aux(nlinemod,file) )

      do file=1,8
      select case (adjustl(adjustr(modtype)))
      case ('MARCSbin')
      if (file.eq.1) then
       write(*,*) 'models are MARCS binary format'
      endif
      call extract_bin(FILE_IN(file),teff(file),logg(file),metal(file),
     & ndepth(file),xlr(file),taus_aux(:,file),tauR_aux(:,file),
     & T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),xit_aux(:,file),
     & rr_aux(:,file),sph(file),xkapref_aux(:,file),
     & Prad_aux(:,file),Pturb_aux(:,file),rhox_aux(:,file),
     & ro_aux(:,file),abu(:,file))
      case ('MARCSweb')
      if (file.eq.1) then
       write(*,*) 'models are MARCS web format'
      endif
      call extract_ascii(FILE_IN(file),teff(file),logg(file),
     & metal(file),ndepth(file),xlr(file),taus_aux(:,file),
     & tauR_aux(:,file),T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),
     & xit_aux(:,file),rr_aux(:,file),sph(file),xkapref_aux(:,file),
     & Prad_aux(:,file),Pturb_aux(:,file),rhox_aux(:,file),
     & ro_aux(:,file),abu(:,file))
      case ('Phoenix')
      if (file.eq.1) then
       write(*,*) 'models are Phoenix format'
      endif
       call extract_Phoenix(FILE_IN(file),teff(file),logg(file),
     & metal(file),ndepth(file),xlr(file),taus_aux(:,file),
     & tauR_aux(:,file),T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),
     & xit_aux(:,file),rr_aux(:,file),sph(file),xkapref_aux(:,file),
     & Prad_aux(:,file),Pturb_aux(:,file),rhox_aux(:,file),
     & ro_aux(:,file),abu(:,file))
       case ('Atlas')
       if (file.eq.1) then
       write(*,*) 'models are Atlas format'
      endif
      call extract_Atlas(FILE_IN(file),teff(file),logg(file),
     & metal(file),ndepth(file),xlr(file),taus_aux(:,file),
     & tauR_aux(:,file),T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),
     & xit_aux(:,file),rr_aux(:,file),sph(file),xkapref_aux(:,file),
     & Prad_aux(:,file),Pturb_aux(:,file),rhox_aux(:,file),
     & ro_aux(:,file),abu(:,file))
      case default
       write(*,*) 'model type unknown:',modtype
       stop
      end select

      if (.not.(((sph(1).and.sph(file))).or.
     &    ((.not.(sph(1))).and.(.not.(sph(file)))))) then
      write(*,*) 'geometry compatibility problem with'
       write(*,78) file,teff(file),logg(file),metal(file)
      stop
      end if
      check=check.and.(xlr(file).eq.xlr(1))
      verif=verif.and.(ndepth(file).eq.ndepth(1))
      write(*,78) file,teff(file),logg(file),metal(file)
      end do
 78   format('model',i2,'  Teff=',f8.0,'  logg=',f5.2,'  z=',f6.2)


      write(*,79) temp_ref,logg_ref,z_ref
 79   format('Interpolation point : Teff=',f8.0,'  logg=',f5.2,
     &                            '  z=',f6.2)
      teff(9) = temp_ref
      logg(9) = logg_ref
      metal(9) = z_ref
************** check if files are length and depth ref compatible *******
      if (.not.(check)) then
         write(*,*) 'All the models do not have the same'
          write(*,*) 'lambda ref'
          write(*,*) 'no interpolation done'
          stop

       else
      if (.not.(verif)) then
         write(*,*) 'WARNING : All the models do not have the same'
         write(*,*) 'number of layers, resampling to',
     &                 ndepth(1),'layers'
      end if
      ndepth_ref=ndepth(1)
      lambda_ref=xlr(1)

********* calculation of the interpolation point(x,y,z) in {Teff,logg,z} space******************

       allocate(taus(ndepth_ref,nfile),tauR(ndepth_ref,nfile),
     & T(ndepth_ref,nfile),Pe(ndepth_ref,nfile),Pg(ndepth_ref,nfile),
     & xit(ndepth_ref,nfile),rr(ndepth_ref,nfile),
     &  xkapref(ndepth_ref,nfile),Prad(ndepth_ref,nfile),
     &   Pturb(ndepth_ref,file),rhox(ndepth_ref,file),
     &          ro(ndepth_ref,file))
       taus=taus_aux(1:ndepth_ref,:)
       tauR=tauR_aux(1:ndepth_ref,:)
       T=T_aux(1:ndepth_ref,:)
       Pe=Pe_aux(1:ndepth_ref,:)
       Pg=Pg_aux(1:ndepth_ref,:)
       xit=xit_aux(1:ndepth_ref,:)
       rr=rr_aux(1:ndepth_ref,:)
       xkapref=xkapref_aux(1:ndepth_ref,:)
       Prad=Prad_aux(1:ndepth_ref,:)
       Pturb=Pturb_aux(1:ndepth_ref,:)
       rhox=rhox_aux(1:ndepth_ref,:)
       ro=ro_aux(1:ndepth_ref,:)

         xinf=inf(teff)
         yinf=inf(logg)
         zinf=inf(metal)
         xsup=sup(teff)
         ysup=sup(logg)
         zsup=sup(metal)
         if (xsup.eq.xinf) then
            teffpoint=0.0
         else
          teffpoint=(temp_ref-xinf)/(xsup-xinf)
         end if
         if (ysup.eq.yinf) then
            loggpoint=0.0
         else
          loggpoint=(logg_ref-yinf)/(ysup-yinf)
         end if
         if (zsup.eq.zinf) then
            metpoint=0.0
         else
          metpoint=(z_ref-zinf)/(zsup-zinf)
         end if
         extrapol=((teffpoint.lt.0).or.(teffpoint.gt.1)
     &         .or.(loggpoint.lt.0).or.(loggpoint.gt.1)
     &         .or.(metpoint.lt.0).or.(metpoint.gt.1))

         if (teffpoint.lt.0) then
           signteffpoint = -1
           teffpoint = -1*teffpoint
         else
            signteffpoint = +1
         end if
         if (loggpoint.lt.0) then
           signloggpoint = -1
           loggpoint = -1 * loggpoint
         else
            signloggpoint = +1
         end if
         if (metpoint.lt.0) then
           signmetpoint = -1
           metpoint = -1 *  metpoint
           else
            signmetpoint = +1
         end if


         extension=INDEX(FILE_IN(out),'.',
     &                        .true.)-1
         if (extension <=0) then
          extension=len_trim(FILE_IN(out))
         endif
         if (extrapol) then
         write(0,*) '!!!  WARNING : model extrapolation  !!!'
          FILE_IN(out)=FILE_IN(out)(1:extension)//'.ext'
         else
          FILE_IN(out)=FILE_IN(out)(1:extension)//'.int'
         end if

*
*******resample each layer of each input model on a common depth basis(tau5000 or tauRoss, see resample routine)*****************
!if you don't want to resample all the model to the same depth scale, just comment the following line
        call resample(taus,tauR,T,Pe,Pg,Prad,Pturb,xit,rr,xkapref,rhox,
     &                                                               ro)

****** initialisation of empirical constants for optimized interpolation (see TM thesis)*************

         lin_dif(1,1)=0                          !tau5000 vs Teff
         lin_dif(1,2)=0                          ! ...    vs logg
         lin_dif(1,3)=0                          ! ...    vs z
         lin_dif(2,1)=0                          !tauross vs Teff
         lin_dif(2,2)=0                          ! ...    vs logg
         lin_dif(2,3)=0                          ! ...    vs z
         lin_dif(3,1)=0.15                       !T       vs Teff
         lin_dif(3,2)=0.3                        ! ...    vs logg
         lin_dif(3,3)=1-(temp_ref/4000)**2.0     ! ...    vs z
         lin_dif(4,1)=0.15                       !logPe   vs Teff
         lin_dif(4,2)=0.06                       ! ...    vs logg
         lin_dif(4,3)=1-(temp_ref/3500)**2.5     ! ...    vs z
         lin_dif(5,1)=-0.4                       !logPg   vs Teff
         lin_dif(5,2)=0.06                       ! ...    vs logg
         lin_dif(5,3)=1-(temp_ref/4100)**4       ! ...    vs z
         lin_dif(6,1)=0                          !xit     vs Teff
         lin_dif(6,2)=0                          ! ...    vs logg
         lin_dif(6,3)=0                          ! ...    vs z
         lin_dif(7,1)=0                          !rr      vs Teff
         lin_dif(7,2)=0                          ! ...    vs logg
         lin_dif(7,3)=0                          ! ...    vs z
         lin_dif(8,1)=-0.15                      !logxkapref vs Teff
         lin_dif(8,2)=-0.12                      ! ...    vs logg
         lin_dif(8,3)=1-(temp_ref/3700)**3.5     ! ...    vs z
         lin_dif(9,1)=0                          !Prad   vs Teff
         lin_dif(9,2)=0                          ! ...    vs logg
         lin_dif(9,3)=0                          ! ...    vs z
         lin_dif(10,1)=0                         !Pturb   vs Teff
         lin_dif(10,2)=0                         ! ...    vs logg
         lin_dif(10,3)=0                         ! ...    vs z
         lin_dif(11,1)=0                         !rhox    vs Teff
         lin_dif(11,2)=0                         ! ...    vs logg
         lin_dif(11,3)=0                         ! ...    vs z
         lin_dif(12,1)=0                         !rho     vs Teff
         lin_dif(12,2)=0                         ! ...    vs logg
         lin_dif(12,3)=0                         ! ...    vs z

         if (optimize) then
          write(*,*) 'optimized interpolation applied for standard compo
     &sition models'
         else
            lin_dif=0.
             write(*,*) 'linear interpolation applied'
         end if
!these constants are calibrated on a broad range of stellar parameters; scale them now to the present one.
            power(:,1)= 1-(lin_dif(:,1)*(abs(xsup-xinf)/(7000-3800)))
            power(:,2)= 1-(lin_dif(:,2)*(abs(ysup-yinf)/(5-0.0)))
            power(:,3)= 1-(lin_dif(:,3)*(abs(zsup-zinf)/(0-(-4))))

****** interpolation of each component of the atmosphere (taus,teff,Pe,Pg,microt,rr) and at each layer *****************

        do k=1,ndepth_ref
          x=signteffpoint*(teffpoint)**power(1,1)
          y=signloggpoint*(loggpoint)**power(1,2)
          z= signmetpoint*(metpoint)**power(1,3)
          call blend_103(x,y,z,taus(k,1),taus(k,2),
     &     taus(k,3),taus(k,4),taus(k,5),taus(k,6),taus(k,7),taus(k,8)
     &     ,taus(k,out))

          x=signteffpoint*(teffpoint)**power(2,1)
          y=signloggpoint*(loggpoint)**power(2,2)
          z= signmetpoint*(metpoint)**power(2,3)
          call blend_103(x,y,z,tauR(k,1),tauR(k,2),
     &     tauR(k,3),tauR(k,4),tauR(k,5),tauR(k,6),tauR(k,7),tauR(k,8)
     &     ,tauR(k,out))

          x=signteffpoint*(teffpoint)**power(3,1)
          y=signloggpoint*(loggpoint)**power(3,2)
          z=signmetpoint*(metpoint)**power(3,3)
          call blend_103(x,y,z,T(k,1),T(k,2),T(k,3),T(k,4)
     &     ,T(k,5),T(k,6),T(k,7),T(k,8),T(k,out))

          x=signteffpoint*(teffpoint)**power(4,1)
          y=signloggpoint*(loggpoint)**power(4,2)
          z= signmetpoint*(metpoint)**power(4,3)
          call blend_103(x,y,z,Pe(k,1),Pe(k,2),Pe(k,3),Pe(k,4)
     &     ,Pe(k,5),Pe(k,6),Pe(k,7),Pe(k,8),Pe(k,out))

          x=signteffpoint*(teffpoint)**power(5,1)
          y=signloggpoint*(loggpoint)**power(5,2)
          z= signmetpoint*(metpoint)**power(5,3)
          call blend_103(x,y,z,Pg(k,1),Pg(k,2),Pg(k,3),Pg(k,4)
     &     ,Pg(k,5),Pg(k,6),Pg(k,7),Pg(k,8),Pg(k,out))

          x=signteffpoint*(teffpoint)**power(6,1)
          y=signloggpoint*(loggpoint)**power(6,2)
          z= signmetpoint*(metpoint)**power(6,3)
          call blend_103(x,y,z,xit(k,1),xit(k,2),
     &     xit(k,3),xit(k,4),xit(k,5),xit(k,6),xit(k,7),xit(k,8)
     &     ,xit(k,out))

          x=signteffpoint*(teffpoint)**power(7,1)
          y=signloggpoint*(loggpoint)**power(7,2)
          z= signmetpoint*(metpoint)**power(7,3)
          call blend_103(x,y,z,rr(k,1),rr(k,2),
     &     rr(k,3),rr(k,4),rr(k,5),rr(k,6),rr(k,7),rr(k,8)
     &     ,rr(k,out))

          x=signteffpoint*(teffpoint)**power(8,1)
          y=signloggpoint*(loggpoint)**power(8,2)
          z= signmetpoint*(metpoint)**power(8,3)
          call blend_103(x,y,z,xkapref(k,1),xkapref(k,2),
     &     xkapref(k,3),xkapref(k,4),xkapref(k,5),xkapref(k,6),
     &     xkapref(k,7),xkapref(k,8),xkapref(k,out))

          x=signteffpoint*(teffpoint)**power(9,1)
          y=signloggpoint*(loggpoint)**power(9,2)
          z= signmetpoint*(metpoint)**power(9,3)
          call blend_103(x,y,z,Prad(k,1),Prad(k,2),
     &     Prad(k,3),Prad(k,4),Prad(k,5),Prad(k,6),
     &     Prad(k,7),Prad(k,8),Prad(k,out))

          x=signteffpoint*(teffpoint)**power(10,1)
          y=signloggpoint*(loggpoint)**power(10,2)
          z= signmetpoint*(metpoint)**power(10,3)
          call blend_103(x,y,z,Pturb(k,1),Pturb(k,2),
     &     Pturb(k,3),Pturb(k,4),Pturb(k,5),Pturb(k,6),
     &     Pturb(k,7),Pturb(k,8),Pturb(k,out))

          x=signteffpoint*(teffpoint)**power(11,1)
          y=signloggpoint*(loggpoint)**power(11,2)
          z= signmetpoint*(metpoint)**power(11,3)
          call blend_103(x,y,z,rhox(k,1),rhox(k,2),
     &     rhox(k,3),rhox(k,4),rhox(k,5),rhox(k,6),
     &     rhox(k,7),rhox(k,8),rhox(k,out))

          x=signteffpoint*(teffpoint)**power(12,1)
          y=signloggpoint*(loggpoint)**power(12,2)
          z= signmetpoint*(metpoint)**power(12,3)
          call blend_103(x,y,z,ro(k,1),ro(k,2),
     &     ro(k,3),ro(k,4),ro(k,5),ro(k,6),
     &     ro(k,7),ro(k,8),ro(k,out))
       end do
       do k=1,92
         abu(k,out)=metpoint*(abu(k,2)-abu(k,1))+abu(k,1)
       enddo
      ndepth(out)=ndepth_ref
      xlr(out)=lambda_ref
      sph(out)=sph(1)

**do not work yet *now calculate rhox*interpolation is done instead********
c     write(*,*) 'now calculate rhox'
c     allocate(rhox(ndepth_ref,nfile))
c     do file=1,9
c     call calcrhox(tauR(:,file),xkapref(:,file),ndepth_ref,
c    &                                                 rhox(:,file))
c     enddo

**********calculate estimated error********
      if (optimize) then
         write(*,*) 'now calculate error'
      call calc_error(xinf,xsup,yinf,ysup,zinf,zsup,temp_ref,
     & logg_ref,z_ref)
      endif



******** write interpolated model in file nber out ***********
      write(*,*) 'now write result in ',trim(outtype),'format:'
      write(*,*) trim(FILE_IN(out))
      open(unit=out,file=FILE_IN(out))

      select case (adjustl(adjustr(outtype)))
c  ----output in turbospectrum format
      case ('TurboSpectrum')
      call write_turbo(out,ndepth_ref,lambda_ref,logg_ref,sph(out),
     &                     taus(:,out),T(:,out),Pe(:,out),Pg(:,out),
     &                 xit(:,out),rr(:,out),rhox(:,out),taur(:,out),
     &                                              ro(:,out),FILE_IN)
c  ----output in turbospectrum format
      case ('MARCSweb')
      mass=0.0
      call write_MARCSweb(out,FILE_IN,mass,xit(:,out),sph(out),extrapol,
     &  temp_ref,logg_ref,z_ref,rr(:,out),ndepth_ref,
     &    tauR(:,out),taus(:,out),T(:,out),Pe(:,out),Pg(:,out),
     &  Prad(:,out),Pturb(:,out),
     &               xkapref(:,out),rhox(:,out),ro(:,out),abu(:,out))

C------output in ATLAS9 format
      case ('Atlas')
      call write_atlas(out,ndepth_ref,abu(:,out),temp_ref,
     &  logg_ref,z_ref,Pe(:,out),Prad(:,out),rhox(:,out),T(:,out),
     &    Pg(:,out),xit(:,out),xkapref(:,out))


c------output in MOOG format
       case ('MOOG')
       call write_MOOG(out,sph(1),temp_ref,logg_ref,z_ref,ndepth_ref,
     &            lambda_ref,taus(:,out),T(:,out),Pe(:,out),Pg(:,out),
     &     rhox(:,out),xit(:,out))


c------output default
       case default
       write(*,*) 'no specific output format selected:
     & default format is used'
       call write_default(out,temp_ref,logg_ref,z_ref,ndepth_ref,
     &        lambda_ref,tauR(:,out),taus(:,out),T(:,out),Pe(:,out),
     &        Pg(:,out),Prad(:,out),xit(:,out),rhox(:,out),rr(:,out),
     &            ro(:,out))
       end select
       close(out)


********** write a file compatible for sm (used for a control plot) ************
      open (unit=24,file='models.sm')
      do file=1,9
       call write_default(24,teff(file),logg(file),metal(file),
     &      ndepth_ref,
     &        lambda_ref,tauR(:,file),taus(:,file),T(:,file),
     &       Pe(:,file),Pg(:,file),Prad(:,file),xit(:,file),
     &    rhox(:,file),rr(:,file),ro(:,file))
      end do
!case of a 10th comparison model
      if (test) then
      file = 10
      read(*,*) FILE_IN(file)
      write(*,*) 'load comparison file', trim(FILE_IN(file))
      read(*,*) modtype
      select case (adjustl(adjustr(modtype)))
      case ('MARCSbin')
      call extract_bin(FILE_IN(file),teff(file),logg(file),metal(file),
     & ndepth(file),xlr(file),taus_aux(:,file),tauR_aux(:,file),
     & T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),xit_aux(:,file),
     & rr_aux(:,file),sph(file),xkapref_aux(:,file),
     & Prad_aux(:,file),Pturb_aux(:,file),rhox_aux(:,file),
     & ro_aux(:,file),abu(:,file))
      case( 'MARCSweb')
      call extract_ascii(FILE_IN(file),teff(file),logg(file),
     & metal(file),ndepth(file),xlr(file),taus_aux(:,file),
     & tauR_aux(:,file),T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),
     & xit_aux(:,file),rr_aux(:,file),sph(file),xkapref_aux(:,file),
     & Prad_aux(:,file),Pturb_aux(:,file),rhox_aux(:,file),
     & ro_aux(:,file),abu(:,file))
      case ('Phoenix')
      call extract_Phoenix(FILE_IN(file),teff(file),logg(file),
     & metal(file),ndepth(file),xlr(file),taus_aux(:,file),
     & tauR_aux(:,file),T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),
     & xit_aux(:,file),rr_aux(:,file),sph(file),xkapref_aux(:,file),
     & Prad_aux(:,file),Pturb_aux(:,file),rhox_aux(:,file),
     & ro_aux(:,file),abu(:,file))
      case ('Atlas')
      call extract_Atlas(FILE_IN(file),teff(file),logg(file),
     & metal(file),ndepth(file),xlr(file),taus_aux(:,file),
     & tauR_aux(:,file),T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),
     & xit_aux(:,file),rr_aux(:,file),sph(file),xkapref_aux(:,file),
     & Prad_aux(:,file),Pturb_aux(:,file),rhox_aux(:,file),
     & ro_aux(:,file),abu(:,file))
      case default
       write(*,*) 'model type unknown:',outtype
      end select
c      call calcrhox(tauR_aux(:,file),xkapref_aux(:,file),ndepth(file),
c     &                                                 rhox(:,file))

       call write_default(24,teff(file),logg(file),metal(file),
     &    ndepth(file),xlr(file),tauR_aux(:,file),taus_aux(:,file),
     &   T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),Prad_aux(:,file),
     & xit_aux(:,file),rhox(:,file),rr_aux(:,file),ro_aux(:,file))
      end if
      close(24)

      if (extrapol) then
          write (*,*) 'extrapolation done'
          else
      write (*,*) 'interpolation done'
      end if
      end if

      deallocate(taus,tauR,T,Pe,Pg,xit,taus_aux,tauR_aux,T_aux,Pe_aux,
     & Pg_aux,xit_aux,rr_aux,rr,rhox,rhox_aux,Prad,ro,ro_aux)

      end


c----------------------------------------------------------------------------------------



c--------------------------------------------------------------------------------
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
c--------------------------------------------------------------------------------














c---------------------------------------------------------------------------------
      subroutine extract_bin(FILE,TEFF,grav,metal,ndepth,xlr_ref,tau5,
     &                tauR,temp,prese,presg,xit,rad,sph,kappa,presrad,
     &                 presturb,rhox,rho,abu)

!extracted from osplot.f 07/2003, to get tau,T,Pe,Pg,microturb from a model+ rhox 11/2018
      implicit none
      integer :: ndp,ndepth,k,nlp,nlb
      parameter(ndp=200)
      CHARACTER*117 ADUM
      CHARACTER*256 FILE,file2
      CHARACTER*30 COMMENT
      real :: abu(92)
      real metal,grav,xlr_ref,TEFF,GG,radius,
     &         mass_tot,mass_H,mass_He,mass_Z
      logical :: sph
      real :: tau(ndp),t(ndp),z(ndp),ptot(ndp),prad(ndp),
     &  pg(ndp),pturb(ndp),pe(ndp),ro(ndp),xmass(ndp),xkapr(ndp)
      real :: gradp(ndp),gravity(ndp),pcheck(ndp),rr(ndp),
     &   xit(ndp),geff(ndp),gradptur(ndp),dp(ndp),taus(ndp),xlr(30),
     &   coldens(ndp),rhoxx(ndp)
      real :: tau5(ndp),tauR(ndp),temp(ndp),prese(ndp),rhox(ndp),
     &   presg(ndp),presrad(ndp),presturb(ndp),rad(ndp),kappa(ndp),
     &   rho(ndp)
      real :: xlb(155000),w(155000),fluxme(155000)
      real :: presmo(30,ndp),ptio(ndp)
      real :: bPPR(NDP),bPPT(NDP),bPP(NDP),bGG(NDP),
     & bZZ(NDP),bDD(NDP),
     & bVV(NDP),bFFC(NDP),bPPE(NDP),bTT(NDP),
     & bTAULN(NDP),erad(ndp)
      integer :: NbTAU,IbTER
      common /struct/ tau,t,z,ptot,prad,pg,pturb,pe,ro,rr,taus,xlr,
     &                nlp,xkapr,rhoxx
      common /spectr/ nlb,xlb,w,fluxme
      common /pressure/ presmo,ptio
      common /binstruc/bPPR,bPPT,bPP,bGG,
     & bZZ,bDD,
     & bVV,bFFC,bPPE,bTT,
     & bTAULN,NbTAU,IbTER,erad
      common /radius/ radius

      OPEN(UNIT=10,FILE=FILE,STATUS='OLD',FORM='UNFORMATTED')
c     &     convert='big_endian')
ccc     &     RECL=152600)
*
      CALL READMO(10,NDEPTH,TEFF,GG,metal,abu,sph)
c         open(21,file=file2,status='unknown')

      do k=1,ndepth
        tau5(k)=log10(taus(k))
        tauR(k)=log10(tau(k))
         temp(k)=t(k)
         prese(k)=log10(pe(k))
         presg(k)= log10(pg(k))
         presrad(k)=log10(prad(k))
         presturb(k)= log10(pturb(k))
         kappa(k)=log10(xkapr(k))
         rhox(k)=rhoxx(k)
         rho(k)=ro(k)
      end do
         xit=2.0
         xlr_ref=xlr(nlp)
         grav=log10(GG)
*        rad=rr
         if(.not.sph) radius=0.0
         rad=radius-z
         mass_Z=0
      close(10)
      END
C
      SUBROUTINE READMO(IARCH,JTAU,TEFF,G,metal,abu,spherical)
C        THIS ROUTINE READS ONE MODEL, TO GET INFO ON PRAD
C             ( All features taken from listmo )
      PARAMETER (NDP=200)
C
      DIMENSION ABUND(16),TKORRM(NDP),FCORR(NDP),TAU(NDP),TAUS(NDP),
     *T(NDP),PE(NDP),PG(NDP),PRAD(NDP),PTURB(NDP),XKAPR(NDP),RO(NDP),
     *CP(NDP),CV(NDP),AGRAD(NDP),Q(NDP),U(NDP),V(NDP),ANCONV(NDP),
     *PRESMO(30,NDP),FCONV(NDP),RR(NDP),Z(NDP),EMU(NDP),HNIC(NDP)
     *,NJ(16),XLR(30),IEL(16),ANJON(16,5),PART(16,5),PROV(50,20+1),
     *ABSKA(50),SPRIDA(50),XLB(155000),FLUXME(155000),FLUMAG(155000),
     & PEP(16),FLUXMEC(155000),RHOXX(NDP),
     * ABNAME(50),SOURCE(50),PTOT(NDP)
      DIMENSION W(155000),UW(12),BW(21),VW(25)
      CHARACTER*10 DAG,NAME,NAMEP,KLOCK
      CHARACTER*8 ABNAME,SOURCE
      DIMENSION WAVFLX(10)
      INTEGER readpb
      dimension PTIO(NDP)
      real*8 dluminosity
      real abSc,abTi,abV,abMn,abCo,metal
      real inputabund(92),abu(92)
      character*2 inputaname(92)
      logical spherical
      common /binstruc/ dummy(11*ndp+2),erad(ndp)
      common /struct/ tau,t,z,ptot,prad,pg,pturb,pe,ro,rr,taus,xlr,
     &                nlp,xkapr,rhoxx
      common /spectr/ nlb,xlb,w,fluxme
      common /pressure/ presmo,ptio
      common /radius/ radius
      DATA UW/0.145,0.436,0.910,1.385,1.843,2.126,2.305,2.241,1.270,
     *0.360,0.128,0.028/,BW/0.003,0.026,0.179,0.612,1.903,2.615,2.912,
     *3.005,2.990,2.876,2.681,2.388,2.058,1.725,1.416,1.135,0.840,0.568,
     *0.318,0.126,0.019/,VW/0.006,0.077,0.434,1.455,2.207,2.703,2.872,
     *2.738,2.505,2.219,1.890,1.567,1.233,0.918,0.680,0.474,0.312,0.200,
     *0.132,0.096,0.069,0.053,0.037,0.022,0.012/
      DATA NAME/'LOCAL'/,NAMEP/'PARSONS'/
      DATA A,B/.34785485,.65214515/
      IREAD=5
C
      do i=1,2
        read(iarch)
      enddo
      read(iarch)
     &     teff,flux,g,palfa,pny,py,pbeta,iline,istral,mihal,
     &            idrab1,idrab2,idrab3,idrab4,idrab5,idrab6,
     &            itmax,nel
      do i=1,2
        read(iarch) jtau
      enddo
      ntpo=0
      do k=1,jtau
        read(iarch) kr,tau(k)
        tauk=alog10(tau(k))+10.01
        ktau=tauk
        if(abs(tauk-ktau).gt.0.02) go to 131
        if(ktau.eq.10) k0=k
        ntpo=ntpo+1
  131   continue
      enddo
      read(iarch)(nj(i),i=1,nel),nlp
      do k=1,ntpo
        do j=1,nlp
          if (j.eq.1) then
            do i=1,nel
              read(iarch)
            enddo
          endif
          read(iarch)
        enddo
      enddo
      read(iarch,err=191) nlb,(xlb(j),fluxme(j),j=1,nlb),(w(j),j=1,nlb)
      read(iarch,end=191,err=191,iostat=readpb) (fluxmec(j),j=1,nlb)
      read(iarch,err=193,end=193,iostat=readpb) inputabund,inputaname
      if (readpb.ne.0) then
 191  write(*,*) 'error during binary read'
      stop
      endif
 193  metal=inputabund(26)-7.45
      alpha=inputabund(22)-4.90+metal
      abu=inputabund

      REWIND IARCH
      READ(IARCH)
      erad=-1.e30
      READ(IARCH) INORD,DAG,KLOCK
      READ(IARCH) TEFF,FLUX,G,PALFA,PNY,PY,PBETA,ILINE,ISTRAL,
     &                MIHAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &            ITMAX,NEL,(ABUND(I),I=1,NEL),abSc,abTi,abV,abMn,abCo
      GLOG=ALOG10(G)
      FNORD=0.1*INORD
C        CONVERT TO 'PHYSICAL FLUX'
      FLUX=3.14159*FLUX
      DO 2 I=1,NEL
    2 ABU(I)=ALOG10(ABUND(I))+12.
      READ(IARCH)JTAU,NCORE,DIFLOG,TAUM,RADIUS,(RR(K),K=1,JTAU)
      if (jtau.gt.ndp) then
        print*, 'ERROR !!! Jtau (number of depths of model) = ',jtau
        print*, ' is larger than ndp!! Increase NDP.'
        stop
      endif
      if (radius.le.2.) then
         spherical = .false.
      else
         spherical = .true.
         rr=radius-rr
      endif

      READ(IARCH)JTAU,(TKORRM(I),I=1,JTAU),(FCORR(K),K=1,JTAU)
      NTPO=0
      DO 3 K=1,JTAU
        READ(IARCH) KR,TAU(K),TAUS(K),Z(K),T(K),PE(K),PG(K),PRAD(K),
     &              PTURB(K),XKAPR(K),RO(K),EMU(K),CP(K),CV(K),
     &              AGRAD(K),Q(K),U(K),V(K),ANCONV(K),HNIC(K),NMOL,
     &              (PRESMO(J,K),J=1,NMOL),ptio(k)
        TAUK=ALOG10(TAU(K))+10.01
        KTAU=TAUK
        IF(ABS(TAUK-KTAU).GT.0.02) GO TO 31
        IF(KTAU.EQ.10) K0=K
        NTPO=NTPO+1
   31   CONTINUE
    3 CONTINUE
c      Z0=Z(K0)
c      DO 5 I=1,JTAU
c        Z(I)=Z(I)-Z0
c        i1=min(i+1,jtau)
        PTOT(I)=PG(I)+PRAD(I)+0.5*(pturb(i)+pturb(i1))
    5 CONTINUE
      do I=1,JTAU
       RHOXX(I)=PTOT(I)/G
      enddo
***
      READ(IARCH)(NJ(I),I=1,NEL),NLP,(XLR(I),I=1,NLP)
     & ,NPROV,NPROVA,NPROVS,(ABNAME(KP),SOURCE(KP),KP=1,NPROV)
c      DO 22 KTAU=1,NTPO
c      DO 20 IE=1,NEL
c      NJP=NJ(IE)
c      READ(IARCH) KR,TAUI,TI,PEI,IEL(IE),ABUND(IE),
c     &            (ANJON(IE,JJ),JJ=1,NJP),(PART(IE,JJ),JJ=1,NJP)
c   20 CONTINUE
c      DO 21 KLAM=1,NLP
c      READ(IARCH) KR,TAUIL,(PROV(J,KLAM),J=1,NPROV),
c     &            ABSKA(KLAM),SPRIDA(KLAM)
   21 CONTINUE
   22 continue
c      READ(IARCH) NLB,(XLB(J),FLUXME(J),J=1,NLB),(W(J),J=1,NLB)
C CONVERT TO 'PHYSICAL' FLUXES
c      DO 24 J=1,NLB
c   24 FLUXME(J)=3.14159*FLUXME(J)

c      dluminosity=0.
c      do 25 j=1,nlb
c       dluminosity=dluminosity+fluxme(j)*w(j)
25    continue
c      dluminosity=dluminosity*4.*3.14159*radius**2/3.82d33
c      ddddd=real(dluminosity)
c      print*,'luminosity: ',dluminosity*3.82d33,' erg/s  = ',ddddd,
c     &  ' solar luminosities'
***
      RETURN
         END
*****************************************************************************

c---------------------------------------------------------------------------------
      subroutine extract_ascii(FILE,TEFF,grav,metal,ndepth,xlr_ref,tau5,
     &                tauR,temp,prese,presg,xit,rad,sph,xkapr,pres_rad,
     &                 pres_turb,rhox,ro,abu)
c     adapted from P. DeLaverny
      implicit none
      integer :: ndp,k
      parameter(ndp=200)
      integer :: imod,idum,ndepth
      CHARACTER*117 ADUM
      CHARACTER*256 FILE
      CHARACTER*30 COMMENT
      CHARACTER*50 MOCODE
      real metal,radius,mass,grav,xlr_ref,TEFF,GG,xic
      logical :: sph
      real :: tau(ndp),t(ndp),z(ndp),ptot(ndp),prad(ndp),
     &  pg(ndp),pturb(ndp),pe(ndp),ro(ndp),xmass(ndp),xkapr(ndp)
      real :: gradp(ndp),gravity(ndp),pcheck(ndp),rr(ndp),
     &  xit(ndp),geff(ndp),gradptur(ndp),dp(ndp),taus(ndp),xlr(30),
     &   coldens(ndp)
      real :: tau5(ndp),tauR(ndp),temp(ndp),prese(ndp),
     &   presg(ndp),rad(ndp),emu(ndp),vconv(ndp),fconv(ndp),
     &   pres_rad(ndp),Rhox(ndp),pres_turb(ndp)
      real :: dimension xlb(155000),w(155000),fluxme(155000)
      real :: abu(92),abu2(6)
      real :: mass_He,mass_H,mass_Z




          imod =10
          OPEN(UNIT=imod,FILE=FILE,STATUS='OLD')
          read(imod,'(a)') mocode
c          print*,mocode,' = mocode'
          if (mocode(1:1).eq.'p' .or. mocode(1:3).eq.'sun') then
            print*,' this model is PLANE PARALLEL'
          else if (mocode(1:1).eq.'s') then
            print*,' this model is SPHERICALLY SYMMETRIC'
          else
            print*,' This model may not be a NewMARCS model!'
          endif
          sph=(mocode(1:1).eq.'s')
          xlr_ref=5000
          read(imod,*)TEFF
          read(imod,*)
          read(imod,*)grav
          grav=log10(grav)
          read(imod,*)xic
          read(imod,*)mass
          read(imod,*)metal
          read(imod,*)radius
          read(imod,*)
          read(imod,*)
	  read(imod,*) mass_H,mass_He,mass_Z
	  read(imod,*)
	  read(imod,'(10f7.2)') abu
          read(imod,*) ndepth
          read(imod,*)
          read(imod,*)
          do k=1,ndepth
            read(imod,*) idum,tauR(k),tau5(k),rad(k),temp(k),
     &                   Pe(k),Pg(k),Prad(k),Pturb(k)
c            write(*,*)  idum,tauR(k),tau5(k),rad(k),temp(k),
c     &                   Pe(k),Pg(k),Prad(k),Pturb(k)
            prese(k)  = log10(Pe(k))
            presg(k)  = log10(Pg(k))
            pres_rad(k) = log10(Prad(k))
            pres_turb(k) = log10(Pturb(k))
            xit(k) = xic
            rad(k)=radius-rad(k)

          enddo
          read(imod,*)
          do k=1,ndepth
            read(imod,870) idum,tauR(k),xkapr(k),ro(k),emu(k),
     &                   Vconv(k),Fconv(k),rhox(k)
 870        format(i3,x,f5.2,e12.4,e12.4,f6.3,e11.3,f8.5,e14.6)
            xkapr(k) = log10(xkapr(k))
          enddo
         close(imod)
      END

CPhoenix type models ----------------------------------------
      subroutine extract_Phoenix(FILE,TEFF,grav,metal,ndepth,xlr_ref,
     &        tau5,tauR,temp,prese,presg,xit,rad,sph,xkapr,pres_rad,
     &                 pres_turb,rhox,rho,abu)

      implicit none
      integer :: ndp,ndepth,n,j,k,name(9)
      parameter(ndp=200)
      CHARACTER*256 FILE
      character*46 :: string
      character*5:: dum
      real :: metal,radius,mass,grav,xlr_ref,TEFF,microt
      logical :: sph
      real :: tau5(ndp),tauR(ndp),taubis,temp(ndp),prese(ndp),
     &   presg(ndp),rad(ndp),rho(ndp),xkapr(ndp)
      real :: tau(ndp),tauross(ndp),t(ndp),pe(ndp),xkappar(ndp),
     &   pg(ndp),r(ndp),pres_rad(ndp),xit(ndp),rr(ndp),
     &   prad(ndp),Rhox(ndp),Pturb(ndp),pres_turb(ndp),mu(ndp)
      real :: abu(92),massfrac(92),abu2(6)
      real :: mass_He,mass_H,mass_Z

      write(*,*) 'YOU NEED TO IMPLEMENT READING OF RHOX AND XKAPR!!!!!!'
      stop

      OPEN(UNIT=10,FILE=FILE,STATUS='OLD')
      sph=.true.
      do while (trim(adjustl(adjustr(string)))
     &   .ne. 'Element abundances :')
      read(10,'(a32)') string
      enddo
      do
      read(10,*,err=1147)
      read(10,*,err=1147)
      read(10,*,err=1147) dum,(name(j),j=1,9)
      read(10,*,err=1147) dum,(abu2(j),j=1,9)
      read(10,*,err=1147)
      read(10,*,err=1147) dum,(massfrac(j),j=1,9)
      do j=1,9
      k=int(name(j)/100)
      if (k>=1.and.k<=92) then
      abu(k)=abu2(j)
      endif
      enddo
      enddo

 1147 metal=abu(26)-7.52

      rewind(10)
      do while (trim(adjustl(adjustr(string)))
     &   .ne. 'statistical velocity (km/sec):')
      read(10,'(a32)') string
      enddo
      backspace(10)
      read(10,'(a46)') string
      string=adjustl(adjustr(string(index(string,':',.true.)+1:)))
      read(string,'(f5.2)') microt
      rewind(10)
      do while (trim(adjustl(adjustr(string)))
     &   .ne. 'Mie mostly converged: j,err*10^6,w,a,scat,abs')
      read(10,'(a46)') string
      enddo
      backspace(10)
      read(10,*) dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,xlr_ref
      if (xlr_ref < 1000) then
      xlr_ref=xlr_ref*10000
      endif
      rewind(10)
      do while (trim(adjustl(adjustr(string)))
     &   .ne. 'teff')
      read(10,*) dum,string
      enddo
      backspace(10)
      read(10,*) dum,dum,dum,TEFF
      do while (trim(adjustl(adjustr(string)))
     &   .ne. 'log(g):')
      read(10,*) string
      enddo
      backspace(10)
      read(10,*) dum,grav
      rewind(10)
      do while (trim(adjustl(adjustr(string)))
     &   .ne. 'phoenix: saving structure')
      read(10,'(a40)') string
      enddo
      read(10,*)
      tau(1)=0.0
      do while (tau(1)<= 1.D-10)
 1191    read(10,*,err=1191) dum,tau(1)
      enddo
      backspace(10)
      n=1
      do
       read(10,*,err=1196) dum,tau(n),T(n),Pg(n),Pe(n),rho(n),mu(n),
     &    rr(n)
       n=n+1
      enddo
 1196 n=n-1


      do while (trim(adjustl(adjustr(string))).ne.
     & 'Pressure and temperature structure:')
      read(10,'(a40)') string
      enddo
      do while (taubis.ne.tau(1))
 1205    read(10,*,err=1205) dum,taubis
      enddo
      backspace(10)
      k=1
      do
      read(10,*,err=1218) dum,taubis,dum,dum,dum,Prad(k),Pturb(k)
      if (taubis.ne.tau(k)) then
      write(6,*) 'we are no reading the same tau scale in Pres and temp'
      stop
      endif
      k=k+1
      enddo
 1218 continue

      do while (trim(adjustl(adjustr(string))).ne.
     & 'opacity averages in the comoving frame:')
      read(10,'(a40)') string
      enddo
      do while (taubis.ne.tau(1))
 1226    read(10,*,err=1226) dum,taubis
      enddo
      backspace(10)
      k=1
      do
      read(10,*,err=1231) dum,taubis,dum,dum,dum,dum,dum,dum,dum,
     &    xkappar(k)
      if (taubis.ne.tau(k)) then
      write(6,*) 'we are no reading the same tau scale in opac'
      stop
      endif
      k=k+1
      enddo

 1231 backspace(10)
      do while (trim(adjustl(adjustr(string))).ne.
     & 'tau averages in the comoving frame:')
      read(10,'(a40)') string
      enddo

      do while (taubis.ne.tau(1))
 1222    read(10,*,err=1222) dum,taubis
      enddo
      backspace(10)
      k=1
      do
      read(10,*,err=1229) dum,taubis,dum,dum,dum,tauross(k)
      if (taubis.ne.tau(k)) then
      write(6,*) 'we are no reading the same tau scale in frame'
      stop
      endif
      k=k+1
      enddo
 1229 k=k-1
      if (k<n) then
      n=k
      else
      if  (k>n) then
      k=n
      endif
      endif
      ndepth=n
       do k=1,ndepth
       n=k
       temp(k)=T(n)
       Presg(k)=log10(Pg(n))
       Prese(k)=log10(Pe(n))
       Pres_rad(k)=log10(Prad(n))
       Pres_turb(k)= log10(Pturb(n))
       tau5(k)=log10(tau(n))
       tauR(k)=log10(tauross(n))
       xit(k)=microt
       xkapr(k) = log10(xkappar(n))
       rad(k)=rr(n)
       enddo
      close(10)
      end

c Atlas model
       subroutine extract_Atlas(FILE,TEFF,grav,metal,ndepth,xlr_ref,
     &        tau5,tauR,temp,prese,presg,xit,rad,sph,xkapr,pres_rad,
     &                 pres_turb,rhox,ro,abu)

* Kurucz models. Reading + rinteg taken from moog. 06/04-2001 BPz+ST (Sivarani)
       implicit none
      integer :: k,ndp,ndepth,i,neletstop,nelet(6)
      parameter(ndp=200)
      CHARACTER*256 FILE
      character*46 :: string
      character*1:: dum
      real metal,radius,mass,grav,xlr_ref,TEFF
      logical :: sph
      real :: tau5(ndp),tauR(ndp),taubis(ndp),temp(ndp),prese(ndp),
     &   presg(ndp),rad(ndp),kaprefmass(ndp),xit(ndp)
      real :: tau(ndp),tauross(ndp),t(ndp),pe(ndp),
     &   pg(ndp),r(ndp),pres_rad(ndp),xkapr(ndp),
     &   prad(ndp),Rhox(ndp),Pturb(ndp),pres_turb(ndp),
     &    flxconv(ndp),vconv(ndp),ro
      real :: abu(92),massfrac(92),abu2(6)
      real :: mass_He,mass_H,mass_Z
      real :: tottau,accrad,first
      real,external :: rinteg
c!!!!!!!!!!!!!!!!!!!!!!!!!!
c NOT READY -- TO BE CHECKED!!!!!!!!!!!!1
c!!!!!!!!!!!!!!!!!!!!!!!!!!111
      write(*,*) 'KURUCZ MODEL NOT READY '
      stop
       open(unit=10,file=FILE,status='old')
cTEFF   3500.  GRAVITY 0.00000 LTE
cTITLE  [0.0a] VTURB=2  L/H=1.25 NOVER NEW ODF
c OPACITY IFOP 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0
c CONVECTION ON   1.25 TURBULENCE OFF  0.00  0.00  0.00  0.00
cABUNDANCE SCALE   1.00000 ABUNDANCE CHANGE 1 0.91930 2 0.07824
c ABUNDANCE CHANGE  3 -10.94  4 -10.64  5  -9.49  6  -3.52  7  -4.12  8  -2.81
c ABUNDANCE CHANGE  9  -7.48 10  -3.56 11  -5.71 12  -4.06 13  -5.57 14  -4.09
c ABUNDANCE CHANGE 15  -6.59 16  -4.31 17  -6.54 18  -5.24 19  -6.92 20  -5.28
c ABUNDANCE CHANGE 21  -8.87 22  -6.62 23  -8.04 24  -6.37 25  -6.65 26  -4.54
c ABUNDANCE CHANGE 27  -7.12 28  -5.79 29  -7.83 30  -7.44 31  -9.16 32  -8.63
c ABUNDANCE CHANGE 33  -9.67 34  -8.63 35  -9.41 36  -8.73 37  -9.44 38  -9.07
c ABUNDANCE CHANGE 39  -9.80 40  -9.44 41 -10.62 42 -10.12 43 -20.00 44 -10.20
c ABUNDANCE CHANGE 45 -10.92 46 -10.35 47 -11.10 48 -10.27 49 -10.38 50 -10.04
c ABUNDANCE CHANGE 51 -11.04 52  -9.80 53 -10.53 54  -9.87 55 -10.91 56  -9.91
c ABUNDANCE CHANGE 57 -10.87 58 -10.46 59 -11.33 60 -10.54 61 -20.00 62 -11.03
c ABUNDANCE CHANGE 63 -11.53 64 -10.92 65 -11.69 66 -10.90 67 -11.78 68 -11.11
c ABUNDANCE CHANGE 69 -12.04 70 -10.96 71 -11.98 72 -11.16 73 -12.17 74 -10.93
c ABUNDANCE CHANGE 75 -11.76 76 -10.59 77 -10.69 78 -10.24 79 -11.03 80 -10.91
c ABUNDANCE CHANGE 81 -11.14 82 -10.09 83 -11.33 84 -20.00 85 -20.00 86 -20.00
c ABUNDANCE CHANGE 87 -20.00 88 -20.00 89 -20.00 90 -11.95 91 -20.00 92 -12.54
c ABUNDANCE CHANGE 93 -20.00 94 -20.00 95 -20.00 96 -20.00 97 -20.00 98 -20.00
c ABUNDANCE CHANGE 99 -20.00
cREAD DECK6 72 RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB, FLXCNV,VCONV,VELSND
        sph=.false.
        read(10,*) dum,teff,dum,grav
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)  dum,dum,metal
        metal=log10(metal)
        do while (neletstop.ne.99)
        read(10,*) dum,dum,(nelet(i),abu2(i),i=1,6)
        do i=1,6
        abu(nelet(i))=abu2(i)+12.04+metal
        enddo
        neletstop=nelet(1)
        enddo
        read(10,*) dum,dum,ndepth

        do k=1,ndepth
          read (10,*,err=1298) rhox(k),t(k),pg(k),
     &    pe(k),kaprefmass(k),accrad,xit(k),flxconv(k),vconv(k)
          pe(k)=pe(k)*T(k)*1.38054e-16
        enddo
 1298   first = rhox(1)*kaprefmass(1)
        tottau = rinteg(rhox,kaprefmass,tau,ndepth,first)

        do k=1,ndepth
           if (k.eq.1) then
           tauR(k)= rhox(1)*kaprefmass(1)
           else
            tauR(k) = tau(k-1) + tau(k)
            tau5(k) = 0.0  !! not defined in Atlas,need to compute it
            prese(k)  = log10(Pe(k))
            presg(k)  = log10(Pg(k))
            pres_rad(k) = 0.0        !! to be defined
            pres_turb(k) = 0.0      !to be defined withh fluxconv,vconv
            xit(k) = xit(k)
            rad(k)=0.0
           endif
         enddo
         close(10)
        end




c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine resample(taus,tauR,T,Pe,Pg,Prad,xit,rr,xkapref,rhox,ro)
      implicit  none
      integer :: nlinemod,file,nfile,k,i
      real,dimension(:,:) :: taus,tauR,T,Pe,Pg,Prad,xit,rr,xkapref,rhox
     &                                                           ,ro
      real,dimension(:,:),allocatable :: taubas
      real,dimension(:),allocatable :: tauresample,taustemp,tauRtemp,
     & Ttemp,Petemp,Pgtemp,Pradtemp,xittemp,rrtemp,xkapreftemp,rhoxtemp,
     & rotemp

      INTERFACE
      function SevalSingle(u,x,y)
      REAL,INTENT(IN) :: u  ! abscissa at which the spline is to be evaluated
      REAL,INTENT(IN),DIMENSION(:) :: x ! abscissas of knots
      REAL,INTENT(IN),DIMENSION(:):: y ! ordinates of knots
      real SevalSingle
      end
      END INTERFACE

      nlinemod=size(taus,1)
      nfile=size(taus,2)-3

      allocate(tauresample(nlinemod),taustemp(nlinemod),
     &  tauRtemp(nlinemod),Ttemp(nlinemod)
     & ,Petemp(nlinemod),Pgtemp(nlinemod),Pradtemp(nlinemod),
     &   xittemp(nlinemod),rhoxtemp(nlinemod),rotemp(nlinemod),
     & rrtemp(nlinemod),xkapreftemp(nlinemod),
     &   taubas(nlinemod,size(taus,2)))

!!!!choose here the depth basis for interpolation (tau5000,tauRoss)!!!!!
      taubas=tauR
      write(*,*) 'resample models on common depth basis: tauRoss'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call common_depth_basis(tauresample,taubas,nlinemod,nfile)
c now do the resampling with the common tau
       do file=1,nfile
         do k=1,nlinemod
       taustemp(k)=SevalSingle(tauresample(k),taubas(:,file)
     &                                                  ,taus(:,file))

       tauRtemp(k)=SevalSingle(tauresample(k),taubas(:,file)
     &                                                  ,tauR(:,file))
       Ttemp(k)=SevalSingle(tauresample(k),taubas(:,file),T(:,file))
       Petemp(k)=SevalSingle(tauresample(k),taubas(:,file),Pe(:,file))
       Pgtemp(k)=SevalSingle(tauresample(k),taubas(:,file),Pg(:,file))
       Pradtemp(k)=SevalSingle(tauresample(k),taubas(:,file)
     &                                                   ,Prad(:,file))
       xittemp(k)=SevalSingle(tauresample(k),taubas(:,file),xit(:,file))
       rrtemp(k)=SevalSingle(tauresample(k),taubas(:,file),rr(:,file))
       xkapreftemp(k)=
     &        SevalSingle(tauresample(k),taubas(:,file),xkapref(:,file))
       rhoxtemp(k)=
     &        SevalSingle(tauresample(k),taubas(:,file),rhox(:,file))
       rotemp(k)=
     &        SevalSingle(tauresample(k),taubas(:,file),ro(:,file))
         end do

         taus(:,file)=taustemp
         tauR(:,file)=tauRtemp
         T(:,file)=Ttemp
         Pe(:,file)=Petemp
         Pg(:,file)=Pgtemp
         Prad(:,file)=Pradtemp
         xit(:,file)=xittemp
         rr(:,file)=rrtemp
         xkapref(:,file)=xkapreftemp
         rhox(:,file)=rhoxtemp
         ro(:,file)=rotemp
       end do

       deallocate(tauresample,taustemp,tauRtemp,Ttemp,Petemp
     &,Pgtemp,Pradtemp,xittemp,rrtemp,xkapreftemp,taubas)

       end


!*******************************************************************************

      subroutine common_depth_basis(tauresample,tau,nlinemod,nfile)
      implicit none
      integer :: file,nlinemod,nfile
      real,dimension(nlinemod,nfile) :: tau
      real,dimension(nlinemod) :: tauresample

      tauresample=0
c initialize the common tau(5000) with min depth = max of the min depth of the models
c                                  and max depth = min of the max depth of the models
c essential for  the resampling with cubic spline
      tauresample(1)=tau(1,1)
      tauresample(nlinemod)=tau(nlinemod,1)
      do file=2,nfile
         if (tauresample(1).lt.tau(1,file)) then
            tauresample(1)=tau(1,file)
          end if
         if (tauresample(nlinemod).gt.tau(nlinemod,file)) then
            tauresample(nlinemod)=tau(nlinemod,file)
         end if
      end do
      call blend_i_0d1 ( tauresample, nlinemod )
      end




!*******************************************************************************

      subroutine blend_i_0d1 ( x, m )
!
!
!! BLEND_I_0D1 extends indexed scalar data at endpoints along a line.
!http://orion.math.iastate.edu/burkardt/f_src/f_src.html
!
!  Diagram:
!
!    ( X1, ..., ..., ..., ..., ..., XM )
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
!    15 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X(M).
!
!    On input, X(1) and X(M) contain scalar values which are to be
!    interpolated through the entries X(2) through X(M).  It is assumed
!    that the dependence of the data is linear in the vector index I.
!
!    On output, X(2) through X(M-1) have been assigned interpolated
!    values.
!
!    Input, integer M, the number of entries in X.

      implicit none
!
      integer m
!
      integer i
      real r
      real x(m)
!
      do i = 2, m - 1

        r = real ( i - 1 ) / real ( m - 1 )

        call blend_101 ( r, x(1), x(m), x(i) )

      end do

      return
      end

      subroutine blend_101 ( r, x0, x1, x )
!
!*******************************************************************************
!
!! BLEND_101 extends scalar endpoint data to a line.
!http://orion.math.iastate.edu/burkardt/f_src/f_src.html
!
!  Diagram:
!
!    0-----r-----1
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
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the coordinate where an interpolated value is desired.
!
!    Input, real X0, X1, the data values at the ends of the line.
!
!    Output, real X, the interpolated data value at (R).
!
      implicit none
!
      real r
      real x
      real x0
      real x1
!
      x = ( 1.0E+00 - r ) * x0 + r * x1

      return
      end
!----------------------------------------------------------------------------

      REAL FUNCTION SevalSingle(u,x,y)
! ---------------------------------------------------------------------------
!http://www.pdas.com/fmm.htm
!  PURPOSE - Evaluate the cubic spline function
!     Seval=y(i)+b(i)!(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
!           where  x(i) <= u < x(i+1)

!  NOTES- if u<x(1), i=1 is used;if u>x(n), i=n is used

      REAL,INTENT(IN) :: u  ! abscissa at which the spline is to be evaluated
      REAL,INTENT(IN),DIMENSION(:) :: x ! abscissas of knots
      REAL,INTENT(IN),DIMENSION(:):: y ! ordinates of knots
      REAL, DIMENSION(:),allocatable :: b,c,d ! linear,quadratic,cubic coeff

      INTEGER, SAVE :: i=1
      INTEGER :: j, k, n
      REAL:: dx

      INTERFACE
         subroutine FMMsplineSingle(x, y, b, c, d)
        REAL,DIMENSION(:), INTENT(IN)  :: x ! abscissas of knots
        REAL,DIMENSION(:), INTENT(IN)  :: y ! ordinates of knots
        REAL,DIMENSION(:), INTENT(OUT) :: b ! linear coeff
        REAL,DIMENSION(:), INTENT(OUT) :: c ! quadratic coeff.
        REAL,DIMENSION(:), INTENT(OUT) :: d ! cubic coeff.
         end
      END INTERFACE
!----------------------------------------------------------------------------
      n=SIZE(x)
      allocate(b(n),c(n),d(n))
      call FMMsplineSingle(x, y, b, c, d)
!.....First check if u is in the same interval found on the
!        last call to Seval.............................................
       IF (  (i<1) .OR. (i >= n) ) i=1
       IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
         i=1   ! binary search
         j=n+1

         DO
           k=(i+j)/2
           IF (u < x(k)) THEN
             j=k
           ELSE
        i=k
           END IF
           IF (j <= i+1) EXIT
         END DO
       END IF

        dx=u-x(i)   ! evaluate the spline
        SevalSingle=y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))

        RETURN
        deallocate(b,c,d)
      END Function SevalSingle  ! -------------------------------------------------------


      SUBROUTINE FMMsplineSingle(x, y, b, c, d)
! ---------------------------------------------------------------------------
!http://www.pdas.com/fmm.htm
! PURPOSE - Compute the coefficients b,c,d for a cubic interpolating spline
!  so that the interpolated value is given by
!    s(x) = y(k) + b(k)*(x-x(k)) + c(k)*(x-x(k))**2 + d(k)*(x-x(k))**3
!      when x(k) <= x <= x(k+1)
!  The end conditions match the third derivatives of the interpolated curve to
!  the third derivatives of the unique polynomials thru the first four and
!  last four points.
!  Use Seval or Seval3 to evaluate the spline.
        REAL,DIMENSION(:), INTENT(IN)  :: x ! abscissas of knots
        REAL,DIMENSION(:), INTENT(IN)  :: y ! ordinates of knots
        REAL,DIMENSION(:), INTENT(OUT) :: b ! linear coeff
        REAL,DIMENSION(:), INTENT(OUT) :: c ! quadratic coeff.
        REAL,DIMENSION(:), INTENT(OUT) :: d ! cubic coeff.

        INTEGER:: k,n
        REAL:: t,aux
        REAL,PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
!----------------------------------------------------------------------------
       n=SIZE(x)

       IF (n < 3) THEN   ! Straight line - special case for n < 3
         b(1)=ZERO
         IF (n == 2) b(1)=(y(2)-y(1))/(x(2)-x(1))
         c(1)=ZERO
         d(1)=ZERO
         IF (n < 2) RETURN
         b(2)=b(1)
         c(2)=ZERO
         d(2)=ZERO
         RETURN
       END IF

!.....Set up tridiagonal system.........................................
!.    b=diagonal, d=offdiagonal, c=right-hand side
        d(1)=x(2)-x(1)
        c(2)=(y(2)-y(1))/d(1)
       DO k=2,n-1
         d(k)=x(k+1)-x(k)
         b(k)=TWO*(d(k-1)+d(k))
         c(k+1)=(y(k+1)-y(k))/d(k)
         c(k)=c(k+1)-c(k)
       END DO

!.....End conditions.  third derivatives at x(1) and x(n) obtained
!.       from divided differences.......................................
       b(1)=-d(1)
       b(n)=-d(n-1)
       c(1)=ZERO
       c(n)=ZERO
       IF (n > 3) THEN
         c(1)=c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
         c(n)=c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
         c(1)=c(1)*d(1)*d(1)/(x(4)-x(1))
         c(n)=-c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))
       END IF

       DO k=2,n    ! forward elimination
         t=d(k-1)/b(k-1)
         b(k)=b(k)-t*d(k-1)
         c(k)=c(k)-t*c(k-1)
       END DO

       c(n)=c(n)/b(n)   ! back substitution ( makes c the sigma of text)
       DO k=n-1,1,-1
         c(k)=(c(k)-d(k)*c(k+1))/b(k)
       END DO

!.....Compute polynomial coefficients...................................
       b(n)=(y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
       DO k=1,n-1
         b(k)=(y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
         d(k)=(c(k+1)-c(k))/d(k)
         c(k)=THREE*c(k)
       END DO
       c(n)=THREE*c(n)
       d(n)=d(n-1)

       RETURN
       END Subroutine FMMsplineSingle ! ---------------------------------------------------

      SUBROUTINE NaturalSplineSingle(x,y,b,c,d)
! ---------------------------------------------------------------------------
! PURPOSE - Construct the natural spline thru a set of points
! NOTES - A natural spline has zero second derivative at both endpoints.

       REAL,INTENT(IN),DIMENSION(:):: x,y   ! coordinates of knots
       REAL,INTENT(OUT),DIMENSION(:):: b,c,d  ! cubic coeff.

       INTEGER:: k,n
       REAL,PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
!-----------------------------------------------------------------------
       n=SIZE(x)

       IF (n < 3) THEN   ! Straight line - special case for n < 3
         b(1)=ZERO
         IF (n == 2) b(1)=(y(2)-y(1))/(x(2)-x(1))
         c(1)=ZERO
         d(1)=ZERO
         b(2)=b(1)
         c(2)=ZERO
         d(2)=ZERO
         RETURN
       END IF

       d(1:n-1) = x(2:n)-x(1:n-1)  ! Put the h-array of the text into array d

!.....Set up the upper triangular system in locations 2 thru n-1 of
!        arrays b and c. B holds the diagonal and c the right hand side.
       b(2)=TWO*(d(1)+d(2))
       c(2)=(y(3)-y(2))/d(2)-(y(2)-y(1))/d(1)
       DO  k=3,n-1
         b(k)=TWO*(d(k-1)+d(k))-d(k-1)*d(k-1)/b(k-1)
       c(k)=(y(k+1)-y(k))/d(k)-(y(k)-y(k-1))/d(k-1)-d(k-1)*c(k-1)/b(k-1)
       END DO

       c(n-1)=c(n-1)/b(n-1)   ! Back substitute to get c-array
       DO  k=n-2,2,-1
         c(k)=(c(k)-d(k)*c(k+1))/b(k)
       END DO
       c(1)=ZERO
       c(n)=ZERO   ! c now holds the sigma array of the text


!.....Compute polynomial coefficients ..................................
       b(n)=(y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
       DO  k=1,n-1
         b(k)=(y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
         d(k)=(c(k+1)-c(k))/d(k)
         c(k)=THREE*c(k)
       END DO
       c(n)=THREE*c(n)
       d(n)=d(n-1)
       RETURN

       END Subroutine NaturalSplineSingle



!----------------------------------------------------------------
      subroutine calcrhox(tau,kappa,ndepth,rhox)
c     2 ways to calculate rhox : int(ro*dx) or int(1/kappa*dtau)
c      A&A 387, 595-604 (2002)
      implicit none
      integer :: i,ndepth
      real :: first
      real, dimension(ndepth) :: tau,kappa,rhox,f,x
      real :: tot
      real, external :: rinteg


      f=(1/10**kappa)
      x=10**(tau)
      first = x(1)*f(1)
      tot=rinteg(x,f,rhox,ndepth,first)
      do i=2,ndepth
           rhox(i) = rhox(i-1) + rhox(i)
      enddo
      end

      real function rinteg(x,f,fint,n,start)
c******************************************************************************
c     This routine is from ATLAS6
c******************************************************************************
      implicit none
      integer :: n,i,n1
      real x(5000), f(5000), fint(5000)
      real a(5000), b(5000), c(5000)
      real :: start

      call parcoe (f,x,a,b,c,n)
      fint(1) = start
      rinteg = start
      n1 = n - 1
      do 10 i=1,n1
         fint(i+1)= (a(i)+b(i)/2.*(x(i+1)+x(i))+
     .     c(i)/3.*((x(i+1)+x(i))*x(i+1)+x(i)*x(i)))*(x(i+1)-x(i))
10    rinteg = rinteg + fint(i+1)

      return
      end




      subroutine parcoe(f,x,a,b,c,n)

      implicit none
      integer :: n,n1,j,j1
      real f(5000), x(5000), a(5000), b(5000), c(5000)
      real :: d,wt

      c(1)=0.
      b(1)=(f(2)-f(1))/(x(2)-x(1))
      a(1)=f(1)-x(1)*b(1)
      n1=n-1
      c(n)=0.
      b(n)=(f(n)-f(n1))/(x(n)-x(n1))
      a(n)=f(n)-x(n)*b(n)
      if(n.eq.2)return
      do 1 j=2,n1
      j1=j-1
      d=(f(j)-f(j1))/(x(j)-x(j1))
      c(j)=f(j+1)/((x(j+1)-x(j))*(x(j+1)-x(j1)))-f(j)/((x(j)-x(j1))*
     1(x(j+1)-x(j)))+f(j1)/((x(j)-x(j1))*(x(j+1)-x(j1)))
      b(j)=d-(x(j)+x(j1))*c(j)
    1 a(j)=f(j1)-x(j1)*d+x(j)*x(j1)*c(j)
      c(2)=0.
      b(2)=(f(3)-f(2))/(x(3)-x(2))
      a(2)=f(2)-x(2)*b(2)
      c(3)=0.
      b(3)=(f(4)-f(3))/(x(4)-x(3))
      a(3)=f(3)-x(3)*b(3)
      if(c(j).eq.0.)go to 2
      j1=j+1
      wt=abs(c(j1))/(abs(c(j1))+abs(c(j)))
      a(j)=a(j1)+wt*(a(j)-a(j1))
      b(j)=b(j1)+wt*(b(j)-b(j1))
      c(j)=c(j1)+wt*(c(j)-c(j1))
    2 continue
      a(n1)=a(n)
      b(n1)=b(n)
      c(n1)=c(n)
      return
      end

***********************************************************************************
      subroutine calc_error(xinf,xsup,yinf,ysup,zinf,zsup,teff_ref,
     &   logg_ref,z_ref)
      implicit none
      real :: xinf,xsup,yinf,ysup,zinf,zsup,teff_ref,logg_ref,z_ref
      real ::  error_T,error_Pe,error_Pg,error_kappa
      real :: errorTeffT,errorloggT,errorzT,
     &        errorTeffPe,errorloggPe,errorzPe,
     &         errorTeffPg,errorloggPg,errorzPg,
     &         errorTeffkappa,errorloggkappa,errorzkappa

! values read out of the figures of the manual and scaled down o the according step
              errorTeffT=0.055/32
              errorloggT=0.008/5
              errorzT=0.015/4
              errorTeffPe=0.65/32
              errorloggPe=0.4/5
              errorzPe=0.38/4
              errorTeffPg=0.25/32
              errorloggPg=0.23/5
              errorzPg=0.35/4
              errorTeffkappa=0.8/32
              errorloggkappa=0.36/5
              errorzkappa=0.38/4


      error_T=min(abs(xsup-teff_ref),abs(teff_ref-xinf))/100*errorTeffT+
     &         min(abs(ysup-logg_ref),abs(logg_ref-yinf))*errorloggT +
     &         min(abs(zsup-z_ref),abs(z_ref-zinf))*errorzT
      write(*,1409) 'estimated max error on T =',error_T*100,'%'


      error_Pe=min(abs(xsup-teff_ref),abs(teff_ref-xinf))/100
     &                                                   *errorTeffPe +
     &       + min(abs(ysup-logg_ref),abs(logg_ref-yinf))*errorloggPe +
     &         min(abs(zsup-z_ref),abs(z_ref-zinf))*errorzPe
       write(*,1409) 'estimated max error on Pe =',error_Pe*100,'%'


      error_Pg=min(abs(xsup-teff_ref),abs(teff_ref-xinf))/100
     &                                                   *errorTeffPg +
     &         min(abs(ysup-logg_ref),abs(logg_ref-yinf))*errorloggPg +
     &         min(abs(zsup-z_ref),abs(z_ref-zinf))*errorzPg
       write(*,1409) 'estimated max error on Pg =',error_Pg*100,'%'


      error_kappa=min(abs(xsup-teff_ref),abs(teff_ref-xinf))/100
     &                                                *errorTeffkappa +
     &       min(abs(ysup-logg_ref),abs(logg_ref-yinf))*errorloggkappa +
     &       min(abs(zsup-z_ref),abs(z_ref-zinf))*errorzkappa
      write(*,1409) 'estimated max error on kappa =',error_kappa*100,'%'

 1409 format(a30,f5.1,a2)
      end

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine  write_turbo(out,n,xlr,logg,sph,
     &                     taus,T,Pe,Pg,xit,rr,rhox,taur,ro,FILE)
       implicit none
       logical :: sph
       integer :: k,out,n,fic
       real :: xlr,logg
       real, dimension(n) :: taus,T,Pe,Pg,xit,
     &         rr,taur,rhox,ro
       character*256, dimension(10) :: FILE


       if (sph) then
        write(*,*) 'spherical model'
        write(out,1967) n,xlr,logg
 1967   format('''sphINTERPOL''',1x,i3,f8.0,2x,f5.2,1x,'0 0.00')
         do k=1,n
           write(out,1968) taus(k),T(k),Pe(k),
     &                    Pg(k),xit(k),rr(k),taur(k),rhox(k),ro(k)
         enddo
 1968   format(f8.4,1x,f8.2,3(1x,f8.4),1x,e15.6,1x,f8.4,x,es14.6,x
     &      ,es12.4)

        else
        write(*,*) 'plane parallel model'
        write(out,1966) n,xlr,logg
1966     format('''ppINTERPOL''',1x,i3,f8.0,2x,f5.2,1x,'0 0.00')
         do k=1,n
           write(out,1965) taus(k),T(k),Pe(k),
     &                    Pg(k),xit(k),rr(k),taur(k),rhox(k),ro(k)
1965       format(f8.4,1x,f8.2,3(1x,f8.4),1x,e15.6,1x,f8.4,x,es14.6,x
     &      ,es12.4)
          enddo
       end if
       do fic=1,8
       write(out,*) trim(FILE(fic))
       enddo
       end


      subroutine write_MARCSweb(out,FILE,rmass,xit,sph,extrapol,
     &  temp_ref,logg_ref,z_ref,rr,ndepth_ref,
     &    tauR,taus,T,Pe,Pg,presrad,pturb,xkapref,rhox,rho,abundances)

c Bengt Edvardsson 01/2013
      implicit none
      logical :: sph,extrapol,lwarning
      integer :: ndepth_ref,k,i,out,today(3)
      real :: temp_ref,logg_ref,z_ref,xite,cc,rmass,radius,
     &          Xmass,Ymass,Zmass,alpha
      real :: abundances(92),dluminosity
      real, dimension(ndepth_ref) :: tauR,taus,depth,T,Pe,Pg,presrad,
     &  pturb,rho,Emu,Vconv,Anconv,rhox,xit,xkapref,rr
      real, dimension(32,ndepth_ref) :: partial
      character*256 :: FILE(10),filename,f1,f2,f3,f4,f5,f6,f7,f8
      character*8 :: today_char

      write(*,*) 'MARCS web format output REMAIN TO BE finished'
      stop
      call makename(filename,temp_ref,logg_ref,z_ref)
* filename(out) is the constructed output file name
        write(out,'(a)') filename
***********************
       xite=real(xit(1))
       call idate(today)
       write(today_char,'(a)') today(3),today(2),today(1)
        if(extrapol) then
          write(out,2310) temp_ref,today_char,5.6704e-5*temp_ref**4,
     &                   10.**logg_ref,xite
 2310     format(  f7.0,'      Teff [K]. ',
     &                  '  EXTRAPOLATED   model; yyyymmdd=',a8,
     &                  ''/,
     &           1p,e12.4,0p,' Flux [erg/cm2/s]',/,
     &           1p,e12.4,0p,' Surface gravity [cm/s2]',/,
     &           f5.1,'        Microturbulence parameter [km/s]')
        else
          write(out,2311) temp_ref,today_char,5.6704e-5*temp_ref**4,
     &                   10.**logg_ref,xite
 2311     format(  f7.0,'      Teff [K]. ',
     &                  '  INTERPOLATED   model; yyyymmdd=',a8,
     &                  ''/,
     &           1p,e12.4,0p,' Flux [erg/cm2/s]',/,
     &           1p,e12.4,0p,' Surface gravity [cm/s2]',/,
     &           f5.1,'        Microturbulence parameter [km/s]')
        endif
        if(sph) then
          write(out,2315) rmass
 2315     format(f5.1,'        Mass [Msun]')
        else
          write(out,2317)
 2317     format('  0.0        No mass for plane-parallel models')
        endif
        alpha=abundances(22)-4.90+z_ref
        write(out,2320)  z_ref,alpha
 2320   format(f6.2,f6.2,  ' Metallicity [Fe/H] and [alpha/Fe]')
        if(sph) then
          radius=6.9598e10*sqrt(rmass/10.**(logg_ref-4.44))
        else
          radius=radius
        endif
        dluminosity=rmass/10.**(logg_ref-4.44)*(temp_ref/5777.)**4
        if(sph) then
          write(out,2322)  radius,dluminosity
 2322     format(1p,e12.4,0p,' Radius [cm] at Tau(Rosseland)=1.0',/,
     &           f12.5,      ' Luminosity [Lsun]')
        else
          write(out,2324)  radius,dluminosity
 2324     format(1p,e12.4,0p,' 1 cm radius for plane-parallel models',/,
     &           1p,e12.4,0p,' Luminosity [Lsun] FOR A RADIUS OF 1 cm!')
        endif
        write(out,'(a)') '1.50 8.00 0.076 0.00 are the convection
     & parameters: alpha, nu, y and beta'
        call xyz(92,abundances,Xmass,Ymass,Zmass)
        cc=abundances(6)-8.39+z_ref
        if(cc.eq.-0.13) then
          write(out,2328) Xmass,Ymass,Zmass
 2328     format(f9.5,f8.5,1p,e9.2,0p,' are X, Y and Z, 12C/13C=20')
          print *,'Mildly CN-cycled model'
        else if(cc.eq.-0.38) then
          write(out,2330) Xmass,Ymass,Zmass
 2330     format(f9.5,f8.5,1p,e9.2,0p,' are X, Y and Z, 12C/13C=4')
          print *,'Heavily CN-cycled model'
        else if(cc.eq.+0.00 .or. cc.eq.+0.08) then
          write(out,2332) Xmass,Ymass,Zmass
 2332   format(f9.5,f8.5,1p,e9.2,0p,
     &         ' are X, Y and Z, 12C/13C=89 (=solar)')
        else
          stop 'carbon abundance not expected'
        endif
        write(out,2336) abundances
 2336   format('Logarithmic chemical number abundances, H always 12.00',
     &          10(/,10f7.2))
        if(dluminosity.ge.1.e6) then
          lwarning=.true.
        endif
* structure
        write(out,1230) ndepth_ref
 1230   format(i4,' Number of depth points')
        if(extrapol) then
          write(out,1232)
 1232     format('Model structure,           EXTRAPOLATED')
        else
          write(out,1233)
 1233     format('Model structure,           INTERPOLATED')
        endif
        write(out,1234)
 1234   format(' k lgTauR  lgTau5    Depth     T  ',
     &         '      Pe          Pg         Prad       Pturb')
        do k=1,ndepth_ref
        depth(k)=rr(k)-radius
          write(out,1236) k,tauR(k),taus(k),depth(k),
     &                   T(k),10.**pe(k),10.**pg(k),
     &                   10.**presrad(k),pturb(k)
 1236     format(i3,f6.2,f8.4,1p,e11.3,0p,f8.1,1p,5e12.4,0p)
* Format changed 2012-04-03
        enddo
* thermodynamics
        write(out,2338)
 2338   format(' k lgTauR    KappaRoss   Density   Mu      Vconv',
     &         '   Fconv/F      RHOX')
        do k=1,ndepth_ref
          if(Vconv(k).lt.10.0) Vconv(k)=0.0
          if(Anconv(k).lt.0.0) Anconv(k)=0.0
          write(out,210) k,tauR(k),10.**xkapref(k),
     &                  10.**rho(k),Emu(k),Vconv(k),
     &                  Anconv(k),rhox(k)
  210     format(i3,f6.2,1p,2e12.4,0p,f6.3,1p,e11.3,0p,f8.5,1p,e14.6,0p)
* Format changed 2012-04-03
        enddo
* molecular partial pressures
        if(extrapol) then
          write(out,1240)
 1240     format('Assorted logarithmic partial pressures,       ',
     &           '         EXTRAPOLATED')
        else
          write(out,1241)
 1241     format('Assorted logarithmic partial pressures,       ',
     &           '         INTERPOLATED')
        endif
        write(out,1242)
 1242   format(' k  lgPgas   H I    H-     H2     H2+    H2O',
     &         '    OH     CH     CO     CN     C2  ')
        do k=1,ndepth_ref
          write(out,1244) k,pg(k),(partial(i,k),i=1,10)
 1244     format(i3,f7.3,11f7.2)
        enddo
        write(out,1246)
 1246   format(' k    N2     O2     NO     NH     TiO',
     &         '   C2H2    HCN    C2H    HS     SiH    C3H')
        do k=1,ndepth_ref
          write(out,1248) k,(partial(i,k),i=11,21)
 1248     format(i3,12f7.2)
        enddo
        write(out,1250)
 1250   format(' k    C3     CS     SiC   SiC2    NS  ',
     &         '   SiN    SiO    SO     S2     SiS   Other')
        do k=1,ndepth_ref
          write(out,1248) k,(partial(i,k),i=22,32)
        enddo
* list the input files for later reference and debugging
        f1=FILE(1)
        f2=FILE(2)
        f3=FILE(3)
        f4=FILE(4)
        f5=FILE(5)
        f6=FILE(6)
        f7=FILE(7)
        f8=FILE(8)
        write(out,*) trim(f1)//'   '//trim(f2)//'   '//trim(f3)//'   '//
     &              trim(f4)//'   '//trim(f5)//'   '//trim(f6)//'   '//
     &              trim(f7)//'   '//trim(f8)

      end

*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
      subroutine makename(filename,temp_ref,logg_ref,z_ref)
* construct the output file name
      implicit none
      integer      :: nfile
      parameter (nfile=11)
      character*76 :: filename(nfile),file
      character    :: ctemp*4,clogg*4,cz*5,class*2
      real         :: temp_ref,logg_ref,z_ref,alpha,c
*
      write(ctemp,'(i4)') nint(temp_ref)
      if(logg_ref.ge.0.0) then
        write(clogg,'("+",f3.1)') logg_ref
      else
        write(clogg,'(f4.1)') logg_ref
      endif
      if(z_ref.ge.0.0) then
        write(cz,'("+",f4.2)') z_ref
      else
        write(cz,'(f5.2)') z_ref
      endif
      file=filename(1)
******** composition class:
      read(file,'(32x,f5.2,2x,f5.2)') alpha,c
      if(c.eq.-0.13) then
        class='mc'
      else if(c.eq.-0.38) then
        class='hc'
      else if(alpha.eq.-0.40) then
        class='an'
      else if(z_ref.gt.-1.0 .and. alpha.eq.0.4) then
        class='ae'
      else if(z_ref.lt.0.0 .and. alpha.eq.0.0) then
        class='ap'
      else
        class='st'
      endif
********
*4000_g+1.0_m1.0_t02_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod
      filename=file(1:1)//ctemp//'_g'//clogg//file(12:21)//class//
     &             '_z'//cz//file(31:72)//'.int'
      return
      end
*
*
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
      subroutine xyz(n,abund4,xx,yy,zz)
* compute number fraction and mass fraction of elements
* program solarabund.f does this too
      implicit none
      integer i,n
      character*2 lel(92)
      real abund4(92),xx,yy,zz,norm
      real*8 abund(92),weight(92),sumass,weighthamu
      real*8 w,rf,fract,fmass,weightamu,abumetals,wmetals,x,y,z,sumab
*
      data weight /
     &1.0079,4.0026, 6.941, 9.012, 10.81, 12.01, 14.01, 15.99, 19.00,
     & 20.18, 22.99, 24.30, 26.98, 28.09, 30.97, 32.06, 35.45, 39.95,
     & 39.10, 40.08, 44.96, 47.90, 50.94, 52.00, 54.94, 55.85, 58.93,
     & 58.71, 63.55, 65.37, 69.72, 72.59, 74.92, 78.96, 79.90, 83.80,
     & 85.47, 87.62, 88.91, 91.22, 92.91, 95.94, 98.91,101.07,102.91,
     &106.4 ,107.87,112.40,114.82,118.69,121.75,127.60,126.90,131.30,
     &132.91,137.34,138.91,140.12,140.91,144.24,146.  ,150.4 ,151.96,
     &157.25,158.93,162.50,164.93,167.26,168.93,170.04,174.97,178.49,
     &180.95,183.85,186.2 ,190.2 ,192.2 ,195.09,196.97,200.59,204.37,
     &207.19,208.98,210.  ,210.  ,222.  ,223.  ,226.03,227.  ,232.04,
     &230.04,238.03 /
      data lel /
     & 'H ' , 'He' , 'Li' , 'Be' , 'B ' , 'C ' , 'N ' , 'O ' , 'F ' ,
     & 'Ne' , 'Na' , 'Mg' , 'Al' , 'Si' , 'P ' , 'S ' , 'Cl' , 'Ar' ,
     & 'K ' , 'Ca' , 'Sc' , 'Ti' , 'V ' , 'Cr' , 'Mn' , 'Fe' , 'Co' ,
     & 'Ni' , 'Cu' , 'Zn' , 'Ga' , 'Ge' , 'As' , 'Se' , 'Br' , 'Kr' ,
     & 'Rb' , 'Sr' , 'Y ' , 'Zr' , 'Nb' , 'Mo' , 'Tc' , 'Ru' , 'Rh' ,
     & 'Pd' , 'Ag' , 'Cd' , 'In' , 'Sn' , 'Sb' , 'Te' , 'I ' , 'Xe' ,
     & 'Cs' , 'Ba' , 'La' , 'Ce' , 'Pr' , 'Nd' , 'Pm' , 'Sm' , 'Eu' ,
     & 'Gd' , 'Tb' , 'Dy' , 'Ho' , 'Er' , 'Tm' , 'Yb' , 'Lu' , 'Hf' ,
     & 'Ta' , 'W ' , 'Re' , 'Os' , 'Ir' , 'Pt' , 'Au' , 'Hg' , 'Tl' ,
     & 'Pb' , 'Bi' , 'Po' , 'At' , 'Rn' , 'Fr' , 'Ra' , 'Ac' , 'Th' ,
     & 'Pa' , 'U ' /
* sum all masses * abundances
      if(n.ne.92) stop 'subr XYZ not 92 abundance numbers'
      w=0.0d0
      rf=0.0d0
      if(abs(abund4(1)-12.00e0).gt.1.0e-7) then
        norm=12.0d0-abund4(1)
        do i=1,92
          abund4(i)=abund4(i)+norm
        enddo
      endif
      do i=92,3,-1
        abund(i)=dble(abund4(i))
        if(abund(i).gt.-9.) then
          rf=rf + 10.d0**abund(i)
          w=w + weight(i)*10.d0**abund(i)
        endif
      enddo
      abumetals=dlog10(rf)
*     print *,'abumetals=',abumetals
      wmetals=w
      abund(2)=dble(abund4(2))
      rf=rf + 10.d0**abund(2)
      w=w + weight(2)*10.d0**abund(2)
      abund(1)=dble(abund4(1))
      rf=rf + 10.d0**abund(1)
      w=w + weight(1)*10.d0**abund(1)
*     print *,'Z=',wmetals/w
*     print *,'Y=',weight(2)*10.d0**abund(2)/w
*     print *,'X=',weight(1)*10.d0**abund(1)/w
*
      weightamu=w/rf
* weightamu is the mean weight of an atom (in amu)
*     print 1000
      sumab=0.0
      sumass=0.0
      do 20 i=1,30
        if(abund(i).gt.-3.) then
          fract=10.**abund(i)/rf
          fmass=10.**abund(i)*weight(i)/w
          if(i.eq.1) then
            X=fmass
            weightHamu=weightamu/fract
          else if(i.eq.2) then
            Y=fmass
          endif
        else
          fract=0.0
          fmass=0.0
        endif
        if(i.gt.2) then
          sumab=sumab+fract
          sumass=sumass+fmass
*         print 1010,i,lel(i),weight(i),abund(i),fract,fmass
*1010     format(1x,i2,3x,a2,f10.4,10x,f5.2,8x,1pe10.2,6x,e10.2)
        else
*         print 1020,i,lel(i),weight(i),abund(i),fract,fmass
*1020     format(1x,i2,3x,a2,f10.4,10x,f5.2,8x,f10.6,6x,f10.6)
        endif
*1000   format('Element Weight(amu)  log abundance   Number fraction',
*    &         ' Mass fraction')
   20 continue
      Z=sumass
*     print *,'======================================================'
*     write(*,1030) abumetals
*1030 format(' log. abundance of metal atoms    =   ',f8.5)
*     write(*,1040) sumab
*1040 format(' metal (Z>2) atom number fraction =   ',f8.5)
*     write(*,1050) sumass
*1050 format(' metal (Z>2) atom mass fraction   =   ',f8.5)
*     print *,
*    &'mean mass of a metal atom        =',sumass/sumab*weightamu,' amu'
*     print *,'mean mass per atom               =',weightamu,' amu'
*     print *,'mean mass per hydrogen atom      =',weightHamu,' amu'
*     write(*,1060) x,y,z,x+y+z
*1060 format(' X=',f8.5,'   Y=',f8.5,'   Z=',f8.5,'   X+Y+Z=',f8.5)
      xx=real(x)
      yy=real(y)
      zz=real(z)
      return
      end
*
*
ccccccc
      subroutine write_atlas(out,n,temp_ref,logg_ref,
     &           z_ref,Pe,Prad,abu,rhox,T,Pg,xit,xkapref)

c adapted from Mucciarelli 12/2012
      implicit none
      integer :: out,k,J,n,codex(6)
      real :: yfraction,xfraction,mass_He,mass_H,mass_Z,
     &       temp_ref,logg_ref,z_ref
      real, dimension(n) :: Prad,rhox,T,Pg,Pe,xkapref,accrad,xit
      real :: abu(92),abu2(6)

      call xyz(92,abu,mass_H,mass_He,mass_Z)
       yfraction=mass_He/(4-3*mass_He)
       xfraction=1-yfraction
       write(out,'(A4,2x,F6.0,2x,A7,1x,F7.5,1x,A3)')
     &'TEFF',temp_ref,'GRAVITY',logg_ref,'LTE'
       write(out,'(A18)')'TITLE interpoltoATLAS'
       write(out,'(A53)')
     &' OPACITY IFOP 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0'
       write(out,'(A59)')
     & 'CONVECTION ON   1.50 TURBULENCE OFF  0.00  0.00  0.00  0.00'

       write(out,'(A15,F10.5,A20,F7.5,A3,F7.5)')
     &'ABUNDANCE SCALE',10**z_ref,' ABUNDANCE CHANGE 1 ',xfraction,
     &' 2 ',yfraction

       do J=2,n-1
        ACCRAD(J)=0.5*( (PRAD(J)-PRAD(J-1) ) / (RHOX(J)-
     &  RHOX(J-1))+
     & (PRAD(J+1)-PRAD(J) )/(RHOX(J+1)-RHOX(J) ) )

        ACCRAD(1)=(PRAD(J)-PRAD(J-1))/(RHOX(J)-
     &  RHOX(J-1))
        ACCRAD(n)=(PRAD(J)-PRAD(J-1))/(RHOX(J)-
     &  RHOX(J-1))
       enddo


       do k=3,93,6
         do j=1,6
	   abu2(j)=abu(k-1+j)-z_ref-12.04
	   codex(j)=(k-1)+j
	   if (abu2(j).lt.-20) abu2(j)=-20.00
	 enddo
	 write(out,'(A17,6(I3,F7.2))')
     &     'ABUNDANCE CHANGE',codex(1),abu2(1),
     &     codex(2),abu2(2),codex(3),abu2(3),
     &     codex(4),abu2(4),codex(5),abu2(5),
     &     codex(6),abu2(6)
       enddo
       write(out,'(A27)')'ABUNDANCE CHANGE  99 -20.00'

       write(out,'(A10,I3,A33)')'READ DECK6 ',n,
     &'RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB'

       do k=1,n
        write(out,'(e15.8,2x,f7.1,5(1x,e9.3))')
     &       rhox(k),t(k),10**Pg(k),
     &      (10**Pe(k))/(1.38054e-16*t(k)),10**(xkapref(k)),
     &       ACCRAD(k),xit(k)*1e5
       enddo

       write(out,'(A7,E10.4)')'PRADK0 ',PRAD(1)
       write(out,'(A48)')
     &'BEGIN                    ITERATION  01 COMPLETED'

       end


ccccc
       subroutine write_MOOG(out,sph,teff,logg,z,n,
     &            lambda_ref,taus,T,Pe,Pg,rhox,xit)
Cadapted from J. Sobeck 01/2013

       implicit none
       logical :: sph
       integer :: k,n,out
       real :: teff,logg,z,lambda_ref
       real,dimension(n) :: taus,T,Pe,Pg,rhox,xit

       if (sph) then
*** write interpolated, spherical model to file with MOOG-format
        write(out,19676)
19676   format('WEB2MARC')
        write(out,19672) teff,logg,z
        write(out,19673) n
19672   format('INTERPOL-sph, Teff = ',f6.1,' K ','log g = ',f4.2,
     &         ' dex ','[M/H] = ', f5.2)
19673   format('NTAU        ',i3)
        write(out,19674) lambda_ref
19674   format('       ',f7.1)
          do k=1,n
           write(out,19675) k,taus(k),T(k),Pe(k),Pg(k),
     &     rhox(k)
19675      format(i5,1x,f8.4,1x,f8.2,2(1x,f8.4),1x,e15.6)
          enddo
        write(out,'(a4,f5.3,a4)') '    ', xit(1), 'e+00'
        write(out,'(a14,f5.2)') 'NATOMS     0  ', z
        write(out,'(a12)') 'NMOL      19'
        else
*** MOOG-Format Plane Parallel
        write(out,19686)
19686   format('WEB2MARC')
        write(out,19682) teff,logg,z
        write(out,19685) n
19685   format('NTAU        ',i3)
19682   format('INTERPOL-ppl, Teff = ',f6.1,' K ','log g = ',f4.2,
     &         ' dex',' [M/H] = ', f5.2)
        write(out,19683) lambda_ref
19683   format('       ',f7.1)
          do k=1,n
           write(out,19684) k,taus(k),T(k),Pe(k),Pg(k),
     &     rhox(k)
19684      format(i5,1x,f8.4,1x,f8.2,2(1x,f8.4),1x,e15.6)
          enddo
        write(out,'(a4,f5.3,a4)') '    ', xit(1), 'e+00'
        write(out,'(a14,f5.2)') 'NATOMS     0  ', z
        write(out,'(a12)') 'NMOL      19'
       end if
       write(out,31455)
       write(out,31456)
       write(out,31457)
31455  format('      606.0    106.0    607.0    608.0    107.0    108.0
     &    112.0    707.0')
31456  format('      708.0    808.0     12.1  60808.0  10108.0    101.0
     &      6.1      7.1')
31457  format('        8.1    822.0     22.1')
       end

cc
       subroutine write_default(out,teff,logg,z,n,
     &        lambda_ref,tauR,taus,T,Pe,
     &        Pg,Prad,xit,rhox,rr,ro)

       implicit none
       integer :: out,k,n
       real :: teff,logg,z,lambda_ref
       real, dimension(n) ::tauR,taus,T,Pe,
     &        Pg,Prad,xit,rhox,rr,ro
       character*100 :: string

       write(string,2387) teff,logg,z,lambda_ref,n
 2387  format('Teff= ',f5.0,' logg= ',f6.2,' [m/H]= ',f6.2,
     &    ' lambda= ',f8.1,' nber layer= ',i4)
       write(out,'(a)') string
        write(out,19672)
19672   format('    k log(tau_ross) log(tau)  T    log(Pe)   log(Pg)
     &   log(Prad)  microt rhox  radius ro')
         do k=1,n
           write(out,19681) k,tauR(k),taus(k),T(k),Pe(k),
     &        Pg(k),Prad(k),xit(k),rhox(k),rr(k),ro(k)
19681      format(i5,1x,f8.4,1x,f8.4,8(1x,f10.4))
          enddo
        end
