program tnf_xy_onelevel

! Wave-activity flux based on Eq (38) of Takaya and Nakamura (2001)
!  for horizontal compoents

! Takaya, K., Nakamura, H., 2001: 
!   A formulation of a phase-independent wave-activity flux 
!   for stationary and migratory quasi-geostrophic eddies 
!   on zonally-varying basic flow. J. Atmos. Sci., 58, 608--627. 
 
! Input file format
!  4-byte, real
!  unformatted direct access 
!  lon; 144 (0 - 357.5E)
!  lat;  73 (90S - 90N)
!  1-level 

!  Z(m),U(m/s),V(m/s) for 
!   fields to be subtracted by basic state
!   and basic state.

! Output data
!   x,y component 
!   (144x73)

!   p = p/1000.
!   f=2*7.29E-5 *sin(psi)

! Examples of comilation 
! ifort -o  tnf_xy_onelevel.run tnf_xy_onelevel.f90 -assume byterecl
! ifc -o tnf_xy_onelevel.run tnf_xy_onelevel.f90 -lIEPCF90
! g95 -o tnf_xy_onelevel.run tnf_xy_onelevel.f90 
! gfortran -o tnf_xy_onelevel.run tnf_xy_onelevel.f90
! (latter used for abyss.earthsci.unimelb.edu.au)

! This program works like this; 
!
! ./tnf_xy_onelevel.run << EOF
! 250
! z250.dat
! z250.clim.dat
! u250.dat
! u250.clim.dat
! v250.dat
! v250.clim.dat
! tnfx250.dat
! tnfy250.dat
! EOF

  implicit none

  real,parameter :: unavl=9.999e+20 !data not available
  integer,parameter :: mlon=144,mlat=73
  integer,parameter :: mmax=mlon*mlat,lmax=mlon*mlat
! grid interval in degree
  real,parameter :: dlati=2.5,dloni=2.5

! The limit latitudes of QG approximation
!  ( This may be changed)
! 43-> 15N, 31-> 15S
  integer,parameter :: neqlats=31,neqlatn=43

! Unit number for input file
  integer,parameter :: un_z=10,un_zb=11
  integer,parameter :: un_u=14,un_ub=15
  integer,parameter :: un_v=16,un_vb=17
! Unit number for output file
  integer,parameter :: un_tnfx=30,un_tnfy=31


  real :: umdd(mlon,mlat),ubmdd(mlon,mlat)
  real :: vmdd(mlon,mlat),vbmdd(mlon,mlat)
  real :: zmdd(mlon,mlat),zbmdd(mlon,mlat)
  real :: stabmdd(mlon,mlat)
  real :: tnfxdd(mlon,mlat),tnfydd(mlon,mlat)

  integer irec,irecb,ilev,ilon,ilat

  real :: prem
  real,parameter :: gra = 9.8
  real,parameter :: ubmin=1.,tday=60.*60.*24.,erd=6.371E+06
  real :: pai,dgrd,rloni,rloni2, rlati,rlati2,oum2

! constant
  pai   = acos(-1.)
  dgrd  = pai/180.
  rloni = dloni*dgrd
  rloni2= 2.*rloni
  rlati = dlati*dgrd
  rlati2= 2.*rlati
  oum2 = 2.*2.*pai/tday
  
! Note that block size is multiplied by 4. 
! Some compilers do not need this.

  write(6,*) "input level (hPa)"
  read(*,*) prem

  prem = prem/1000.

  write(6,*) "input file name of z"
  call of_udo (un_z,mmax*4)
  
  write(6,*) "input file name of z of basic state"
  call of_udo (un_zb,mmax*4)

  write(6,*) "input file name of u"
  call of_udo (un_u,mmax*4)

  write(6,*) "input file name of u of basic state"
  call of_udo (un_ub,mmax*4)

  write(6,*) "input file name of v "
  call of_udo (un_v,mmax*4)

  write(6,*) "input file name of v of basic state"
  call of_udo (un_vb,mmax*4)

! for output file
  write(6,*) "The output file name of TN flux x"
  call of_udr (un_tnfx,lmax*4)
  write(6,*) "The output file name of TN flux y"
  call of_udr (un_tnfy,lmax*4)

  irec=0
  irecb=0

! for irec
  do 

! If records of anomaly and basic state fields are different, 
!  modify the increment.

     irec=irec+1

! This is an example that the basic state is climatology for each month
     irecb=irecb+1
     if(irecb .eq. 13)  irecb = 1

! test
!     if(irec .eq. 2) exit

     read(un_z,rec=irec,err=300)  zmdd
     read(un_zb,rec=irecb,err=300) zbmdd

     read(un_u,rec=irec,err=300)  umdd
     read(un_ub,rec=irecb,err=300) ubmdd
     read(un_v,rec=irec,err=300 ) vmdd
     read(un_vb,rec=irecb,err=300) vbmdd

     call tnflux_xy_comp
         
     write(un_tnfx,rec=irec)tnfxdd
     write(un_tnfy,rec=irec)tnfydd

  enddo ! irec

300 write(*,*) "file ended"
  stop

contains

    subroutine tnflux_xy_comp
      implicit none
      integer :: ilatn,ilats
      real :: termxU,termxV,termyU,termyV
      real :: wspeed
      real :: stab,ub,vb,psiam,zam
      real :: uam,vam,dudx,dvdx,dudy,dvdy

      real :: rlat,coriol,coslat
      integer :: ilon,ilone,ilonw
      real :: uamn,uams,uame,uamw,vame,vamw
      real :: fx,fy
! You can input zero into flux where the flux is undefined.      
      do ilat=1,mlat
         do ilon=1,mlon
            tnfxdd(ilon,ilat)=unavl
            tnfydd(ilon,ilat)=unavl
         enddo
      enddo

      do ilat=2,mlat-1

! If the latitude is north to south, change the following 4 lines 
!  and definition of neqlats and neqlatn

         if(ilat .gt. neqlats .and. ilat .lt.neqlatn) cycle

         ilatn=ilat+1
         ilats=ilat-1

         rlat=dgrd*(90.-real((ilat-1))*dlati)

         coriol=oum2*sin(rlat)
         coslat=cos(rlat)

         do ilon=1,mlon

! If westerly wind of the basic state is smaller than ubmin or easterly, 
! flux is unavailable.

            ub=ubmdd(ilon,ilat)
            if(ub.lt.ubmin) cycle

            vb=vbmdd(ilon,ilat)
            wspeed = sqrt(ub*ub+vb*vb)

            ilone=ilon+1
            if(ilone.gt.mlon) ilone=1

            ilonw=ilon-1
            if(ilonw.lt.1) ilonw=mlon

! anomaly
            zam=zmdd(ilon,ilat)-zbmdd(ilon,ilat)
            
! geopotential height -> QG streamfunction

            psiam=gra/coriol*zam

            uam=umdd(ilon,ilat)-ub
            vam=vmdd(ilon,ilat)-vb

            uamn=umdd(ilon,ilatn)-ubmdd(ilon,ilatn)
            uams=umdd(ilon,ilats)-ubmdd(ilon,ilats)

            uame=umdd(ilone,ilat)-ubmdd(ilone,ilat)
            uamw=umdd(ilonw,ilat)-ubmdd(ilonw,ilat)

            vame=vmdd(ilone,ilat)-vbmdd(ilone,ilat)
            vamw=vmdd(ilonw,ilat)-vbmdd(ilonw,ilat)

            dudx = (uame-uamw)/(erd*rloni2*coslat)
            dvdx = (vame-vamw)/(erd*rloni2*coslat)

            dudy = (uamn-uams)/(erd*rlati2)


            termxU=vam*vam-psiam*dvdx
            termxV=psiam*dudx-uam*vam
            termyU=termxV
            termyV=uam*uam+psiam*dudy

! tnf
            fx=ub*termxU+vb*termxV
            fy=ub*termyU+vb*termyV

            tnfxdd(ilon,ilat)=prem*coslat/(2.*wspeed)*fx
            tnfydd(ilon,ilat)=prem*coslat/(2.*wspeed)*fy
         enddo
      enddo

      return
    end subroutine tnflux_xy_comp

  end program tnf_xy_onelevel

! external subroutine

  subroutine of_udo (un,bs)
! unformatted, direct access, old file
  implicit none
  character(160) :: file
  integer :: un,bs ! unit number and block size
  integer :: ier

  read(5,'(A160)') file
  write(6,'(A160)') file
  open(unit=un,status='old',file=file,FORM='UNFORMATTED', &
  & recl=bs,access='direct',iostat=ier)

  select case (ier)
  case(1:)
     write(6,*) " has error in opening."
     stop
  case(0)
     write(6,*) " opened successfully."
     return
  case default
     write(6,*) " file end"
     stop
  end select


end subroutine of_udo

subroutine of_udr (un,bs)
! unformatted, direct access, replace file
  implicit none
  character(160) :: file
  integer :: un,bs ! unit number and block size
  integer :: ier

  read(5,'(A160)') file
  write(6,'(A160)') file
  open(unit=un,status='replace',file=file,FORM='UNFORMATTED', &
  & recl=bs,access='direct',iostat=ier)

  select case (ier)
  case(1:)
     write(6,*) " has error in opening."
     stop
  case(0)
     write(6,*) " opened successfully."
     return
  case default
     write(6,*) " file end"
     stop
  end select

end subroutine of_udr
