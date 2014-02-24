!   This file is part of SIMPUT.
!
!   SIMPUT is free software: you can redistribute it and/or modify it
!   under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   any later version.
!
!   SIMPUT is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!
!   For a copy of the GNU General Public License see
!   <http://www.gnu.org/licenses/>.
!
!
!   Copyright 2007-2014 Christian Schmid, FAU


module simput
  use, intrinsic :: iso_c_binding
  implicit none

  ! Define catalog data structure.
  type, bind(c) :: simctl
     type(c_ptr)          :: fptr
     integer(kind=c_long) :: nentries
     integer(kind=c_int)  :: csrcid
     integer(kind=c_int)  :: csrcname
     integer(kind=c_int)  :: cra
     integer(kind=c_int)  :: cdec
     integer(kind=c_int)  :: cimgrota
     integer(kind=c_int)  :: cimgscal
     integer(kind=c_int)  :: ce_min
     integer(kind=c_int)  :: ce_max
     integer(kind=c_int)  :: cflux
     integer(kind=c_int)  :: cspec
     integer(kind=c_int)  :: cimage
     integer(kind=c_int)  :: ctiming
     real(kind=c_float)   :: fra
     real(kind=c_float)   :: fdec
     real(kind=c_float)   :: fimgrota
     real(kind=c_float)   :: fe_min
     real(kind=c_float)   :: fe_max
     real(kind=c_float)   :: fflux
     type(c_ptr)          :: filename
     type(c_ptr)          :: filepath
     type(c_ptr)          :: srcbuff
     type(c_ptr)          :: midpspecbuff
     type(c_ptr)          :: phlistbuff
     type(c_ptr)          :: lcbuff
     type(c_ptr)          :: psdbuff
     type(c_ptr)          :: imgbuff
     type(c_ptr)          :: specbuff
     type(c_ptr)          :: krlcbuff
     type(c_ptr)          :: arf
  end type simctl

  ! Define source data structure.
  type, bind(c) :: simsrc
     integer(kind=c_long) :: src_id
     type(c_ptr)          :: src_name
     real(kind=c_double)  :: ra
     real(kind=c_double)  :: dec
     real(kind=c_float)   :: imgrota
     real(kind=c_float)   :: imgscal
     real(kind=c_float)   :: e_min
     real(kind=c_float)   :: e_max
     real(kind=c_float)   :: eflux
     real(kind=c_float)   :: phrate
     type(c_ptr)          :: spectrum
     type(c_ptr)          :: image
     type(c_ptr)          :: timing
  end type simsrc

  ! Define spectrum data structure.
  ! C version.
  type, bind(c) :: simspec
     integer(kind=c_long) :: nentries
     type(c_ptr)          :: energy
     type(c_ptr)          :: pflux
     type(c_ptr)          :: refflux
     type(c_ptr)          :: distr
     type(c_ptr)          :: name
     type(c_ptr)          :: fileref
  end type simspec


  ! C function and subroutine interfaces.
  interface
     function opensimputctlg(filename,mode,msname,msspectrum,msimage,mstiming,status) bind(c,name="openSimputCtlg")
       use, intrinsic :: iso_c_binding
       character(kind=c_char),dimension(*) :: filename
       integer(kind=c_int),value           :: mode
       integer(kind=c_int),value           :: msname
       integer(kind=c_int),value           :: msspectrum
       integer(kind=c_int),value           :: msimage
       integer(kind=c_int),value           :: mstiming
       integer(kind=c_int)                 :: status
       type(c_ptr)                         :: opensimputctlg
     end function opensimputctlg
     
     subroutine freesimputctlg(cat,status) bind(c,name="freeSimputCtlg")
       use, intrinsic :: iso_c_binding
       type(c_ptr)          :: cat
       integer(kind=c_int)  :: status
     end subroutine freesimputctlg
     
     function getsimputsrc(cat, row, status) bind(c,name="getSimputSrc")
       use, intrinsic :: iso_c_binding
       type(c_ptr),value          :: cat
       integer(kind=c_long),value :: row
       integer(kind=c_int)        :: status
       type(c_ptr)                :: getsimputsrc
     end function getsimputsrc

     function getsimputsrcmidpspec(cat, src, time, mjdref, status) bind(c,name="getSimputSrcMIdpSpec")
       use, intrinsic :: iso_c_binding
       type(c_ptr),value          :: cat
       type(c_ptr),value          :: src
       real(kind=c_double),value  :: time
       real(kind=c_double),value  :: mjdref
       integer(kind=c_int)        :: status
       type(c_ptr)                :: getsimputsrcmidpspec
     end function getsimputsrcmidpspec
     
     subroutine getsimputmidpspecval(spec, row, energy, pflux, status) bind(c,name="getSimputMIdpSpecVal")
       use, intrinsic             :: iso_c_binding
       type(c_ptr),value          :: spec
       integer(kind=c_long),value :: row
       real(kind=c_float)         :: energy
       real(kind=c_float)         :: pflux
       integer(kind=c_int)        :: status
     end subroutine getsimputmidpspecval

     ! TODO? simputSetARF(struct ARF* const arf);
     ! TODO? simputSetRndGen(double(*rndgen)(void));
     ! TODO? getSimputPhotonEnergy(const SimputSource* const src, const double time, const double mjdref, int* const status);
     ! TODO? getSimputPhotonRate(const SimputSource* const src, const double time, const double mjdref, int* const status);
     ! TODO? getSimputPhotonTime(const SimputSource* const src, double prevtime, const double mjdref, int* const failed, int* const status);
     ! TODO? getSimputPhotonCoord()
  end interface
   
  
  ! Definition of fortran wrappers around c functions.
  contains

    ! This subroutine opens a SIMPUT source catalog. It returns a pointer
    ! to a catalog data structure.
    subroutine simopctl(filename,mode,msname,msspectrum,msimage,mstiming,cat,status)
      use, intrinsic :: iso_c_binding
      character(kind=c_char),dimension(*) :: filename
      integer(kind=c_int),value           :: mode
      integer(kind=c_int),value           :: msname
      integer(kind=c_int),value           :: msspectrum
      integer(kind=c_int),value           :: msimage
      integer(kind=c_int),value           :: mstiming
      type(simctl), pointer               :: cat
      integer(kind=c_int)                 :: status
      type(c_ptr)                         :: ccat
      
      ccat=opensimputctlg(filename,mode,msname,msspectrum,msimage,mstiming,status)
      call c_f_pointer(ccat, cat)
    end subroutine simopctl

    ! This subroutine closes an open SIMPUT source catalog defined by the 
    ! catalog data structure and releases the occupied memory.
    subroutine simfrctl(cat,status)
      use, intrinsic :: iso_c_binding
      type(simctl), pointer :: cat
      integer(kind=c_int)   :: status
      type(c_ptr)           :: ccat

      ccat=c_loc(cat)
      call freesimputctlg(ccat,status)
    end subroutine simfrctl

    ! This function returns a pointer to a SIMPUT source data structure
    ! for the source in the specified line of the catalog.
    subroutine simlcsrc(cat, row, src, status)      
      use, intrinsic :: iso_c_binding
      type(simctl), pointer      :: cat
      integer(kind=c_long)       :: row
      type(simsrc), pointer      :: src
      integer(kind=c_int)        :: status
      type(c_ptr)                :: ccat
      type(c_ptr)                :: csrc

      ccat=c_loc(cat)
      csrc=getsimputsrc(ccat,row,status)
      call c_f_pointer(csrc, src)
    end subroutine simlcsrc

    ! This subroutine returns a pointer to a data structure containing the
    ! SIMPUT mission-independent spectrum of the specified source.
    subroutine simrspec(cat,src,time,mjdref,spec,status)
      use, intrinsic :: iso_c_binding
      type(simctl),pointer       :: cat
      type(simsrc),pointer       :: src
      real(kind=c_double), value :: time
      real(kind=c_double), value :: mjdref
      integer(kind=c_int)        :: status
      type(simspec), pointer     :: spec
      type(c_ptr)                :: ccat
      type(c_ptr)                :: csrc
      type(c_ptr)                :: cspec

      ccat=c_loc(cat)
      csrc=c_loc(src)
      cspec=getsimputsrcmidpspec(ccat,csrc,time,mjdref,status)
      call c_f_pointer(cspec, spec)
    end subroutine simrspec

    ! This function returns the energy and flux value for a particular 
    ! entry in the SIMPUT mission-independent spectrum data structure. 
    ! The energy is given in [keV], the flux in [photons/s/cm**2/keV].
    subroutine simspecv(spec, row, energy, pflux, status)      
      use, intrinsic :: iso_c_binding
      type(simspec), pointer     :: spec
      integer(kind=c_long)       :: row
      real(kind=c_float)         :: energy
      real(kind=c_float)         :: pflux
      integer(kind=c_int)        :: status
      type(c_ptr)                :: cspec

      cspec=c_loc(spec)
      call getsimputmidpspecval(cspec, row, energy, pflux, status)
    end subroutine simspecv
   
end module simput

