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
     integer(kind=c_int)  :: clc
     real(kind=c_float)   :: fra
     real(kind=c_float)   :: fdec
     real(kind=c_float)   :: fimgrota
     real(kind=c_float)   :: fe_min
     real(kind=c_float)   :: fe_max
     real(kind=c_float)   :: fflux
     type(c_ptr)          :: filename
     type(c_ptr)          :: filepath
     type(c_ptr)          :: srcbuff
  end type simctl

  ! Define source data structure.
  type, bind(c) :: simsrc
     integer(kind=c_long) :: src_id
     type(c_ptr)          :: src_name
     real(kind=c_double)  ::ra
     real(kind=c_double)  :: dec
     real(kind=c_float)   :: imgrota
     real(kind=c_float)   :: imgscal
     real(kind=c_float)   :: e_min
     real(kind=c_float)   :: e_max
     real(kind=c_float)   :: eflux
     type(c_ptr)          :: spectrum
     type(c_ptr)          :: image
     type(c_ptr)          :: lightcur
     type(c_ptr)          :: filename
     type(c_ptr)          :: filepath
  end type simsrc

  ! Define spectrum data structure.
  ! C version.
  type, bind(c) :: simspec
     integer(kind=c_long) :: nentries
     type(c_ptr)          :: energy
     type(c_ptr)          :: pflux
     type(c_ptr)          :: distr
     type(c_ptr)          :: name
     type(c_ptr)          :: fileref
  end type simspec


  ! C function and subroutine interfaces.
  interface
     function opensimputcatalog(filename,mode,status) bind(c,name="openSimputCatalog")
       use, intrinsic :: iso_c_binding
       character(kind=c_char),dimension(*) :: filename
       integer(kind=c_int),value           :: mode
       integer(kind=c_int)                 :: status
       type(c_ptr)                         :: opensimputcatalog
     end function opensimputcatalog
     
     subroutine freesimputcatalog(catalog,status) bind(c,name="freeSimputCatalog")
       use, intrinsic :: iso_c_binding
       type(c_ptr)          :: catalog
       integer(kind=c_int)  :: status
     end subroutine freesimputcatalog
     
     function loadcachesimputsource(catalog, row, status) bind(c,name="loadCacheSimputSource")
       use, intrinsic :: iso_c_binding
       type(c_ptr),value          :: catalog
       integer(kind=c_long),value :: row
       integer(kind=c_int)        :: status
       type(c_ptr)                :: loadcachesimputsource
     end function loadcachesimputsource

     function returnsimputsrcspec(src, time, mjdref, status) bind(c,name="returnSimputSrcSpec")
       use, intrinsic :: iso_c_binding
       type(c_ptr),value          :: src
       real(kind=c_double),value  :: time
       real(kind=c_double),value  :: mjdref
       integer(kind=c_int)        :: status
       type(c_ptr)                :: returnsimputsrcspec
     end function returnsimputsrcspec
     
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
    subroutine simopctl(filename,mode,catalog,status)
      use, intrinsic :: iso_c_binding
      character(kind=c_char),dimension(*) :: filename
      integer(kind=c_int),value           :: mode
      type(simctl), pointer :: catalog
      integer(kind=c_int)   :: status
      type(c_ptr)           :: ccatalog
      
      ccatalog=opensimputcatalog(filename,mode,status)
      call c_f_pointer(ccatalog, catalog)      
    end subroutine simopctl

    ! This subroutine closes an open SIMPUT source catalog defined by the 
    ! catalog data structure and releases the occupied memory.
    subroutine simfrctl(catalog,status)
      use, intrinsic :: iso_c_binding
      type(simctl), pointer :: catalog
      integer(kind=c_int)   :: status
      type(c_ptr)           :: ccatalog

      ccatalog=c_loc(catalog)
      call freesimputcatalog(ccatalog,status)
    end subroutine simfrctl


    ! This function returns a pointer to a SIMPUT source data structure
    ! for the source in the specified line of the catalog.
    subroutine simlcsrc(catalog, row, src, status)      
      use, intrinsic :: iso_c_binding
      type(simctl), pointer      :: catalog
      integer(kind=c_long)       :: row
      type(simsrc), pointer      :: src
      integer(kind=c_int)        :: status
      type(c_ptr)                :: ccatalog
      type(c_ptr)                :: csrc

      ccatalog=c_loc(catalog)
      csrc=loadcachesimputsource(ccatalog,row,status)
      call c_f_pointer(csrc, src)
    end subroutine simlcsrc


    ! This subroutine returns a pointer to a data structure containing the
    ! SIMPUT mission-independent spectrum of the specified source.
    subroutine simrspec(src,time,mjdref,spec,status)
      use, intrinsic :: iso_c_binding
      type(simsrc),pointer       :: src
      real(kind=c_double), value :: time
      real(kind=c_double), value :: mjdref
      integer(kind=c_int)        :: status
      type(simspec), pointer     :: spec
      type(c_ptr)                :: csrc
      type(c_ptr)                :: cspec

      csrc=c_loc(src)
      cspec=returnsimputsrcspec(csrc,time,mjdref,status)
      call c_f_pointer(cspec, spec)
    end subroutine simrspec
   
end module simput

