module simput
  implicit none

  ! TODO Define catalog data structure.
  ! TODO Define source data structure.
  ! TODO Define spectrum data structure.

  ! This function opens a SIMPUT source catalog. It returns a pointer
  ! to a catalog data structure.
  interface
     function simopctl(file,mode,status) bind(c,name="openSimputCatalog")
       use, intrinsic :: iso_c_binding
       character(kind=c_char),dimension(*) :: file
       integer(kind=c_int),value           :: mode
       integer(kind=c_int)                 :: status
       type(c_ptr)                         :: simopctl
     end function simopctl
  end interface

  ! This routine closes an open SIMPUT source catalog defined by the 
  ! catalog data structure and releases the occupied memory.
  interface
     subroutine simfrctl(catalog,status) bind(c,name="freeSimputCatalog")
       use, intrinsic :: iso_c_binding
       type(c_ptr)          :: catalog
       integer(kind=c_int)  :: status
     end subroutine simfrctl
  end interface

  ! This function returns a pointer to a SIMPUT source data structure
  ! for the source in the specified line of the catalog.
  interface
     function simlcsrc(catalog, row, status) bind(c,name="loadCacheSimputSource")
       use, intrinsic :: iso_c_binding
       type(c_ptr),value          :: catalog
       integer(kind=c_long),value :: row
       integer(kind=c_int)        :: status
       type(c_ptr)                :: simlcsrc
     end function simlcsrc
  end interface

  ! This function returns a pointer to a data structure containing the
  ! SIMPUT mission-independent spectrum of the specified source.
  interface
     function simrspec(source, time, mjdref, status) bind(c,name="returnSimputSrcSpec")
       use, intrinsic :: iso_c_binding
       type(c_ptr),value          :: source
       real(kind=c_double),value  :: time
       real(kind=c_double)        :: mjdref
       integer(kind=c_int)        :: status
       type(c_ptr)                :: simrspec
     end function simrspec
  end interface

  ! This routine returns the RA and Dec coordinates of the specified source.
  interface
     subroutine simphcrd(source, ra, dec, status) bind(c,name="getSimputPhotonCoord")
       use, intrinsic :: iso_c_binding
       type(c_ptr),value   :: source
       real(kind=c_double) :: ra
       real(kind=c_double) :: dec
       integer(kind=c_int) :: status
     end subroutine simphcrd
  end interface
       
! TODO? simputSetARF(struct ARF* const arf);
! TODO? simputSetRndGen(double(*rndgen)(void));
! TODO? getSimputPhotonEnergy(const SimputSource* const src, const double time, const double mjdref, int* const status);
! TODO? getSimputPhotonRate(const SimputSource* const src, const double time, const double mjdref, int* const status);
! TODO? getSimputPhotonTime(const SimputSource* const src, double prevtime, const double mjdref, int* const failed, int* const status);

end module simput

