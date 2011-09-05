module simput
  implicit none

  interface
     function simopctl(file,mode,status) bind(c,name="openSimputCatalog")
       use, intrinsic :: iso_c_binding
       character(kind=c_char),dimension(*) :: file
       integer(kind=c_int),value           :: mode
       integer(kind=c_int)                 :: status
       type(c_ptr)                         :: simopctl
     end function simopctl
  end interface

  interface
     subroutine simfrctl(catalog,status) bind(c,name="freeSimputCatalog")
       use, intrinsic :: iso_c_binding
       type(c_ptr)          :: catalog
       integer(kind=c_int)  :: status
     end subroutine simfrctl
  end interface

  interface
     function simrtsrc(catalog, row, status) bind(c,name="returnSimputSource")
       use, intrinsic :: iso_c_binding
       type(c_ptr),value          :: catalog
       integer(kind=c_long),value :: row
       integer(kind=c_int)        :: status
       type(c_ptr)                :: simrtsrc
     end function simrtsrc
  end interface

  interface
     subroutine simphcrd(source, ra, dec, status) bind(c,name="getSimputPhotonCoord")
       use, intrinsic :: iso_c_binding
       type(c_ptr),value   :: source
       real(kind=c_double) :: ra
       real(kind=c_double) :: dec
       integer(kind=c_int)        :: status
     end subroutine simphcrd
  end interface
       
end module simput


! TODO Split up in module and sample program.
program test
  use, intrinsic :: iso_c_binding
  use simput
  implicit none

  ! PI.
  real, parameter :: pi = 3.14159265

  ! FITS error status.
  integer       :: status = 0

  ! SIMPUT catalog.
  type(c_ptr)   :: catalog
  ! Single source entry in the catalog.
  type(c_ptr)   :: source
  ! Right ascension and declination of the source [rad].
  real(kind=c_double) :: ra, dec


  ! Open the FITS file containing the SIMPUT catalog.
  ! (mode=0 for read only access)
  write(*,*) "Open SIMPUT catalog ..."
  catalog = simopctl("ASM_Catalog.fits" // c_null_char , 0, status)

  ! Get the first source from the catalog (numbering starts at 1).
  write(*,*) "Get first source from the catalog ..."
  source = simrtsrc(catalog, 1, status)

  ! Determine the coordinates of the first source.
  call simphcrd(source, ra, dec, status)

  ! Print the coordinates of the source.
  write(*,10) ra*180/pi, dec*180/pi
10 format (1x, 'RA=', F6.2, 1x, ' Dec=' F6.2)

  ! Close the SIMPUT catalog FITS file.
  write(*,*) "Close SIMPUT catalog ..."
  call simfrctl(catalog, status)

  ! Print the error status.
  write(*,*) "Error status: ", status
  write(*,*) "Finished."

end program test
