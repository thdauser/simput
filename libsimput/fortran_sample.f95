program fortran_sample
  use, intrinsic :: iso_c_binding
  use simput
  implicit none

  ! PI.
  real, parameter :: pi = 3.14159265

  ! Row number in the catalog.
  integer(kind=c_long) :: row = 1 

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
  source = simlcsrc(catalog, row, status)

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

end program fortran_sample
