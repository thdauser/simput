program fortran_sample
  use, intrinsic :: iso_c_binding
  use simput
  implicit none

  ! PI.
  real, parameter :: pi = 3.14159265
  ! MJDREF.
  real(kind=c_double), parameter :: mjdref = 55000.0
  ! Time.
  real(kind=c_double), parameter :: time = 0.0

  ! SIMPUT catalog.
  type(simctl),pointer   :: catalog
  ! Single source entry in the catalog.
  type(simsrc),pointer   :: src
  ! Source spectrum.
  type(simspec), pointer :: spec

  ! Row number in the catalog (counter variable).
  integer(kind=c_long) :: row = 1 
  ! Bin number in the spectrum (counter variable).
  integer(kind=c_long) :: entry

  ! Spectral data.
  real(kind=c_float) :: energy, flux

  ! FITS error status.
  integer :: status = 0

  ! Open the FITS file containing the SIMPUT catalog
  ! (mode=0 for read only access).
  ! The filename has to be terminated with a null character,
  ! since this is required for C strings.
  write(*,*) "Open SIMPUT catalog ..."
  call simopctl("ASM_Catalog.fits" // c_null_char, 0, catalog, status)

  ! Determine the number of entries in the catalog.
  write(*,*) "Catalog contains", catalog%nentries, "entries"

  ! Loop over all sources.
  do row=1,catalog%nentries,1

     ! Get a source from the catalog (numbering starts at 1).
     write(*,*) " Get source number", row, "from the catalog"
     call simlcsrc(catalog, row, src, status)

     ! Print the coordinates of the source.
     write(*,10) src%ra*180/pi, src%dec*180/pi
10   format (1x, ' RA=', F6.2, 1x, ' Dec=', F6.2)

     ! Determine the spectrum of the source.
     call simrspec(src, time, mjdref, spec, status)

     ! Print the number of data points of the spectrum.
     write(*,*) " Spectrum contains", spec%nentries, "data points"

  end do
  ! End of loop over all catalog entries.



  ! Print out the spectrum of the first source in the catalog.
  row=1
  call simlcsrc(catalog, row, src, status)
  call simrspec(src, time, mjdref, spec, status)

  ! Loop over all entries in the spectrum.
  do entry=0,spec%nentries-1,1
     ! Determine the energy and the flux value of this entry.
     call simspecv(spec, entry, energy, flux, status)

     ! Output
     write(*,20) energy, flux
20   format (1x, ' Energy=', F6.2, 1x, 'keV, Photon Flux=', F8.2, ' photons/s/cm**2/keV')
  end do
  ! END of printing out the spectrum of the first source in the catalog.

       

  ! Close the SIMPUT catalog FITS file.
  write(*,*) "Close SIMPUT catalog ..."
  call simfrctl(catalog, status)
 
  ! Print the error status.
  write(*,*) "Error status: ", status
  write(*,*) "Finished."

end program fortran_sample

