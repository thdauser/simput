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
  type(simctl),pointer   :: cat
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
  call simopctl("simput.fits" // c_null_char, 0, 0, 0, 0, 0,cat, status)

  ! Determine the number of entries in the catalog.
  write(*,*) "Catalog contains", cat%nentries, "entries"

  ! Loop over all sources.
  do row=1,cat%nentries,1

     ! Get a source from the catalog (numbering starts at 1).
     write(*,*) " Get source number", row, "from the catalog"
     call simlcsrc(cat, row, src, status)

     ! Print the coordinates of the source.
     write(*,10) src%ra*180/pi, src%dec*180/pi
10   format (1x, ' RA=', F6.2, 1x, ' Dec=', F6.2)

     ! Determine the spectrum of the source.
     call simrspec(cat, src, time, mjdref, spec, status)

     ! Print the number of data points of the spectrum.
     write(*,*) " Spectrum contains", spec%nentries, "data points"

  end do
  ! End of loop over all catalog entries.



  ! Print out the spectrum of the first source in the catalog.
  write(*,*) "Print the spectrum of the first source ..."
  row=1
  call simlcsrc(cat, row, src, status)
  call simrspec(cat, src, time, mjdref, spec, status)

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
  call simfrctl(cat, status)
 
  ! Print the error status.
  write(*,*) "Error status: ", status
  write(*,*) "Finished."

end program fortran_sample

