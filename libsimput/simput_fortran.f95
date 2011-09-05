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

