version [2.5.0]
  - Update testing framework
  - Fix bug with simputmerge only working on files in
    the current working directory
  - Update external libraries
  - Fix Mac installation issues
  - Fix RMF validity check writing temporary files, which caused
    issues in some batch processing

sssssss [2.4.10]
  - fixes a bug in RMF handling

version [2.4.9]
  - updates internal handling of RMFs

version [2.4.8]
  - fixes an error while writing very long input parameters
  - adds new RMF validity check
  - adds new function to calculate source extends from images
    or photon lists via great circle distance

version [2.4.7]
  - fixes a bug that caused problems in Sixte ero_vis tool
    * The positions of point-like sources were not taken into account correctly
      in ero_vis since v2.4.6 which led to wrong GTIs.

version [2.4.6]
  - fixes visibility of extended sources
    * The centers of extended sources were determined without use of their
      RA/Dec placements in the source catalog, making image extensions not
      placed at 0/0 not visible
  - if an image reference is specified for a source, check for valid IMGSCAL
    values (positive, nonzero)
  - updates "getSimputSrcExt"
    * now calulates RA, Dec of Img Center (necessitates Sixte version > 2.5.11)
  - updates cfitsio library in simput
  - updates simputmultispec
    * now checks the indices in an image reference and
      throws an error if they are out of bounds

version [2.4.5]
  - updates copyright and license notices
  - fixes segfault in simputmultispec

version [2.4.4]
  - updates the Mac dependency solver

version [2.4.3]
  - updates install.txt with information on pgplot

version [2.4.2]
  - fixes bug for time variable light curve

version [2.4.1]
  - fixes spectrum caching
    (fixes segfault when simulating more than ten sources)

version [2.4.0]
  - setting ARF to 0.0 for any bit NOT >=0
    (fixes SIXTE crash if an ARF contained NULL values)
  - increases normalization criterion in loadRMF to 10%
