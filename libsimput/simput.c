/////////////////////////////////////////////////////////////////
// Function Definitions.
/////////////////////////////////////////////////////////////////


float getSimputSrcBandFlux(SimputCatalog* const cat,
			   const SimputSource* const src,
			   const double time, const double mjdref,
			   int* const status)
{
  // Determine the spectrum valid for the specified point of time.
  SimputMIdpSpec* spec=
    returnSimputSrcSpec(cat, src, time, mjdref, status);
  CHECK_STATUS_RET(*status, 0.);

  return(getSimputSpecBandFlux(spec, src->e_min, src->e_max, status));
}


