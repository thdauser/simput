/*
   This file is part of SIMPUT.

   SIMPUT is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIMPUT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "rmf.h"
#include "common.h"


struct RMF* getRMF(int* const status)
{
  struct RMF* rmf=(struct RMF*)malloc(sizeof(struct RMF));
  CHECK_NULL_RET(rmf, *status, "memory allocation for RMF failed", rmf);

  // Initialize.
  rmf->NumberChannels     =0;
  rmf->NumberEnergyBins   =0;
  rmf->NumberTotalGroups  =0;
  rmf->NumberTotalElements=0;
  rmf->FirstChannel       =0;
  rmf->isOrder            =0;
  rmf->NumberGroups       =NULL;
  rmf->FirstGroup         =NULL;
  rmf->FirstChannelGroup  =NULL;
  rmf->NumberChannelGroups=NULL;
  rmf->FirstElement       =NULL;
  rmf->OrderGroup         =NULL;
  rmf->LowEnergy          =NULL;
  rmf->HighEnergy         =NULL;
  rmf->Matrix             =NULL;
  rmf->ChannelLowEnergy   =NULL;
  rmf->ChannelHighEnergy  =NULL;
  rmf->AreaScaling        =0.0;
  rmf->ResponseThreshold  =0.0;
  rmf->ChannelType[0]='\0';
  rmf->RMFVersion[0]='\0';
  rmf->EBDVersion[0]='\0';
  rmf->Telescope[0]='\0';
  rmf->Instrument[0]='\0';
  rmf->Detector[0]='\0';
  rmf->Filter[0]='\0';
  rmf->RMFType[0]='\0';
  rmf->RMFExtensionName[0]='\0';
  rmf->EBDExtensionName[0]='\0';

  return(rmf);
}


struct RMF* loadRMF(char* const filename, int* const status)
{
  // First, check validity of RMF
  checkRMF(filename, status);
  CHECK_STATUS_RET(*status, NULL);

  struct RMF* rmf=getRMF(status);
  CHECK_STATUS_RET(*status, rmf);

  // Load the RMF from the FITS file using the HEAdas access routines
  // (part of libhdsp).
  fitsfile* fptr=NULL;
  fits_open_file(&fptr, filename, READONLY, status);
  CHECK_STATUS_RET(*status, rmf);

  // Read the 'SPECRESP MATRIX' or 'MATRIX' extension:
  *status=ReadRMFMatrix(fptr, 0, rmf);
  CHECK_STATUS_RET(*status, rmf);

  // Print some information:
  headas_chat(5, "RMF loaded with %ld energy bins and %ld channels\n",
	      rmf->NumberEnergyBins, rmf->NumberChannels);

  // Check if the RMF file contains matrix rows with a sum of more than 1.
  // In that case the RSP probably also contains the mirror ARF, what should
  // not be the case for this simulation. Row sums with a value of less than
  // 1 should actually also not be used, but can be handled by the simulation.
  long igrp;
  double min_sum=1.;
  double total_sum = 0.0;

  double sumresp;
  long ie, i, j;
  for (ie=0; ie<rmf->NumberEnergyBins; ie++) {

    /* sum up the response in this energy */

    sumresp = 0.0;

    for (i=0; i<rmf->NumberGroups[ie]; i++) {
       igrp = i + rmf->FirstGroup[ie];
       for (j=0; j<rmf->NumberChannelGroups[igrp]; j++) {
         sumresp += rmf->Matrix[j+rmf->FirstElement[igrp]];
       }
    }
   if (sumresp>1.1) {
      SIMPUT_ERROR("RMF contains rows with a sum > 1.1 (probably contains ARF)");
      *status=EXIT_FAILURE;
      return(rmf);
    }
    if (sumresp<min_sum) {
      min_sum=sumresp;
    }
    total_sum += sumresp/rmf->NumberEnergyBins;
  }

  if (total_sum<0.995) {
    SIMPUT_WARNING("RMF is not normalized");
    printf(" (total average sum of each RMF row %e) \n",total_sum);
  }

  // Check consistency of FirstChannelGroup values.
  // Note: This check is currently ignored in checkRMF because of a bug in the
  // HEASoft ftchkrmf code (see heacore/heasp/rmf.cxx). Here they throw an error
  // when F_CHAN >= NumberChannels, but the correct condition should be
  // F_CHAN > NumberChannels.
  for (long ii = 0; ii < rmf->NumberTotalGroups; ii++) {
    if ( (rmf->FirstChannelGroup[ii] < rmf->FirstChannel) ||
         (rmf->FirstChannelGroup[ii] > rmf->NumberChannels) ) {
      char msg[SIMPUT_MAXSTR];
      snprintf(msg, sizeof(msg), "FirstChannelGroup has invalid value (%li) for group %li. "
        "Should be >= %li and <= %li", rmf->FirstChannelGroup[ii], ii, rmf->FirstChannel,
        rmf->NumberChannels);
      SIMPUT_WARNING(msg);
    }
  }

  // Close the open FITS file.
  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, rmf);

  // Read the EBOUNDS extension.
  loadEbounds(rmf, filename, status);
  CHECK_STATUS_RET(*status, rmf);

  return(rmf);
}


void loadArfRmfFromRsp(char* const filename,
		       struct ARF** arf,
		       struct RMF** rmf,
		       int* const status)
{
  *rmf=getRMF(status);
  CHECK_STATUS_VOID(*status);

  // Load the RSP from the FITS file using the HEAdas access routines
  // (part of libhdsp).
  fitsfile* fptr=NULL;
  fits_open_file(&fptr, filename, READONLY, status);
  CHECK_STATUS_VOID(*status);

  // Read the 'SPECRESP MATRIX' or 'MATRIX' extension:
  *status=ReadRMFMatrix(fptr, 0, *rmf);
  CHECK_STATUS_VOID(*status);

  // Print some information:
  headas_chat(5, "RSP loaded with %ld energy bins and %ld channels\n",
	      (*rmf)->NumberEnergyBins, (*rmf)->NumberChannels);

  // Allocate memory for the ARF.
  *arf=getARF(status);
  CHECK_STATUS_VOID(*status);

  // Produce an ARF from the RSP data.
  (*arf)->NumberEnergyBins=(*rmf)->NumberEnergyBins;
  (*arf)->LowEnergy=(float*)malloc((*arf)->NumberEnergyBins*sizeof(float));
  CHECK_NULL_VOID((*arf)->LowEnergy, *status,
		  "memory allocation for energy bins failed");
  (*arf)->HighEnergy=(float*)malloc((*arf)->NumberEnergyBins*sizeof(float));
  CHECK_NULL_VOID((*arf)->HighEnergy, *status,
		  "memory allocation for energy bins failed");
  (*arf)->EffArea=(float*)malloc((*arf)->NumberEnergyBins*sizeof(float));
  CHECK_NULL_VOID((*arf)->EffArea, *status,
		  "memory allocation for effective area failed");
  strcpy((*arf)->Telescope, (*rmf)->Telescope);
  strcpy((*arf)->Instrument, (*rmf)->Instrument);
  strcpy((*arf)->Detector, (*rmf)->Detector);
  strcpy((*arf)->Filter, (*rmf)->Filter);

  // Calculate the row sums.
  long bincount;
  int notrmf=0;
  for (bincount=0; bincount<(*rmf)->NumberEnergyBins; bincount++) {
    double sum=0.;
    long chancount;
    for (chancount=0; chancount<(*rmf)->NumberChannels; chancount++) {
      sum+=ReturnRMFElement(*rmf, chancount, bincount);
    }

    // Check if this might be an RMF and not an RSP, as it should be.
    if ((0==notrmf)&&(fabs(sum-1.0)<0.001)) {
      SIMPUT_WARNING("response matrix declared as RSP file looks like RMF");
    } else {
      notrmf=1;
    }

    (*arf)->LowEnergy[bincount]=(*rmf)->LowEnergy[bincount];
    (*arf)->HighEnergy[bincount]=(*rmf)->HighEnergy[bincount];
    (*arf)->EffArea[bincount]=(float)sum;
  }

  // Normalize the RSP to obtain an ARF.
  NormalizeRMF(*rmf);

  // Close the open FITS file.
  fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);

  // Read the EBOUNDS extension.
  loadEbounds(*rmf, filename, status);
  CHECK_STATUS_VOID(*status);
}


// Calls the HEASoft FTOOL ftchkrmf in the following way:
// ftchkrmf infile=infile outfile=outfile clobber=yes
static void call_ftchkrmf(char* const infile, char* const outfile,
                          int* const status) {
  // Check for previous error
  CHECK_STATUS_VOID(*status);

  // Piece together the ftchkrmf call
  char cmd[SIMPUT_MAXSTR];
  snprintf(cmd, sizeof(cmd), "ftchkrmf infile=%s outfile=%s clobber=yes",
           infile, outfile);

  // Call ftchkrmf
  FILE *fp = popen(cmd, "r");

  // Cleanup
  pclose(fp);
}


// Prints a warning about failed ftchkrmf validity check on filename
void print_ftchkrmf_warning(char* const filename, int* warning_printed) {
  if (!(*warning_printed)) { // Only print this warning once
    char msg[SIMPUT_MAXSTR];
    snprintf(msg, sizeof(msg), "RMF validity check failed for %s\n"
             "The following problems have been identified", filename);
    SIMPUT_WARNING(msg);
    headas_chat(3, "\n=== Start of RMF validity check report ===\n\n");
    *warning_printed = 1;
  }
}


// Checks the output of ftchkrmf (on filename) for errors.
static void check_ftchkrmf_output(char* const filename, char* const outfile,
                                  int* const status) {
  // Check for previous program errors
  CHECK_STATUS_VOID(*status);

  // Open outfile
  FILE* fp = fopen(outfile, "r");
  CHECK_NULL_VOID(fp, *status, "Failed to read ftchkrmf output");

  // Check outfile line by line (should be empty or only contain optional hints
  // if RMF is valid)
  size_t len1 = SIMPUT_MAXSTR;
  char* err_str1 = malloc(len1 * sizeof(*err_str1));
  size_t len2 = len1;
  char* err_str2 = malloc(len2 * sizeof(*err_str1));;
  ssize_t read;
  int warning_printed = 0;
  while ( (read = getline(&err_str1, &len1, fp)) != -1 ) {
    // Ignore warnings about optional properties
    if (strstr(err_str1, "optional") != NULL) {
      // Skip this line (i.e., the error description) and the
      // following (i.e., the error value)
      read = getline(&err_str1, &len1, fp);
      continue;
    }

    // Ignore warning if it is only about missing HDUVERS1 keyword (obsolete since 1998)
    if (strstr(err_str1, "mandatory keywords") != NULL) {
      // Read next line (i.e., the keyword name)
      read = getline(&err_str2, &len2, fp);
      if (strncmp(err_str2, " HDUVERS1 \n", 12) == 0) {
        continue;
      }

      // Missing keyword is not just HDUVERS1, so print warning and error message
      print_ftchkrmf_warning(filename, &warning_printed);
      headas_chat(3, "%s", err_str1);
      headas_chat(3, "%s", err_str2);
      continue;
    }

    // Ignore warning about invalid FirstChannelGroup values. This is checked
    // in loadRMF.
    if (strstr(err_str1, "FirstChannelGroup has invalid value") != NULL) {
      // Just skip this line
      continue;
    }

    // Ignore strings that are just the newline character
    if (err_str1[0] == '\n') continue;

    // The error is not optional or about missing HDUVERS1 keyword, so print it.
    print_ftchkrmf_warning(filename, &warning_printed);
    headas_chat(3, "%s", err_str1);
  }

  // Cleanup
  if (warning_printed) {
    headas_chat(3, "\n=== End of RMF validity check report ===\n\n");
  }
  fclose(fp);
  free(err_str1);
  free(err_str2);
}


void checkRMF(char* const filename, int* const status) {
  // Check for previous error
  CHECK_STATUS_VOID(*status);

  // Check if ftchkrmf from the HEASoft FTOOLS if available on the system
  if (system("which ftchkrmf > /dev/null 2>&1")) {
    headas_chat(5, "\n*** Warning in %s: Failed to run ftchkrmf (should "
                   "be part of ftools) for %s! ***\n", __func__, filename);
    return;
  }

  // Call ftchkrmf to check validity of RMF
  char* outfile = "rmf.log";
  call_ftchkrmf(filename, outfile, status);

  // Check ftchkrmf output for errors
  check_ftchkrmf_output(filename, outfile, status);

  // Remove logfile created by ftchkrmf
  if (access(outfile, F_OK) == 0) {
    remove(outfile);
  }
}


void freeRMF(struct RMF* const rmf)
{
  if (NULL!=rmf) {
    if (NULL!=rmf->NumberGroups) {
      free(rmf->NumberGroups);
    }
    if (NULL!=rmf->FirstGroup) {
      free(rmf->FirstGroup);
    }
    if (NULL!=rmf->FirstChannelGroup) {
      free(rmf->FirstChannelGroup);
    }
    if (NULL!=rmf->NumberChannelGroups) {
      free(rmf->NumberChannelGroups);
    }
    if (NULL!=rmf->FirstElement) {
      free(rmf->FirstElement);
    }
    if (NULL!=rmf->OrderGroup) {
      free(rmf->OrderGroup);
    }
    if (NULL!=rmf->LowEnergy) {
      free(rmf->LowEnergy);
    }
    if (NULL!=rmf->HighEnergy) {
      free(rmf->HighEnergy);
    }
    if (NULL!=rmf->Matrix) {
      free(rmf->Matrix);
    }
    if (NULL!=rmf->ChannelLowEnergy) {
      free(rmf->ChannelLowEnergy);
    }
    if (NULL!=rmf->ChannelHighEnergy) {
      free(rmf->ChannelHighEnergy);
    }

    free(rmf);
  }
}

/** this routine is a clone of the HEASP routine "ReturnChannel", with the only exception that
 *  the sixte RNG is used
 */
static void ReturnChannelSixte(struct RMF *rmf, float energy, int NumberPhoton, long *channel)
{
  long upper, lower, middle, energybin;
  int i, j, igrp, ichan;
  float *sumresponse;
  float *RandomNumber;

  /* initialize the output array to -1s in the event that either the input energy is
     outside the response range or that the response does not sum to unity and events
     can all off the end of the channels. */

  for (i=0; i<NumberPhoton; i++) channel[i] = -1;

  lower = 0;
  upper = rmf->NumberEnergyBins-1;

  /* trap the case of the energy being outside the response range */

  if ( energy < rmf->LowEnergy[lower] || energy > rmf->HighEnergy[upper] ) return;

  /* find the energy bin associated with the input energy - assumes the energies are in increasing order */

  while ( upper - lower > 1 ) {
    middle = (upper + lower)/2;
    if ( energy < rmf->HighEnergy[middle] ) {
      upper = middle;
    } else {
      lower = middle;
    }
  }
  if ( energy > rmf->HighEnergy[lower] ) {
    energybin = upper;
  } else {
    energybin = lower;
  }

  /* generate an array of size channel each element of which is the integrated response up to and
     including that channel */

  sumresponse = (float *) malloc(rmf->NumberChannels*sizeof(float));
  for (i=0; i<rmf->NumberChannels; i++) sumresponse[i] = 0.;
  for (i=0; i<rmf->NumberGroups[energybin]; i++) {
    igrp = i + rmf->FirstGroup[energybin];
    for (j=0; j<rmf->NumberChannelGroups[igrp]; j++) {
      ichan = j + rmf->FirstChannelGroup[igrp]-rmf->FirstChannel;
      sumresponse[ichan] = rmf->Matrix[j+rmf->FirstElement[igrp]];
    }
  }
  for (i=1; i<rmf->NumberChannels; i++) sumresponse[i] += sumresponse[i-1];

  /* generate random numbers between 0 and 1 */

  int status = EXIT_SUCCESS;
  RandomNumber = (float *) malloc(NumberPhoton*sizeof(float));
  for (i=0; i<NumberPhoton; i++) RandomNumber[i] = (float) getRndNum(&status);

  CHECK_STATUS_VOID(status);

  /* loop round the photons */

  for (i=0; i<NumberPhoton; i++) {

  /* find the array element containing this random number. note that we do
     not assume that the total response sums to 1 - if the random number
     exceeds the total response then we assume that the event fell off the
     end of the channel array and return a -1 */

    lower = 0;
    upper = rmf->NumberChannels-1;

    if ( RandomNumber[i] <= sumresponse[upper] ) {
      while ( upper - lower > 1 ) {
	middle = (upper + lower)/2;
	if ( RandomNumber[i] < sumresponse[middle] ) {
	  upper = middle;
	} else {
	  lower = middle;
	}
      }
      if ( RandomNumber[i] > sumresponse[lower] ) {
	channel[i] = upper;
      } else {
	channel[i] = lower;
      }

  /* correct the channel number for the first channel number in use in t
he response matrix */

      channel[i] += rmf->FirstChannel;

    }

  }

  /* memory tidy-up */

  free(sumresponse);
  free(RandomNumber);

  return;
}


void returnRMFChannel(struct RMF *rmf,
		      const float energy,
		      long* const channel)
{
  ReturnChannelSixte(rmf, energy, 1, channel);
}




void loadEbounds(struct RMF* rmf, char* const filename, int* const status)
{
  // Check if the rmf data structure has been initialized.
  if (NULL==rmf) {
    *status=EXIT_FAILURE;
    SIMPUT_ERROR("memory for RMF data structure not allocated");
    return;
  }

  // Load the EBOUNDS from the FITS file using the HEAdas access routines
  // (part of libhdsp).
  fitsfile* fptr=NULL;
  fits_open_file(&fptr, filename, READONLY, status);
  CHECK_STATUS_VOID(*status);

  // Read the EBOUNDS extension:
  *status=ReadRMFEbounds(fptr, 0, rmf);
  CHECK_STATUS_VOID(*status);

  // Close the open FITS file.
  fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);
}


long getEBOUNDSChannel(const float energy, const struct RMF* const rmf)
{
  // In case there is no RMF, just return a negativ number (-1).
  if (NULL==rmf) return(-1);

  // Check if the charge is outside the range of the energy bins defined
  // in the EBOUNDS table. In that case the return value of this function
  // is '-1'.
  if (rmf->ChannelLowEnergy[0] > energy) {
    return(-1);
  } else if (rmf->ChannelHighEnergy[rmf->NumberChannels-1] < energy) {
    return(-1);
    //return(rmf->NumberChannels - 1 + rmf->FirstChannel);
  }

  // Perform a binary search to obtain the detector PHA channel
  // that corresponds to the given detector charge.
  long min, max, mid;
  min=0;
  max=rmf->NumberChannels-1;
  while (max>min) {
    mid=(min+max)/2;
    if (rmf->ChannelHighEnergy[mid]<energy) {
      min=mid+1;
    } else {
      max=mid;
    }
  }

  // Return the PHA channel.
  return(min+rmf->FirstChannel);
}


void getEBOUNDSEnergyLoHi(const long channel,
			  const struct RMF* const rmf,
			  float* const lo,
			  float* const hi,
			  int* const status)
{
  // Subtract the channel offset (EBOUNDS may either start at 0 or at 1).
  long mchannel=channel-rmf->FirstChannel;
  if ((mchannel<0) || (mchannel>=rmf->NumberChannels)) {
    char msg[SIMPUT_MAXSTR];
    sprintf(msg, "channel %ld is outside allowed range (%ld-%ld) ",
	    channel, rmf->FirstChannel, rmf->FirstChannel+rmf->NumberChannels-1);
    SIMPUT_ERROR(msg);
    *status=EXIT_FAILURE;
    return;
  }

  // Return the lower and upper energy value of the specified PHA
  // channel according to the EBOUNDS table.
  *lo=rmf->ChannelLowEnergy[mchannel];
  *hi=rmf->ChannelHighEnergy[mchannel];
}
