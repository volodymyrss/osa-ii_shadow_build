/************************************************************************
 * FILE:        ii_shadow_build_main.cpp
 * VERSION:     1.14.12
 * COMPONENT:   ii_shadow_build
 * AUTHOR:      Christophe. Couvreur,     SAp-CEA,  ccouvreu@cea.fr
 * DESCRIPTION: component subroutine source code
 * HISTORY:
 * 
 *   NB: previous versions are made by S.Chazalmartin
 *   CC, 1.4.11,             Selection of events based on corrected 
 *                           rise time =[16-116]
 *   CC, 1.4.12, 12/07/2007, SPR 4666 (take into account the drift of 
 *                           energy for Lowthreshold)
 *                           - The correction is based on ibis_isgr_energy-6.0
 *                           - new parameter "protonDOL", see ibis_isgri_energy
 *   NP 1.5      5/03/2009  - SCREW 2117
 *   JZ,IC,PL 1.6 2/02/2012 New energy correction + 
 *                          evolution of the LT energy resolution
 *   PL, 1.9     02/04/2012  modify PAR1_..._corrPH1 parameters to correct <50 keV behavior
 *   CF  2.0                 SPR 05084 implemented with Paris agreement for OSA 10.1
 *   CF  2.!     05/11/2015  SCREW 2624 implemented with Paris agreement for OSA 10.2
 *   CF  2.2     09/12/2016  SCREW XXXX Piotr Lubinski's patch for high rates
 *************************************************************************/

#include "ii_shadow_build.h"

int main(int argc, char *argv[]) 
{
  int              
    status=    ISDC_OK,                      // ISDC_OK means no without mishaps                              
    RILstatus= ISDC_OK, 
    HowMany=   5,                            // Max EFF maps to write at once: for memory considerations      
    NumImaBin= 1,                            // number of energy channels                                   
    NoisyDetFlag= 0;                         // Spectral Noisy Pixels Detection Flag
  unsigned char   
    detailedOutput= 0;                       // Detailed output Yes = 1 No = 0 = default                     
  char            
    UserRowFilter[DAL_FILE_NAME_STRING]= "", // User-defined ROW filter on ISGRI events
    InGTIName[DAL_FILE_NAME_STRING]= "",
    outputLevel[PIL_LINESIZE]= "";                               
  dal_byte        
    MinRiseTime= 0,
    MaxRiseTime= 0x7F;
    
  dal_float       
    **EnergyBounds= NULL;                    // Energy bounds                               
  dal_double
    SCWduration,                             // Size of SCW
    TimeLen;                                 // Size of time bins
  dal_element
    *idxREVcontext= NULL,                    // DOL to the index of ISGRI REVOLUTION Contexts
    *idxHK3maps=    NULL,                    // DOL to the index of ISGRI HK3 noisy Maps
    *NewGRP=        NULL;                    // DOL to the SWG


  //switchOn time variables
  //double my_switchOnTime;     /* ijd of switch on time */
  int chatter=4;

  //hk1 variables
  dal_element *isgrHK1_Ptr = NULL;
  double meanT,meanBias[8];               /* mean of the 8 mce temperatures */
  dal_dataType DALtype_hk1;
  double tstart_grp;
  //int revol_scw;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //                              Initialize the common library stuff
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  if((status= CommonInit(COMPONENT_NAME, COMPONENT_VERSION, argc, argv))!=ISDC_SINGLE_MODE) 
    {
      RILstatus= RILlogMessage(NULL, Log_1,   "CommonInit status = %d", status);
      RILstatus= RILlogMessage(NULL, Log_1,   "number of command line arguments = %d", argc);
      RILstatus= RILlogMessage(NULL, Log_1,   "program name : %s", argv[0]);
      RILstatus= RILlogMessage(NULL, Error_2, "Program aborted : could not initialize.");
      CommonExit(status);
    }
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //               Get the necessary INPUT PARAMETERS given in the ii_shadow_build.par file.           
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  RILstatus= RILlogMessage(NULL,Log_1,"---------- Retrieve program PARAMETERS ----------");
  status= GetPars(&NewGRP, &idxREVcontext, &idxHK3maps, InGTIName, &EnergyBounds, &NumImaBin, &TimeLen, 
		  UserRowFilter, &MinRiseTime, &MaxRiseTime,
		  &NoisyDetFlag, outputLevel, &detailedOutput);
  if(status  != ISDC_OK ) Abort(NewGRP, "Program aborted with status %d : could not retrieve parameters.", status);
  
  // Get the science window ID stuff    
  dal_dataType 
    DALtype= DAL_CHAR;   
  T_OSM_ATTRIBUTE
    OSMattribute;

  RILstatus= RILlogMessage(NULL,Log_1,"-------- Get the Science Window (SCW) IDs --------"); 
  if((status= DALattributeGet(NewGRP, ERTFIRST, DALtype, &OSMattribute.ertfirst, NULL, NULL, status))!=ISDC_OK) 
    Abort(NewGRP, "Program aborted with status %d: could not retrieve attribute (ERTFIRST).", status);
  if((status= DALattributeGet(NewGRP, ERTLAST, DALtype, &OSMattribute.ertlast, NULL, NULL, status))!=ISDC_OK) 
    Abort(NewGRP, "Program aborted with status %d: could not retrieve attribute (ERTLAST).", status);
  if((status= DALattributeGet(NewGRP, SWID, DALtype, &OSMattribute.swid, NULL, NULL, status))!=ISDC_OK) 
    Abort(NewGRP, "Program aborted with status %d: could not retrieve attribute (SWID).", status);
  if((status= DALattributeGet(NewGRP, SWTYPE, DALtype, &OSMattribute.sw_type, NULL, NULL, status))!=ISDC_OK) 
    Abort(NewGRP, "Program aborted with status %d: could not retrieve attribute (SWTYPE).", status);
  if((status= DALattributeGet(NewGRP, SWBOUND, DALtype, &OSMattribute.swbound, NULL, NULL, status))!=ISDC_OK) 
   Abort(NewGRP, "Program aborted with status %d: could not retrieve attribute (SWBOUND).", status);
  if((status= DALattributeGet(NewGRP, BCPPID, DALtype, &OSMattribute.bcppid, NULL, NULL, status))!=ISDC_OK) 
    Abort(NewGRP, "Program aborted with status %d: could not retrieve attribute (BCPPID).", status);
  if((status= DALattributeGet(NewGRP, PREVSWID, DALtype, &OSMattribute.prevswid, NULL, NULL, status))!=ISDC_OK) 
    Abort(NewGRP, "Program aborted with status %d: could not retrieve attribute %s.", status);
  DALtype= DAL_INT;
  if((status= DALattributeGet(NewGRP, REVOL, DALtype, &OSMattribute.revol, NULL, NULL, status))!=ISDC_OK) 
    Abort(NewGRP, "Program aborted with status %d: could not retrieve attribute (REVOL).", status);
  strcpy(OSMattribute.outputlevel, outputLevel);

  // Inform user
  RILstatus= RILlogMessage(NULL, Log_1, "ERTFIRST     : %s", OSMattribute.ertfirst);
  RILstatus= RILlogMessage(NULL, Log_1, "ERTLAST      : %s", OSMattribute.ertlast);
  RILstatus= RILlogMessage(NULL, Log_1, "SWID         : %s", OSMattribute.swid);
  if(detailedOutput) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "SWTYPE       : %s", OSMattribute.sw_type);
      RILstatus= RILlogMessage(NULL, Log_1, "SWBOUND      : %s", OSMattribute.swbound);
      RILstatus= RILlogMessage(NULL, Log_1, "BCPPID       : %s", OSMattribute.bcppid);
      RILstatus= RILlogMessage(NULL, Log_1, "PREVSWID     : %s", OSMattribute.prevswid);
      RILstatus= RILlogMessage(NULL, Log_1, "REVOL        : %d", OSMattribute.revol);
      RILstatus= RILlogMessage(NULL, Log_1, "BINNING TYPE : %s", OSMattribute.outputlevel);
    }
  OSMattribute.rise_min= (int)MinRiseTime;
  OSMattribute.rise_max= (int)MaxRiseTime;
  


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //                             CHECK for INPUT and OUTPUT FILES 
  //      isgri_detector_shadowgram.fits      --> for the images            --> 0001 if not existant   
  //      isgri_detector_spectra.fits         --> for the spectra           --> 0010 if not existant   
  //      isgri_detector_lightcurve.fits      --> for the light curves      --> 0100 if not existant   
  //      isgri_detector_effi_shadowgram.fits --> for the efficiency images --> 1000 if not existant   
  // If the ISDC structures already exist then this function will return a pointer to them.            
  // If there are no events in the data: exit without creating any new data structures.                
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  if(detailedOutput)
    RILstatus= RILlogMessage(NULL,Log_1,"-------- Check INPUTS and initialize OUTPUTS --------");

  int
    IndexAction= 0; 
  long                      
    NumDeadTimes= 0; 
  OBTime                  
    StartTime,
    EndTime;
  dal_element
    *DTtable=    NULL, // DOL to Table containing the dead time values 
    *idxRAWshad= NULL, // DOL to Index for RAW images                                          
    *idxEFFshad= NULL, // DOL to Index for EFF images            
    *REVcontext= NULL; // DOL to REV context                              

  // Check files and retrieve SCW start/end
  status= ChkFilesExist(NewGRP, &StartTime, &EndTime, 
			&DTtable, &NumDeadTimes, 
			&idxRAWshad, &idxEFFshad, 
			idxREVcontext, &REVcontext,
			&IndexAction, detailedOutput);
  if(status != ISDC_OK){
    switch(status) {
    case ERR_ISGR_OSM_FILE_NECESSARY_NOTFOUND :
      RILstatus= RILlogMessage(NULL, Log_1, "Some index(es) do not exist.");
      break;
    case ERR_ISGR_OSM_FILE_NOTFOUND :
      RILstatus= RILlogMessage(NULL, Log_1, "Some file(s) do not exist.");
      status= ISDC_OK;
      Abort(NewGRP, "Error : input structures not found : status = %d", status);
      break;
    default :
      Abort(NewGRP, "Error : output structures not found : status = %d", status);
      break;
    }
  }

  if(status==ERR_ISGR_OSM_FILE_NECESSARY_NOTFOUND)
    {
      char *TmpString= (char *)calloc(8*sizeof(long)+1, sizeof(char));
      lDecToBin((long)IndexAction, TmpString);
      RILstatus= RILlogMessage(NULL, Log_1, "Some index(es) do not exist --> will NOT create them now %s", 
			       TmpString+strlen(TmpString)-4);
      free(TmpString);
      CommonExit(status);
    }

 

  /*#################################################################*/
  /* Locate bintable for ISGRI temperature and bias */
  /*#################################################################*/
  status=DALobjectFindElement(NewGRP, DS_ISGR_HK, &isgrHK1_Ptr, status);
  if (status != ISDC_OK)
    RILlogMessage(NULL, Warning_2, "%13s bintable NOT found.", DS_ISGR_HK);
  else if (chatter > 2)
    RILlogMessage(NULL, Log_0, "%13s bintable found.", DS_ISGR_HK);
  
  meanT=KEY_DEF_TEMP;
  for (int i=0; i<8; i++) meanBias[i]=KEY_DEF_BIAS;

  /*#################################################################*/
  /* Calculation of Temperature*/
  /*#################################################################*/
  status=ibis_energyIsgrHkCal(NewGRP,StartTime, EndTime,
			      &meanT,meanBias, chatter, status);
  if (status != ISDC_OK) {
    RILlogMessage(NULL, Warning_2, "Reverting from status=%d to ISDC_OK",
		  status);
      RILlogMessage(NULL, Warning_2, "Using constant ISGRI temperature and bias (%+6.2f %+6.1f)", 
                                   meanT, KEY_DEF_BIAS);
  }
    //else if (chatter > 3) {
    //RILlogMessage(NULL, Log_0, "Mean ISGRI module bias (V):");
    //strcpy(logString, "");
    //for (i=0; i<8; i++) {
    //sprintf(hkName, " %+6.1f", meanBias[i]);
    //  strcat(logString, hkName);
    //}
    //RILlogMessage(NULL, Log_0, logString);
    //}
  status=ISDC_OK;


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //                         Output Information about SCW TIME boundaries.                                        
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  // SCW Boundaries
  double
    JulianStart= 0.,
    JulianEnd=   0.;
  OBTime
    OBTStart= StartTime,
    OBTEnd= EndTime;
  if(detailedOutput)
    {
      if((status= DAL3AUXconvertOBT2IJD(NewGRP, TCOR_ANY, 1, &StartTime, &JulianStart, status))!=ISDC_OK) 
	{
	  RILstatus= RILlogMessage(NULL, Log_1,"main : Error in converting OBT times to ISDC Julian days, status = %d.", status);
	  CommonExit(status);
	}
      if((status= DAL3AUXconvertOBT2IJD(NewGRP, TCOR_ANY, 1, &EndTime, &JulianEnd, status))!=ISDC_OK) 
	{
	  RILstatus= RILlogMessage(NULL, Log_1, 
				  "main : Error in converting OBT times to ISDC Julian days, status = %d.", status);
	  CommonExit(status);
	}
      RILstatus= RILlogMessage(NULL,Log_1,"The full SCW (OBTimes)= [ %llu, %llu ].", OBTStart, OBTEnd);
      RILstatus= RILlogMessage(NULL,Log_1,"    (ISDC Julian days)= [ %f, %f ]", JulianStart, JulianEnd);
    }

  if((status= DAL3AUXconvertOBT2IJD(NewGRP, TCOR_ANY, 1, &StartTime, &OSMattribute.tstart, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "Error in converting OBT times to ISDC Julian days, status = %d.", status);
      CommonExit(status);
    }
  if((status= DAL3AUXconvertOBT2IJD(NewGRP, TCOR_ANY, 1, &EndTime, &OSMattribute.tstop, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "Error in converting OBT times to ISDC Julian days, status = %d.", status);
      CommonExit(status);
    }
  if((status= DAL3AUXconvertOBT2IJD(NewGRP, TCOR_ANY, 1, &OBTStart, &OSMattribute.tfirst, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "Error in converting OBT times to ISDC Julian days, status = %d.", status);
      CommonExit(status);
    }
  if((status= DAL3AUXconvertOBT2IJD(NewGRP, TCOR_ANY, 1, &OBTEnd, &OSMattribute.tlast, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "Error in converting OBT times to ISDC Julian days, status = %d.", status);
      CommonExit(status);
    }

  // SCW Duration
  OBTime 
    deltaOBT= 0;
  double
    TotalScwDt= 0.;
  if((status= DAL3GENelapsedOBT(EndTime, StartTime, &deltaOBT, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "Error in effecting OBT time differences, status = %d.", status);
      CommonExit(status);
    }
  if((status= DeltaOBT(EndTime, StartTime, &TotalScwDt, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Error_1, "Fatal error in OBT time difference, status = %d.", status);
      return status;
    }
  if(detailedOutput) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "Time interval of SCW : %llu (OBT)", deltaOBT);
      RILstatus= RILlogMessage(NULL, Log_1, "Time interval of SCW : %.3f (sec)", TotalScwDt);
    }


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //                      Read the REVOLUTION CONTEXT and the HK3 status                                             
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  RILstatus= RILlogMessage(NULL,Log_1,"-------- Get REVOL Context and HK3 status --------");

  // Retrieve Pixel Revolution Status
  dal_double **LowThreshMap= NULL; 
  if((LowThreshMap=(dal_double **)calloc(ISGRI_SIZE, sizeof(dal_double *)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for LowThreshold image  map.", status);
  for(int i=0; i<ISGRI_SIZE; i++) 
    if((LowThreshMap[i]=(dal_double *)calloc(ISGRI_SIZE, sizeof(dal_double )))==NULL) 
      Abort(NewGRP, "Error in allocating memory for LowThreshold image map : i = %d.", i);
  dal_int **ONpixelsREVmap= NULL;
  if((ONpixelsREVmap= (dal_int**)calloc(ISGRI_SIZE, sizeof(dal_int*)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for pixels on / off map.", status);
  for(int i=0; i<ISGRI_SIZE; i++) 
    if((ONpixelsREVmap[i]= (dal_int*)calloc(ISGRI_SIZE, sizeof(dal_int)))==NULL) 
      Abort(NewGRP, "Error in allocating memory for  pixels on / off map.", status);
  if((status= GetREVcontext(REVcontext, OSMattribute.revol, EndTime,
			    LowThreshMap, ONpixelsREVmap, detailedOutput)) !=ISDC_OK )
    {
      RILstatus= RILlogMessage(NULL,Error_1,"Retrieving Revolution Context failed, status %d.",status);
      CommonExit(status);
    }

  ///////////////////////////////////////////////////////////////////////////////////
  //                        Correction of energy drift
  ///////////////////////////////////////////////////////////////////////////////////
  
  status= DAL3AUXconvertOBT2IJD(NewGRP, TCOR_ANY, 1, &StartTime, &JulianStart, status); 
  
  if (status != ISDC_OK){
    RILstatus= RILlogMessage(NULL, Warning_1,
			     "main : Error in converting OBT times to ISDC Julian days, status = %d."
			     , status);
    //CommonExit(status);
    status=ISDC_OK;
  }
  status= EnerDriftCorr(LowThreshMap,OSMattribute.revol);
  if(status != ISDC_OK) {
    RILlogMessage(NULL,Error_1,"no energy drift correction for Low Threshold Map");
    return status;
  }

  // Free memory here for proton tables : to limit necessary memory space
  int i;


  // Retrieve Pixel HK3 Status
  dal_int **ONpixelsHK3map= NULL;
  if((ONpixelsHK3map= (dal_int**)calloc(ISGRI_SIZE,sizeof(dal_int*)))==NULL) 
    Abort(NewGRP,"Error in allocating memory for ONpixelsHK3 map[].",ERR_ISGR_OSM_MEMORY_ALLOC);
  for(int y=0;y<ISGRI_SIZE;y++) 
    if((ONpixelsHK3map[y]= (dal_int*)calloc(ISGRI_SIZE,sizeof(dal_int)))==NULL)
      Abort(NewGRP,"Error in allocating memory for ONpixelsHK3 map[][]",ERR_ISGR_OSM_MEMORY_ALLOC);
  if((status= GetHK3status(idxHK3maps, StartTime, EndTime, ONpixelsHK3map, detailedOutput)) != ISDC_OK)
    {
      RILstatus= RILlogMessage(NULL,Error_1,"Error while retrieving HK3 pixels live time, status= %d.", status);
      CommonExit(status);
    }

  // Output stats
  long NumPixON=0;
  double MeanLT=0.;
  for(int y=0;y<ISGRI_SIZE;y++)
    for(int z=0;z<ISGRI_SIZE;z++)
      if( (LowThreshMap[y][z]>0) && ONpixelsREVmap[y][z] && ONpixelsHK3map[y][z] )
	{
	  MeanLT+= LowThreshMap[y][z];
	  NumPixON++;
	}
  RILstatus= RILlogMessage(NULL,Log_1,"Number of Pixels ON (in REV and SCW) with good LowThresh = %ld.",NumPixON);
  if(NumPixON)
    {
      MeanLT/= NumPixON;
      RILstatus= RILlogMessage(NULL,Log_1,"Mean LowThreshold for these pixels= %.2fkeV.",MeanLT);
    }
  

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                    Get HOUSEKKEPING DATA :
  //                          Good Times Intervals and Module DeadTimes.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  RILstatus= RILlogMessage(NULL,Log_1,"--------  Get  SCW  DATA  --------");

  int
    NumGTI= 1,           
    NumGTIXtra= 0,                                                
    NumIntervals= 0,                      
    *NumGTIperMCE= NULL;             
  long                      
    *NumEventsPerMCE= NULL;                               
  dal_float       
    **DeadTimes= NULL; 
  OBTime                                   
    *GTIstart=        NULL,                                      
    *GTIstop=         NULL,
    **GTIstartPerMCE= NULL,                            
    **GTIstopPerMCE=  NULL, 
    *OBTdeadTimes=    NULL;

  // Allocations                           
  if((NumGTIperMCE=(int *)calloc(IBIS_NUM_BLOCK, sizeof(int)))==NULL)
    Abort(NewGRP, "Error in allocating memory for NumGTIperMCE files : status = %d", status);
  if((GTIstopPerMCE=(OBTime **)calloc(IBIS_NUM_BLOCK, sizeof(OBTime)))==NULL)
    Abort(NewGRP, "Error in allocating memory for NumGTIperMCE files : status = %d", status);
  if((GTIstartPerMCE=(OBTime **)calloc(IBIS_NUM_BLOCK, sizeof(OBTime)))==NULL)
    Abort(NewGRP, "Error in allocating memory for NumGTIperMCE files : status = %d", status);
  if((NumEventsPerMCE=(long *)calloc(IBIS_NUM_BLOCK, sizeof(long)))==NULL)
    Abort(NewGRP, "Error in allocating memory for NumGTIperMCE files : status = %d", status);
     
  // Retrieve HK data
  if((status= GetHKdata(NewGRP, DTtable, InGTIName, &GTIstart, &GTIstop, &NumIntervals, 
			GTIstartPerMCE, GTIstopPerMCE, NumGTIperMCE, &NumGTIXtra, &DeadTimes, &OBTdeadTimes,
			&NumDeadTimes, StartTime, EndTime, detailedOutput))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "Error GetHKdata status = %d", status);
      CommonExit(status);
    }

  // Match GTIs boundaries to SCW limits
  NumGTI= NumIntervals;
  OBTStart= (GTIstart[0]>StartTime) ? GTIstart[0] : StartTime;
  OBTEnd= (GTIstop[NumGTI-1]<EndTime) ? GTIstop[NumGTI-1] : EndTime;
  GTIstart[0]= OBTStart;
  GTIstop[NumGTI-1]= OBTEnd;

  // SPR 3777: Check GTIs are time coherent and not overlapping (=> outnumbering events...) 
  for(int i=0;i<NumGTI;i++)
    {
      if(GTIstop[i]<GTIstart[i])
	{ 
	  status= ERR_ISGR_OSM_DATA_INCONSISTENCY;
	  RILstatus= RILlogMessage(NULL,Warning_1,"%s: GTI_%d is not time coherent!",InGTIName,i+1);
	}
      if(i>0 && GTIstop[i-1]>GTIstart[i])
	{
	  status= ERR_ISGR_OSM_DATA_INCONSISTENCY;
	  RILstatus= RILlogMessage(NULL,Warning_1,"%s: GTI_%d overlaps GTI_%d",InGTIName,i,i+1);
	}
    }
  if(status)
    {
      RILstatus= RILlogMessage(NULL,Error_1,"Terminating process with status %d",status);
      Abort(NewGRP,"Check the GTIs and rerun.",status);
    }
  
  // Resulting Boundaries of ISGRI GTIs   
  if(detailedOutput)
    RILstatus= RILlogMessage(NULL, Log_1, "=> ISGRI GTIs'Up&Low boundaries (OBT) = [ %llu, %llu ]", 
			     GTIstart[0], GTIstop[NumGTI-1]);


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //                              Define SCW interval of analysis
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  RILstatus= RILlogMessage(NULL, Log_1, "-------- Define SCW interval of analysis --------");

  // Interval boundaries
  if((status= DAL3AUXconvertOBT2IJD(NewGRP, TCOR_ANY, 1, &OBTStart, &JulianStart, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1,"main : Error in converting OBT times to ISDC Julian days, status = %d.", status);
      CommonExit(status);
    }
  if((status= DAL3AUXconvertOBT2IJD(NewGRP, TCOR_ANY, 1, &OBTEnd, &JulianEnd, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, 
			       "main : Error in converting OBT times to ISDC Julian days, status = %d.", status);
      CommonExit(status);
    }
  RILstatus= RILlogMessage(NULL,Log_1,"SCW interval of analysis (IJD)= [ %f, %f ]", JulianStart, JulianEnd);
  if(detailedOutput)
    RILstatus= RILlogMessage(NULL,Log_1,"SCW interval of analysis (OBT)= [ %llu, %llu ]", OBTStart, OBTEnd);

  // SCW of analysis duration
  if((status= DAL3GENelapsedOBT(OBTEnd, OBTStart, &deltaOBT, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Error_1, "Error in effecting OBT time difference, status= %d.", status);
      CommonExit(status);
    }
  if((status= DeltaOBT(OBTEnd, OBTStart, &SCWduration, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Error_1, "Error in effecting OBT time difference, status= %d.", status);
      return status;
    }
  RILstatus= RILlogMessage(NULL,Log_1,"=> SCW duration= %.3f (sec)", SCWduration);
  if(detailedOutput)
    RILstatus= RILlogMessage(NULL,Log_1,"=> SCW duration= %llu (OBT)", deltaOBT);
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                          Get EVENT DATA
  //       This section might require iterating over good time intervals if memory is a problem. 
  //            If No events are present, abort here before creating any output structures.            
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  RILstatus= RILlogMessage(NULL,Log_1,"---------- Select ISGRI EVENTS ----------");
 
  char 
    RowFilter[200]= "";
    //RowFilterCorr[200]= "";
  long 
    NumEvents= 0;
  dal_byte        
    *IsgriY=          NULL,                                       
    *IsgriZ=          NULL,                                     
    *IsgriRiseTime=   NULL,
    *IsgriSelectFlag= NULL;
  dal_double
    *IsgriEnergy=     NULL,
    *ChannelMin=      NULL,                                       
    *ChannelMax=      NULL; 
  OBTime 
    *IsgriTime=       NULL;
  
  // Apply User-defined Filtering
  //sprintf(RowFilter,"RISE_TIME>=%d && RISE_TIME<=%d",MinRiseTime,MaxRiseTime);//selection on raw risetime
  sprintf(RowFilter,"ISGRI_PI>=%d && ISGRI_PI<=%d",MinRiseTime,MaxRiseTime);//selection on corr rise time
  
  if(UserRowFilter[0]!='\0')
    {
      RILstatus= RILlogMessage(NULL,Log_1, "Applying User Events_Filter: %s",UserRowFilter);
      sprintf(RowFilter,"%s && %s",RowFilter,UserRowFilter);
    }
  
  // This routine also selects only ISGRI GTIs compliant events
  if((status= DAL3IBISselectEvents(NewGRP, ISGRI_EVTS, ALL_PAR, NumGTI, GTIstart, GTIstop, RowFilter, status))!=ISDC_OK) 
    switch(status) 
      {
      case DAL3IBIS_NO_IBIS_EVENTS :
	RILstatus= RILlogMessage(NULL, Warning_1, "No IBIS events found: nothing to be done.");
	RILstatus= RILlogMessage(NULL, Warning_1, "Terminating this process.");
	status= ISDC_OK;
	CommonExit(status);
	break;
      default :
	RILstatus= RILlogMessage(NULL, Warning_1, "Cannot select IBIS events!");
	break;
      }

  if((status= DAL3IBISgetNumEvents(&NumEvents, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Warning_1, "Cannot determine number of IBIS events!");
      CommonExit(status);
    }
  if(!NumEvents) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "No ISGRI events were found. Nothing to do");
      CommonExit(ISDC_OK);
    }
  if(detailedOutput) 
    RILstatus= RILlogMessage(NULL, Log_1, "Initial number of ISGRI events= %ld,", NumEvents);

  // Allocations 
  // Nota: We are working in energy for these shadowgrams thus it is                      
  // IsgriEnergy that is read. In all cases the data will be used as floating point values.
  if((IsgriY= (dal_byte *)calloc(NumEvents, sizeof(dal_byte)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for IsgriY.", ERR_ISGR_OSM_MEMORY_ALLOC);
  if((IsgriZ= (dal_byte *)calloc(NumEvents, sizeof(dal_byte)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for IsgriZ.", ERR_ISGR_OSM_MEMORY_ALLOC);
  if((IsgriEnergy= (dal_double *)calloc(NumEvents, sizeof(dal_double)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for IsgriEnergy.", ERR_ISGR_OSM_MEMORY_ALLOC);
  if((IsgriTime= (OBTime *)calloc(NumEvents, sizeof(OBTime)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for IsgriTime.", ERR_ISGR_OSM_MEMORY_ALLOC);
  if((IsgriRiseTime= (dal_byte *)calloc(NumEvents, sizeof(dal_byte)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for IsgriTime.", ERR_ISGR_OSM_MEMORY_ALLOC);
  if((IsgriSelectFlag= (dal_byte *)calloc(NumEvents, sizeof(dal_byte)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for IsgriTime.", ERR_ISGR_OSM_MEMORY_ALLOC);  
  if((ChannelMin= (dal_double *)calloc(NUM_PHA, sizeof(dal_double)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for IsgriEnergy.", ERR_ISGR_OSM_MEMORY_ALLOC);
  if((ChannelMax= (dal_double *)calloc(NUM_PHA, sizeof(dal_double)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for IsgriEnergy.", ERR_ISGR_OSM_MEMORY_ALLOC); 

  // Retrieve the detector coordinates 
  DALtype= DAL_BYTE;
  if((status= DAL3IBISgetEvents(ISGRI_Y, &DALtype, IsgriY, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "(1) DAL3IBISgetEvents status %d", status);
      CommonExit(status);
    }
  if((status= DAL3IBISgetEvents(ISGRI_Z, &DALtype, IsgriZ, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "(2) DAL3IBISgetEvents status %d", status);
      CommonExit(status);
    } 
  // Retrieve rise times for future selection 
  if((status= DAL3IBISgetEvents(RISE_TIME, &DALtype, IsgriRiseTime, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "(3) DAL3IBISgetEvents status %d", status);
      CommonExit(status);
    }
  if((status= DAL3IBISgetEvents(SELECT_FLAG, &DALtype, IsgriSelectFlag, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "(4) DAL3IBISgetEvents status %d", status);
      CommonExit(status);
    }
  // Retrieve the energy amplitude 
  DALtype=DAL_DOUBLE;
  if((status= DAL3IBISgetEvents(ISGRI_ENERGY, &DALtype, IsgriEnergy, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "(5) DAL3IBISgetEvents status %d", status);
      CommonExit(status);
    }
  // Retrieve event OBTimes 
  // (at this point we can construct the Time Index array) 
  DALtype= DAL3_OBT;
  if((status= DAL3IBISgetEvents(OB_TIME, &DALtype, IsgriTime, status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "(6) DAL3IBISgetEvents status %d", status);
      CommonExit(status);
    }
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                        SELECT  Events 
  //
  // We Exclude events in case:  a/ Time[event] is out of ISGRI GTIs, 
  //                             b/ Time[event] is out of MCE GTIs, 
  //                             c/ E[event]<LowThreshold criterium, 
  //                             d/ event is SELECT_FLAG Tagged.
  //
  // BUT We keep events from DEAD (LT=63) and OFF pixels (in REV or HK3): only their Efficiency is set 
  // to 0 (for diagnostic reasons). This is the same for later NOISY analysis (it only affects EFFI).   
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  int 
    IsgrY=0,
    IsgrZ=0,
    mce=  0,
    NumPixTAGGED= 0;
  long 
    NumSFexcluded=     0,
    NumLTexcluded=     0,
    NumMceGTIexcluded= 0,
    NumSelected=       0;
  double
    LTkeV= 0.,
    Alpha= 3.,
    Sigma= (FWHM0+FWHM1*OSMattribute.revol)/2.36;
    
  // Allocations
  long *TmpLongArray= NULL;
  if((TmpLongArray=(long *)calloc(NumEvents, sizeof(long)))==NULL) 
    Abort(NewGRP, "Error in allocating memory TmpLongArray.", ERR_ISGR_OSM_MEMORY_ALLOC);
  dal_int **SelectFlagMap= NULL;
  if( (SelectFlagMap= (dal_int**)calloc(ISGRI_SIZE,sizeof(dal_int*)))==NULL ) 
    Abort(NewGRP, "Error in allocating memory for SelectFlagMap.", ERR_ISGR_OSM_MEMORY_ALLOC);
  for(int i=0; i<ISGRI_SIZE; i++)
    if( (SelectFlagMap[i]= (dal_int*)calloc(ISGRI_SIZE,sizeof(dal_int)))==NULL ) 
      Abort(NewGRP, "Error in allocating memory for SelectFlagMap[].", ERR_ISGR_OSM_MEMORY_ALLOC);

  // Loop on all events
  for(long i=0;i<NumEvents;i++) 
    {
      IsgrY= IsgriY[i];
      IsgrZ= IsgriZ[i];
      mce=   YZtoModN(IsgrY,IsgrZ);
      LTkeV= LowThreshMap[IsgrY][IsgrZ];

      // If Event was not tagged noisy
      if(!IsgriSelectFlag[i]) 
	{
	  // IF EVENT has ( ENERGY>LowThreshold[pixel]-Alpha*Sigma )
	  if(IsgriEnergy[i]>=LTkeV-Alpha*Sigma)
	    for(long j=0;j<NumGTIperMCE[mce];j++)
	      {
		// IF EVENT arrived during a module GTI
		if( (GTIstartPerMCE[mce][j]<=IsgriTime[i]) && (IsgriTime[i]<=GTIstopPerMCE[mce][j]) ) 
		  {
		    TmpLongArray[NumSelected]=i;
		    NumSelected++;
		    NumEventsPerMCE[mce]++;
		    break;
		  }
		else if(j==NumGTIperMCE[mce]-1)
		  NumMceGTIexcluded++;
	      }
	  else
	    NumLTexcluded++;
	}
      // If Event was tagged noisy=> tag pixel NOISY!
      else
	{
	  if(!SelectFlagMap[IsgrY][IsgrZ])
	    NumPixTAGGED++;
	  SelectFlagMap[IsgrY][IsgrZ]= 1;
	  NumSFexcluded++;
	}	
    }
	
  // Selection results 
  RILstatus= RILlogMessage(NULL, Log_1, "MCE Good Times Interval Selection=> excluded %ld events.", NumMceGTIexcluded);
  RILstatus= RILlogMessage(NULL, Log_1, "Low Thresholds criteria Selection=> excluded %ld events.", NumLTexcluded);
  RILstatus= RILlogMessage(NULL, Log_1, "SELECT_FLAG criterium Selection  => excluded %ld events.", NumSFexcluded);
  RILstatus= RILlogMessage(NULL, Log_1, "                                 => excluded %d pixels.", NumPixTAGGED);
  RILstatus= RILlogMessage(NULL, Log_1, "Resulting Number of ISGRI events= %ld events.", NumSelected);
  if(detailedOutput) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "Events are distributed accordingly: ");
      for(int j=0;j<IBIS_NUM_BLOCK;j++)
	RILstatus= RILlogMessage(NULL, Log_1, "MCE%d : %ld viz. %.2f in this Module.",j,NumEventsPerMCE[j], 
				 100.*(float)NumEventsPerMCE[j]/(float)NumSelected);
    }
  
  // Re-allocations due to selections
  dal_byte        
    *TmpByte=   NULL;
  dal_double
    *TmpDouble= NULL;
  OBTime 
    *TmpOBTime= NULL;
  
  if((TmpByte=(unsigned char *)calloc(NumSelected, sizeof(unsigned char)))==NULL) 
    Abort(NewGRP, "Error in allocating memory TmpByte.", ERR_ISGR_OSM_MEMORY_ALLOC);
  if((TmpDouble=(double *)calloc(NumSelected, sizeof(double)))==NULL) 
    Abort(NewGRP, "Error in allocating memory TmpDouble.", ERR_ISGR_OSM_MEMORY_ALLOC);
  if((TmpOBTime=(OBTime *)calloc(NumSelected, sizeof(OBTime)))==NULL) 
    Abort(NewGRP, "Error in allocating memory TmpOBTime.", ERR_ISGR_OSM_MEMORY_ALLOC);
  
  for(long i=0; i<NumSelected; i++) 
    TmpByte[i]=IsgriY[TmpLongArray[i]];
  memcpy((void *)IsgriY, (const void *)TmpByte, NumSelected*sizeof(unsigned char));
  for(long i=0; i<NumSelected; i++) 
    TmpByte[i]=IsgriZ[TmpLongArray[i]];
  memcpy((void *)IsgriZ, (const void *)TmpByte, NumSelected*sizeof(unsigned char));
  for(long i=0; i<NumSelected; i++) 
    TmpByte[i]=IsgriRiseTime[TmpLongArray[i]];;
  memcpy((void *)IsgriRiseTime, (const void *)TmpByte, NumSelected*sizeof(unsigned char));
  for(long i=0; i<NumSelected; i++) 
    TmpOBTime[i]=IsgriTime[TmpLongArray[i]];;
  memcpy((void *)IsgriTime, (const void *)TmpOBTime, NumSelected*sizeof(OBTime));
  for(long i=0; i<NumSelected; i++) 
    TmpDouble[i]=IsgriEnergy[TmpLongArray[i]];;
  memcpy((void *)IsgriEnergy, (const void *)TmpDouble, NumSelected*sizeof(double));
  
  if(TmpByte) 
    { free(TmpByte); TmpByte= NULL; }
  if(TmpOBTime) 
    { free(TmpOBTime); TmpOBTime= NULL; }
  if(TmpDouble) 
    { free(TmpDouble); TmpDouble= NULL; }
  
  if(!NumSelected) 
    {
      RILstatus= RILlogMessage(NULL, Warning_1, "Selection criteria have reduced the number of ISGRI events");
      RILstatus= RILlogMessage(NULL, Warning_1, "from %ld to zero.", NumEvents);
      RILstatus= RILlogMessage(NULL, Warning_1, "Terminating process with ISDC_OK.");
      CommonExit(ISDC_OK);
    }
  NumEvents= NumSelected;
  if((IsgriY= (dal_byte *)realloc(IsgriY, NumEvents*sizeof(dal_byte)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for IsgriY.", ERR_ISGR_OSM_MEMORY_ALLOC);
  if((IsgriZ= (dal_byte *)realloc(IsgriZ, NumEvents*sizeof(dal_byte)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for IsgriZ.", ERR_ISGR_OSM_MEMORY_ALLOC);
  if((IsgriEnergy= (dal_double *)realloc(IsgriEnergy, NumEvents*sizeof(dal_double)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for IsgriEnergy.", ERR_ISGR_OSM_MEMORY_ALLOC);
  if((IsgriTime= (OBTime *)realloc(IsgriTime, NumEvents*sizeof(OBTime)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for IsgriTime.", ERR_ISGR_OSM_MEMORY_ALLOC);
  
  if(IsgriRiseTime)
    { free(IsgriRiseTime); IsgriRiseTime= NULL; }
  if(IsgriSelectFlag) 
    { free(IsgriSelectFlag); IsgriSelectFlag= NULL; }
  if(TmpLongArray) 
    { free(TmpLongArray); TmpLongArray=NULL; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                   TIME AND ENERGY ALLOCATIONS
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  RILstatus= RILlogMessage(NULL,Log_1,"-------- TIME BINNING information --------");

  int 
    NumTime= 1;
  dal_element
    **DOLimages=        NULL,                         
    **DOLeffiImages=    NULL; 

  // Time binning parameters 
  if(TimeLen<0) 
    TimeLen=(fabs(TimeLen)<SCWduration) ? SCWduration : fabs(TimeLen);
  if(TimeLen<ZERO) 
    Abort(NewGRP, "Time bin value is NULL!", status);
  else
    {
      RILstatus= RILlogMessage(NULL, Log_1, "Data range(OBT) = [ %llu, %llu ]", IsgriTime[0], IsgriTime[NumEvents-1]);
      RILstatus= RILlogMessage(NULL, Log_1, "SCW duration    = %.3fsec", SCWduration);
      RILstatus= RILlogMessage(NULL, Log_1, "Time Bin Length = %.3fsec", TimeLen);
      NumTime= (int)floor(SCWduration/TimeLen) + 1;
      if(!NumTime)
	{ 
	  RILstatus= RILlogMessage(NULL, Warning_1, "Number of time bins should not be NULL!");
	  RILstatus= RILlogMessage(NULL, Warning_1, "Reverting to TimeBin= SCW and NumBins= 1.");
	  TimeLen= SCWduration;
	  NumTime= 1;
	}
    }

  // Shadowgrams elements allocations             
  if((DOLimages=(dal_element **)calloc(NumImaBin*NumTime, sizeof(dal_element *)))==NULL)
    Abort(NewGRP, "Error in allocating memory for DOLimages files : status = %d", status);
  if((DOLeffiImages=(dal_element **)calloc(NumImaBin*NumTime, sizeof(dal_element *)))==NULL)
    Abort(NewGRP, "Error in allocating memory for DOLlightcurves files : status = %d", status);

  // Raw Shadowgrams file creation 
  char TemplateName[DAL_FILE_NAME_STRING]= ""; 
  strcpy(TemplateName, ISGR_SHD_TPL);
  if((status= MkOutputFiles(idxRAWshad, DOLimages, NumImaBin, NumTime, TemplateName, status))!=ISDC_OK)
    Abort(NewGRP, "Error: unable to create output SHAD structure (status %d)", status); 

  // Eff Shadowgrams file creation 
  strcpy(TemplateName, ISGR_SHD_EFFI_TPL);
  if((status= MkOutputFiles(idxEFFshad, DOLeffiImages, NumImaBin, NumTime, TemplateName, status))!=ISDC_OK)
    Abort(NewGRP, "Error: unable to create output EFFI structure (status %d)", status);

  // Time intervals Allocations
  RILstatus= RILlogMessage(NULL, Log_1, "=> Divides SCW into %d Time Bins:", NumTime);
  long *TimeIdx= NULL; 
  if((TimeIdx=(long *)calloc(NumTime, sizeof(long)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for time indices.", status);
  OBTime *ImageOBTStart= NULL;
  if((ImageOBTStart=(OBTime *)calloc(NumTime, sizeof(OBTime)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for ImageOBTStart.", status);
  OBTime *ImageOBTEnd= NULL;
  if((ImageOBTEnd=(OBTime *)calloc(NumTime, sizeof(OBTime)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for ImageOBTEnd.", status);
  bool *BTIflag= NULL;                           
  if((BTIflag= (bool*)calloc(NumTime,sizeof(bool)))==NULL)
    Abort(NewGRP, "Error in allocating memory for BTIflag.", status);

  // Define the NumTime initial OBT_intervals [ start[BinT]; end[BinT] ]     
  long j= 0;
  OBTEnd= OBTStart;
  for(int BinT=0;BinT<NumTime;BinT++) 
    {
      TimeIdx[BinT]= -1;
      // At the start of each iteration, start[BinT]=end[BinT-1] and end[BinT]=start[BinT]+TimeLen 
      ImageOBTStart[BinT]= OBTEnd;
      if((status= DAL3GENskipOBT(OBTEnd, (OBTime)TimeLen*DAL3_OBT_SECOND, &OBTEnd, status))!=ISDC_OK)
	{
	  RILstatus= RILlogMessage(NULL, Error_1, " Fatal error in computing OBTime + dt, status = %d.", status);
	  CommonExit(status);
	}
       ImageOBTEnd[BinT]= OBTEnd;
     
      // if there are events left to process (if all NumEvents have not yet been considered before) 
      if(j<NumEvents) 
	{
	  // if end[BinT] exceeds end[scw] 
	  if(OBTEnd>EndTime) 
	    OBTEnd= EndTime;
	  
	  // find TimeIdx[BinT] = rank of last event of OBT_Interval BinT; its date is IsgriTime[TimeIdx[BinT]] 
	  BinarySearch(&OBTEnd, IsgriTime, j, NumEvents-1, TimeIdx+BinT, TREE_LLONG);
	  
	  // if rank of last event found is lower or equal than the one of first event of OBT_interval 
	  if(TimeIdx[BinT]<=j) 
	    {
	      BTIflag[BinT]= TRUE;
	      RILstatus= RILlogMessage(NULL, Log_1, " No event found in Time_Bin[%d]= [ %llu, %llu ] (OBT)", 
				       BinT+1, ImageOBTStart[BinT], ImageOBTEnd[BinT]);
	    }
	  
	  // if there are events in the OBT_Interval 
	  else 
	    {
	      ImageOBTStart[BinT]= IsgriTime[j];
	      ImageOBTEnd[BinT]=   IsgriTime[TimeIdx[BinT]]; 
	      j= TimeIdx[BinT]+1;
	      RILstatus= RILlogMessage(NULL, Log_1, "Time_Bin[%d]= [ %llu, %llu ] (OBT)", 
				       BinT+1, ImageOBTStart[BinT], ImageOBTEnd[BinT]);
	    }
	}
      
      // if first event of OBT_Interval belongs to another SCW, don't  
      else 
	{
	  RILstatus= RILlogMessage(NULL, Warning_1, "[ Start, End ] = [ %llu, %llu ] ==> [ %ld, %llu ]", 
				   ImageOBTStart[BinT], ImageOBTEnd[BinT], 
				   TimeIdx[BinT], IsgriTime[TimeIdx[BinT]]);
	  RILstatus= RILlogMessage(NULL, Warning_1, "These times are NOT within the science window.");
	  BTIflag[BinT]= TRUE;
	}
    }
	  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //                                   FIND REMAINING NOISY PIXELS                                 //
  /////////////////////////////////////////////////////////////////////////////////////////////////// 
      
  // Allocate Empty Spec Noisy Map
  dal_double **SpecNoisyMap= NULL;
  if((SpecNoisyMap= (dal_double **)calloc(ISGRI_SIZE,sizeof(dal_double*)))==NULL) 
    Abort(NewGRP,"Error in allocating memory for Noisy Flag map[].",ERR_ISGR_OSM_MEMORY_ALLOC);
  for(int y=0;y<ISGRI_SIZE;y++) 
    if((SpecNoisyMap[y]= (dal_double*)calloc(ISGRI_SIZE,sizeof(dal_double)))==NULL)
      Abort(NewGRP,"Error in allocating memory for Noisy Flag map[][]",ERR_ISGR_OSM_MEMORY_ALLOC);
  
  // We here use a fitting algorithm based on pixels'spectra over the whole SCW 	   
  if(NoisyDetFlag)
    {
      RILstatus= RILlogMessage(NULL,Log_1,"-------- Find Remaining NOISY PIXELS --------");
      int NumSpecNoisy= 0;
      if((status= SpecNoisyPixels(IsgriY, IsgriZ, IsgriEnergy, NumEvents, StartTime, EndTime,
				  GTIstart, GTIstop, NumGTI, GTIstartPerMCE, GTIstopPerMCE, NumGTIperMCE, 
				  DeadTimes, OBTdeadTimes, NumDeadTimes, 
				  ONpixelsREVmap, ONpixelsHK3map, SelectFlagMap, LowThreshMap, OSMattribute.revol,
				  SpecNoisyMap, &NumSpecNoisy, detailedOutput)) != ISDC_OK)
	{
	  RILstatus= RILlogMessage(NULL,Error_1,"Error while extracting noisy pixels, status= %d.", status);
	  CommonExit(status);
	}
      RILstatus= RILlogMessage(NULL,Log_1,"=> Found %d remaining Noisy Pixels in this SCW.",NumSpecNoisy);
    }
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //                                   BUILD RAW/EFF SHADOWGRAMS                                   //
  ///////////////////////////////////////////////////////////////////////////////////////////////////   
  
  RILstatus= RILlogMessage(NULL,Log_1,"-------- BUILD RAW/EFF SHADOWGRAMS --------");

  // Allocate generic ISGRI Time Efficiency map  
  dal_double **TimeEffMap= NULL;
  if((TimeEffMap= (dal_double **)calloc(ISGRI_SIZE,sizeof(dal_double*)))==NULL) 
    Abort(NewGRP,"Error in allocating memory for Time Efficiency map[].",ERR_ISGR_OSM_MEMORY_ALLOC);
  for(int y=0;y<ISGRI_SIZE;y++) 
    if((TimeEffMap[y]= (dal_double*)calloc(ISGRI_SIZE,sizeof(dal_double)))==NULL)
      Abort(NewGRP,"Error in allocating memory for Time Efficiency map[][]",ERR_ISGR_OSM_MEMORY_ALLOC);

  // RAW and EFF Shadowgrams allocations
  dal_double  ***IsgriSHD= NULL;
  if((IsgriSHD=(dal_double ***)calloc(HowMany, sizeof(dal_double **)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for image", status);
  for(int i=0; i<HowMany; i++) 
    if((IsgriSHD[i]=(dal_double **)calloc(ISGRI_SIZE, sizeof(dal_double *)))==NULL) 
      Abort(NewGRP, "Error in allocating memory for image : i = %d", i);
  for(int i=0; i<HowMany; i++)
    for(int j=0; j<ISGRI_SIZE; j++) 
      if((IsgriSHD[i][j]=(dal_double *)calloc(ISGRI_SIZE, sizeof(dal_double)))==NULL) 
	Abort(NewGRP, "Error in allocating memory for image : i = %d", i);
  dal_double ***IsgriEffSHD= NULL;                                         
  if((IsgriEffSHD=(dal_double ***)calloc(HowMany, sizeof(dal_double **)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for image efficiency map.", status);
  for(int i=0; i<HowMany; i++) 
    if((IsgriEffSHD[i]=(dal_double **)calloc(ISGRI_SIZE, sizeof(dal_double *)))==NULL) 
      Abort(NewGRP, "Error in allocating memory for image efficiency map : i = %d.", i);
  for(int i=0; i<HowMany; i++)
    for(int j=0; j<ISGRI_SIZE; j++) 
      if((IsgriEffSHD[i][j]=(dal_double *)calloc(ISGRI_SIZE, sizeof(dal_double)))==NULL) 
	Abort(NewGRP, "Error in allocating memory for image efficiency map : i = %d", i);

  // Shadowgrams details allocations                                        
  long *NumOuts= NULL;                                                 
  if((NumOuts= (long*)calloc(HowMany, sizeof(long)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for NumOuts.", status);
  dal_double *MeanEff= NULL;
  if((MeanEff= (dal_double*)calloc(HowMany, sizeof(dal_double)))==NULL) 
    Abort(NewGRP, "Error in allocating memory for MeanEff.", status);
  
  // Loop on each Time Bin   
  long 
    IdxStart= 0; 
  dal_double 
    DeadC= 0.,
    OnTime= 0.,
    MeanTimeEff= 0.; 
  for(int BinT=0;BinT<NumTime;BinT++)
    {
      RILstatus= RILlogMessage(NULL,Log_1," ");
      RILstatus= RILlogMessage(NULL,Log_1,"WORKING in TimeBin[%d]= [ %llu, %llu ]",
			       BinT+1, ImageOBTStart[BinT], ImageOBTEnd[BinT]);
      
      // If Time Intervall is OK + if there are events inside
      if(!BTIflag[BinT]) 
	{
	  // ============================== Number of processed events  ============================== 
	  long NumEntries= *(TimeIdx+BinT)-IdxStart;
	  RILstatus= RILlogMessage(NULL, Log_1, "Number of processed events= %ld (indices= %ld to %ld).", 
				   NumEntries, IdxStart, TimeIdx[BinT]);

	  // ================= Compute ISGRI ONTIME, DEADC and Time Efficiencies =====================   
	  if((status= ComputeTimeEffMap(IsgriTime[IdxStart], IsgriTime[TimeIdx[BinT]],
					GTIstart, GTIstop, NumGTI, GTIstartPerMCE, GTIstopPerMCE, NumGTIperMCE,
					DeadTimes, OBTdeadTimes, NumDeadTimes,
					ONpixelsREVmap, ONpixelsHK3map, SelectFlagMap, TimeEffMap, 
					&MeanTimeEff, &OnTime, &DeadC, detailedOutput)) != ISDC_OK)
	    {
	      RILstatus= RILlogMessage(NULL,Error_1,"Error while Computing Time Efficiencies, status= %d.", status);
	      CommonExit(status);
	    }
	  RILstatus= RILlogMessage(NULL, Log_1, "=> Average Time Efficiency over the TimeBin= %.3f%c.", 
				   100*MeanTimeEff,'%');

	  // ======================================== IMAGES =========================================
	  // Compute the image from the event data and the efficiency map: one image/Eband/TimeBin.       
	  int 
	    m=0,
	    Iteration=0;	      
	  while(Iteration<NumImaBin) 
	    {

	      // EnergyBin info 
	      m= (NumImaBin-Iteration)>HowMany ? HowMany : NumImaBin-Iteration;
	      RILstatus= RILlogMessage(NULL,Log_1," ");
	      RILstatus= RILlogMessage(NULL,Log_1,"WORKING in EnergyRange [%.2f,%.2f[keV (includes %d Ebins)",
				       EnergyBounds[Iteration][0],EnergyBounds[Iteration+m-1][1],m);
	      
	      // Compute IsgriSHD 
	      if((status= MkImage(&IsgriY[IdxStart], &IsgriZ[IdxStart], &IsgriEnergy[IdxStart], NumEntries, 
				  &EnergyBounds[Iteration], m, IsgriSHD, NumOuts)) != ISDC_OK)
		{
		  RILstatus= RILlogMessage(NULL,Error_1,"Error Filling RAW shadowgram, status= %d.", status);
		  CommonExit(status);
		}	     
	      
	      // Compute IsgriEffSHD
	      if((status= MkeffImage(TimeEffMap, LowThreshMap, SpecNoisyMap, 
				     &EnergyBounds[Iteration], m, OSMattribute.revol, IsgriEffSHD, MeanEff)) != ISDC_OK)
		{
		  RILstatus= RILlogMessage(NULL,Error_1,"Error Filling EFFI shadowgram, status= %d.", status);
		  CommonExit(status);
		}
	      
	      // Shadowgram details
	      for(int bin=0;bin<m;bin++)
		RILstatus= RILlogMessage(NULL, Log_1, "Eband [%.2f,%.2f[: NumEvents= %6ld, MeanEfficiency= %.3f%c.",
					 EnergyBounds[Iteration+bin][0], EnergyBounds[Iteration+bin][1],
					 NumOuts[bin],100*MeanEff[bin],'%');
	      
	      // Write SHADOWGRAMS and related KEYWORDS to output files
	      if((status= WriteImage(NewGRP, &DOLimages[Iteration+BinT*NumImaBin], &DOLeffiImages[Iteration+BinT*NumImaBin], 
				     IsgriSHD, IsgriEffSHD, StartTime, EndTime, IsgriTime[IdxStart], IsgriTime[TimeIdx[BinT]], 
				     &EnergyBounds[Iteration], m,(char*)"ENERGY",(char*)"", OnTime, DeadC, MeanTimeEff)) != ISDC_OK)
		{
		  RILstatus= RILlogMessage(NULL,Error_1,"Error Writing shadowgrams, status= %d.", status);
		  CommonExit(status);
		}
	      
	      // Clean matrices 
	      for(int l=0; l<m; l++) 
		for(int y=0; y<ISGRI_SIZE; y++)
		  for(int z=0; z<ISGRI_SIZE; z++) 
		    {
		      IsgriSHD[l][y][z]=0.0;
		      IsgriEffSHD[l][y][z]=0.0;
		    }
	      
	      // Increment EnergyRange 
	      Iteration+=m;
	      
	    }// end Energy block (while LOOP)	
	  
	} // If TimeBin is OK
      
      // If Time interval is invalid, write empty shadowgrams
      else 
	{
	  RILstatus= RILlogMessage(NULL, Log_1, "No event in this Time interval.");	  
	  int Iteration=0;	      
	  while(Iteration<NumImaBin) 
	    {	 
	      int m= (NumImaBin-Iteration)>HowMany ? HowMany : NumImaBin-Iteration;
	      if((status= WriteImage(NewGRP, &DOLimages[Iteration+BinT*NumImaBin], &DOLeffiImages[Iteration+BinT*NumImaBin], 
				     IsgriSHD, IsgriEffSHD, StartTime, EndTime, ImageOBTStart[BinT], ImageOBTEnd[BinT], 
				     &EnergyBounds[Iteration], m,(char*)"ENERGY",(char*)"", OnTime, DeadC, MeanTimeEff)) != ISDC_OK)
		{
		  RILstatus= RILlogMessage(NULL,Error_1,"Error Writing shadowgrams, status= %d.", status);
		  CommonExit(status);
		}
	      Iteration+=m;
	    }
	}
      
      // Increment TimeBin 
      IdxStart= (TimeIdx[BinT]+1) <= NumEvents ? (TimeIdx[BinT]+1) : NumEvents;
    }
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////
  //                            UPDATE shadowgrams and indexes
  /////////////////////////////////////////////////////////////////////////////////////////////
 
  RILstatus= RILlogMessage(NULL,Log_1," ");
  RILstatus= RILlogMessage(NULL,Log_1,"-------- UPDATE shadowgrams and EXIT --------");
  
  // Fill Images Attributes
  unsigned long AttributeCode= 0;
  AttributeCode= REVOL_NUM|SWID_NUM|SWTYPE_NUM|SWBOUND_NUM|OUTPUT_NUM|RISE_MIN_NUM|RISE_MAX_NUM;
  for(int k=0; k<NumImaBin*NumTime; k++)
    {
      status= WriteAttributes(DOLimages[k], AttributeCode, OSMattribute, status);
      status= WriteAttributes(DOLeffiImages[k], AttributeCode, OSMattribute, status);
    }

  // Fill Index Attributes
  // SPR-4146: don't fill OUTPUT_NUM|TFIRST_NUM|TLAST_NUM|RISE_MIN_NUM|RISE_MAX_NUM (non constant values)
  AttributeCode= 0;
  AttributeCode= REVOL_NUM|SWID_NUM|SWTYPE_NUM|SWBOUND_NUM|TSTART_NUM|TSTOP_NUM;
  status= WriteAttributes(idxRAWshad, AttributeCode, OSMattribute, status);
  status= WriteAttributes(idxEFFshad, AttributeCode, OSMattribute, status);
  status= DAL3GENattributeCopy(NewGRP,idxRAWshad,"OBTSTART,OBTEND",status);// SPR 2806 
  status= DAL3GENattributeCopy(NewGRP,idxEFFshad,"OBTSTART,OBTEND",status);// SPR 2806 
  
  // Update Indexes
  if((status= UpdateIndex(idxRAWshad, DOLimages, NumImaBin*NumTime, detailedOutput, status))!=ISDC_OK)
    Abort(NewGRP, "Error in updating index.", ERR_ISGR_OSM_OUTPUT_INDEX_CREATION);
  if((status= UpdateIndex(idxEFFshad, DOLeffiImages, NumImaBin*NumTime, 
			  detailedOutput, status))!=ISDC_OK) 
    Abort(NewGRP, "Error in updating index.", ERR_ISGR_OSM_OUTPUT_INDEX_CREATION);
  RILstatus= RILlogMessage(NULL, Log_1, "Finished producing RAW/EFF shadowgrams (status %d).", status);
  
  if((status= CommonStampObject(idxRAWshad, "ISGRI Detector Images.", status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Error_1, "Writing in time stamping index : status = %d.", status);
      CommonExit(status);
    }
  if((status= CommonStampObject(idxEFFshad, "ISGRI Detector Images Efficiencies.", status))!=ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Error_1, "Writing in time stamping index : status = %d.", status);
      CommonExit(status);
    }
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////
  //                                     DE-ALLOCATIONS
  /////////////////////////////////////////////////////////////////////////////////////////////
  
  if(detailedOutput) 
    RILstatus= RILlogMessage(NULL, Log_1, "Freeing arrays.");
  
  // Free log(CHI2)_method obtained Noisy Maps
  for(int y=0;y<ISGRI_SIZE;y++) 
    if(SpecNoisyMap[y])
      { free(SpecNoisyMap[y]); SpecNoisyMap[y]= NULL; }
  if(SpecNoisyMap)
    { free(SpecNoisyMap); SpecNoisyMap= NULL; }
  
  // Free SELECT_FLAG_method obtained Noisy Maps
  for(int i=0; i<ISGRI_SIZE; i++)
    if(SelectFlagMap[i])
      { free(SelectFlagMap[i]); SelectFlagMap[i]= NULL; }
  if(SelectFlagMap)
    { free(SelectFlagMap); SelectFlagMap= NULL; }
  
  // Free time indices 
  if(ImageOBTEnd)
    { free(ImageOBTEnd); ImageOBTEnd= NULL; }
  if(ImageOBTStart)
    { free(ImageOBTStart); ImageOBTStart= NULL; }
  if(TimeIdx)
    { free(TimeIdx); TimeIdx= NULL; }
  if(BTIflag)
    { free(BTIflag); BTIflag= NULL; }
  
  // Free allocated arrays for the channels / energy bounds 
  if(ChannelMax)
    { free(ChannelMax); ChannelMax= NULL; }
  if(ChannelMin)
    { free(ChannelMin); ChannelMin= NULL; }
  
  // Free allocated arrays for the energies 
  if(IsgriEnergy)
    { free(IsgriEnergy); IsgriEnergy= NULL; }
  for(int i=0; i<NumImaBin; i++) 
    if(EnergyBounds[i])
      { free(EnergyBounds[i]); EnergyBounds[i]= NULL; }
  if(EnergyBounds)
    { free(EnergyBounds); EnergyBounds= NULL; }
  
  // Free arrays for detector coordinates 
  if(IsgriZ)
    { free(IsgriZ); IsgriZ= NULL; }
  if(IsgriY)
    { free(IsgriY); IsgriY= NULL; }
  if(IsgriRiseTime)
    { free(IsgriRiseTime); IsgriRiseTime= NULL; }
  
  // Deallocate pixel efficiency, ON/OFF pixels, low threshold images  
  for(int i=0; i<ISGRI_SIZE; i++)
    {
      if(TimeEffMap[i])
	{ free(TimeEffMap[i]); TimeEffMap[i]= NULL; }
      if(ONpixelsREVmap[i])
	{ free(ONpixelsREVmap[i]); ONpixelsREVmap[i]= NULL; }
      if(ONpixelsHK3map[i])
	{ free(ONpixelsHK3map[i]); ONpixelsHK3map[i]= NULL; }
      if(LowThreshMap[i])
	{ free(LowThreshMap[i]); LowThreshMap[i]= NULL; }
    }
  if(TimeEffMap)
    { free(TimeEffMap); TimeEffMap= NULL; }
  if(ONpixelsREVmap)
    { free(ONpixelsREVmap); ONpixelsREVmap= NULL; }
  if(ONpixelsHK3map)
    { free(ONpixelsHK3map); ONpixelsHK3map= NULL; }
  if(LowThreshMap)
    { free(LowThreshMap); LowThreshMap= NULL; }
  
  // Deallocate arrays for Isgri shadowgrams 
  for(int i=0; i<HowMany; i++)
    for(int j=0; j<ISGRI_SIZE; j++)
      {
	if(IsgriEffSHD[i][j])
	  { free(IsgriEffSHD[i][j]); IsgriEffSHD[i][j]= NULL; }
	if(IsgriSHD[i][j])
	  { free(IsgriSHD[i][j]); IsgriSHD[i][j]= NULL; }
      }
  for(int i=0; i<HowMany; i++)
    {
      if(IsgriEffSHD[i])
	{ free(IsgriEffSHD[i]); IsgriEffSHD[i]= NULL; }
      if(IsgriSHD[i])
	{ free(IsgriSHD[i]); IsgriSHD[i]= NULL; }
    }
  if(IsgriEffSHD)
    { free(IsgriEffSHD); IsgriEffSHD= NULL; }
  if(IsgriSHD)
    { free(IsgriSHD); IsgriSHD= NULL; }
  fflush(NULL);                            

  // Deallocate the deadtime arrays 
  for(int i=0; i<NumDeadTimes; i++) 
    if(DeadTimes[i])
      { free(DeadTimes[i]); DeadTimes[i]= NULL; }
  if(DeadTimes)
    { free(DeadTimes); DeadTimes= NULL; }
  if(OBTdeadTimes)
    { free(OBTdeadTimes); OBTdeadTimes= NULL; }
  
  // Deallocate the good time interval arrays
  for(int i=0;i<IBIS_NUM_BLOCK;i++)
    {
      if(GTIstopPerMCE[i]) 
	{ free(GTIstopPerMCE[i]); GTIstopPerMCE[i]= NULL; }
      if(GTIstartPerMCE[i]) 
	{ free(GTIstartPerMCE[i]); GTIstartPerMCE[i]= NULL; }
    }
  if(GTIstopPerMCE) 
    { free(GTIstopPerMCE); GTIstopPerMCE= NULL; }
  if(GTIstartPerMCE) 
    { free(GTIstartPerMCE); GTIstartPerMCE= NULL; }
  if(NumGTIperMCE) 
    { free(NumGTIperMCE); NumGTIperMCE= NULL; }
  if(NumEventsPerMCE) 
    { free(NumEventsPerMCE); NumEventsPerMCE= NULL; } 
  if(GTIstop)
    { free(GTIstop); GTIstop= NULL; }
  if(GTIstart)
    { free(GTIstart); GTIstart= NULL; }
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                           COMMON  EXIT
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  
  Abort(NewGRP, "Program terminated normally", status);
}
