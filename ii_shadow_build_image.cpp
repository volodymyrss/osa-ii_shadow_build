#include "ii_shadow_build_f.h"


/////////////////////////////////////////////////////////////////////////////////////////////
//                    Algorithm to determine remaining noisy pixels.
//                         The flag is 0 in case pixel is OK; 
//       In case it is NOISY, it is its distance to "mean state", in sigma unit.
/////////////////////////////////////////////////////////////////////////////////////////////

int SpecNoisyPixels(dal_byte      *IsgrY,           // Table of Y coordinate of events
		    dal_byte      *IsgrZ,           // Table of Z coordinate of events 
		    dal_double    *IsgrE,           // Table of Energy of events
		    long          NumEvents,        // Number of ISGRI events			  
		    OBTime        SCWstart,         // SCW start
		    OBTime        SCWend,           // SCW end
		    OBTime        *GTIstart,        // Good time intervals, starting times               
		    OBTime        *GTIstop,         // Good time intervals, finishing times              
		    int           NumGTI,           // Number of good time intervals                     
		    OBTime        **GTIstartPerMCE, // Good time intervals per module, starting times    
		    OBTime        **GTIstopPerMCE,  // Good time intervals per module, stopping times    
		    int           *NumGTIperMCE,    // Number of good time intervals  per module 
		    dal_float     **DeadTimes,      // DeadTimes for each module, over the SCW
		    OBTime        *OBTdeadTimes,    // ISGRI DeadTimes dates over the SCW
		    long          NumDeadTimes,     // Number ISGRI DeadTimes over the SCW
		    dal_int       **ONpixelsREVmap, // Map of Pixels Status for this REV (REV context)
		    dal_int       **ONpixelsHK3map, // Map of Pixels Status for this SCW (HK3 data)
		    dal_int       **SelectFlagMap,  // Map of Noisy Pixels, as detected by SELECT_FLAG method
		    dal_double    **LowThreshMap,   // Map of Low Thresholds (keV)
		    int           revol_scw,        // Revolution number
		    dal_double    **SpecNoisyMap,   // Output: Remaining Noisy pixels map
		    int           *NumSpecNoisy,    // Output: Number of found noisy pixels
		    unsigned char detailedOutput,    // Detailed output
            ISGRI_efficiency_struct *ptr_ISGRI_efficiency
            )  
{
  int 
    status=     ISDC_OK,
    RILstatus=  ISDC_OK,
    NumBinEner= 45,
    NumBinChi2= 1000;
  dal_double 
    **TimeEffMap=  NULL,
    ***EffMap=     NULL,
    *BinEner=      NULL,
    *SpecMean=     NULL,
    *BinChi2=      NULL,
    *DistChi2=     NULL,
    **Chi2Map=     NULL,  
    ***SpecPix=    NULL,
    criteria= 2;

  // ==================================================================
  //                   Allocate generic modules arrays 
  // ==================================================================

  if((TimeEffMap= (dal_double**)calloc(ISGRI_SIZE,sizeof(dal_double*)))==NULL) 
    {
      RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : TimeEffMap.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if((EffMap= (dal_double***)calloc(ISGRI_SIZE,sizeof(dal_double**)))==NULL) 
    {
      RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : EffMap.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if((BinEner= (dal_double*)calloc(NumBinEner+1, sizeof(dal_double)))==NULL) 
    {
      RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : BinEner.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if((SpecMean= (dal_double*)calloc(NumBinEner, sizeof(dal_double)))==NULL) 
    {
      RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : SpecMean.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if((BinChi2= (dal_double*)calloc(NumBinChi2+1, sizeof(dal_double)))==NULL) 
    {
      RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : BinChi2.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if((DistChi2= (dal_double*)calloc(NumBinChi2, sizeof(dal_double)))==NULL) 
    {
      RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : DistChi2.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if((Chi2Map= (dal_double**)calloc(ISGRI_SIZE, sizeof(dal_double*)))==NULL) 
    {
      RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : Chi2Map.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if((SpecPix= (dal_double***)calloc(ISGRI_SIZE, sizeof(dal_double**)))==NULL) 
    {
      RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : SpecPix.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }

  for(int y=0;y<ISGRI_SIZE;y++)
    {
      if((TimeEffMap[y]= (dal_double*)calloc(ISGRI_SIZE,sizeof(dal_double)))==NULL) 
	{
	  RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : TimeEffMap[].");
	  return ERR_ISGR_OSM_MEMORY_ALLOC;
	}
      if((EffMap[y]= (dal_double**)calloc(ISGRI_SIZE,sizeof(dal_double*)))==NULL) 
	{
	  RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : EffMap[].");
	  return ERR_ISGR_OSM_MEMORY_ALLOC;
	}
      if((Chi2Map[y]= (dal_double*)calloc(ISGRI_SIZE,sizeof(dal_double)))==NULL) 
	{
	  RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : Chi2Map[].");
	  return ERR_ISGR_OSM_MEMORY_ALLOC;
	}
      if((SpecPix[y]= (dal_double**)calloc(ISGRI_SIZE,sizeof(dal_double*)))==NULL) 
	{
	  RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : SpecPix[].");
	  return ERR_ISGR_OSM_MEMORY_ALLOC;
	}
      for(int z=0;z<ISGRI_SIZE;z++)
	{
	  if((EffMap[y][z]= (dal_double*)calloc(NumBinEner,sizeof(dal_double)))==NULL) 
	    {
	      RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : EffMap[][].");
	      return ERR_ISGR_OSM_MEMORY_ALLOC;
	    }
	  if((SpecPix[y][z]= (dal_double*)calloc(NumBinEner,sizeof(dal_double)))==NULL) 
	    {
	      RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : SpecPix[][].");
	      return ERR_ISGR_OSM_MEMORY_ALLOC;
	    }
	}
    }

  
  // ==================================================================
  //                       Binnings definitions
  // ================================================================== 
      
  // FWHM Binning definition: compute full energy range
  dal_double DeltaE=0.;
  for(int bin=0;bin<NumBinEner+1;bin++)
    {
      if((ceil(0.054*bin)-0.054*bin)>=0.5) 
	DeltaE += 0.5*(1+floor(0.054*bin));
      else 
	DeltaE += 0.5*(1+ceil(0.054*bin));
      BinEner[bin] = 12.0+DeltaE;
    }

  // Chi2 values binning: range [1;4]
  for(int bin=0;bin<=NumBinChi2;bin++)
    BinChi2[bin] = 1.0 + 3.0*bin/NumBinChi2;

  
  // ==================================================================
  //      Compute Efficiency of ISGRI pixels over the whole SCW
  // ==================================================================

  RILlogMessage(NULL,Log_1,"Computing pixels' Efficiencies over the SCW and on 256 bands...");
  
  // Compute Time Efficiency Map
  dal_double tmpD1=0., tmpD2=0., tmpD3=0.;
  if((status= ComputeTimeEffMap(SCWstart, SCWend, 
				GTIstart, GTIstop, NumGTI, GTIstartPerMCE, GTIstopPerMCE, NumGTIperMCE, 
				DeadTimes, OBTdeadTimes, NumDeadTimes, 
				ONpixelsREVmap, ONpixelsHK3map, SelectFlagMap, 
				TimeEffMap, &tmpD1, &tmpD2, &tmpD3, 0)) != ISDC_OK)
    {
      RILstatus= RILlogMessage(NULL,Error_1,"Error while Computing Time Efficiency, status= %d.", status);
      CommonExit(status);
    }
  
  // Compute Total Efficiency= Energetic Efficiency * Time Efficiency:
  // Use analytical approximation of Energetic Efficiency
  // ( Nota: bins are narrow => don't need to weight with ARF(E)*pow(E,-2) )
  for(int y=0;y<ISGRI_SIZE;y++) 
    for(int z=0;z<ISGRI_SIZE;z++)
      if(TimeEffMap[y][z]>ZERO)
	for(int bin=0;bin<NumBinEner;bin++) {

	  EffMap[y][z][bin]= TimeEffMap[y][z] * LTfunction((BinEner[bin+1]+BinEner[bin])/2,y,z,ptr_ISGRI_efficiency);
    }
  
  // ==================================================================
  //                Compute mean and pixels Spectra  
  // ==================================================================

  RILlogMessage(NULL,Log_1,"Computing pixels' Spectra...");

  for(long j=0;j<NumEvents;j++)
    // If pixel is valid
    if(TimeEffMap[IsgrY[j]][IsgrZ[j]]>ZERO)
      for(int bin=0;bin<NumBinEner;bin++)
	if(IsgrE[j]>=BinEner[bin] && IsgrE[j]<BinEner[bin+1])
	  {
	    SpecPix[IsgrY[j]][IsgrZ[j]][bin]++;
	    break; // We can use it because Ebands are adjacent
	  }
      
  // Normalize spectra
  int NumPixOn;
  for(int bin=0;bin<NumBinEner;bin++)
    {
      NumPixOn=0;
      for(int y=0;y<ISGRI_SIZE;y++) 
	for(int z=0;z<ISGRI_SIZE;z++) 
	  if(EffMap[y][z][bin]>ZERO)
	    {
	      SpecPix[y][z][bin]/= EffMap[y][z][bin];
	      SpecMean[bin]+= SpecPix[y][z][bin];
	      NumPixOn++; 
	    }
      if(NumPixOn>0)
	SpecMean[bin]/= NumPixOn;
    }

      
  // ==================================================================
  //                       Compute Chi2Map 
  // ==================================================================

  if(detailedOutput)
    RILlogMessage(NULL,Log_1,"Computing log10(Chi2(spectra))...");

  for(int y=0;y<ISGRI_SIZE;y++) 
    for(int z=0;z<ISGRI_SIZE;z++)   
      // If pixel is valid
      if(TimeEffMap[y][z]>ZERO)
	{
	  for(int bin=0;bin<NumBinEner;bin++)
	    if(SpecMean[bin]>ZERO)
	      Chi2Map[y][z]+= pow((SpecPix[y][z][bin]-SpecMean[bin]),2)/SpecMean[bin]; 
	  if(Chi2Map[y][z]>ZERO)
	    {
	      Chi2Map[y][z]= log(Chi2Map[y][z])/log(10.);
	      for(int m=0;m<NumBinChi2;m++)
		if( (Chi2Map[y][z]>=BinChi2[m]) && (Chi2Map[y][z]<BinChi2[m+1]) )
		  {
		    DistChi2[m]++;
		    break;
		  }
	      if(Chi2Map[y][z]>=BinChi2[NumBinChi2])
		DistChi2[NumBinChi2-1]++;
	    }
	}

      
  // ==================================================================
  //        Estimate Mean and Std Deviation for this Distribution
  // ==================================================================

  // Sum probabilities pi
  double 
    Tot=0.,
    Mean=0.,
    StdDev=0.;
  for(int m=0;m<NumBinChi2;m++)
    Tot+= DistChi2[m];
  if(Tot>ZERO)
    {
      // Mean~= Sum[xi.pi]/Sum[pi]
      for(int m=0;m<NumBinChi2;m++)
	Mean+= DistChi2[m]*BinChi2[m];
      Mean/= Tot;
      
      // Var~= Sum[(xi-Mean)².pi]/Sum[pi]
      for(int m=0;m<NumBinChi2;m++)
	StdDev+= DistChi2[m] * pow(BinChi2[m]-Mean,2);
      StdDev= sqrt(StdDev/Tot);
    }

  // ==================================================================
  //                     Determine Noisy pixels
  // ==================================================================
  double NoisyThresh= Mean+criteria*StdDev;
  if(Mean<NoisyThresh)
    {
      if(detailedOutput)
	RILlogMessage(NULL,Log_1,"=> Mean= %.3f, Sigma= %.3f",Mean,StdDev);
      int N=0;
      for(int y=0;y<ISGRI_SIZE;y++)
	for(int z=0;z<ISGRI_SIZE;z++)
	  if(Chi2Map[y][z]>NoisyThresh) 
	    {
	      SpecNoisyMap[y][z]= (Chi2Map[y][z]-Mean)/StdDev;
	      N++;
	    }
      *NumSpecNoisy= N;
    }
  else
    {
      RILlogMessage(NULL,Warning_1,"Found Std Deviation is null!");
      RILlogMessage(NULL,Warning_1,"Assuming no noisy pixels.");
    }
  
  // ==================================================================
  //                        Free memory and exit
  // ==================================================================

  for(int y=0;y<ISGRI_SIZE;y++)
    {
      for(int z=0;z<ISGRI_SIZE;z++)
	{
	  if(SpecPix[y][z])
	    { 
	      free(SpecPix[y][z]);
	      SpecPix[y][z]= NULL;
	    }
	  if(EffMap[y][z])
	    { 
	      free(EffMap[y][z]);
	      EffMap[y][z]= NULL;
	    }
	}
      if(TimeEffMap[y])
	{ 
	  free(TimeEffMap[y]);
	  TimeEffMap[y]= NULL;
	}
      if(SpecPix[y])
	{ 
	  free(SpecPix[y]);
	  SpecPix[y]= NULL;
	}
      if(Chi2Map[y])
	{ 
	  free(Chi2Map[y]);
	  Chi2Map[y]= NULL;
	}
      if(EffMap[y])
	{ 
	  free(EffMap[y]);
	  EffMap[y]= NULL;
	}
    }
  if(TimeEffMap)
    { 
      free(TimeEffMap);
      TimeEffMap= NULL;
    }
  if(EffMap)
    { 
      free(EffMap);
      EffMap= NULL;
    }
  if(SpecPix)
    { 
      free(SpecPix);
      SpecPix= NULL;
    }
  if(Chi2Map)
    { 
      free(Chi2Map);
      Chi2Map= NULL;
    }
  if(BinEner)
    {
      free(BinEner);
      BinEner= NULL;
    }
  if(BinChi2)
    {
      free(BinChi2);
      BinChi2= NULL;
    }
  if(SpecMean)
    {
      free(SpecMean); 
      SpecMean= NULL;
    }
  if(DistChi2)
    {
      free(DistChi2); 
      DistChi2= NULL;
    }
  return status;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//               Compute ISGRI and Modules' ONTIMEs for this Time BIN, by handling relevant GTIs                      
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int ComputeTimeEffMap(OBTime        OBTstart,         // Time Bin starting time
		      OBTime        OBTend,           // Time Bin finishing time 
		      OBTime        *GTIstart,        // Good time intervals, starting times               
		      OBTime        *GTIstop,         // Good time intervals, finishing times              
		      int           NumGTI,           // Number of good time intervals                     
		      OBTime        **GTIstartPerMCE, // Good time intervals per module, starting times    
		      OBTime        **GTIstopPerMCE,  // Good time intervals per module, stopping times    
		      int           *NumGTIperMCE,    // Number of good time intervals  per module  
		      dal_float     **DeadTimes,      // DeadTimes for each module, over the SCW
		      OBTime        *OBTdeadTimes,    // ISGRI DeadTimes dates over the SCW
		      long          NumDeadTimes,     // Number ISGRI DeadTimes over the SCW
		      dal_int       **ONpixelsREVmap, // Map of Pixels Status for this REV (REV context)
		      dal_int       **ONpixelsHK3map, // Map of Pixels Status for this SCW (HK3 data)
		      dal_int       **SelectFlagMap,  // Map of Noisy Pixels, as detected by SELECT_FLAG method
		      dal_double    **TimeEffMap,     // Output: ISGRI Time Efficiency map
		      dal_double    *MeanTimeEff,     // Output: ISGRI Average Time Efficiency
		      dal_double    *OnTime,          // Output: ISGRI ONTIME (sec)
		      dal_double    *DeadC,           // Output: ISGRI DeadTime Correction factor (%)
		      unsigned char detailedOutput)   // Detailed output
{
  int 
    status=    ISDC_OK,
    RILstatus= ISDC_OK,
    NumMergedGTI= 0,
    *NumMCEGTI= NULL;
  dal_double 
    *ModOnTime=   NULL,
    *ModDTfactor= NULL;
  OBTime        
    *MergedGTIstart= NULL,
    *MergedGTIstop=  NULL,
    *MergedStart=    NULL,
    *MergedStop=     NULL,
    **MCEGTIstart=   NULL,
    **MCEGTIstop=    NULL;

  // ==================================================================
  //                 Check TimeBin interval coherence
  //  Nota: [OBTstart,OBTend] are times of events inside the Time Bin 
  //       BEFORE it is merged with GTIs (which is done here)
  // ==================================================================

  dal_double dt;
  if((status= DeltaOBT(OBTend, OBTstart, &dt, status))<0)
    {
      RILstatus=RILlogMessage(NULL,Error_1,"OBTime difference error for TimeBin [ %llu, %llu ]!",
			      OBTstart, OBTend);
      return status;
    }
  *MeanTimeEff= 0.; 
  *OnTime= 0.;
  *DeadC= 0.;

  // ==================================================================
  //                     Modules allocations
  // ================================================================== 

  if( (ModOnTime= (dal_double*)calloc(IBIS_NUM_BLOCK,sizeof(dal_double))) == NULL)
    {
      RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation Error for ModOnTime dal_double*.");  
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    } 
  if( (ModDTfactor= (dal_double*)calloc(IBIS_NUM_BLOCK,sizeof(dal_double))) == NULL)
    {
      RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation Error for ModDTfactor dal_double*.");  
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    } 

  // ==================================================================
  //           Merge ISGRI GTIs with the TimeBin interval
  // ==================================================================

  if(detailedOutput)
    RILstatus= RILlogMessage(NULL, Log_1, "Merging Time Interval with ISGRI GTIs."); 
  if( (MergedGTIstart= (OBTime*)calloc(NumGTI,sizeof(OBTime))) == NULL ) 
    {
      RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation Error while merging GTIs with time interval.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if( (MergedGTIstop= (OBTime*)calloc(NumGTI,sizeof(OBTime))) == NULL ) 
    {
      RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation Error while merging GTIs with time interval.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if( (status= DAL3HKmergeGTI(NumGTI,GTIstart,GTIstop, 1,&OBTstart,&OBTend,
			      &NumMergedGTI,MergedGTIstart,MergedGTIstop, status)) != ISDC_OK )
    {
      RILstatus= RILlogMessage(NULL, Error_1, "Error while merging GTIs with time interval: DAL3HKmergeGTI Error %d.",status);
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  // Check Number of resulting merged intervals
  if(!NumMergedGTI)
    {
      RILstatus= RILlogMessage(NULL, Warning_1, "Merging ISGRI GTIs with the Time BIN gave a NULL interval!");
      RILstatus= RILlogMessage(NULL, Warning_1, "Setting ISGRI and modules OnTimes to 0.");
      *OnTime= 0.;
      for(int i=0; i<IBIS_NUM_BLOCK; i++)      
	ModOnTime[i]= 0.;
      RILstatus= RILlogMessage(NULL, Warning_1, "Returning ISDC_OK");
      return ISDC_OK;
    }
  // Re-allocation for (merged) global GTIs 
  if((MergedGTIstop= (OBTime *)realloc(MergedGTIstop ,NumMergedGTI*sizeof(OBTime)))==NULL)
    {
      RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation Error while merging GTIs with time interval.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if((MergedGTIstart=(OBTime *)realloc(MergedGTIstart,NumMergedGTI*sizeof(OBTime)))==NULL)
    {
      RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation Error while merging GTIs with time interval.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }

  // ==================================================================
  //                      Merge MCEGTIs with GTIs
  // ==================================================================

  if(detailedOutput)
    RILstatus= RILlogMessage(NULL, Log_1, "Merging ISGRI GTIs with MCE_GTIs."); 
  if( (NumMCEGTI= (int*)calloc(IBIS_NUM_BLOCK,sizeof(int))) == NULL)
    {
      RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation Error while building MCE_GTIs int*.");  
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }  
  if( (MCEGTIstart= (OBTime**)calloc(IBIS_NUM_BLOCK,sizeof(OBTime*))) == NULL)
    {
      RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation Error while building MCE_GTIs OBTime**.");   
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if( (MCEGTIstop=  (OBTime**)calloc(IBIS_NUM_BLOCK,sizeof(OBTime*))) == NULL)
    {
      RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation Error while building MCE_GTIs OBTime**.");    
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }  
  for(int i=0; i<IBIS_NUM_BLOCK; i++) 
    {
      if( (MCEGTIstart[i]= (OBTime*)calloc(NumMergedGTI+NumGTIperMCE[i],sizeof(OBTime))) == NULL)
	{
	  RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation Error while building MCE_GTIs OBTime*."); 
	  return ERR_ISGR_OSM_MEMORY_ALLOC;
	}
      if( (MCEGTIstop[i]= (OBTime*)calloc(NumMergedGTI+NumGTIperMCE[i],sizeof(OBTime))) == NULL)
	{
	  RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation Error while building MCE_GTIs OBTime*.");  
	  return ERR_ISGR_OSM_MEMORY_ALLOC;
	} 
      if( (status= DAL3HKmergeGTI(NumMergedGTI,MergedGTIstart,MergedGTIstop,NumGTIperMCE[i],GTIstartPerMCE[i],GTIstopPerMCE[i],
				  NumMCEGTI+i,MCEGTIstart[i],MCEGTIstop[i],status)) != ISDC_OK )
	{
	  RILstatus= RILlogMessage(NULL, Error_1, "Error while merging ISGRI GTIs with MCE_GTIs: DAL3HKmergeGTI Error %d",status);
	  return ERR_ISGR_OSM_MEMORY_ALLOC;
	} 
      // Check Number of resulting merged intervals
      if(!NumMCEGTI[i])
	{
	  RILstatus= RILlogMessage(NULL, Warning_1, "Merging Time Bin GTIs with the MCE GTIs of module[%d] gave a NULL interval!",i);
	  RILstatus= RILlogMessage(NULL, Warning_1, "Setting ModOnTime[mce%d] to 0.",i);    
	}
      else
	{
	  // Re-Allocation for (merged) MCE_GTI[i]     
	  if( (MCEGTIstart[i]= (OBTime*)realloc(MCEGTIstart[i],NumMCEGTI[i]*sizeof(OBTime))) == NULL)
	    {
	      RILstatus= RILlogMessage(NULL, Error_2, "Memory reallocation Error while building MCE_GTIs OBTime*.");   
	      return ERR_ISGR_OSM_MEMORY_ALLOC;
	    }
	  if( (MCEGTIstop[i]= (OBTime*)realloc(MCEGTIstop[i],NumMCEGTI[i]*sizeof(OBTime))) == NULL)
	    {
	      RILstatus= RILlogMessage(NULL, Error_2, "Memory reallocation Error while building MCE_GTIs OBTime*.");   
	      return ERR_ISGR_OSM_MEMORY_ALLOC;
	    } 
	}
    }

  // ====================================================================================
  //  MCE'Effective observation time and DeadTime correction factor are based on MCEGTIs
  // ====================================================================================

  long Nmod=0, NumDTmod;
  dal_double tmp;
  dal_double det;
  OBTime OBTdead1;
  OBTime OBTdead2;
  OBTime OBTdif;
  OBTime OBTround;
  /* 4s period in OBT units */
  dal_double de4s= 4194300.616;
  dal_double delti;
  dal_double sumti;
  dal_double sumco;
  dal_double sumtime;
  dal_double wedefactor;
  dal_double wedeco;
  int ireset;
  int iseconds;
  int ifracsec;
  
  wedeco = 0.;
  for(int mce=0;mce<IBIS_NUM_BLOCK;mce++) 
    {
      // Compute ModOnTime
      ModOnTime[mce]= 0.;
      if(NumMCEGTI[mce])
	for(int i=0;i<NumMCEGTI[mce];i++) 
	  {
	    if( (status= DeltaOBT(MCEGTIstop[mce][i], MCEGTIstart[mce][i], &tmp, status)) <0 ) 
	      {
		RILstatus= RILlogMessage(NULL,Error_1,"MCE%d, MCE_GTI_%d: OBTime difference fatal error (status= %d).",mce,i,status);
		return status;
	      }
	    ModOnTime[mce]+= tmp;
	  }

      // Compute ModDTfactor and DEADC  
      NumDTmod= 0;
      sumco = 0.;
      sumtime = 0.;
      wedefactor = 0.;
      for(long dt=0;dt<NumDeadTimes;dt++) {

     // Dead time bin limits
        OBTdead1 = OBTdeadTimes[dt]-de4s;
        OBTdead2 = OBTdeadTimes[dt]+de4s;		

        sumti = 0.;

	for(long i=0;i<NumMCEGTI[mce];i++) {

       // Dead time limit within GTI (normal situation)
	  if( (OBTdead1 > MCEGTIstart[mce][i]) && (OBTdead2 < MCEGTIstop[mce][i]) )
	    { 
            sumti = sumti + 8.0;
            if ((mce == 0) && (dt > 59 && dt < 70)) 
              { 
              RILstatus= RILlogMessage(NULL, Log_1, "dt= %ld, OBTdeadTimes= %lld, OBTdead1= %lld, OBTdead2= %lld", 
                           dt,i,OBTdeadTimes[dt],OBTdead1,OBTdead2); 
              RILstatus= RILlogMessage(NULL, Log_1, "i= %ld, MCEGTIstart= %lld, MCEGTIstop= %lld", 
                           i,MCEGTIstart[mce][i],MCEGTIstop[mce][i]); 
	      RILstatus= RILlogMessage(NULL, Log_1, "sumti= %.4f", sumti); 
              }				 
            }

       // GTI within the dead time limit
	  if( (OBTdead2 >= MCEGTIstart[mce][i]) && (OBTdead1 <= MCEGTIstop[mce][i]) )
	    { 

            iseconds = 0;
            ifracsec = 0;

       // Entire GTI within the dead time limit
            if( (OBTdead1 <= MCEGTIstart[mce][i]) && (OBTdead2 >= MCEGTIstop[mce][i]) )
              {
              status = DAL3GENroundOBT(MCEGTIstop[mce][i]-MCEGTIstart[mce][i], 0.0001, &OBTround, status);		  
	      status = DAL3GENsplitOBT(OBTround,&ireset,&iseconds,&ifracsec,status);
              }
              
       // GTI on the lower border
            if( (OBTdead1 > MCEGTIstart[mce][i]) && (OBTdead2 >= MCEGTIstop[mce][i]) )
              {
              status = DAL3GENroundOBT(MCEGTIstop[mce][i]-OBTdead1, 0.0001, &OBTround, status);		  
	      status = DAL3GENsplitOBT(OBTround,&ireset,&iseconds,&ifracsec,status);
              }        
              
       // GTI on the upper border
            if( (OBTdead1 <= MCEGTIstart[mce][i]) && (OBTdead2 < MCEGTIstop[mce][i]) )
              {
              status = DAL3GENroundOBT(OBTdead2-MCEGTIstart[mce][i], 0.0001, &OBTround, status);		  
	      status = DAL3GENsplitOBT(OBTround,&ireset,&iseconds,&ifracsec,status);
              }        
              
            delti = iseconds + 4.0*ifracsec/de4s;
            sumti = sumti + delti;
          
       // Print control sequence
            if ((mce == 0) && (dt > 59 && dt < 70)) 
              { 
              RILstatus= RILlogMessage(NULL, Log_1, "dt= %ld, OBTdeadTimes= %lld, OBTdead1= %lld, OBTdead2= %lld", 
                           dt,i,OBTdeadTimes[dt],OBTdead1,OBTdead2); 
              RILstatus= RILlogMessage(NULL, Log_1, "i= %ld, MCEGTIstart= %lld, MCEGTIstop= %lld", 
                           i,MCEGTIstart[mce][i],MCEGTIstop[mce][i]); 
	      RILstatus= RILlogMessage(NULL, Log_1, "OBTround= %lld, iseconds= %d, ifracsec= %d", 
	  		 OBTround,iseconds,ifracsec); 
	      RILstatus= RILlogMessage(NULL, Log_1, "delti= %.4f, sumti= %.4f", delti, sumti); 
              }				 
	    }

       // Standard method
	  if( (OBTdeadTimes[dt]>=MCEGTIstart[mce][i]) && (OBTdeadTimes[dt]<=MCEGTIstop[mce][i]) )
	    {
	      NumDTmod++;
	      ModDTfactor[mce]+= DeadTimes[dt][mce];
	      break;
	    }
	    
       // End of GTI loop   
	  }

     // Summed correction for a given dead time interval
        sumtime = sumtime + sumti;
        sumco = sumco + sumti*DeadTimes[dt][mce];

     // End of dead time loop
	}   

   // Weighted dead time correction for a given module
      if (sumtime > 0.0001)
        {
        wedefactor = 1.-sumco/sumtime;
        wedeco = wedeco + wedefactor;
        }   
      
   // Print the result    	
      RILstatus= RILlogMessage(NULL, Log_1, "MCE%d, sumco= %.4f, sumtime= %.4f, wedefactor= %.4f, wedeco= %.3f", 
				 mce,sumco,sumtime,wedefactor,wedeco); 

   // Standard method	
      if(NumDTmod)
	{
	  ModDTfactor[mce]= 1.-ModDTfactor[mce]/NumDTmod;
	  *DeadC+= ModDTfactor[mce];
	  Nmod++;
	}

      // Info
    //  if(detailedOutput)
	RILstatus= RILlogMessage(NULL, Log_1, "MCE%d: OnTime= %.3fsec, NumDTmod= %d, DTfactor= %.3f%c", 
				 mce,ModOnTime[mce],NumDTmod,100*ModDTfactor[mce],'%'); 
				 
   // Replace standard method with the weighted mean
      if (sumtime > 0.) {
        ModDTfactor[mce] = wedefactor;
      }

 // End of MCE loop				 
    }
  
  // Compute keyword DEADC for this TimeBin
  if(Nmod)
//    *DeadC/= Nmod;
    *DeadC= wedeco/Nmod;
  else
    {
      *DeadC=1.;
      RILstatus= RILlogMessage(NULL, Warning_1, "No ISGRI DeadTimes => DT Correction Factor is set to 1.");
    }

  // Compute keyword ONTIME for this TimeBin (out of Merged ISGRI GTIs)
  *OnTime= 0.;
  for(int i=0;i<NumMergedGTI;i++) 
    {
      if( (status= DeltaOBT(MergedGTIstop[i],MergedGTIstart[i],&tmp,status)) <0 ) 
	{
	  RILstatus= RILlogMessage(NULL,Error_1,"MCE%d: OBTime difference fatal error (status= %d). Please check input GTIs.",i,status);
	  return status ;
	}
      *OnTime+= tmp;
    }
  if(*OnTime<ZERO)
    {
      RILstatus= RILlogMessage(NULL, Warning_1, "ISGRI OnTime is invalid (%gsec)",*OnTime);
      return ISDC_OK;
    }
  else
    RILstatus= RILlogMessage(NULL, Log_1, "ISGRI: ONTIME= %.3fsec, DEADC= %.3f%c", *OnTime, 100.*(*DeadC),'%');

  // ====================================================================================
  //                          Compute ISGRI Time Efficiency Map
  //  Note: apply threshold on it to prevent low efficiency pixels from being used later
  // ====================================================================================

  int mce;
  *MeanTimeEff= 0.;
  for(int y=0;y<ISGRI_SIZE;y++) 
    for(int z=0;z<ISGRI_SIZE;z++)
      {
	if(ONpixelsREVmap[y][z] && ONpixelsHK3map[y][z] && !SelectFlagMap[y][z])
	  {
	    mce= YZtoModN(y,z);
	    tmp= ModDTfactor[mce]*ModOnTime[mce]/(*OnTime);
	    if(tmp>=TIME_EFF_THRESH)
	      {
		TimeEffMap[y][z]= tmp;
		*MeanTimeEff+= TimeEffMap[y][z];
	      }
	    else
	      TimeEffMap[y][z]= 0.;
	  }
	else
	  TimeEffMap[y][z]= 0.;
      }
  *MeanTimeEff/= ISGRI_SIZE*ISGRI_SIZE;

  // ====================================================================================
  //                                Free memory and exit  
  // ====================================================================================

  if(MergedGTIstart) 
    {
      free(MergedGTIstart);
      MergedGTIstart= NULL;
    }
  if(MergedGTIstop) 
    {
      free(MergedGTIstop);
      MergedGTIstop= NULL;
    } 
  if(NumMCEGTI) 
    {
      free(NumMCEGTI);
      NumMCEGTI= NULL;
    } 
  if(ModOnTime) 
    {
      free(ModOnTime);
      ModOnTime= NULL;
    }  
  if(ModDTfactor) 
    {
      free(ModDTfactor);
      ModDTfactor= NULL;
    }
  if(MCEGTIstart) 
    {
      for(int i=0;i<IBIS_NUM_BLOCK;i++)
	if(MCEGTIstart[i])
	  {
	    free(MCEGTIstart[i]);
	    MCEGTIstart[i]= NULL;
	  }
      free(MCEGTIstart);
      MCEGTIstart= NULL;
    }
  if(MCEGTIstop) 
    {
      for(int i=0;i<IBIS_NUM_BLOCK;i++)
	if(MCEGTIstop[i])
	  {
	    free(MCEGTIstop[i]);
	    MCEGTIstop[i]= NULL;
	  }
      free(MCEGTIstop);
      MCEGTIstop= NULL;
    }
  if(MergedStart) 
    {
      free(MergedStart);
      MergedStart= NULL;
    }
  if(MergedStop) 
    {
      free(MergedStop);
      MergedStop= NULL;
    }
  return ISDC_OK;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//                                           Create IMAGE
//                      CAUTION: Energy band must be [INF;SUP[, not [INF;SUP]!
//////////////////////////////////////////////////////////////////////////////////////////////////////////

int MkImage(dal_byte      *IsgriY,        // Detector coordinate on Y                                 
            dal_byte      *IsgriZ,        // Detector coordinate on Z                                 
            dal_double    *IsgriEnergy,   // Phase height amplitude raw / corrected or deposed energy 
            long          NumRow,         // Number of events                                         
            dal_float     **EnergyBounds, // Energy/channel intervals                                 
            int           NumImaBin,      // Number of energy intervals
            dal_double    ***IsgriSHD,    // Output: resultant image 
	    long          *NumEvents)     // Output: Number of events in these Ebands
{
  for(int k=0;k<NumImaBin;k++)
    NumEvents[k]=0;

  for(long i=0;i<NumRow;i++)
    for(int k=0; k<NumImaBin; k++) 
      if( (IsgriEnergy[i]>=(double)EnergyBounds[k][0]) && (IsgriEnergy[i]<(double)EnergyBounds[k][1]) ) 
	{
	  IsgriSHD[k][IsgriY[i]][IsgriZ[i]]++;
	  NumEvents[k]++;
	}

  return ISDC_OK;
}
	    

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                       The actual image efficiency is computed in this function.   
//                   NOTE: we here use the analytical ERF function, via a call to LTfunction                        
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MkeffImage(dal_double    **TimeEffMap,    // Pixels Time Efficiency map                                        
	       dal_double    **LowThresImage, // Low Threshold map (keV)
	       dal_double    **SpecNoisyMap,  // Map of Noisy Pixels, as detected by spectral method
	       dal_float     **EnergyBounds,  // Energy/channel intervals                                    
	       int           NumImaBin,       // Number of energy bands                                      
	       int           revol_scw,       // Revolution number
	       dal_double    ***IsgriEffSHD,  // Output: Image efficiency map, one per energy band                   
	       dal_double    *MeanEff,
           ISGRI_efficiency_struct *ptr_ISGRI_efficiency)        // Output: Means of the efficiency maps to be calculated
{
  int
    RILstatus= ISDC_OK,
    NumSteps;
  double
    //Step value < min(bin_of_RMF_full_resolution)
    Step= 0.478,
    NewStep, Eband, E, tmp, Int, Tot;
  dal_double 
    **EngEffMap= NULL;


  // ==================================================================
  //                    Allocate generic  EffEng map
  // ==================================================================

  if((EngEffMap= (dal_double**)calloc(ISGRI_SIZE,sizeof(dal_double*)))==NULL) 
    {
      RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : EngEffMap.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  for(int IsgrY=0;IsgrY<ISGRI_SIZE;IsgrY++)
    if((EngEffMap[IsgrY]= (dal_double*)calloc(ISGRI_SIZE,sizeof(dal_double)))==NULL) 
      {
	RILstatus=RILlogMessage(NULL,Warning_2,"Error in memory allocation : EngEffMap[].");
	return ERR_ISGR_OSM_MEMORY_ALLOC;
      }

  // ==================================================================
  //                        LOOP on Energy bands
  // ==================================================================

  for(int k=0;k<NumImaBin;k++)
    {
      MeanEff[k]= 0;
      Eband= EnergyBounds[k][1]-EnergyBounds[k][0];

      // 1- If Energy Band wider than Step, integer Energetic Efficiency function
      if( (Eband>ZERO) && (Eband>Step) )
	{
	  // Determine smooth integration step (for Energetic efficiency computation)
	  NumSteps= (int)rint(Eband/Step);
	  NewStep= Eband/NumSteps;
	  
	  // LOOP on pixels
	  for(int IsgrY=0;IsgrY<ISGRI_SIZE;IsgrY++)
	    for(int IsgrZ=0;IsgrZ<ISGRI_SIZE;IsgrZ++)
	      {
		// Make computation only for VALID and NON-NOISY pixels = 
		// those who have a non_null time efficiency
		// and who were NOT detected NOISY by log(Chi2) method
		if( (SpecNoisyMap[IsgrY][IsgrZ]<ZERO) && (TimeEffMap[IsgrY][IsgrZ]>ZERO) )
		  {
		    Int= 0; 
		    Tot= 0;
		    for(int i=0;i<NumSteps;i++)
		      {
			// Increment Energy, take E at the middle of bin:
			E= EnergyBounds[k][0] + (i+0.5)*NewStep;
			// Use LTshape analytical approximation + Add "ARF*power_law" 
			tmp= ARF(E) * pow(E,-2);
			Int+= tmp*LTfunction(E,IsgrY,IsgrZ,ptr_ISGRI_efficiency);
			Tot+= tmp;
		      }
		    if(Tot>ZERO)
		      EngEffMap[IsgrY][IsgrZ]= Int/Tot;
		    else
		      EngEffMap[IsgrY][IsgrZ]= 0;
		  }
		else
		  EngEffMap[IsgrY][IsgrZ]= 0;
	      }// pixel LOOP
	}

      // 2- If Energy Band narrower than NewStep, get Energetic Efficiency function middle value
      else if( (Eband>ZERO) && (Eband<=Step) )
	{
	  for(int IsgrY=0;IsgrY<ISGRI_SIZE;IsgrY++)
	    for(int IsgrZ=0;IsgrZ<ISGRI_SIZE;IsgrZ++)
	      {
		if( (SpecNoisyMap[IsgrY][IsgrZ]<ZERO) && (TimeEffMap[IsgrY][IsgrZ]>ZERO) )
		  {
		    E= (EnergyBounds[k][0]+EnergyBounds[k][1])/2;
		    EngEffMap[IsgrY][IsgrZ]= LTfunction(E,IsgrY,IsgrZ,ptr_ISGRI_efficiency);
		  }
		else
		  EngEffMap[IsgrY][IsgrZ]= 0;
	      }
	}

      // 3- If Energy Band is NULL, assume Energetic Efficiency Map is NULL
      else
	{
	  RILstatus=RILlogMessage(NULL,Warning_2,"Energy Range [%.2f,%.2f] is NULL!",EnergyBounds[k][0],EnergyBounds[k][1]);
	  RILstatus=RILlogMessage(NULL,Warning_2,"Program will assume Efficiency Map is NULL.");
	  for(int l=0;l<8;l++) 
	    SetImageMCEtoZero(EngEffMap,l);
	}

      // 4- Compute TOTAL EFFICIENCY
      for(int IsgrY=0;IsgrY<ISGRI_SIZE;IsgrY++)
	for(int IsgrZ=0;IsgrZ<ISGRI_SIZE;IsgrZ++)
	  {
	    // Select VALID pixels: check ENG and TEMP efficiencies
	    if( (EngEffMap[IsgrY][IsgrZ]>ZERO) && (TimeEffMap[IsgrY][IsgrZ]>ZERO) )
	      {
		IsgriEffSHD[k][IsgrY][IsgrZ]= EngEffMap[IsgrY][IsgrZ] * TimeEffMap[IsgrY][IsgrZ];	   
		MeanEff[k]+= IsgriEffSHD[k][IsgrY][IsgrZ];
	      }
	    // else, kill pixel
	    else
	      IsgriEffSHD[k][IsgrY][IsgrZ]= 0;
	  }
      MeanEff[k]/= ISGRI_SIZE*ISGRI_SIZE;

    } // Eband loop

	  
  // ==================================================================
  //                              Exit
  // ==================================================================

  for(int IsgrY=0;IsgrY<ISGRI_SIZE;IsgrY++)
    if(EngEffMap[IsgrY])
      { 
	free(EngEffMap[IsgrY]); 
	EngEffMap[IsgrY]= NULL;
      }
  if(EngEffMap)
    { 
      free(EngEffMap); 
      EngEffMap= NULL;
    }
  return ISDC_OK;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                  WRITE IMAGES AND ITS KEYWORDS
//////////////////////////////////////////////////////////////////////////////////////////////////////////

int WriteImage(dal_element   *NewGRP,        // Pointer to the Science window                           
	       dal_element   **image,        // Pointer to the object in which image is written         
	       dal_element   **efficiency,   // Pointer to the object in which efficiencies are written 
	       dal_double    ***IsgriSHD,    // Image matrices                                            
	       dal_double    ***IsgriEffSHD, // Efficiency matrices
	       OBTime        StartTime,      // Start time of science window                            
	       OBTime        EndTime,        // End Time of science window                              
	       OBTime        t0,             // Lower boundary of time interval                         
	       OBTime        t1,             // Upper boundary of time interval                         
	       dal_float     **EnergyBounds, // energy/channel intervals                                
	       int           NumImaBin,      // number of energy intervals                              
	       char          *BndTyp,        // band type                                               
	       char          *ChTyp,         // channel type                    
	       double        OnTime,         // Effective time computed with global GTIs  
	       double        DeadC,          // Dead Time Correction factor
	       double        MeanTimeEff)    // Mean of the Time Efficiency map of the current Time Bin
{
  int              
    status=    ISDC_OK, 
    RILstatus= ISDC_OK,
    NumAxis=   2;
  char
    *name= NULL;
  unsigned long 
    AttributeCode= 0;
  long
    *StartIdx= NULL, 
    *EndIdx=   NULL;
  double  
    dt=0.0,            // A time difference.           
    *TmpArray= NULL;   // A temporary array just to write the images to the fits file  
  OBTime
    TmpOBTime= 0;
  dal_dataType   
    type;
  T_OSM_ATTRIBUTE 
    OSMattribute;

  // Allocations
  if((name= (char*)malloc(128*sizeof(char)))==NULL) 
    {
      RILstatus=RILlogMessage(NULL, Error_2, "Memory allocation Error for name.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if((StartIdx= (long *)calloc(NumAxis, sizeof(long)))==NULL) 
    {
      RILstatus=RILlogMessage(NULL, Error_2, "Memory allocation Error for StartIdx.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if((EndIdx= (long *)calloc(NumAxis, sizeof(long)))==NULL) 
    {
      RILstatus=RILlogMessage(NULL, Error_2, "Memory allocation Error for EndIdx.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  if((TmpArray= (double *)malloc(ISGRI_SIZE*ISGRI_SIZE*sizeof(double)))==NULL)
    {
      RILstatus=RILlogMessage(NULL, Error_2, "Memory allocation Error for TmpArray.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }

  // Compute non energy related Keywords
  for(int i=0;i<NumAxis;i++) 
    {
      StartIdx[i]=1;
      EndIdx[i]=ISGRI_SIZE;
    }
  OSMattribute.obtstart=StartTime;
  OSMattribute.obtstop=EndTime;
  TmpOBTime=StartTime;
  if((status=DAL3AUXconvertOBT2IJD(NewGRP, TCOR_ANY, 1, &TmpOBTime, &OSMattribute.tstart, status))!=ISDC_OK) 
    {
      RILstatus=RILlogMessage(NULL, Log_1, 
			      "11/ writeImage : Error in converting OBT times to ISDC Julian days, status = %d.", status);
      return status;
  }
  TmpOBTime=EndTime;
  if((status=DAL3AUXconvertOBT2IJD(NewGRP, TCOR_ANY, 1, &TmpOBTime, &OSMattribute.tstop, status))!=ISDC_OK) 
    { 
      RILstatus=RILlogMessage(NULL, Log_1, 
			      "22/ writeImage : Error in converting OBT times to ISDC Julian days, status = %d.", status);
      return status;
    }
  TmpOBTime=t0;
  if((status=DAL3AUXconvertOBT2IJD(NewGRP, TCOR_ANY, 1, &TmpOBTime, &OSMattribute.tfirst, status))!=ISDC_OK) 
    {
      RILstatus=RILlogMessage(NULL, Log_1, 
			      "33/ writeImage : Error in converting OBT times to ISDC Julian days, status = %d.", status);
      return status;
    }
  TmpOBTime=t1;
  if((status=DAL3AUXconvertOBT2IJD(NewGRP, TCOR_ANY, 1, &TmpOBTime, &OSMattribute.tlast, status))!=ISDC_OK) 
    {
      RILstatus=RILlogMessage(NULL, Log_1, 
			      "44/ writeImage : Error in converting OBT times to ISDC Julian days, status = %d.", status);
      return status;
    }
  if((status=DeltaOBT(t1, t0, &dt, status))!=ISDC_OK) 
    {
      RILstatus=RILlogMessage(NULL, Error_1, 
			      "writeImage : fatal error in OBT time difference, status = %d.", status);
      return status;
    }
  OSMattribute.telapse=  (dal_int)dt;
  strcpy(OSMattribute.chantype, ChTyp);
  strcpy(OSMattribute.bandtype, BndTyp);
  OSMattribute.obt_acq=  0;
  OSMattribute.bintime=  0;
  OSMattribute.numrange= 0;
  OSMattribute.int_time= 0;
  OSMattribute.ontime=   OnTime;
  OSMattribute.deadc=    DeadC;
  OSMattribute.livetime= OnTime * DeadC;

  // Loop on each Energy Bin
  AttributeCode= OBTSTART_NUM|OBTEND_NUM|TSTART_NUM|TSTOP_NUM|TFIRST_NUM|TLAST_NUM|TELAPSE_NUM|
    BANDTYPE_NUM|RISE_MIN_NUM|RISE_MAX_NUM|ONTIME_NUM|LIVETIME_NUM|DEADC_NUM|
    CHANTYPE_NUM|CHANMIN_NUM|CHANMAX_NUM|E_MIN_NUM|E_MAX_NUM|EXPOSURE_NUM; 
 
  for(int k=0; k<NumImaBin; k++) 
    {
      // Compute energy related Keywords (EXPOSURE must take pixels efficiency into account)
      OSMattribute.chanmin=  (int)EnergyBounds[k][0];
      OSMattribute.chanmax=  (int)EnergyBounds[k][1];
      OSMattribute.e_min=    EnergyBounds[k][0];
      OSMattribute.e_max=    EnergyBounds[k][1];
      OSMattribute.exposure= OnTime * MeanTimeEff;

      // Write Keywords attributes
      status=WriteAttributes(image[k], AttributeCode, OSMattribute, status);

      // Write the shadowgram image to the output file. 
      type= DAL_DOUBLE;
      for(int i=0;i<ISGRI_SIZE;i++) 
        for(int j=0;j<ISGRI_SIZE;j++) 
	  TmpArray[ISGRI_SIZE*j+i]= IsgriSHD[k][i][j];
      if((status=DALarrayPutSection(image[k], NumAxis, StartIdx, EndIdx, type, 
				    (void *)TmpArray, status))!=ISDC_OK) 
	{
	  RILstatus= RILlogMessage(NULL, Warning_1, 
				   "writeImage : error writing image (ISGRI SHADOWGRAM).");
	  RILstatus= RILlogMessage(NULL, Warning_1, 
				   "writeImage : reverting from %d to ISDC_OK to continue.", status);
	  status= ISDC_OK;
	}
      status= CommonStampObject(image[k], "ISGRI Detector Images", status);
    }
  
  // Write the efficiency image to the output file. 
  for(int k=0; k<NumImaBin; k++) 
    {
      // Compute energy related Keywords (EXPOSURE must take pixels efficiency into account)
      OSMattribute.chanmin=  (int)EnergyBounds[k][0];
      OSMattribute.chanmax=  (int)EnergyBounds[k][1];
      OSMattribute.e_min=    EnergyBounds[k][0];
      OSMattribute.e_max=    EnergyBounds[k][1];
      OSMattribute.exposure= OnTime * MeanTimeEff;

      // Write Keywords attributes
      status=WriteAttributes(efficiency[k], AttributeCode, OSMattribute, status);

      // Write the efficiency image to the output file. 
      type= DAL_DOUBLE;      
      for(int i=0; i<ISGRI_SIZE; i++)
        for(int j=0; j<ISGRI_SIZE; j++) 
	  TmpArray[ISGRI_SIZE*j+i]=IsgriEffSHD[k][i][j]; 
      if((status=DALarrayPutSection(efficiency[k], NumAxis, StartIdx, EndIdx, type, (void *)TmpArray, status))!=ISDC_OK) 
	{
	  RILstatus=RILlogMessage(NULL, Warning_1, "IbisIsgrOSMwriteImage : error writing image efficiency (ISGRI EFFICIENCY).");
	  RILstatus=RILlogMessage(NULL, Warning_1, "IbisIsgrOSMwriteImage : reverting from %d to ISDC_OK to continue.", status);
	  status=ISDC_OK;
	}
      status=CommonStampObject(efficiency[k], "ISGRI Detector Efficiency Images", status);
    }

  // Free memory and exit
  free(TmpArray);
  free(EndIdx);
  free(StartIdx);
  free(name);
  return status;
}


/////////////////////////////////////////////////////////////////////////////////////////////
//                     Compute the ISGRI Low Threshold function.
//
// Purpose: we compute, at the given Energy, the value of the Low Threshold function (0:1) 
//          using an analytical approximation of the theoretical ISGRI Low THreshold function,
//          naming:
//                  CONV( Heaviside(LTkeV), Gauss(sigma=FWHM/2.36) )
//   
// Work   : The used approximation is an ERF function re-centered on LTkeV and dilated  by a 
//          fitted factor of dilatation. We use the "TMATH" library approximation of ERF.
//          Added function: power law in "-2.3" to match ISGRI background decrease.
//
// Nota   : Return 0 for pixels with dummy or 63 step (LTkeV is null)
//          so as to kill them ( LTfunction=0 => Energetic_Efficiency=0 => Total_Efficiency=0 ). 
//
/////////////////////////////////////////////////////////////////////////////////////////////


double LTfunction(double energy,
          int y,
          int z,
          ISGRI_efficiency_struct *ptr_ISGRI_efficiency
		  ) {
    double f;
    int chatter=0;

    if (y==10 && z==10) chatter=10;

    DAL3IBIS_get_ISGRI_efficiency(energy,y,z,ptr_ISGRI_efficiency,&f,chatter,0);

    if (chatter>9)
        printf("%.5lg\n",f);

    return f;
}

 
double LTfunction_anal(double Energy,
		  double LTkeV,
		  int revol_scw)
{
  // Preliminary check
  if(LTkeV<1)
    return 0.;
  
  // Parameters of the Chebyshev fit
  const dal_double
    a1 = -1.26551223,   a2 = 1.00002368,
    a3 =  0.37409196,   a4 = 0.09678418,
    a5 = -0.18628806,   a6 = 0.27886807,
    a7 = -1.13520398,   a8 = 1.48851587,
    a9 = -0.82215223,  a10 = 0.17087277;  
  dal_double 
    Sigma= (FWHM0+FWHM1*revol_scw)/2.36,
    InflexionPoint= LTkeV,
    factor= 1.41,  // factor of dilatation to fit ISGRI resolution
    dilatation= factor*Sigma,
    x= (Energy-InflexionPoint)/dilatation,
    z= fabs(x),
    t= 1/(1+0.5*z);

  if(z<=0) 
    return 0.5; // Correction: erfc(0)= 1, LTfunction=0.5 and no 0.
  
  double v= t * exp((-z*z)+a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*(a9+t*a10)))))))));
  if(x<0) 
    v= 2-v;   // erfc(-x)= 2-erfc(x)  
  return 1-v/2;
}


/////////////////////////////////////////////////////////////////////////////////////////////
//            Compute ISGRI detection efficiency (cm²) at a given incident Energy.
//                             Use a fitted analytical expression.
/////////////////////////////////////////////////////////////////////////////////////////////
 
double ARF(double Energy)
{
  dal_double
    Emin= 13.44,
    a0= 84.3, a1= 13.44, a2= 0.85, a3= 25900, a4= 27200, a5= 2.47;

  if(Energy<=Emin)
    return 0.;
  else
    return ( a0 * pow((Energy-a1),a2) * (a3+Energy) / (a4+pow(Energy,a5)) );
}


/////////////////////////////////////////////////////////////////////////////////////////////
//                                   Set IMAGE to ZERO
/////////////////////////////////////////////////////////////////////////////////////////////

void SetImageMCEtoZero(double **Image, 
		       int    MCEid)
{
  int 
    StartY= (MCEid<4) ? 64 : 0,
    StartZ= (3-MCEid%4)*32;

  for(int i=StartY;i<StartY+64;i++)
    for(int j=StartZ;j<StartZ+32;j++)
      Image[i][j]= 0.;
}


/************************************************************************
 * FUNCTION:  ibis_energyIsgrHkCal
 * DESCRIPTION:
 *  Reads the converted HK1 needed for ISGRI LUT1
 *  Returns ISDC_OK if everything is fine, else returns an error code.
 * ERROR CODES:
 *  DAL error codes
 *
 * PARAMETERS:
 *  workGRP       dal_element *      in  DOL of the working group
 *  obtStart           OBTime        in  start of time range
 *  obtEnd             OBTime        in  end   of time range
 *  meanT              double *     out  calculated ISGRI mean temperature
 *  meanBias[]         double    in/out  calculated mean bias per MCE
 *  chatter               int        in  verbosity level
 *  status                int        in  input status
 * RETURN:            int     current status
 ************************************************************************/
int ibis_energyIsgrHkCal(dal_element *workGRP,
                         OBTime       obtStart,
                         OBTime       obtEnd,
                         double      *meanT,
                         double       meanBias[8],
                         int          chatter,
                         int          status)
{
  int     j, totMCE;
  char    hkName[20], num[3];
  long    i, nValues,
          totVal[8];
  double  myMean, myTot,
          meanTscw[8]={-51.0,-51.0,-51.0 -51.0, -51.0,-51.0,-51.0,-51.0},
         *hkBuff=NULL;
  OBTime *obtime,
          startTime=DAL3_NO_OBTIME,
          endTime = DAL3_NO_OBTIME;
  dal_dataType dataType;
  double DtempH1[8] = {0.43, -0.39, -0.77, 0.84, -0.78, 1.09, -0.08, -0.31};
             /* delta from ISGRI mean temperature, to check probe temp2 is OK */

  if (status != ISDC_OK) return status;
  /* SPR 3686: if OBT limits are valid, use S1 PRP OBT limits */
  if (obtStart != DAL3_NO_OBTIME) {

    if ( (obtime=(OBTime *)calloc(1, sizeof(OBTime))) == NULL)
      return( ERR_ISGR_OSM_MEMORY_ALLOC);
      //return(I_ISGR_ERR_MEMORY);
    endTime=obtEnd; 
    /* check if OBT are not too close */
    do {
    status=DAL3GENelapsedOBT(obtEnd, obtStart, &startTime, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Warning_1, "Error calculating elapsed OBT");
      RILlogMessage(NULL, Warning_1, "Reverting from status=%d to ISDC_OK",
                                    status);
      status=ISDC_OK;
      break;
    }
    *obtime=(DAL3_OBT_SECOND) * (SEC_DELTA_MIN);
    if (startTime < *obtime) {
      RILlogMessage(NULL, Warning_1, "Last OBT (%020lld) too close from first OBT (%020lld)",
                                    obtEnd, obtStart);
      status=DAL3GENskipOBT(obtStart, *obtime, &startTime, status);
      if (status != ISDC_OK) {
        RILlogMessage(NULL, Warning_1, "Error adding %d seconds to first OBT", SEC_DELTA_MIN);
        RILlogMessage(NULL, Warning_1, "Reverting from status=%d to ISDC_OK",
                                      status);
        status=ISDC_OK;
        break;
      }
      endTime=startTime;
      RILlogMessage(NULL, Warning_1, "Adding %d seconds to first OBT, last OBT = %020lld",
                                    SEC_DELTA_MIN, endTime);
    }
    } while (0);
    free(obtime);
    startTime=obtStart;

  }
  obtime=NULL;
  totMCE=0;
  for (j=0; j<8; j++) {

    strcpy(hkName, KEY_MCE_TEMP);
    sprintf(num, "%d", j);
    strcat(hkName, num);
    /* get the information for the requested data */
    status=DAL3HKgetValueInfo(workGRP, hkName, IBIS, DAL3HK_CONVERTED,
                              startTime, endTime,  &nValues, &dataType, status);
    if (status != ISDC_OK) break;
    if (nValues < 1) {
      RILlogMessage(NULL, Warning_1, "%13s has NO valid row.", DS_ISGR_HK);
      status=ERR_ISGR_OSM_MEMORY_ALLOC;//I_ISGR_ERR_BAD_INPUT;
      break;
    }
    /* now allocate the required memory and get the data values */
    status=DALallocateDataBuffer((void **)&hkBuff, nValues*sizeof(double), status);
    status=DALallocateDataBuffer((void **)&obtime, nValues*sizeof(OBTime), status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_1, "Cannot allocate buffers for HK data: %s",
                                    hkName);
      break;
    }
    /* retrieve the housekeeping data from the fits file */
    dataType=DAL_DOUBLE;
    status=DAL3HKgetValues(workGRP, hkName, IBIS, DAL3HK_CONVERTED,
                           startTime, endTime, obtime, hkBuff, dataType, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_1, "Cannot get buffers for HK data: %s",
                                    hkName);
      break;
    }
    totVal[j]=0;
    myMean=0.0;
    for (i=0; i<nValues; i++) {
      if (hkBuff[i]  >  KEY_MIN_TEMP) { myMean+=hkBuff[i]; totVal[j]++; }
    }
    if (totVal[j] < 1) {
      RILlogMessage(NULL, Warning_1, "Column %19s has NO valid data.", hkName);
      RILlogMessage(NULL, Warning_1, "Continue to next HK data");
    }
    else {
      meanTscw[j]=myMean/totVal[j];
      totMCE++;
    }

  }
  if (status != ISDC_OK) {
    /* big problem: don't even try to get next HK data */
    DAL3HKfreeData(status);
    if (hkBuff != NULL) DALfreeDataBuffer((void *)hkBuff, ISDC_OK);
    if (obtime != NULL) DALfreeDataBuffer((void *)obtime, ISDC_OK);
    return status;
  }
  if (totMCE == 0) {
    myMean=KEY_DEF_TEMP;
    RILlogMessage(NULL, Warning_1, "Using default ISGRI mean temperature: %+6.2f degC",
                                  myMean);
  }
  else {
    myMean=0.0; myTot=0.0;
    for (j=0; j<8; j++)
      if (totVal[j] > 0) { myMean+=meanTscw[j]; myTot+=totVal[j]; }
    myMean/=totMCE;
    if (chatter > 1)
      RILlogMessage(NULL, Log_1, "Mean temp. (%05.1f values) on %d MCEs: %+6.2f degC",
                                myTot/totMCE, totMCE, myMean);
    /* Check probe is OK, otherwise re-computes the mean with valid values */
    i=totMCE;
    totMCE=0;
    for (j=0; j<8; j++)
      if (totVal[j] > 0) {
        if (fabs(meanTscw[j]-myMean-DtempH1[j]) > KEY_RMS_TEMP) {
          RILlogMessage(NULL, Warning_2,
                             "REJECTING mean temp. on MDU%d: %+6.2f degC",
                             j, meanTscw[j]);
          totVal[j]=0;
        }
        else totMCE++;
      }
    if (i != totMCE) {
      if (totMCE == 0) {
        myMean=KEY_DEF_TEMP;
        RILlogMessage(NULL, Warning_2,
                           "NO mean temp. on MDU fit, Cconvert_obt2ijd(time0,gradient,obt0bis,constant,obt_bias0)HANGE DtempH1 array");
      }
      else {
        myMean=0.0; myTot=0.0;
        for (j=0; j<8; j++)
          if (totVal[j] > 0) { myMean+=meanTscw[j]; myTot+=totVal[j]; }
        myMean/=totMCE;
        if (chatter > 1)
          RILlogMessage(NULL, Log_1, "NEW mean temp. (%05.1f values) on %d MCEs: %+6.2f degC",
                                    myTot/totMCE, totMCE, myMean);
      }
    }
  }
  *meanT=myMean;

  totMCE=0;
  for (j=0; j<8; j++) {

    strcpy(hkName, KEY_MCE_BIAS);
    sprintf(num, "%d", j);
    strcat(hkName, num);
    /* get the information for the requested data */
    status=DAL3HKgetValueInfo(workGRP, hkName, IBIS, DAL3HK_CONVERTED,
                              startTime, endTime,  &nValues, &dataType, status);
    if (status != ISDC_OK) break;
    if (nValues < 1) {
      RILlogMessage(NULL, Warning_1, "%13s has NO valid row.", DS_ISGR_HK);
      status= ERR_ISGR_OSM_MEMORY_ALLOC;//I_ISGR_ERR_BAD_INPUT;
      break;
    }
    /* now allocate the required memory and get the data values */
    status=DALallocateDataBuffer((void **)&hkBuff, nValues*sizeof(double), status);
    status=DALallocateDataBuffer((void **)&obtime, nValues*sizeof(OBTime), status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_1, "Cannot allocate buffers for HK data: %s",
                                    hkName);
      break;
    }
    /* retrieve the housekeeping data from the fits file */
    dataType=DAL_DOUBLE;
    status=DAL3HKgetValues(workGRP, hkName, IBIS, DAL3HK_CONVERTED,
                           startTime, endTime, obtime, hkBuff, dataType, status);
    if (status != ISDC_OK) {
      RILlogMessage(NULL, Error_1, "Cannot get buffers for HK data: %s",
                                    hkName);
      break;
    }
    totVal[j]=0;
    myMean=0.0;
    for (i=0; i<nValues; i++) {
      if (hkBuff[i]  <  KEY_MAX_BIAS) { myMean+=hkBuff[i]; totVal[j]++; }
    }
    if (totVal[j] < 1) {
      RILlogMessage(NULL, Warning_1, "Column %19s has NO valid data.", hkName);
      RILlogMessage(NULL, Warning_1, "Continue to next HK data");
      /* meanBias[j] not changed, already contains default KEY_DEF_BIAS */
    }
    else {
      meanBias[j]=myMean/totVal[j];
      totMCE++;
    }

  }
  if (totMCE) {
    myMean=0.0; myTot=0.0;
    for (j=0; j<8; j++)
      if (totVal[j] > 0) { myMean+=meanBias[j]; myTot+=totVal[j]; }
    myMean/=totMCE;
    if (chatter > 1)
      RILlogMessage(NULL, Log_1, "Mean bias (%05.1f values) on %d MCEs: %+6.1f V",
                                myTot/totMCE, totMCE, myMean);
  }
  else
    RILlogMessage(NULL, Warning_1, "Using default ISGRI mean bias: %+6.1f V",
                                  KEY_DEF_BIAS);
  /* in any case must call this function to free internal buffers */
  DAL3HKfreeData(status);
  if (hkBuff != NULL) DALfreeDataBuffer((void *)hkBuff, ISDC_OK);
  if (obtime != NULL) DALfreeDataBuffer((void *)obtime, ISDC_OK);
  return status;
}

