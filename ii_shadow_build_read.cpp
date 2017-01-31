#include "ii_shadow_build_f.h"


/////////////////////////////////////////////////////////////////////////////////////////////
//                                 RETRIEVE PROGRAM PARAMETERS
/////////////////////////////////////////////////////////////////////////////////////////////

int GetPars(dal_element   **newGRP,          // DOL to the SCW group
            dal_element   **idxREVcontext,   // DOL to the index of REVOLUTION contexts
            dal_element   **idxHK3maps,      // DOL to the index of HK3 noisy maps
            char          *InGTIName,        // Name of the GTIs to be used
            char          *InEFFCDOL,       // DOL of the EFFC structure
            dal_float     ***EnergyBounds,   // Energy binning                
            int           *NumImaBin,        // Number of energy channels                     
            dal_double    *TimeLen,          // Time bin length  
	    char          *UserRowFilter,    // User-defined ROW filter on ISGRI events
            unsigned char *MinRiseTime,      // Minimum value of the Corrected RiseTime selection       
            unsigned char *MaxRiseTime,      // Maximum value of the Corrected RiseTime selection    
	    int           *NoisyDetFlag,     // Spectral Noisy Pixels Detection Flag
            char          *outputLevel,      // OutputLevel for shadow build                  
            unsigned char *detailedOutput   // Detailed Output   
	    )

{
  int   
    status=    ISDC_OK,
    RILstatus= ISDC_OK,
    makeUnique= 1, 
    Clobber;
  char
    *inGRPpar=  (char*)"inSWGGRP",
    *inDOLpar=  (char*)"inRawEvts,inPrwEvts,inSrwEvts,inPrpEvts,inCorEvts,inDead,inGTI,inNoisList", 
    *outGRPpar= (char*)"outSWGGRP",
    *outDOLpar= (char*)"outRawShadow,outEffShadow",
    *TmpEnergyValues= NULL,
    *parName=         NULL,
    *fileName=        NULL,
    *details=         NULL,
    colname[32];  
  float
    *tmpFloat= NULL;
  dal_dataType 
    DALtype= DAL_FLOAT;  
  dal_element   
    *EnergyList= NULL; 

  //char *pname=NULL;
  
  // Allocations
  if((parName=(char *)malloc(PIL_LINESIZE*sizeof(char)))==NULL)
    {
      RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation error: parName"); 
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  
  if((fileName=(char *)malloc(PIL_LINESIZE*sizeof(char)))==NULL)
    {
      RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation error: fileName"); 
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    } 
  if((details=(char *)malloc(PIL_LINESIZE*sizeof(char)))==NULL)
    {
      RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation error: details"); 
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    } 
  if((TmpEnergyValues=(char *)malloc(PIL_LINESIZE*sizeof(char)))==NULL) 
    {
      RILstatus= RILlogMessage(NULL, Error_2, "Memory allocation error: TmpEnergyValues"); 
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  
  // PARAMETER "detailed output or not" 
  strcpy(parName, "details");
  if((status=PILGetString(parName, details))!=ISDC_OK) 
    if(status==PIL_NOT_FOUND) 
      {
	*detailedOutput=0;
	strcpy(details, "NO");
	RILstatus=RILlogMessage(NULL, Log_1, "Parameter <<details>> not listed, using NO.");
	status=ISDC_OK;
      }
  if(!strcmp(details,"YES")) 
    *detailedOutput=1;
  else               
    *detailedOutput=0;

  // PARAMETER "proton table name"
  strcpy(parName, "protonDOL");
 
  // Open the science window 
  if((status=CommonPreparePARsStrings(inGRPpar, inDOLpar, outGRPpar, outDOLpar, makeUnique, newGRP, 
				      &Clobber, status))!=ISDC_OK) 
    {
      RILstatus=RILlogMessage(NULL, Error_1, "CommonPreparePARsStrings Fatal Error : status = %d", status);
      return status;
    }
  
  // PARAMETER "DOL of index of REV contexts""
  strcpy(parName, "idxLowThre");
  if((status=PILGetString(parName, fileName))!=ISDC_OK) 
    {
      if(status==PIL_NOT_FOUND) 
	{
	  RILstatus=RILlogMessage(NULL, Warning_1, "Parameter <<idxLowThre>> not listed.");
	  RILstatus=RILlogMessage(NULL, Warning_1, "This parameter is mandatory must be a valid file.");
	  RILstatus=RILlogMessage(NULL, Warning_1, "Please check the parameter file.");
	  return status;
	}
    }
  else
    {
      RILstatus=RILlogMessage(NULL, Log_1, "REVOL CONTEXT File            : %s",fileName);
      if((status=DALobjectOpen(fileName, idxREVcontext, status))!=ISDC_OK) 
	{
	  RILstatus=RILlogMessage(NULL, Error_1, "Error opening %s.", fileName);
	  RILstatus=RILlogMessage(NULL, Error_1, "Revering from status = %d  to ISDC_OK to continue.", status);
	  idxREVcontext=NULL; 
	  status=ISDC_OK; 
	}
    }
  
  // PARAMETER "DOL of index of HK3 noisy maps"
  strcpy(parName, "idxNoisy");
  if((status=PILGetString(parName, fileName))!=ISDC_OK) {
    if(status==PIL_NOT_FOUND) 
      {
	RILstatus=RILlogMessage(NULL, Warning_1, "Parameter <<idxNoisy>> not listed.");
	RILstatus=RILlogMessage(NULL, Warning_1, "This parameter is mandatory must be a valid file.");
	RILstatus=RILlogMessage(NULL, Warning_1, "Please check the parameter file.");
	return(status);
      }
  }
  else 
    {
      if(detailedOutput) 
	RILstatus=RILlogMessage(NULL, Log_1, "Noisy Index File= %s", fileName);
      if((status=DALobjectOpen(fileName, idxHK3maps, status))!=ISDC_OK) 
	{
	  RILstatus=RILlogMessage(NULL, Error_1, "Error opening %s.", fileName);
	  RILstatus=RILlogMessage(NULL, Error_1, "Revering from status = %d  to ISDC_OK to continue.", status);
	  idxHK3maps=NULL; 
	  status=ISDC_OK; 
	}
      else 
	RILstatus=RILlogMessage(NULL, Log_1, "HK3 Noisy index successfully open."); 
    }
  
  // PARAMETER "Energy bands file"
  TmpEnergyValues[0]='\0';
  strcpy(parName, "inEnergyValues");
  if((status=PILGetString(parName, TmpEnergyValues))!=ISDC_OK) 
    {
      if(status==PIL_NOT_FOUND) 
	{
	  RILstatus=RILlogMessage(NULL, Warning_1, "Parameter <<inEnergyValues>> not listed.");
	  RILstatus=RILlogMessage(NULL, Warning_1, 
				  "This parameter is mandatory must be a valid file or NULL.");
	  RILstatus=RILlogMessage(NULL, Warning_1, "Please check the parameter file.");
	  return status;
	}
    }
  else if(TmpEnergyValues[0]!='\0')
    {
      RILstatus=RILlogMessage(NULL, Log_1, "Energies FILE                 : %s", TmpEnergyValues);
      if(*detailedOutput)
	RILstatus=RILlogMessage(NULL, Log_1, "... to be read as a list of energy intervals.");
    }
  
  // PARAMETER "GTIs' name"
  InGTIName[0]='\0';
  strcpy(parName, "gti_name");
  if((status=PILGetString(parName, InGTIName))!=ISDC_OK)
    {
      if(status==PIL_NOT_FOUND) 
	{
	  RILstatus=RILlogMessage(NULL, Warning_1, "Parameter <<InGTIName>> not listed.");
	  RILstatus=RILlogMessage(NULL, Warning_1, 
				  "This parameter is mandatory must be a valid file or NULL.");
	  RILstatus=RILlogMessage(NULL, Warning_1, "Please check the parameter file.");
	  return status;
	}
    }
  
  InEFFCDOL[0]='\0';
  strcpy(parName, "inEFFC");
  if((status=PILGetString(parName, InEFFCDOL))!=ISDC_OK)
    {
      if(status==PIL_NOT_FOUND) 
	{
	  RILstatus=RILlogMessage(NULL, Warning_1, "Parameter <<InEFFCDOL>> not listed.");
	  RILstatus=RILlogMessage(NULL, Warning_1, 
				  "This parameter is mandatory must be a valid file or NULL.");
	  RILstatus=RILlogMessage(NULL, Warning_1, "Please check the parameter file.");
	  return status;
	}
    }
    else
    {
        RILstatus=RILlogMessage(NULL, Log_1, "DOL for the ISGR-EFFC-MOD",InEFFCDOL);
    }
  
  // PARAMETER "User-defined Row Filter"
  UserRowFilter[0]='\0';
  strcpy(parName, "isgri_row_filter");
  if((status=PILGetString(parName, UserRowFilter))!=ISDC_OK)
    {
      if(status==PIL_NOT_FOUND) 
	{
	  RILstatus=RILlogMessage(NULL, Warning_1, "Parameter <<UserRowFilter>> not listed.");
	  RILstatus=RILlogMessage(NULL, Warning_1, 
				  "This parameter is mandatory must be a valid file or NULL.");
	  RILstatus=RILlogMessage(NULL, Warning_1, "Please check the parameter file.");
	  return status;
	}
    }
  
  // PARAMETER "Number of Energy Bands"
  strcpy(parName, "isgri_e_num");
  if((status=PILGetInt(parName, NumImaBin))!=ISDC_OK)
    {
      if(status==PIL_NOT_FOUND) 
	{
	  *NumImaBin=0;
	  RILstatus=RILlogMessage(NULL, Warning_1, "Parameter <<isgri_e_num>> not listed.");
	  RILstatus=RILlogMessage(NULL, Warning_1, "This parameter is mandatory >= 1.");
	  RILstatus=RILlogMessage(NULL, Warning_1, "Please check the parameter file.");
	  return status;
	}
    }

  // Get the energy bin boundaries but first allocate the necessary memory. 
  // If the number of energy intervals does not correspond to the number    
  // of boundaries, an error WILL MOST LIKELY RESULT. On the otherhand, if  
  // NumImaBin=0 then use the IBIS nominal energy range. finally if the     
  // requested number of energy bins is les than zero (0) read the values   
  // from a file.                                                           
  if(*NumImaBin>0 && *NumImaBin<11) 
    {
      if((*EnergyBounds=(float **)calloc(*NumImaBin, sizeof(float *)))==NULL) 
	return ERR_ISGR_OSM_MEMORY_ALLOC;
      for(int i=0; i<*NumImaBin; i++) 
	if((*(*EnergyBounds+i)=(float *)calloc(2, sizeof(float)))==NULL) 
	  return ERR_ISGR_OSM_MEMORY_ALLOC;

      if((tmpFloat=(float *)calloc(*NumImaBin, sizeof(float)))==NULL)
	return ERR_ISGR_OSM_MEMORY_ALLOC;
      
      strcpy(parName, "isgri_e_min");
      if((status=PILGetReal4Vector(parName, *NumImaBin, tmpFloat))!=ISDC_OK) 
	return status;
      for(int i=0; i<*NumImaBin; i++) 
	*(*(*EnergyBounds+i))=*(tmpFloat+i);
      
      strcpy(parName, "isgri_e_max");
      if((status=PILGetReal4Vector(parName, *NumImaBin, tmpFloat))!=ISDC_OK)
	return status;
      for(int i=0; i<*NumImaBin; i++)
	*(*(*EnergyBounds+i)+1)=*(tmpFloat+i);

      free(tmpFloat);
    }
  else if(*NumImaBin==0) 
    {
      *NumImaBin=1;
      if((*EnergyBounds=(float **)calloc(*NumImaBin, sizeof(float *)))==NULL) 
	return ERR_ISGR_OSM_MEMORY_ALLOC;
      for(int i=0; i<*NumImaBin; i++) 
	if((*(*EnergyBounds+i)=(float *)calloc(2, sizeof(float)))==NULL) 
	  return ERR_ISGR_OSM_MEMORY_ALLOC;
      **(*EnergyBounds)=1.;
      *(*(*EnergyBounds)+1)=2048.;
    }
  else 
    {
      long tmpL= 0;
      if((status=DALobjectOpen(TmpEnergyValues, &EnergyList, status))!=ISDC_OK) 
	{
	  RILstatus=RILlogMessage(NULL, Error_1, "Error opening %s. Aborting execution.", TmpEnergyValues);
	  return status;
	}
      if((status=DALtableGetNumRows(EnergyList, &tmpL, status))!=ISDC_OK) 
	{
	  status=ERR_ISGR_OSM_FILE_NOTFOUND;
	  RILstatus=RILlogMessage(NULL, Warning_2, "Number of rows in table is zero, status = %d.", status);
	}
      *NumImaBin=(int)tmpL;
      if((*EnergyBounds=(float **)calloc(*NumImaBin, sizeof(float *)))==NULL) 
	return ERR_ISGR_OSM_MEMORY_ALLOC;
      for(int i=0; i<*NumImaBin; i++) 
	if((*(*EnergyBounds+i)=(float *)calloc(2, sizeof(float)))==NULL) 
	  return ERR_ISGR_OSM_MEMORY_ALLOC;
      if((tmpFloat=(float *)calloc(*NumImaBin, sizeof(float)))==NULL) 
	return ERR_ISGR_OSM_MEMORY_ALLOC;
      
      strcpy(colname, E_MIN);
      if((status=DALtableGetCol(EnergyList, colname, 0, &DALtype, &tmpL, tmpFloat, status))!=ISDC_OK) 
	{
	  RILstatus=RILlogMessage(NULL, Error_1, "Fatal : could not retrieve column %s. status = %d.", 
				  colname, status);
	  return status;
	}
      for(int i=0; i<*NumImaBin; i++) 
	*(*(*EnergyBounds+i))=*(tmpFloat+i);
      
      strcpy(colname, E_MAX);
      if((status=DALtableGetCol(EnergyList, colname, 0, &DALtype, &tmpL, tmpFloat, status))!=ISDC_OK)
	{
	  RILstatus=RILlogMessage(NULL, Error_1, "Fatal : could not retrieve column %s. status = %d.", 
				  colname, status);
	  return status;
	}
      for(int i=0; i<*NumImaBin; i++)
	*(*(*EnergyBounds+i)+1)=*(tmpFloat+i);
      free(tmpFloat);
    }
  
  RILstatus=RILlogMessage(NULL, Log_1, "Number of ENERGY BINS         : %d", *NumImaBin);
  if(*NumImaBin>=0)
    for(int bin=0;bin<*NumImaBin;bin++)
       RILstatus=RILlogMessage(NULL, Log_1, "            => ENERGY BIN [%d] = [%.2f,%.2f]keV",
			       bin+1, (*EnergyBounds)[bin][0], (*EnergyBounds)[bin][1]);

  // PARAMETER "time interval's length"  
  strcpy(parName, "isgri_t_len");
  if((status=PILGetReal(parName, TimeLen))!=ISDC_OK) 
    if(status==PIL_NOT_FOUND) 
      {
	*TimeLen=1;
	RILstatus=RILlogMessage(NULL, Log_1, "Parameter <<isgri_t_len>> not listed.");
	RILstatus=RILlogMessage(NULL, Warning_1, "This parameter is mandatory >= 1.");
	RILstatus=RILlogMessage(NULL, Warning_1, "Please check the parameter file.");
	return status;
      }
  if(*TimeLen<ZERO)
    {
      RILstatus=RILlogMessage(NULL, Warning_1, "Parameter isgri_t_len shouldn't be NULL!");
      RILstatus=RILlogMessage(NULL, Warning_1, "Setting it to default value 100000.0 sec.");
      *TimeLen= 100000.;
    }
  RILstatus=RILlogMessage(NULL, Log_1, "TIME BIN Length               : %.3fsec", *TimeLen);
  
  // PARAMETER "Min Rise Time" 
  int i=0;
  strcpy(parName, "isgri_min_rise");
  if((status=PILGetInt(parName, &i))!=ISDC_OK) 
    {
      if(status==PIL_NOT_FOUND) 
	{
	  *MinRiseTime=0;
	  RILstatus=RILlogMessage(NULL, Warning_1, "Parameter <<isgri_min_rise>> not listed.");
	  RILstatus=RILlogMessage(NULL, Warning_1, "This parameter is mandatory.");
	  RILstatus=RILlogMessage(NULL, Warning_1, "Please check the parameter file.");
	  return status;
	}
    }
  else 
    *MinRiseTime = (i>=0 && i<256) ? (unsigned char)i : 0 ;
  
  // PARAMETER "Max Rise Time" 
  strcpy(parName, "isgri_max_rise");
  if((status=PILGetInt(parName, &i))!=ISDC_OK) 
    {
      if(status==PIL_NOT_FOUND) 
	{
	  *NumImaBin=0;
	  RILstatus=RILlogMessage(NULL, Warning_1, "Parameter <<isgri_max_rise>> not listed.");
	  RILstatus=RILlogMessage(NULL, Warning_1, "This parameter is mandatory.");
	  RILstatus=RILlogMessage(NULL, Warning_1, "Please check the parameter file.");
	  return status;
	}
    }
  else 
    *MaxRiseTime = (i>=0 && i<256) ? (unsigned char)i : 255 ;
  RILstatus=RILlogMessage(NULL, Log_1, "Selection on corrected rise times");
  RILstatus=RILlogMessage(NULL, Log_1, "CORRECTED RISE TIMES range           : [%d,%d]",*MinRiseTime,*MaxRiseTime);

 
  // PARAMETER "Spectral Noisy Pixels Detection Flag"
  strcpy(parName, "NoisyDetFlag");
  if((status=PILGetInt(parName, NoisyDetFlag))!=ISDC_OK) 
    if(status==PIL_NOT_FOUND) 
      {
	*detailedOutput=0;
	strcpy(details, "NO");
	RILstatus=RILlogMessage(NULL, Log_1, "Parameter <<NoisyDetFlag>> not listed but required.");
	RILstatus=RILlogMessage(NULL, Log_1, "Program will abort.");
	return status;
      }
  RILstatus=RILlogMessage(NULL, Log_1, "Spectral Noisy Detection FLAG : %d", *NoisyDetFlag);
  
  // PARAMETER "outputlevel" 
  strcpy(parName, "outputLevel");
  if((status=PILGetString(parName, outputLevel))!=ISDC_OK) 
    if(status==PIL_NOT_FOUND) 
      {
	*detailedOutput=0;
	strcpy(details, "NO");
	RILstatus=RILlogMessage(NULL, Log_1, "Parameter <<outputLevel>> not listed but required.");
	RILstatus=RILlogMessage(NULL, Log_1, "Program will abort.");
	return status;
      }
  RILstatus=RILlogMessage(NULL, Log_1, "Output LEVEL                  : %s", outputLevel);
  
  // free temporary arrays 
  free(TmpEnergyValues);
  free(details);
  free(parName);
  free(fileName);
  return status;
}


/////////////////////////////////////////////////////////////////////////////////////////////
//                         CHECK INPUT and OUTPUT REQUIRED FILES.
/////////////////////////////////////////////////////////////////////////////////////////////

int ChkFilesExist(dal_element   *NewGRP,         // DOL to the SWG
		  OBTime        *StartTime,      // Start time of the SCW (OBT)
 		  OBTime        *EndTime,        // End time of the SCW (OBT)
                  dal_element   **DTtable,       // Output: DOL to the table of Dead Times
                  long          *NumDeadTimes,   // Output: Number of Dead Times in this table
                  dal_element   **idxRAWshad,    // Output: DOL to the index of the RAW shadowgrams
                  dal_element   **idxEFFshad,    // Output: DOL to the index of the EFF shadowgrams
                  dal_element   *idxREVcontext,  // DOL to the index of the REV Contexts                                  
                  dal_element   **REVcontext,    // Output: DOL to the REV Context to be used
                  int           *IndexAction,
		  unsigned char detailedOutput) //output
{
  int 
    status=    ISDC_OK, 
    RILstatus= ISDC_OK,
    NumMembers= 0, 
    revolution= 0;
  char 
    *name=      NULL, 
    *SelectRow= NULL;
  long 
    revol= 0;
  dal_element  
    *TmpDOL= NULL;
  dal_dataType 
    DALtype= DAL_CHAR;
  OBTime 
    SCWend= 0;
  // variable for proton table
  long numRow=0;
  char   keyVal[DAL_MAX_STRING]; 

  // Generic allocations
  if((name= (char*)malloc(PIL_LINESIZE*sizeof(char)))==NULL) 
    return ERR_ISGR_OSM_MEMORY_ALLOC;
  if((SelectRow= (char*)malloc(PIL_LINESIZE*sizeof(char)))==NULL) 
    return ERR_ISGR_OSM_MEMORY_ALLOC;

  // check to see if dead times have been previously computed 
  strcpy(name, ISGR_DEAD_STA);
  if((status= DALobjectFindElement(NewGRP, name, DTtable, status))!=ISDC_OK) 
    {
      status= ERR_ISGR_OSM_FILE_NOTFOUND;
      RILstatus= RILlogMessage(NULL, Warning_2, "Error finding table %s error = %d.", name, status);
      RILstatus= RILlogMessage(NULL, Warning_2, "Reverting to status = ISDC_OK.");
      RILstatus= RILlogMessage(NULL, Warning_2, "Dead Times will be set to 0%c.",'%');
      status= ISDC_OK;
    }

  // Number of deadtime values                                    
  // Just make sure they are there, if not, then the dtimes will  
  // be simulated. In this programe, DO NOT tamper with the       
  // dead time values, use ibis_isgr_deadtime to calculate them.  
  if((status= DALtableGetNumRows(*DTtable, NumDeadTimes, status))!=ISDC_OK) 
    {
      status= ERR_ISGR_OSM_FILE_NOTFOUND;
      RILstatus= RILlogMessage(NULL, Warning_2, "Number of rows in deadtime table is zero, status = %d.", status);
      status= ISDC_OK;
    }

  // The index files should normally all exist. If they  
  // do not then they can be created by the program. To  
  // be discussed with the ISDC.                         
  // the index for the shadowgrams should normally exist 
  strcpy(name, ISGR_SHD_IDX_TPL);
  if((status= DALobjectFindElement(NewGRP, name, idxRAWshad, status))!=ISDC_OK) 
    {
      *IndexAction+=1;
      RILstatus= RILlogMessage(NULL, Warning_2, "Error finding image index %s error = %d.", 
			       name, ERR_ISGR_OSM_FILE_NECESSARY_NOTFOUND);
      RILstatus= RILlogMessage(NULL, Warning_2, "Reverting to from status %d = ISDC_OK to continue.", status);
      status= ISDC_OK;
    }
  // the index for the efficiency images is required 
  strcpy(name, ISGR_SHD_EFFI_IDX_TPL);
  if((status=DALobjectFindElement(NewGRP, name, idxEFFshad, status))!=ISDC_OK) 
    {
      *IndexAction+=8;
      RILstatus= RILlogMessage(NULL, Warning_2, "Error finding efficiency image index %s error = %d.", 
			       name, ERR_ISGR_OSM_FILE_NECESSARY_NOTFOUND);
      RILstatus= RILlogMessage(NULL, Warning_2, "Reverting to status = ISDC_OK to continue.");
      status= ISDC_OK;
    }
  if(*IndexAction) 
    status= ERR_ISGR_OSM_FILE_NECESSARY_NOTFOUND;

  // Need the context pixel images for the low threshold values. 
  // SPR 4262: from all the contexts that are indexed, use the most recent who
  // occurred before the END of the SCW.                                                     
  strcpy(name, ISGR_CTXT_GRP);
  if((status= DAL3GENindexGetNumMembers(idxREVcontext, name, &NumMembers, status))!=ISDC_OK) 
    {
     //SPR 05084 -- Let the index be considered as a group, when it is read from the IC tree instead of the rev directory
     //	         -- a possible erroneous file will be treated aferwards
      RILstatus= RILlogMessage(NULL, Warning_1, "No index members for low initial context found. Treating it as a Context group."); 
      RILstatus= RILlogMessage(NULL, Warning_1, "Reverting from %d to ISDC_OK to continue.", status);
      *REVcontext = idxREVcontext;
      status= ISDC_OK;
      NumMembers= 0;
    }

  if(!NumMembers) 
    RILstatus= RILlogMessage(NULL, Warning_1, "Index is empty.");
  else 
    {
      if(detailedOutput)
	RILstatus= RILlogMessage(NULL, Log_1, "Index of contexts contains %d items.", NumMembers);

      // Retrieve SCW parameters necessary to chose the relevant context
      if((status= DALattributeGetInt(NewGRP, "REVOL", &revol, NULL, NULL, status))!=ISDC_OK) 
	{
	  RILstatus= RILlogMessage(NULL, Warning_1, "Error getting revolution number."); 
	  RILstatus= RILlogMessage(NULL, Warning_1, "Setting Revolution number to -1 to use first context found."); 
	  RILstatus= RILlogMessage(NULL, Warning_1, "Reverting from %d to ISDC_OK to continue.", status);
	  revol= -1;
	  status= ISDC_OK;
	}
      if((status= DAL3GENattributeGetOBT(NewGRP, "OBTEND", &SCWend, NULL, status))!=ISDC_OK)
	{
	  RILstatus= RILlogMessage(NULL, Warning_1, "Error getting OBTEND.");
	  RILstatus= RILlogMessage(NULL, Warning_1, "Setting to zero (0)."); 
	  RILstatus= RILlogMessage(NULL, Warning_1, "Reverting from %d to ISDC_OK to continue.", status);
	  SCWend= 0;
	  status= ISDC_OK;
	}

      // Determine relevant context
      // Note: members are supposed to be ranked by time => last read is last in time (REV and OBT)
      int 
	Ctxt,
	GoodCtxt=0;
      long
	LastREV;
      OBTime 
	CTXTobt= 0;
      for(Ctxt=0;Ctxt<NumMembers;Ctxt++)
	{
	  // Check revolution number for this context
	  long k=1;
	  DALtype= DAL_INT;
	  if((status=DALtableGetColBins(idxREVcontext, "REVOL", 0, &DALtype, Ctxt+1, Ctxt+1, &k, &revolution, status))!=ISDC_OK) 
	    {
	      RILstatus= RILlogMessage(NULL, Warning_1, "Error getting revolution number from indexed files."); 
	      RILstatus= RILlogMessage(NULL, Warning_1, "Setting Revolution number to -1 to use first context found."); 
	      RILstatus= RILlogMessage(NULL, Warning_1, "Reverting from %d to ISDC_OK to continue.", status);
	      revolution= -1;
	      status= ISDC_OK;
	    }
	  if(revol<(long)revolution) break;

	  // Check time-coherence with the SCW
	  if((status= DAL3GENtableGetOBTBins(idxREVcontext, "CTXT_OBT", 0, Ctxt+1, Ctxt+1, &k, &CTXTobt, status))!=ISDC_OK) 
	    {
	      RILstatus= RILlogMessage(NULL, Warning_1, "Error getting CTXTstart of member %d/%d",Ctxt,NumMembers); 
	      RILstatus= RILlogMessage(NULL, Warning_1, "Setting Revolution to zero (0)."); 
	      RILstatus= RILlogMessage(NULL, Warning_1, "Reverting from %d to ISDC_OK to continue.", status);
	      CTXTobt= 0;
	      status= ISDC_OK;
	    }
	  if(SCWend>=CTXTobt)
	    {
	      GoodCtxt= Ctxt+1;
	      LastREV= revolution;
	    }
	}

      // Check if a good context has been found for this SCW
      if(GoodCtxt==0) 
	{
	  RILstatus= RILlogMessage(NULL, Warning_1, "Given SCW has an OBTEND before launch."); 
	  RILstatus= RILlogMessage(NULL, Warning_1, "Can be due to NRT problems."); 
	  RILstatus= RILlogMessage(NULL, Warning_1, "Will use last context in index.");
	  GoodCtxt= NumMembers;
	}
      else
	{
	  // SPR4262 special case: warn if use a member from another revolution
	  if(LastREV<revol)
	    if(detailedOutput)
	      {
		RILstatus= RILlogMessage(NULL, Warning_1, "Couldn't find any context for REV %04ld.", revol);
		RILstatus= RILlogMessage(NULL, Warning_1, "Will use last found (from REV %04ld).", LastREV);
	      }
	  // SPR4262 special case: warn if index contains unused members
	  if(Ctxt>GoodCtxt)
	    {
	      if(detailedOutput)
		RILstatus= RILlogMessage(NULL, Warning_1, "There are %d unused contexts before orbit end.", Ctxt-GoodCtxt);

	      // ---- SCREW 1793 
	      // SPR4262 special case: hard coding of special REV problem
	      //if(revol==FAULTY_REV_CTXT)
	      //{
	      //GoodCtxt = Ctxt;
	      //if(detailedOutput)
	      //RILstatus= RILlogMessage(NULL, Warning_1, "Will use last context of REV %04ld", FAULTY_REV_CTXT);
	      //}
	      // --------------
	    }
	}

      // Retrieve this relevant context (member in the index)
      Ctxt= GoodCtxt;
      sprintf(SelectRow, "#row==%d", GoodCtxt);
      if((status=DAL3GENindexFindMember(idxREVcontext, name, SelectRow, &Ctxt, &TmpDOL, status))!=ISDC_OK) 
	{
	  *REVcontext= NULL;
	  RILstatus= RILlogMessage(NULL, Warning_1, "Error getting item %d/%d, status = %d.", GoodCtxt, NumMembers, status); 
	  RILstatus= RILlogMessage(NULL, Warning_1, "Reverting from %d to ISDC_OK to continue.", status);
	  status= ISDC_OK;
	}
      else
	{
	  *REVcontext= TmpDOL;
	  if(detailedOutput)
	    RILstatus= RILlogMessage(NULL, Log_1, "Retrieved item %d/%d", GoodCtxt, NumMembers);
	}
	
    } // If there are members in the index
  
  // Read the time interval from the header of the FITS file   
  if((status= DAL3GENattributeGetOBT(NewGRP, OBTSTART, StartTime, NULL, ISDC_OK)) != ISDC_OK ) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "OBTSTART DAL3GENattributeGetOBT status = %d", status);
      return status;
    }
  if((status= DAL3GENattributeGetOBT(NewGRP, OBTEND, EndTime, NULL, ISDC_OK)) != ISDC_OK) 
    {
      RILstatus= RILlogMessage(NULL, Log_1, "OBTEND DAL3GENattributeGetOBT status = %d", status);
      return status;
    }
  // Free memory and exit
  free(SelectRow);
  free(name);
  return status;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                   Retrieve REVOLUTION context: 
//    Read PIXEL Initial_Switch_Status and call DAL3IBIS routine to read and convert Low Thresholds      
//////////////////////////////////////////////////////////////////////////////////////////////////////////

int GetREVcontext(dal_element   *REVcontext,       // DOL to the REV context
		  int           Revol,             // Revolution number of the SCW
		  OBTime        OBTend,            // End Time of the SCW
		  dal_double    **LowThreshMap,    // Output: Map of Low Thresholds (keV)
		  dal_int       **ONpixelsREVmap,  // Output: Map of Pixels Status for this REV
		  unsigned char detailedOutput)
{
  int	    
    status= ISDC_OK,
    RILstatus= ISDC_OK;

 
  // ===================================================
  //      Retrieve and convert Low Thresholds map
  // ===================================================

  // Allocation
  float *BufferF= NULL;
  if((BufferF= (float*)calloc(ISGRI_SIZE*ISGRI_SIZE, sizeof(float)))==NULL) 
    {
      RILstatus= RILlogMessage(NULL,Error_1,"Error in allocating memory for LowThreshold buffer.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }

  // Use Dal3Ibis library: LT(keV) are set to 0 if LT(step) is dummy or 63 (noisy pixel)
  if((status= DAL3IBISGetlowthresholdKev(REVcontext,OBTend,BufferF,status))!=ISDC_OK)
    {
      RILstatus= RILlogMessage(NULL,Error_1,"Getting Revolution LowThresholds failed, status %d.",status);
      // SCREW 1746: if reading failed, set to average value (depends on Revol)
      double MeanLT= 18.2;
      if(Revol>55) MeanLT= 17.2;
      if(Revol>256) MeanLT= 15.3;
      RILstatus= RILlogMessage(NULL,Warning_1,"All Pixels LowThresholds set to %2.1fkeV.",MeanLT);
      for(int y=0;y<ISGRI_SIZE;y++)
	for(int z=0;z<ISGRI_SIZE;z++) 
	  LowThreshMap[y][z]= MeanLT;
      RILstatus= RILlogMessage(NULL,Warning_1,"Reverting to ISDC_OK.");
      status= ISDC_OK;
    }
  else 
    {
      // Rearrange buffer into matrix
      for(int y=0;y<ISGRI_SIZE;y++)
	{
	for(int z=0;z<ISGRI_SIZE;z++) 
	  {	    
	    LowThreshMap[y][z]= BufferF[z*ISGRI_SIZE+y];
	  }
	}
    }

  // De-allocation
  if(BufferF)
    { free(BufferF); BufferF= NULL; }


  // ===================================================
  //     Retrieve Initial Pixel status on this REV
  // ===================================================

  // Retrieve map
  DAL3_Byte BufferB[ISGRI_SIZE][ISGRI_SIZE];
  if((status= DAL3IBISctxtGetImaPar(REVcontext, &OBTend, ISGRI_PIX_STA, BufferB, status)) !=ISDC_OK )
    {
      RILstatus= RILlogMessage(NULL,Warning_1,"Error finding Pixels Initial Status, error= %d.", 
			       ERR_ISGR_OSM_FILE_NECESSARY_NOTFOUND);
      RILstatus= RILlogMessage(NULL,Warning_1,"All Pixels status initialized to ON.");
      for(int y=0;y<ISGRI_SIZE;y++)
	for(int z=0;z<ISGRI_SIZE;z++) 
	  ONpixelsREVmap[y][z]= 1;
      RILstatus= RILlogMessage(NULL,Warning_1,"Reverting to ISDC_OK to continue.");
      status= ISDC_OK;
    }
  else
    {
      // Rearrange matrix
      for(int y=0;y<ISGRI_SIZE;y++)
	for(int z=0;z<ISGRI_SIZE;z++) 
	  ONpixelsREVmap[y][z]= BufferB[z][y];
    }

  // Info
  if(detailedOutput)
    {
      int
	NumON=0,
	NumBadLT=0;
      for(int y=0;y<ISGRI_SIZE;y++)
	for(int z=0;z<ISGRI_SIZE;z++)
	  {
	    // Pixels ON at initial REV status
	    if(ONpixelsREVmap[y][z])
	      NumON++;
	    // Pixels that switched during the SCW
	    if(LowThreshMap[y][z]<ZERO)
	      NumBadLT++;
	  }
      RILstatus= RILlogMessage(NULL, Log_1, "REV stats: NumPixelsON= %d, NumBadLT= %d.", NumON,NumBadLT);
    }

  // Exit
  return status;
}


/////////////////////////////////////////////////////////////////////////////////////////////
//                                 
/////////////////////////////////////////////////////////////////////////////////////////////

int GetHK3status(dal_element   *idxHK3maps,      // DOL to the HK3 noisy maps
		 OBTime        SCWstart,         // SCW Starting time
		 OBTime        SCWend,           // SCW Finishing tine
		 dal_int       **ONpixelsHK3map, // Output: Map[y][z] of percent of Time ON
                 unsigned char detailedOutput)
{
  int	    
    status= ISDC_OK,
    RILstatus= ISDC_OK;
    
  // Allocations
  dal_int *NumMaps= NULL; 
  if((NumMaps= (dal_int*)calloc(IBIS_NUM_BLOCK, sizeof(dal_int)))==NULL)
    {
      RILstatus= RILlogMessage(NULL,Error_1,"Error in allocating memory for NumMaps.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  dal_double **PixelLive= NULL; 
  if((PixelLive=(dal_double **)calloc(ISGRI_SIZE, sizeof(dal_double*)))==NULL) 
    {
      RILstatus= RILlogMessage(NULL,Error_1,"Error in allocating memory for PixelLive.");
      return ERR_ISGR_OSM_MEMORY_ALLOC;
    }
  for(int i=0; i<ISGRI_SIZE; i++) 
    if((PixelLive[i]=(dal_double *)calloc(ISGRI_SIZE, sizeof(dal_double)))==NULL) 
      {
	RILstatus= RILlogMessage(NULL,Error_1,"Error in allocating memory for PixelLive[].");
	return ERR_ISGR_OSM_MEMORY_ALLOC;
      }

  // Use HK3stuff special function
  if((status= getPixelLive(idxHK3maps, SCWstart, SCWend, 
			   PixelLive, NumMaps, detailedOutput+1, status)) != ISDC_OK)
    {
      RILstatus= RILlogMessage(NULL,Error_1,"Error while retrieving HK3 pixels live time ratios, status= %d.", status);
      return status;
    }

  // Fill ONpixelsHK3map:
  // Note: Only Pixels with LiveTime ratio of 1 are kept alive
  int
    NumON=0,
    NumSwitched=0,
    NumOFF=0;
  for(int y=0;y<ISGRI_SIZE;y++)
    for(int z=0;z<ISGRI_SIZE;z++)
      {
	// Pixels allways ON during the SCW
	if(PixelLive[y][z]>=1.)
	  {
	    NumON++;
	    ONpixelsHK3map[y][z]= 1;
	  }
	// Pixels that switched during the SCW
	else if ( (PixelLive[y][z]>ZERO) && (PixelLive[y][z]<1.) )
	  {
	    NumSwitched++;
	    ONpixelsHK3map[y][z]= 0;
	  }
	// Pixels allways OFF during the SCW
	else
	  {
	    NumOFF++;
	    ONpixelsHK3map[y][z]= 0;
	  }
      }
  if(detailedOutput)
    RILstatus= RILlogMessage(NULL, Log_1, "HK3 stats: NumON= %d, NumSwitched= %d, NumOFF= %d.", 
			     NumON, NumSwitched, NumOFF);

  // Free memory and exit
  for(int y=0;y<ISGRI_SIZE;y++)
    if(PixelLive[y])
      { free(PixelLive[y]); PixelLive[y]= NULL;}
  if(PixelLive)
    { free(PixelLive); PixelLive= NULL;}
  if(NumMaps)
    { free(NumMaps); NumMaps= NULL;}
  return status;
}


/////////////////////////////////////////////////////////////////////////////////////////////
//                                 RETRIEVE HOUSEKEEPING DATA
/////////////////////////////////////////////////////////////////////////////////////////////

int GetHKdata(dal_element   *NewGRP, 
              dal_element   *DTtable,
              char          *InGTIName,
              OBTime        **GTIstart, 
              OBTime        **GTIstop, 
              int           *NumIntervals,
              OBTime        **GTIStartPerMCE,
              OBTime        **GTIStopPerMCE,
              int           *NumGTIperMCE,
              int           *NumGTIXtra,
              dal_float     ***DeadTimes,
              OBTime        **OBTdeadTimes, 
              long          *NumDeadTimes,
              OBTime        StartTime, 
              OBTime        EndTime,
              unsigned char detailedOutput)
{
  int          
    status=    ISDC_OK, 
    RILstatus= ISDC_OK,          
    DTflag= 0;  
  char         
    ChNum[4], 
    colname[32],
    *GTINames= NULL; 
  long         
    IBISnum[5];  
  dal_float    
    *TmpFloatArr= NULL;  
  dal_dataType
    DALtype= DAL_FLOAT; 
  OBTime       
    deltaOBT= 0,
    *GTIStartTmp= NULL,
    *GTIStopTmp=  NULL;
  
  // ===================================================================================  
  //                                    Get Events                                        
  // ===================================================================================  
  
  if((status=DAL3IBISshowAllEvents(NewGRP, IBISnum, status))!=ISDC_OK)
    return status;
  
  // Output info 
  RILstatus=RILlogMessage(NULL, Log_1, "ISGRI Events            : %ld",IBISnum[ISGRI_EVTS]);
  if(detailedOutput) 
    {
      RILstatus=RILlogMessage(NULL, Log_1, "PICsIT Single Events    : %ld",IBISnum[PICSIT_SGLE]);
      RILstatus=RILlogMessage(NULL, Log_1, "PICsIT Multiple Events  : %ld",IBISnum[PICSIT_MULE]);
      RILstatus=RILlogMessage(NULL, Log_1, "COMPTON Single Events   : %ld",IBISnum[COMPTON_SGLE]);
      RILstatus=RILlogMessage(NULL, Log_1, "COMPTON Multiple Events : %ld",IBISnum[COMPTON_MULE]);
    }
  
  // Exit if there are no events 
  if(!IBISnum[ISGRI_EVTS]) {
    RILstatus=RILlogMessage(NULL, Warning_1, "No IBIS events found.");
    RILstatus=RILlogMessage(NULL, Warning_1, "Nothing to be done.");
    RILstatus=RILlogMessage(NULL, Warning_1, "Terminating this process.");
    CommonExit(status);
  }
  
  // ===================================================================================  
  //                      Get DeadTimes between EndTime and StartTime                                     
  // ===================================================================================  
  
  if(*NumDeadTimes<=0) 
    {
      DTflag=1;
      if((status=DAL3GENelapsedOBT(EndTime, StartTime, &deltaOBT, status))!=ISDC_OK) 
	{
	  RILstatus=RILlogMessage(NULL, Log_1, "Error in effecting OBT time differences, status = %d.", status);
	  return status;
	}
      *NumDeadTimes=(deltaOBT%(8*DAL3_OBT_SECOND)) ? Quot(deltaOBT, (8*DAL3_OBT_SECOND))+1 : Quot(deltaOBT, (8*DAL3_OBT_SECOND));
    }
  
  if((*DeadTimes=(dal_float **)calloc(*NumDeadTimes, sizeof(dal_float *)))==NULL) 
    return ERR_ISGR_OSM_MEMORY_ALLOC;
  for(long i=0; i<*NumDeadTimes; i++) 
    if((*(*DeadTimes+i)=(dal_float *)calloc(IBIS_NUM_BLOCK, sizeof(dal_float)))==NULL) 
      return ERR_ISGR_OSM_MEMORY_ALLOC;
  if((*OBTdeadTimes=(OBTime *)calloc(*NumDeadTimes, sizeof(OBTime)))==NULL) 
    return ERR_ISGR_OSM_MEMORY_ALLOC;
  if((TmpFloatArr=(dal_float *)calloc(*NumDeadTimes, sizeof(dal_float)))==NULL) 
    return ERR_ISGR_OSM_MEMORY_ALLOC;
  
  switch(DTflag)
    {
    case 0 : 
      {
	for(int i=0; i<IBIS_NUM_BLOCK; i++) 
	  {
	    sprintf(ChNum, "%d", i);
	    strcpy(colname, "II_DEADTIME_");
	    strcat(colname, ChNum);
	    if((status=DALtableGetCol(DTtable, colname, i+2, &DALtype, NumDeadTimes, TmpFloatArr,status))!=ISDC_OK) 
	      return status;
	    for(long k=0; k<*NumDeadTimes; k++)
	      *(*(*DeadTimes+k)+i)=TmpFloatArr[k]; 
	  }
	if((status=DAL3GENtableGetOBT(DTtable, "OB_TIME", 1, NumDeadTimes, *OBTdeadTimes,status))!=ISDC_OK)
	  return status;
      }
      break;
      // This is the deadtime simulation: Has NOT been tested!             
    case 1 : 
      {
	RILstatus=RILlogMessage(NULL, Log_2, "Simulating deadtimes. DTflag = %d", DTflag);
	for(int i=0; i<IBIS_NUM_BLOCK; i++) 
	  for(long k=0; k<*NumDeadTimes; k++) 
	    *(*(*DeadTimes+k)+i)=0.; 
	for(long k=0; k<*NumDeadTimes; k++) 
	  {
	    if((status=DAL3GENskipOBT(StartTime, (OBTime)((long long)k*DAL3_OBT_SECOND), &deltaOBT, 
				      status))!=ISDC_OK) 
	      {
		RILstatus=RILlogMessage(NULL, Error_1, "Fatal error in computing OBTime + dt, status = %d.", 
					status);
		return status;
	      }
	    *(*OBTdeadTimes+k)=deltaOBT;
	  }
      }
      break;
    } 
  // Free memory
  if(TmpFloatArr)
    {
      free(TmpFloatArr);
      TmpFloatArr= NULL;
    }

  // ===================================================================================  
  //                                    Get GTIs                                        
  // ===================================================================================  

  // Allocations  
  *GTIstart= NULL;
  *GTIstop=  NULL;
  if((GTINames=(char *)calloc(32, sizeof(char)))==NULL) 
    return ERR_ISGR_OSM_MEMORY_ALLOC;

  // Loop on modules
  double
    TmpD= 0.,
    AverageDOBT= 0.;
  for(int i=0; i<IBIS_NUM_BLOCK; i++) 
    {
      // Retrieve number of GTIs for this MCE
      sprintf(GTINames, "ISGRI_MCE%d", i);
      if((status=DAL3HKgetNumGTI(NewGRP, ISGRI, &GTINames, 1, NumGTIperMCE+i, status))!=ISDC_OK) 
	{
	  RILstatus=RILlogMessage(NULL, Warning_1, "Could not retrieve number of nominal GTIs : status %d", status);
	  RILstatus=RILlogMessage(NULL, Warning_1, "Reverting to ISDC_OK ...");
	  status=ISDC_OK;
	}
      
      // Retrieve GTIs for this MCE
      if(*(NumGTIperMCE+i)) 
	{
	  if((*(GTIStartPerMCE+i)=(OBTime *)calloc(*(NumGTIperMCE+i), sizeof(OBTime)))==NULL) 
	    return ERR_ISGR_OSM_MEMORY_ALLOC;
	  if((*(GTIStopPerMCE+i)=(OBTime *)calloc(*(NumGTIperMCE+i), sizeof(OBTime)))==NULL) 
	    return ERR_ISGR_OSM_MEMORY_ALLOC;
	  if((status=DAL3HKgetGTI(NewGRP, ISGRI, &GTINames, 1, StartTime, EndTime, NumGTIperMCE+i, 
				  *(GTIStartPerMCE+i), *(GTIStopPerMCE+i), status))!=ISDC_OK) 
	    {
	      RILstatus=RILlogMessage(NULL, Warning_1, "Could not read nominal GTIs : status %d", status);
	      RILstatus=RILlogMessage(NULL, Warning_1, "Reverting to ISDC_OK ...");
	      status=ISDC_OK;
	      if(*(GTIStartPerMCE+i)) 
		{
		  free(*(GTIStartPerMCE+i));
		  *(GTIStartPerMCE+i)=NULL;
		  RILstatus=RILlogMessage(NULL, Warning_1, "Setting number of GTIs for MCE%d to zero.", i);
		  *(NumGTIperMCE+i)=0;
		}
	    }
	}

      // Details about GTIs for this MCE
      if(detailedOutput) 
	{
	  RILstatus= RILlogMessage(NULL, Log_1, "Number of GTIs for MCE%d= %d:",i,NumGTIperMCE[i]);	
	  AverageDOBT=0.;
	  for(int k=0; k<NumGTIperMCE[i]; k++) 
	    {
	      RILstatus= RILlogMessage(NULL, Log_1, "[ %llu, %llu ]", GTIStartPerMCE[i][k], GTIStopPerMCE[i][k]);
	      if((status= DeltaOBT(GTIStopPerMCE[i][k], GTIStartPerMCE[i][k], &TmpD, status))!=ISDC_OK) 
		{
		  RILstatus= RILlogMessage(NULL, Error_1, "Fatal error in OBT time difference, status = %d.", 
					   status);
		  return status;
		}
	      AverageDOBT+=TmpD;
	    }
	  RILstatus= RILlogMessage(NULL, Log_1, "MCE%d: Total time over GTIs= %.3fsec",i,AverageDOBT);
	}
    }
    
  // Get Number of User-defined GTIs 
  sprintf(GTINames, "ISGRI_XTRA_GTI");
  if((status=DAL3HKgetNumGTI(NewGRP, ISGRI, &GTINames, 1, NumGTIXtra, status))!=ISDC_OK) 
    {
      RILstatus=RILlogMessage(NULL, Warning_1, "Could not retrieve number of user defined GTIs : status %d", status);
      RILstatus=RILlogMessage(NULL, Warning_1, "Reverting to ISDC_OK ...");
      status=ISDC_OK;
    }
  if(*NumGTIXtra) 
    RILstatus=RILlogMessage(NULL, Log_1, "Number of supplementary (user defined) GTIs= %d.", *NumGTIXtra);
  
  // Get Number of inGTIs 
  strcpy(GTINames, InGTIName);
  if((status=DAL3HKgetNumGTI(NewGRP, ISGRI, &GTINames, 1, NumIntervals, status))!=ISDC_OK)
    {
      RILstatus=RILlogMessage(NULL, Warning_1, "Could not retrieve number of %s GTIs : status %d", InGTIName, status);
      RILstatus=RILlogMessage(NULL, Warning_1, "Reverting to ISDC_OK ...");
      *NumIntervals=0;
      status=ISDC_OK;
    }
  else if(detailedOutput) 
    RILstatus=RILlogMessage(NULL, Log_1, "Number of %s GTIs= %d", InGTIName, *NumIntervals);
  if(!*NumIntervals) 
    {
      *NumIntervals=1;
      RILstatus=RILlogMessage(NULL, Warning_1, "No Merged GTIs were defined.");
      RILstatus=RILlogMessage(NULL, Warning_1, "Check the GTIs and rerun.");
      RILstatus=RILlogMessage(NULL, Warning_1, "Terminating process with status = ISDC_OK.");
      CommonExit(ISDC_OK);
    }
  
  // Read all inGTIs 
  if((GTIStartTmp= (OBTime *)realloc(GTIStartTmp, *NumIntervals*sizeof(OBTime)))==NULL) 
    return ERR_ISGR_OSM_MEMORY_ALLOC;
  if((GTIStopTmp= (OBTime *)realloc(GTIStopTmp, *NumIntervals*sizeof(OBTime)))==NULL) 
    return ERR_ISGR_OSM_MEMORY_ALLOC;
  
  if((status= DAL3HKgetGTI(NewGRP, ISGRI, &GTINames, 1, StartTime, EndTime, NumIntervals, 
			   GTIStartTmp, GTIStopTmp, status))!=ISDC_OK) 
    {
      // SPR 3068 : don't use DAL3IBISfindEventGaps 
      // In case reading of inGTIs failed, return bad status to main 
      RILstatus=RILlogMessage(NULL, Warning_1,"Could not get %s GTIs : reading failed with status %d",
			      InGTIName,status);
      return status;
    }
  else 
    {
      free(*GTIstart);
      free(*GTIstop);
      if((*GTIstart=(OBTime *)calloc((*NumIntervals), sizeof(OBTime)))==NULL) 
	return ERR_ISGR_OSM_MEMORY_ALLOC;
      if((*GTIstop=(OBTime *)calloc((*NumIntervals), sizeof(OBTime)))==NULL) 
	return ERR_ISGR_OSM_MEMORY_ALLOC;
      for(long i=0; i<*NumIntervals; i++) {
	*(*GTIstart+i)=GTIStartTmp[i];
	*(*GTIstop+i)=GTIStopTmp[i];
      }
    }
  // Free memory
  if(GTINames) 
    {
      free(GTINames);
      GTINames=NULL;
    }
  if(GTIStartTmp) 
    {
      free(GTIStartTmp);
      GTIStartTmp=NULL;
    }
  if(GTIStopTmp) 
    {
      free(GTIStopTmp);
      GTIStopTmp=NULL;
    }
  
  // Output information 
  RILstatus= RILlogMessage(NULL, Log_1, "Number of ISGRI GTIs= %d:", *NumIntervals);
  for(long i=0; i<*NumIntervals; i++) 
    RILstatus= RILlogMessage(NULL, Log_1, "[ %llu, %llu ]", *(*GTIstart+i), *(*GTIstop+i));
    
  // Exit function 
  return status;
}


/////////////////////////////////////////////////////////////////////////////////////////////
//                                Creator of INDEX members.
/////////////////////////////////////////////////////////////////////////////////////////////

int MkOutputFiles(dal_element   *OSMindex,
                  dal_element   **DOLs,
                  int           NumImaBin,
                  int           NumTime,
                  char          *TemplateName,
                  int           status)
// Creates the required output files with the stipulated naming convention. 
// Attaches the output files to the Science Window Group.                   
// Updates the columns in the index file                                    
{
  if(status==ISDC_OK) 
    for(int k=0; k<NumTime; k++) 
      for(int i=0; i<NumImaBin; i++) 
	if((status=DAL3GENindexCreateMember(OSMindex, TemplateName, NULL, DOLs+i+k*NumImaBin, status))!=ISDC_OK) 
	  return status;

  return status;
}


/////////////////////////////////////////////////////////////////////////////////////////////
//                                    Update the INDEX.
/////////////////////////////////////////////////////////////////////////////////////////////

int UpdateIndex(dal_element   *OSMindex,
                dal_element   **DOLs,
                int           NumImaBin,
                unsigned char detailedOutput,
                int           status)
{
  int RILstatus= ISDC_OK;

  for(int i=0; i<NumImaBin; i++)
    // test sur l'erreur= INDEX_NOT_FOUND (ou MULTIPLE_ ou _DOESNT_MATCH ou _KEY_NOT_FOUND) ou MEMBER_NOT_ATTACHED  
    if((status=DAL3GENindexUpdate(DOLs[i], OSMindex, status))!=ISDC_OK) 
      switch(status) 
	{
	  case DAL3GEN_INDEX_KEY_NOT_FOUND :
	    RILstatus=RILlogMessage(NULL,Warning_1,"Column not found status = %d, reverting to status = ISDC_OK.",status);
	    return(ISDC_OK);
	    break;

	  default :
	    if(detailedOutput) 
	      RILstatus=RILlogMessage(NULL, Error_1,"Error update index : status = %d, aborting \t i = %d.",status,i);
	    // CORRECTION for SPR 3196 : return error code generated by DAL3GENindexUpdate to main.c and exit from main 
	    return status;	    
	    break;
	}
  // Execution error on INTEL => need to return good status 
  return status;
}


/////////////////////////////////////////////////////////////////////////////////////////////
//                                 WRITE keywords to fitsfiles.
/////////////////////////////////////////////////////////////////////////////////////////////

int WriteAttributes(dal_element           *Element,
                    const unsigned long   code,
                    const T_OSM_ATTRIBUTE values,
                    int                   InStatus)
{
 int           
   status=    ISDC_OK,
   RILstatus= ISDC_OK;

 if(InStatus!=ISDC_OK) 
   return(InStatus);
 for(unsigned long i=0; i<8*sizeof(unsigned long); i++) 
   {
     switch(code&(1<<i)) 
       {
       case 0 :
	 // do nothing if zero  
	 continue;
	 break;
       case OBTSTART_NUM :
	 if((status=DAL3GENattributePutOBT(Element, OBTSTART, values.obtstart, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (OBTSTART).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case OBTEND_NUM   :  
	 if((status=DAL3GENattributePutOBT(Element, OBTEND, values.obtstop, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (OBTEND).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case TSTART_NUM   : 
	 if((status=DALattributePutReal(Element, TSTART, values.tstart, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (TSTART).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case TSTOP_NUM    :
	 if((status=DALattributePutReal(Element, TSTOP, values.tstop, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (TSTOP).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case TELAPSE_NUM  :
	 if((status=DALattributePutReal(Element, TELAPSE, values.telapse, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (TELAPSE).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case CHANTYPE_NUM :
	 if((status=DALattributePutChar(Element, CHANTYPE, values.chantype, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing attribute (CHANTYPE).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case CHANMIN_NUM  :
	 if((status=DALattributePutInt(Element, CHANMIN, values.chanmin, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (CHANMIN).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case CHANMAX_NUM  :
	 if((status=DALattributePutInt(Element, CHANMAX, values.chanmax, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (CHANMAX).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case E_MIN_NUM    :
	 if((status=DALattributePutReal(Element, E_MIN, values.e_min, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     " : error writing (E_MIN).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     " : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case E_MAX_NUM    :
	 if((status=DALattributePutReal(Element, E_MAX, values.e_max, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     " : error writing (E_MAX).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     " : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case BANDTYPE_NUM :
	 if((status=DALattributePutChar(Element, BANDTYPE, values.bandtype, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing attribute (BANDTYPE).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case OBT_ACQ_NUM  :
	 if((status=DAL3GENattributePutOBT(Element, OBT_ACQ, values.obt_acq, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing attribute (OBT_ACQ).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case BINTIME_NUM  :
	 if((status=DALattributePutInt(Element, BINTIME, values.bintime, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing attribute (BINTIME).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case NUMRANGE_NUM : 
	 if((status=DALattributePutInt(Element, NUMRANGE, values.numrange, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing attribute (NUMRANGE).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case INT_TIME_NUM :
	 if((status=DALattributePutInt(Element, INT_TIME, values.int_time, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (INT_TIME).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case ERTFIRST_NUM  :
       if((status=DALattributePutChar(Element, ERTFIRST, values.ertfirst, NULL, NULL, status))!=ISDC_OK) 
	 {
	   RILstatus=RILlogMessage(NULL, Warning_1, 
				   "IbisOsmIsgrWriteAttributes : error writing attribute (BANDTYPE).");
	   RILstatus=RILlogMessage(NULL, Warning_1, 
				   "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	   status=ISDC_OK;
	 }
       break;
       case ERTLAST_NUM   :
	 if((status=DALattributePutChar(Element, ERTLAST, values.ertlast, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing attribute (BANDTYPE).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case REVOL_NUM    :
	 if((status=DALattributePutInt(Element, REVOL, values.revol, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (CHANMAX).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case SWID_NUM     :
	 if((status=DALattributePutChar(Element, SWID, values.swid, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing attribute (BANDTYPE).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case SWTYPE_NUM   :
	 if((status=DALattributePutChar(Element, SWTYPE, values.sw_type, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing attribute (BANDTYPE).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case SWBOUND_NUM  :
	 if((status=DALattributePutChar(Element, SWBOUND, values.swbound, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing attribute (BANDTYPE).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case BCPPID_NUM   :
	 if((status=DALattributePutChar(Element, BCPPID , values.bcppid, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing attribute (BANDTYPE).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case PREVSWID_NUM :
	 if((status=DALattributePutChar(Element, PREVSWID, values.prevswid, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing attribute (BANDTYPE).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case OUTPUT_NUM :
	 if((status=DALattributePutChar(Element, OUTPUTLEV, values.outputlevel, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing attribute (OUTPUT LEVEL).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case TFIRST_NUM   : 
	 if((status=DALattributePutReal(Element, TFIRST, values.tfirst, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (TFIRST).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case TLAST_NUM    :
	 if((status=DALattributePutReal(Element, TLAST, values.tlast, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (TLAST).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case ONTIME_NUM    :
	 if((status=DALattributePutReal(Element, ONTIME, values.ontime, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (ONTIME).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case LIVETIME_NUM :
	 if((status=DALattributePutReal(Element, LIVETIME, values.livetime, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (LIVETIME).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case EXPOSURE_NUM :
	 if((status=DALattributePutReal(Element, EXPOSURE, values.exposure, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (EXPOSURE).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case DEADC_NUM   :
	 if((status=DALattributePutReal(Element, DEADC, values.deadc, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (DEADC).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case RISE_MIN_NUM   :
	 if((status=DALattributePutInt(Element, RISE_MIN, values.rise_min, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (RISE_MIN).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
       case RISE_MAX_NUM   :
	 if((status=DALattributePutInt(Element, RISE_MAX, values.rise_max, NULL, NULL, status))!=ISDC_OK) 
	   {
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : error writing (RISE_MAX).");
	     RILstatus=RILlogMessage(NULL, Warning_1, 
				     "IbisOsmIsgrWriteAttributes : reverting from %d  to ISDC_OK to continue.", status);
	     status=ISDC_OK;
	   }
	 break;
	 
       default : 
	 RILstatus=RILlogMessage(NULL, Warning_1, 
				 "IbisOsmIsgrWriteAttributes : error writing this value (%lx).", code&(1<<i));
	 continue;
	 break;
       }
   }
 return status;
}


/////////////////////////////////////////////////////////////////////////////////////////////
//                              Abort the program if necessary 
/////////////////////////////////////////////////////////////////////////////////////////////

void Abort(dal_element *NewGRP,
           const char  *msg,
           int         status)
{
  int RILstatus= ISDC_OK;
  RILstatus= RILlogMessage(NULL, Log_1, msg, status);
  // Deallocate the objects now 
  if(NewGRP) 
    status= CommonCloseSWG(NewGRP,status);
  // Always exit with CommonExit 
  CommonExit(status);
}



///////////////////////////////////////////////////////////////////////////////////////////
//                          Reading proton table
/************************************************************************
 * FUNCTION: ReadProTable
 * DESCRIPTION:
 *  Reads the proton data for at t= my_tStart
 *  Returns ISDC_OK if everything is fine, else returns an error code.
 * ERROR CODES:
 *  DAL error codes
 *
 * PARAMETERS:
 *  pTabPtr  dal_element *    in  ISGR-OFFS-MOD
 *  my_tStart      double*    in  time in ijd
 *  C11_90         double*    in  cumulative proton counter (in const*count) for 11-90 Mev
 *  C40_50         double*    in  cumulative proton counter (in const*count) for 40-50 Mev
 *  C50_70         double*    in  cumulative proton counter (in const*count) for 50-70 Mev
 *  Csup130        double*    in  cumulative proton counter (in const*count) for >130 Mev
 *  Csup39         double*    in  cumulative proton counter (in const*count) for >39 Mev
 *  C20_550        double*    in  cumulative proton counter (in const*count) for 20-550 Mev
 *  Cions          double*    in  cumulative ion counter (in const*count) for 150-185 Mev 
 * RETURN:            int     current status
 ************************************************************************/
////////////////////////////////////////////////////////////////////////////////////////////
int ReadProTable(dal_element *pTabPtr,
		 double      *my_tStart,
		 double      *C11_90,
		 double      *C20_27,
		 double	     *C40_50, 
		 double	     *C50_70, 
		 double	     *Csup130, 
		 double	     *Csup39, 
		 double	     *C20_550, 
		 double	     *Cions)
{

  int    status = ISDC_OK;
  dal_dataType type;

  type=DAL_DOUBLE;
  do {
    /*warning: reading begins at row,col=1*/
    status=DALtableGetCol(pTabPtr,"TIME",0,&type,NULL,(void *)my_tStart,status);
    if (status != ISDC_OK) break;
    status=DALtableGetCol(pTabPtr,"C11_90MEV",0,&type,NULL,(void *)C11_90,status);
    if (status != ISDC_OK) break;
    status=DALtableGetCol(pTabPtr,"C20_27MEV",0,&type,NULL,(void *)C20_27,status);
    if (status != ISDC_OK) break;
    status=DALtableGetCol(pTabPtr,"C40_50MEV",0,&type,NULL,(void *)C40_50,status);
    if (status != ISDC_OK) break;
    status=DALtableGetCol(pTabPtr,"C50_70MEV",0,&type,NULL,(void *)C50_70,status);
    if (status != ISDC_OK) break;
    status=DALtableGetCol(pTabPtr,"CSUP130MEV",0,&type,NULL,(void *)Csup130,status);
    if (status != ISDC_OK) break;
    status=DALtableGetCol(pTabPtr,"CSUP39MEV",0,&type,NULL,(void *)Csup39,status);
    if (status != ISDC_OK) break;
    status=DALtableGetCol(pTabPtr,"C20_550MEV",0,&type,NULL,(void *)C20_550,status);
    if (status != ISDC_OK) break;
    status=DALtableGetCol(pTabPtr,"CION150_185MEV",0,&type,NULL,(void *)Cions,status);
    if (status != ISDC_OK) break;
  }while(0);

  if (status != ISDC_OK) {
    RILlogMessage(NULL, Error_2, "Can not read columns from %s bintable NOT found. Status=%d",
		    DS_IREM_CAL,status);
  }
  return status;
}
