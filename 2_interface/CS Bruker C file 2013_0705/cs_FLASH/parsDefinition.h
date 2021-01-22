/****************************************************************
 *
 * $Source: /pv/CvsTree/pv/gen/src/prg/methods/FLASH/parsDefinition.h,v $
 *
 * Copyright (c) 1999-2003
 * Bruker BioSpin MRI GmbH
 * D-76275 Ettlingen, Germany
 *
 * All Rights Reserved
 *
 * $Id: parsDefinition.h,v 1.14.2.3 2009/08/14 14:25:05 sako Exp $
 *
 ****************************************************************/



/****************************************************************/
/* INCLUDE FILES						*/
/****************************************************************/
double parameter ReadGradRatio;
double parameter SliceGradRatio;
double parameter Phase3dInteg;
double parameter Phase2dInteg;
double parameter EvolutionDuration;
double parameter OneRepTime;
double parameter
{
  display_name "Read Spoiler Duration";
  format "%.2f";
  units "ms";
  relations backbone;
} ReadSpoilerDuration;

double parameter
{
  display_name "Read Spoiler Strength";
  format "%.1f";
  units "%";
  relations backbone;
} ReadSpoilerStrength;

double parameter
{
  display_name "Slice Spoiler Duration";
  format "%.2f";
  units "ms";
  relations backbone;
} SliceSpoilerDuration;

double parameter
{
  display_name "Slice Spoiler Strength";
  format "%.1f";
  units "%";
  relations backbone;
} SliceSpoilerStrength;




PV_PULSE_LIST parameter
{
  display_name "Excitation Pulse Shape";
  relations    ExcPulseEnumRelation;
}ExcPulseEnum;


PVM_RF_PULSE_TYPE parameter
{
  display_name "Excitation Pulse";
  relations    ExcPulseRelation;
}ExcPulse;

int parameter 
{
  display_name "Number of Dummy Scans";
  relations dsRelations;
} NDummyScans;

TE_MODE parameter
{
  display_name "TE Optimisation Mode";
  relations EchoTimeModeRels;
} EchoTimeMode;

double parameter
{
  display_name "Inter Slice Delay";
  relations backbone;
  units "ms";
  format "%.2f";
}SliceSegDur;

double parameter SliceSegDelay;
double parameter MinSliceSegDur;

double parameter SliceSegEndDelay;

double parameter 
{
 display_name "Time for Movie:";
 units "ms";
 format "%.2f";
 relations backbone;
}TimeForMovieFrames; 

/* new parameters for SWI Reconstruction */
RecoMeth_MODE parameter
{
  display_name "Reconstruction Mode";
  relations RecoMethModeRel;
}RecoMethMode;

MASK_MODE parameter
{
  display_name "Weighting Mode";
  relations MaskModeRel;
}WeightingMode;

double parameter
{
  display_name "Mask Weighting";
  relations MaskWeightRange;
  format "%.2f";
}MaskWeighting;

double parameter
{
  display_name "Gauss Broadening";
  relations GaussBroadRange;
  format "%.2f";
  units "mm";
}GaussBroadening;

/* new parameters for Compressive Encoding */
Compressive_Mode parameter
{
	display_name "Compressive Mode";
	relations CompressiveModeRel;
}CompressiveMode;

/****************************************************************/
/*	E N D   O F   F I L E					*/
/****************************************************************/

