//------------------------------------------------------------------------------
// File : ppm.cpp
//------------------------------------------------------------------------------
// GLVU : Copyright 1997 - 2002 
//        The University of North Carolina at Chapel Hill
//------------------------------------------------------------------------------
// Permission to use, copy, modify, distribute and sell this software and its 
// documentation for any purpose is hereby granted without fee, provided that 
// the above copyright notice appear in all copies and that both that copyright 
// notice and this permission notice appear in supporting documentation. 
// Binaries may be compiled with this software without any royalties or 
// restrictions. 
//
// The University of North Carolina at Chapel Hill makes no representations 
// about the suitability of this software for any purpose. It is provided 
// "as is" without express or implied warranty.

//============================================================================
// ppm.cpp : Portable Pixel Map image format module
//============================================================================

#include <stdio.h>

#define PPM_VERBOSE 0

//----------------------------------------------------------------------------
// READS AN IMAGE IN FROM A PPM FILE. RETURNS THE COLOR RGB ARRAY AND DIMENSIONS
// PERFORMS AUTO-ALLOCATION OF Color ARRAY IF SET TO NULL BEFORE CALLING; OTHERWISE
// ASSUMES THAT COLOR HAS BEEN PRE-ALLOCED.
//----------------------------------------------------------------------------
void LoadPPM(const char *FileName, unsigned char* &Color, int &Width, int &Height)
{
  FILE* fp = fopen(FileName, "rb");
  if (fp==NULL) 
    { printf("PPM ERROR (ReadPPM) : unable to open %s!\n",FileName);
      Color=NULL; Width=0; Height=0; return; }
  int c,s;
  do{ do { s=fgetc(fp); } while (s!='\n'); } while ((c=fgetc(fp))=='#');
  ungetc(c,fp);
  fscanf(fp, "%d %d\n255\n", &Width, &Height);
#if PPM_VERBOSE
  printf("Reading %dx%d Texture [%s]. . .\n", Width, Height, FileName);
#endif
  int NumComponents = Width*Height*3;
  if (Color==NULL) Color = new unsigned char[NumComponents];
  fread(Color,NumComponents,1,fp);
  fclose(fp);
}

//----------------------------------------------------------------------------
// Writes an unsigned byte RGB color array out to a PPM file.
//----------------------------------------------------------------------------
void WritePPM(const char *FileName, unsigned char* Color, int Width, int Height)
{
  FILE* fp = fopen(FileName, "wb");
  if (fp==NULL) { printf("PPM ERROR (WritePPM) : unable to open %s!\n",FileName); return; }
  fprintf(fp, "P6\n%d %d\n255\n", Width, Height);
  fwrite(Color,1,Width*Height*3,fp);
  fclose(fp);
}
