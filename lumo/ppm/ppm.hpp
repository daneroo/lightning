//------------------------------------------------------------------------------
// File : ppm.hpp
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
// ppm.hpp : Portable Pixel Map image format module
//============================================================================

void LoadPPM(const char *FileName, unsigned char* &Color, int &Width, int &Height);
void WritePPM(const char *FileName, unsigned char* Color, int Width, int Height);
