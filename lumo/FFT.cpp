///////////////////////////////////////////////////////////////////////////////////
// File : FFT.cpp
///////////////////////////////////////////////////////////////////////////////////
//
// LumosQuad - A Lightning Generator
// Copyright 2007
// The University of North Carolina at Chapel Hill
// 
///////////////////////////////////////////////////////////////////////////////////
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  The University of North Carolina at Chapel Hill makes no representations 
//  about the suitability of this software for any purpose. It is provided 
//  "as is" without express or implied warranty.
//
//  Permission to use, copy, modify and distribute this software and its
//  documentation for educational, research and non-profit purposes, without
//  fee, and without a written agreement is hereby granted, provided that the
//  above copyright notice and the following three paragraphs appear in all
//  copies.
//
//  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY WARRANTIES,
//  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
//  FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN
//  "AS IS" BASIS, AND THE UNIVERSITY OF NORTH CAROLINA HAS NO OBLIGATION TO
//  PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
//
//  Please send questions and comments about LumosQuad to kim@cs.unc.edu.
//
///////////////////////////////////////////////////////////////////////////////////
//
//  This program uses OpenEXR, which has the following restrictions:
// 
//  Copyright (c) 2002, Industrial Light & Magic, a division of Lucas
//  Digital Ltd. LLC
// 
//  All rights reserved.
// 
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are
//  met:
//  *       Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//  *       Redistributions in binary form must reproduce the above
//  copyright notice, this list of conditions and the following disclaimer
//  in the documentation and/or other materials provided with the
//  distribution.
//  *       Neither the name of Industrial Light & Magic nor the names of
//  its contributors may be used to endorse or promote products derived
//  from this software without specific prior written permission. 
// 

#include "FFT.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FFT::FFT()
{

}

FFT::~FFT()
{

}

bool FFT::convolve(float* source, float* kernel, int xSource, int ySource, int xKernel, int yKernel)
{
  int x, y, index;
  
  // get normalization params
  float maxCurrent = 0.0f;
  for (x = 0; x < xSource * ySource; x++)
    maxCurrent = (maxCurrent < source[x]) ? source[x] : maxCurrent;
  float maxKernel = 0.0f;
  for (x = 0; x < xKernel * yKernel; x++)
    maxKernel = (maxKernel < kernel[x]) ? kernel[x] : maxKernel;
  float maxProduct = maxCurrent * maxKernel;
 
  // retrieve dimensions
  int xHalf = xKernel / 2;
  int yHalf = yKernel / 2;
  int xResPadded = xSource + xKernel;
  int yResPadded = ySource + yKernel;

  if (xResPadded != yResPadded)
    (xResPadded > yResPadded) ? yResPadded = xResPadded : xResPadded = yResPadded;

  // create padded field
  fftw_complex* padded = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * xResPadded * yResPadded);
  if (!padded)
  {
    cout << " IMAGE: Not enough memory! Try a smaller final image size." << endl;
    return false;
  }
 
  // init padded field
  for (index = 0; index < xResPadded * yResPadded; index++)
    padded[index][0] = padded[index][1] = 0.0f;
  index = 0;
  for (y = 0; y < ySource; y++)
    for (x = 0; x < xSource; x++, index++)
    {
      int paddedIndex = (x + xKernel / 2) + (y + yKernel / 2) * xResPadded;
      padded[paddedIndex][0] = source[index];
    }

  // create padded filter
  fftw_complex* filter = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * xResPadded * yResPadded);
  if (!filter)
  {
    cout << " FILTER: Not enough memory! Try a smaller final image size." << endl;
    return false;
  }

  // init padded filter
  for (index = 0; index < xResPadded * yResPadded; index++)
    filter[index][0] = filter[index][1] = 0.0f;
  
  // quadrant IV
  for (y = 0; y < (yHalf + 1); y++)
    for (x = 0; x < (xHalf + 1); x++)
    {
      int filterIndex = x + xHalf + y * xKernel;
      int fieldIndex = x + (y + yResPadded - (yHalf + 1)) * xResPadded;
      filter[fieldIndex][0] = kernel[filterIndex];
    }

  // quadrant I
  for (y = 0; y < yHalf; y++)
    for (x = 0; x < (xHalf + 1); x++)
    {
      int filterIndex = (x + xHalf) + (y + yHalf + 1) * xKernel;
      int fieldIndex = x + y * xResPadded;
      filter[fieldIndex][0] = filter[fieldIndex][1] = kernel[filterIndex];
    }
  
  // quadrant III
  for (y = 0; y < (yHalf + 1); y++)
    for (x = 0; x < xHalf; x++)
    {
      int filterIndex = x + y * xKernel;
      int fieldIndex = (x + xResPadded - xHalf) + (y + yResPadded - (yHalf + 1)) * xResPadded;
      filter[fieldIndex][0] = filter[fieldIndex][1] = kernel[filterIndex];
    }

  // quadrant II
  for (y = 0; y < yHalf; y++)
    for (x = 0; x < xHalf; x++)
    {
      int filterIndex = x + (y + yHalf + 1) * xKernel;
      int fieldIndex = (x + xResPadded - xHalf) + y * xResPadded;
      filter[fieldIndex][0] = filter[fieldIndex][1] = kernel[filterIndex];
    }
  
  // perform forward FFT on field
  fftw_complex* paddedTransformed = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * xResPadded * yResPadded);
  if (!paddedTransformed)
  {
    cout << " T-IMAGE: Not enough memory! Try a smaller final image size." << endl;
    return false;
  }
  fftw_plan forwardField = fftw_plan_dft_2d(xResPadded, yResPadded, padded, paddedTransformed, FFTW_FORWARD, FFTW_ESTIMATE);  
  fftw_execute(forwardField);

  // perform forward FFT on filter
  fftw_complex* filterTransformed = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * xResPadded * yResPadded);
  if (!filterTransformed)
  {
    cout << " T-FILTER: Not enough memory! Try a smaller final image size." << endl;
    return false;
  }
  fftw_plan forwardFilter = fftw_plan_dft_2d(xResPadded, yResPadded, filter, filterTransformed, FFTW_FORWARD, FFTW_ESTIMATE);  
  fftw_execute(forwardFilter);

  // apply frequency space filter
  for (index = 0; index < xResPadded * yResPadded; index++)
  {
    float newReal = paddedTransformed[index][0] * filterTransformed[index][0] - 
                    paddedTransformed[index][1] * filterTransformed[index][1];
    float newIm   = paddedTransformed[index][0] * filterTransformed[index][1] + 
                    paddedTransformed[index][1] * filterTransformed[index][0];
    paddedTransformed[index][0] = newReal;
    paddedTransformed[index][1] = newIm;
  }
   
  // transform back
  fftw_plan backwardField = fftw_plan_dft_2d(xResPadded, yResPadded, paddedTransformed, padded, FFTW_BACKWARD, FFTW_ESTIMATE);  
  fftw_execute(backwardField);
  
  // copy back into padded
  index = 0;
  for (y = 0; y < ySource; y++)
    for (x = 0; x < xSource; x++, index++)
    {
      int paddedIndex = (x + xKernel / 2) + (y + yKernel / 2) * xResPadded;
      source[index] = padded[paddedIndex][0];
    }

  // clean up
  fftw_free(padded);
  fftw_free(paddedTransformed);
  fftw_free(filter);
  fftw_free(filterTransformed);

  // if normalization is exceeded, renormalize
  float newMax = 0.0f;
  for (x = 0; x < xSource * ySource; x++)
    newMax = (newMax < source[x]) ? source[x] : newMax;
  if (newMax > maxProduct)
  {
    float scale = maxProduct / newMax;
    for (x = 0; x < xSource * ySource; x++)
      source[x] *= scale;
  }

  return true;
}
