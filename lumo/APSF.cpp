///////////////////////////////////////////////////////////////////////////////////
// File : APSF.cpp
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

#include "APSF.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

APSF::APSF(int res) :
  _res(res)
{
  // make sure kernel is odd dimensions
  if (!(_res % 2)) _res++;
  _kernel = new float[_res * _res];

  _q = 0.999;
  _R = 400.0f;
  _D = 2000.0f;
  _T = 1.001f;
  _sigma = 0.5f;
  
  _maxTerms = 600;
  _I0 = 1.0f;
  _retinaSize = 0.01f;
  _eyeSize = 0.025f;
}

APSF::~APSF()
{
  delete[] _kernel;
}

//////////////////////////////////////////////////////////////////////
// Legendre polymonial
//////////////////////////////////////////////////////////////////////
float APSF::legendreM(int m, float mu)
{
  vector<float> memoized;
  memoized.push_back(1.0f);
  memoized.push_back(mu);
 
  for (int x = 2; x <= m; x++)
  {
    float newMemo = ((2.0f * (float)x - 1.0f) * mu * memoized[x - 1] - 
                             ((float)x - 1.0f) * memoized[x - 2]) / (float)x;
    memoized.push_back(newMemo);
  }

  return memoized[m];
}

//////////////////////////////////////////////////////////////////////
// scattering function at a point
//////////////////////////////////////////////////////////////////////
float APSF::pointAPSF(float mu)
{
  float total = 0.0f;

  for (int m = 0; m < _maxTerms; m++)
    total += (gM(_I0, m) + gM(_I0, m + 1)) * legendreM(m, mu);
    
  return total;
}

//////////////////////////////////////////////////////////////////////
// generate a convolution kernel
//////////////////////////////////////////////////////////////////////
void APSF::generateKernelFast()
{
  float dx = _retinaSize / (float)_res;
  float dy = _retinaSize / (float)_res;
  int halfRes = _res / 2;
  float* oneD = new float[_res];
  
  float max = 0.0f;
  float min = 1000.0f;
  int x,y = halfRes;
  for (x = 0; x < _res; x++)
  {
    // calc angle
    float diffX = (x - halfRes) * dx;
    float diffY = (y - halfRes) * dy;
    float distance = sqrt(diffX * diffX + diffY * diffY);
    if ((distance / _eyeSize) > (_R / _D))
      oneD[x] = 0.0f;
    else
    {
      float i = -distance * distance * _D * _D + _eyeSize * _eyeSize * _R * _R + distance * distance * _R * _R;
      i = _eyeSize * _eyeSize * _D - _eyeSize * sqrt(i);
      i /= _eyeSize * _eyeSize + distance * distance;
      float mu = M_PI - atan(_retinaSize / distance) - asin((_D - i) / _R);
      oneD[x] = pointAPSF(cos(mu));
    
      min = (oneD[x] < min) ? oneD[x] : min;
    }
    max = (oneD[x] > max) ? oneD[x] : max;
  }
  
  // floor 
  if (min > 0.0f)
  {
    for (int i = 0; i < _res; i++)
      if (oneD[i] > 0.0f)
        oneD[i] -= min;
    max -= min;
  }
  
  // normalize
  if (max > 1.0f)
  {
    float maxInv = 1.0f / max;
    for (int i = 0; i < _res; i++)
      oneD[i] *= maxInv;
  }

  // interpolate the kernel
  int index = 0;
  for (y = 0; y < _res; y++)
    for (x = 0; x < _res; x++, index++)
    {
      float dx = fabs((float)(x - halfRes));      
      float dy = fabs((float)(y - halfRes));
      float magnitude = sqrtf(dx * dx + dy * dy);
      
      int lower = floor(magnitude);
      if (lower > halfRes - 1)
      {
        _kernel[index] = 0.0f;
        continue;
      }
      float lerp = magnitude - lower;
      _kernel[index] = (1.0f - lerp) * oneD[halfRes + lower] +
                               lerp  * oneD[halfRes + lower + 1];
      
    }

  delete[] oneD;
}

//////////////////////////////////////////////////////////////////////
// save the kernel in binary
//////////////////////////////////////////////////////////////////////
void APSF::write(const char* filename)
{
  // open file
  FILE* file;
  file = fopen(filename, "wb");

  fwrite((void*)&_res,        sizeof(int), 1, file);
  fwrite((void*)&_q,          sizeof(float), 1, file);
  fwrite((void*)&_T,          sizeof(float), 1, file);
  fwrite((void*)&_I0,         sizeof(float), 1, file);
  fwrite((void*)&_sigma,      sizeof(float), 1, file);
  fwrite((void*)&_R,          sizeof(float), 1, file);
  fwrite((void*)&_D,          sizeof(float), 1, file);
  fwrite((void*)&_retinaSize, sizeof(float), 1, file);
  fwrite((void*)&_eyeSize,    sizeof(float), 1, file);
  fwrite((void*)&_maxTerms,   sizeof(int), 1, file);
  fwrite((void*)&_kernel,     sizeof(float) * _res * _res, 1, file);

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// load a binary kernel
//////////////////////////////////////////////////////////////////////
void APSF::read(const char* filename)
{
  // open file
  FILE* file;
  file = fopen(filename, "rb");

  if (_kernel) delete[] _kernel;

  fread((void*)&_res,        sizeof(int), 1, file);
  fread((void*)&_q,          sizeof(float), 1, file);
  fread((void*)&_T,          sizeof(float), 1, file);
  fread((void*)&_I0,         sizeof(float), 1, file);
  fread((void*)&_sigma,      sizeof(float), 1, file);
  fread((void*)&_R,          sizeof(float), 1, file);
  fread((void*)&_D,          sizeof(float), 1, file);
  fread((void*)&_retinaSize, sizeof(float), 1, file);
  fread((void*)&_eyeSize,    sizeof(float), 1, file);
  fread((void*)&_maxTerms,   sizeof(int), 1, file);

  _kernel = new float[_res * _res];
  fread((void*)&_kernel,     sizeof(float) * _res * _res, 1, file);

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// write the kernel to a PPM file
//////////////////////////////////////////////////////////////////////
void APSF::writePPM(const char* filename)
{
  unsigned char* ppm = new unsigned char[3 * _res * _res];

  for (int x = 0; x < _res * _res; x++)
  {
    ppm[3 * x]     = 255 * _kernel[x];
    ppm[3 * x + 1] = 255 * _kernel[x];
    ppm[3 * x + 2] = 255 * _kernel[x];
  }
  WritePPM(filename, ppm, _res, _res);
  
  delete[] ppm;
}
