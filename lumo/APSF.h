///////////////////////////////////////////////////////////////////////////////////
// File : APSF.h
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

#ifndef APSF_H
#define APSF_H

#include <cmath>
#include <gl/glut.h>
#include <vector>
#include "ppm\ppm.hpp"

using namespace std;

#ifndef M_PI
#define M_PI 3.1415926535897931f
#endif

////////////////////////////////////////////////////////////////////
/// \brief Generates the rendering filter
////////////////////////////////////////////////////////////////////
class APSF  
{
public:
  //! constructor
	APSF(int res = 512);
  //! destructor
	virtual ~APSF();

  //! read in an APSF file
  void read(const char* filename);
  //! write out an APSF file
  void write(const char* filename);
  //! write out the current kernel to a PPM file
  void writePPM(const char* filename);
  
  //! generate one line of the kernel and spin it radially
  void generateKernelFast();

  //! resolution of current kernel
  int res() { return _res; };
  
  //! returns the float array for the kernel
  float* kernel() { return _kernel; };
  
private:
  //! kernel resolution
  int     _res;
  //! convolution kernel
  float*  _kernel;

  ////////////////////////////////////////////////////////////////////
  // APSF components
  ////////////////////////////////////////////////////////////////////

  // scattering parameters
  float _q;
  float _T;
  float _I0;
  float _sigma;
  float _R;
  float _D;

  float _retinaSize;
  float _eyeSize;
  
  //! number of coefficients
  int _maxTerms;
  
  //! function value at a point
  float pointAPSF(float mu);
  
  ////////////////////////////////////////////////////////////////////
  // auxiliary functions
  ////////////////////////////////////////////////////////////////////
  float legendreM(int m, float mu);
  float gM(float I0, int m) {
    return (m == 0) ? 0.0f : exp(-(betaM(m, _q) * _T + alphaM(m) * log(_T)));
  };
  float alphaM(float m) { 
    return m + 1.0f;
  };
  float betaM(float m, float q) {
    return ((2.0f * m + 1.0f) / m) * (1.0f - pow(q, (int)m - 1));
  };
  float factorial(float x) {
    return (x <= 1.0f) ? 1.0f : x * factorial(x - 1.0f);
  };
  float choose(float x, float y) {
    return factorial(x) / (factorial(y) * factorial(x - y));
  };
};

#endif
