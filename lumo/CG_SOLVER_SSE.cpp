///////////////////////////////////////////////////////////////////////////////////
// File : CG_SOLVER_SSE.cpp
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

#include "CG_SOLVER_SSE.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CG_SOLVER_SSE::CG_SOLVER_SSE(int maxDepth, int iterations, int digits) :
  CG_SOLVER(maxDepth, iterations, digits)
{
}

CG_SOLVER_SSE::~CG_SOLVER_SSE()
{
  if (_direction) _aligned_free(_direction);
  if (_potential) _aligned_free(_potential);
  if (_residual)  _aligned_free(_residual);
  if (_q)         _aligned_free(_q);

  _direction = NULL;
  _residual  = NULL;
  _q         = NULL;
  _potential = NULL;
}

//////////////////////////////////////////////////////////////////////
// reallocate the sse arrays
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::reallocate()
{
  if (_arraySize >= _listSize) return;
  _arraySize = _listSize * 2;

  if (_arraySize % 4)
    _arraySize += 4 - _arraySize % 4;
  
  if (_direction) _aligned_free(_direction);
  if (_potential) _aligned_free(_potential);
  if (_residual)  _aligned_free(_residual);
  if (_q)         _aligned_free(_q);

  _direction = (float*) _aligned_malloc(_arraySize * sizeof(float), 16); 
  _potential = (float*) _aligned_malloc(_arraySize * sizeof(float), 16); 
  _residual  = (float*) _aligned_malloc(_arraySize * sizeof(float), 16); 
  _q         = (float*) _aligned_malloc(_arraySize * sizeof(float), 16); 

  return;
}

//////////////////////////////////////////////////////////////////////
// solve the linear system
//////////////////////////////////////////////////////////////////////
int CG_SOLVER_SSE::solve(list<CELL*> cells)
{
  // counters
  int x, y, index;
  list<CELL*>::iterator cellIterator;

  // i = 0
  int i = 0;

  // precalculate stencils
  calcStencils(cells);
 
  // reallocate scratch arrays if necessary
  _listSize = cells.size();
  reallocate();
  wipeSSE(_potential);
  wipeSSE(_direction);
  wipeSSE(_residual);
  wipeSSE(_q);
  
  // compute a new lexicographical order
  cellIterator = cells.begin();
  for (x = 0; x < _listSize; x++, cellIterator++)
  {
    CELL* cell = *cellIterator;
    cell->index = x;
    _potential[x] = cell->potential;
  }

  // r = b - Ax
  calcResidual(cells);

  // d = r
  copySSE(_direction, _residual);
  
  // deltaNew = r^T r
  float deltaNew = dotSSE(_residual, _residual);

  // delta0 = deltaNew
  float delta0 = deltaNew;

  // While deltaNew > (eps^2) * delta0
  float eps  = pow(10.0f, (float)-_digits);
  float maxR = 2.0f * eps;
  while ((i < _iterations) && (maxR > eps))
  {
    // q = Ad
    cellIterator = cells.begin();
    for (y = 0; y < _listSize; y++, cellIterator++)
    {
      CELL* currentCell = *cellIterator;
      CELL** neighbors = currentCell->neighbors;
      float* stencil = currentCell->stencil;
     
      float neighborSum = 0.0f;
      for (int x = 0; x < 8; x++)
      {
        if (neighbors[x])
          neighborSum += _direction[neighbors[x]->index] * stencil[x];
      }
      _q[y] = -neighborSum + _direction[y] * currentCell->stencil[8];
    }

    // alpha = deltaNew / (transpose(d) * q)
    float alpha = dotSSE(_q, _direction);
    if (fabs(alpha) > 0.0f)
      alpha = deltaNew / alpha;

    // x = x + alpha * d
    saxpySSE(alpha, _direction, _potential);

    // r = r - alpha * q
    saxpySSE(-alpha, _q, _residual);
    maxR = maxSSE(_residual);
    
    // deltaOld = deltaNew
    float deltaOld = deltaNew;

    // deltaNew = transpose(r) * r
    deltaNew = dotSSE(_residual, _residual);

    // beta = deltaNew / deltaOld
    float beta = deltaNew / deltaOld;

    // d = r + beta * d
    saypxSSE(beta, _residual, _direction);

    // i = i + 1
    i++;
  }

  // copy back into the tree
  cellIterator = cells.begin();
  for (x = 0; x < _listSize; x++, cellIterator++)
    (*cellIterator)->potential = _potential[x];

  return i;
}

//////////////////////////////////////////////////////////////////////
// dot product of two vectors
//////////////////////////////////////////////////////////////////////
float CG_SOLVER_SSE::dotSSE(float* x, float* y)
{
  __m128 sum = _mm_set_ps1(0.0f);
  __m128* xSSE = (__m128*)x;
  __m128* ySSE = (__m128*)y;
  __m128 temp;
  for (int index = 0; index < _arraySize / 4; index++)
  {
    temp = _mm_mul_ps(*xSSE, *ySSE);
    sum  = _mm_add_ps(sum, temp);
    xSSE++;
    ySSE++;
  }
  union u {
    __m128 m;
    float f[4];
  } extract;
  extract.m = sum;
  return extract.f[0] + extract.f[1] + extract.f[2] + extract.f[3];
}

//////////////////////////////////////////////////////////////////////
// scalar 'a' x + y
// Y = aX + Y
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::saxpySSE(float s, float* x, float* y)
{
  __m128* ySSE = (__m128*)y;
  __m128* xSSE = (__m128*)x;
  __m128 sSSE = _mm_set_ps1(s);
  __m128 temp;
  for (int index = 0; index < _arraySize / 4; index++)
  {
    temp = _mm_mul_ps(*xSSE, sSSE);
    *ySSE = _mm_add_ps(*ySSE, temp);

    xSSE++;
    ySSE++;
  }
}

//////////////////////////////////////////////////////////////////////
// scalar 'a' y + x
// Y = aY + X
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::saypxSSE(float s, float* x, float* y)
{
  __m128* ySSE = (__m128*)y;
  __m128* xSSE = (__m128*)x;
  __m128 sSSE = _mm_set_ps1(s);
  __m128 temp;
  for (int index = 0; index < _arraySize / 4; index++)
  {
    temp = _mm_mul_ps(*ySSE, sSSE);
    *ySSE = _mm_add_ps(*xSSE, temp);

    xSSE++;
    ySSE++;
  }
}

//////////////////////////////////////////////////////////////////////
// scalar 'a' y + x
// Y = aY + X
//////////////////////////////////////////////////////////////////////
float CG_SOLVER_SSE::maxSSE(float* x)
{
  __m128 maxFoundSSE= _mm_set_ps1(0.0f);
  __m128* xSSE = (__m128*)x;
  for (int index = 0; index < _arraySize / 4; index++)
  {
    maxFoundSSE = _mm_max_ps(*xSSE, maxFoundSSE);
    xSSE++;
  }
  union u {
    __m128 m;
    float f[4];
  } extract;
  extract.m = maxFoundSSE;
  float maxFound = extract.f[0] > extract.f[1] ? extract.f[0] : extract.f[1];
  maxFound = maxFound > extract.f[2] ? maxFound : extract.f[2];
  return maxFound > extract.f[3] ? maxFound : extract.f[3];
}

//////////////////////////////////////////////////////////////////////
// SSE add
// Y = X + Y
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::addSSE(float* x, float* y)
{
  __m128* ySSE = (__m128*)y;
  __m128* xSSE = (__m128*)x;
  for (int index = 0; index < _arraySize / 4; index++)
  {
    *ySSE = _mm_add_ps(*ySSE, *xSSE);
    xSSE++;
    ySSE++;
  }
}

//////////////////////////////////////////////////////////////////////
// SSE multiply
// Y = X * Y
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::multiplySSE(float* x, float* y)
{
  __m128* ySSE = (__m128*)y;
  __m128* xSSE = (__m128*)x;
  for (int index = 0; index < _arraySize / 4; index++)
  {
    *ySSE = _mm_mul_ps(*ySSE, *xSSE);
    xSSE++;
    ySSE++;
  }
}

//////////////////////////////////////////////////////////////////////
// SSE multiply
// Z = X * Y
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::multiplySSE(float* x, float* y, float* z)
{
  __m128* zSSE = (__m128*)z;
  __m128* ySSE = (__m128*)y;
  __m128* xSSE = (__m128*)x;
  for (int index = 0; index < _arraySize / 4; index++)
  {
    *zSSE = _mm_mul_ps(*ySSE, *xSSE);
    xSSE++;
    ySSE++;
    zSSE++;
  }
}

//////////////////////////////////////////////////////////////////////
// SSE multiply
// Z = W - X * Y
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::multiplySubtractSSE(float* w, float* x, float* y, float* z)
{
  __m128* zSSE = (__m128*)z;
  __m128* ySSE = (__m128*)y;
  __m128* xSSE = (__m128*)x;
  __m128* wSSE = (__m128*)w;
  for (int index = 0; index < _arraySize / 4; index++)
  {
    *zSSE = _mm_mul_ps(*ySSE, *xSSE);
    *zSSE = _mm_sub_ps(*wSSE, *zSSE);
    
    xSSE++;
    ySSE++;
    zSSE++;
    wSSE++;
  }
}

//////////////////////////////////////////////////////////////////////
// SSE set
// X = val
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::setSSE(float* x, float val)
{
  __m128* xSSE = (__m128*)x;
  for (int index = 0; index < _arraySize / 4; index++)
  {
    *xSSE = _mm_set_ps1(val);
    xSSE++;
  }
}

//////////////////////////////////////////////////////////////////////
// SSE set
// X = 0
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::wipeSSE(float* x)
{
  __m128* xSSE = (__m128*)x;
  for (int index = 0; index < _arraySize / 4; index++)
  {
    *xSSE = _mm_setzero_ps();
    xSSE++;
  }
}

//////////////////////////////////////////////////////////////////////
// SSE set
// X = Y
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::copySSE(float* x, float* y)
{
  __m128* ySSE = (__m128*)y;
  for (int index = 0; index < _arraySize / 4; index++)
  {
    _mm_store_ps(x,*ySSE);    
    x += 4;
    ySSE++;
  }
}
