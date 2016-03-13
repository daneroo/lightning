///////////////////////////////////////////////////////////////////////////////////
// File : CG_SOLVER.cpp
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

#include "CG_SOLVER.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CG_SOLVER::CG_SOLVER(int maxDepth, int iterations, int digits) :
  _iterations(iterations), 
  _arraySize(0), _listSize(0), _digits(digits),
  _direction(NULL), _residual(NULL), _q(NULL), _potential(NULL)
{
  // compute the physical size of various grid cells
  _dx = new float[maxDepth + 1];
  _dx[0] = 1.0f;
  for (int x = 1; x <= maxDepth; x++)
    _dx[x] = _dx[x - 1] * 0.5f;
}

CG_SOLVER::~CG_SOLVER()
{
  if (_direction) delete[] _direction;
  if (_residual) delete[] _residual;
  if (_potential) delete[] _potential;
  if (_q) delete[] _q;
}

//////////////////////////////////////////////////////////////////////
// reallocate scratch arrays if necessary
//////////////////////////////////////////////////////////////////////
void CG_SOLVER::reallocate()
{
  // if we have enough size already, return
  if (_arraySize >= _listSize) return;

  // made sure it SSE aligns okay
  _arraySize = _listSize * 2;
  if (_arraySize % 4)
    _arraySize += 4 - _arraySize % 4;
 
  // delete the old ones
  if (_direction) delete[] _direction;
  if (_residual) delete[] _residual;
  if (_q) delete[] _q;

  // allocate the new ones
  _direction = new float[_arraySize];
  _residual = new float[_arraySize];
  _q = new float[_arraySize];

  // wipe the new ones
  for (int x = 0; x < _arraySize; x++)
    _direction[x] = _residual[x] = _q[x] = 0.0f;
}

//////////////////////////////////////////////////////////////////////
// conjugate gradient solver
//////////////////////////////////////////////////////////////////////
int CG_SOLVER::solve(list<CELL*> cells)
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
 
  // compute a new lexicographical order
  cellIterator = cells.begin();
  for (x = 0; x < _listSize; x++, cellIterator++)
    (*cellIterator)->index = x;
  
  // r = b - Ax
  calcResidual(cells);

  // copy residual into easy array
  // d = r
  cellIterator = cells.begin();
  float deltaNew = 0.0f;
  for (x = 0; x < _listSize; x++, cellIterator++)
  {
    _direction[x] = _residual[x];
    deltaNew += _residual[x] * _residual[x];
  }

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

      float neighborSum = 0.0f;
      for (int x = 0; x < 4; x++)
      {
        int j = x * 2;
        neighborSum += _direction[neighbors[j]->index] * currentCell->stencil[j];
        if (neighbors[j+1])
          neighborSum += _direction[neighbors[j+1]->index] * currentCell->stencil[j+1];
      }
      _q[y] = -neighborSum + _direction[y] * currentCell->stencil[8];
    }
    // alpha = deltaNew / (transpose(d) * q)
    float alpha = 0.0f;
    for (x = 0; x < _listSize; x++)
      alpha += _direction[x] * _q[x];
    if (fabs(alpha) > 0.0f)
      alpha = deltaNew / alpha;

    // x = x + alpha * d
    cellIterator = cells.begin();
    for (x = 0; x < _listSize; x++, cellIterator++)
      (*cellIterator)->potential += alpha * _direction[x];

    // r = r - alpha * q
    maxR = 0.0f;
    for (x = 0; x < _listSize; x++)
    {
      _residual[x] -= _q[x] * alpha;
      maxR = (_residual[x] > maxR) ? _residual[x] : maxR;
    }
    
    // deltaOld = deltaNew
    float deltaOld = deltaNew;

    // deltaNew = transpose(r) * r
    deltaNew = 0.0f;
    for (x = 0; x < _listSize; x++)
      deltaNew += _residual[x] * _residual[x];

    // beta = deltaNew / deltaOld
    float beta = deltaNew / deltaOld;

    // d = r + beta * d
    for (x = 0; x < _listSize; x++)
      _direction[x] = _residual[x] + beta * _direction[x];

    // i = i + 1
    i++;
  }

  return i;
}

//////////////////////////////////////////////////////////////////////
// calculate the residuals
//////////////////////////////////////////////////////////////////////
float CG_SOLVER::calcResidual(list<CELL*> cells)
{
  float maxResidual = 0.0f;
  
  list<CELL*>::iterator cellIterator = cells.begin();
  for (int i = 0; i < _listSize; i++, cellIterator++)
  {
    CELL* currentCell = *cellIterator;
    float dx = _dx[currentCell->depth];
    float neighborSum = 0.0f;
    
    for (int x = 0; x < 4; x++)
    {
      int i = x * 2;
      neighborSum += currentCell->neighbors[i]->potential * currentCell->stencil[i];
      if (currentCell->neighbors[i+1])
        neighborSum += currentCell->neighbors[i+1]->potential * currentCell->stencil[i+1];
    }
    _residual[i] = currentCell->b - (-neighborSum + currentCell->potential * currentCell->stencil[8]);
    
    if (fabs(_residual[i]) > maxResidual)
      maxResidual = fabs(_residual[i]);
  }
  return maxResidual;
}

//////////////////////////////////////////////////////////////////////
// compute stencils once and store
//////////////////////////////////////////////////////////////////////
void CG_SOLVER::calcStencils(list<CELL*> cells)
{
  list<CELL*>::iterator cellIterator = cells.begin();
  for (cellIterator = cells.begin(); cellIterator != cells.end(); cellIterator++)
  {
    CELL* currentCell = *cellIterator;
    float invDx = 1.0f / _dx[currentCell->depth];

    // sum over faces
    float deltaSum = 0.0f;
    float bSum = 0.0f;

    for (int x = 0; x < 4; x++)
    {
      int i = x * 2;
      currentCell->stencil[i] = 0.0f;
      currentCell->stencil[i+1] = 0.0f;
      
      if (currentCell->neighbors[i + 1] == NULL) {
        // if it is the same refinement level (case 1)
        if (currentCell->depth == currentCell->neighbors[i]->depth) {
          deltaSum += invDx;
          if (!currentCell->neighbors[i]->boundary)
            currentCell->stencil[i] = invDx;
          else
            bSum += (currentCell->neighbors[i]->potential) * invDx;
        }
        // else it is less refined (case 3)
        else {
          deltaSum += 0.5f * invDx;
          if (!currentCell->neighbors[i]->boundary)
            currentCell->stencil[i] = 0.5f * invDx;
          else
            bSum += currentCell->neighbors[i]->potential * 0.5f * invDx;
        }
      }
      // if the neighbor is at a lower level (case 2)
      else {
        deltaSum += 2.0f * invDx;
        if (!currentCell->neighbors[i]->boundary)
          currentCell->stencil[i] = invDx;
        else
          bSum += currentCell->neighbors[i]->potential * invDx;
        if (!currentCell->neighbors[i+1]->boundary)
          currentCell->stencil[i+1] = invDx;
        else
          bSum += currentCell->neighbors[i+1]->potential * invDx;
      }
    }

    currentCell->stencil[8] = deltaSum;
    currentCell->b = bSum;
  }
}
