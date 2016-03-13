///////////////////////////////////////////////////////////////////////////////////
// File : QUAD_DBM_2D.cpp
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

#include "QUAD_DBM_2D.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

QUAD_DBM_2D::QUAD_DBM_2D(int xRes, int yRes, int iterations) :
  _xRes(xRes),
  _yRes(yRes),
  _bottomHit(0),
  _iterations(iterations),
  _quadPoisson(NULL),
  _dag(NULL),
  _skips(10),
  _twister(123456)
{
  allocate();
  _dag = new DAG(_xRes, _yRes);
  
  // calculate dimensions
  _dx = 1.0f / (float)_xRes;
  _dy = 1.0f / (float)_yRes;
  if (_dx < _dy) _dy = _dx;
  else _dx = _dy;
  
  _maxRes = _xRes * _yRes;
}

QUAD_DBM_2D::~QUAD_DBM_2D()
{
  deallocate();
}

void QUAD_DBM_2D::allocate()
{
  _quadPoisson = new QUAD_POISSON(_xRes, _yRes, _iterations);
  _xRes = _yRes = _quadPoisson->maxRes();
}

void QUAD_DBM_2D::deallocate()
{
  if (_dag)         delete _dag;
  if (_quadPoisson) delete _quadPoisson;
}

//////////////////////////////////////////////////////////////////////
// check neighbors for any candidate nodes
//////////////////////////////////////////////////////////////////////
void QUAD_DBM_2D::checkForCandidates(CELL* cell)
{
  int maxDepth = _quadPoisson->maxDepth();
  
  CELL* north = cell->northNeighbor();
  if (north) {
    if (north->depth == maxDepth) {
      if (!north->candidate) {
        _candidates.push_back(north);
        north->candidate = true;
      }

      CELL* northeast = north->eastNeighbor();
      if (northeast && !northeast->candidate) {
        _candidates.push_back(northeast);
        northeast->candidate = true;
      }
      CELL* northwest = north->westNeighbor();
      if (northwest && !northwest->candidate) {
        _candidates.push_back(northwest);
        northwest->candidate = true;
      }
    }
  }

  CELL* east = cell->eastNeighbor();
  if (east && !east->candidate) {
    _candidates.push_back(east);
    east->candidate = true;
  }
  
  CELL* south = cell->southNeighbor();
  if (south) {
    if (!south->candidate) {
      _candidates.push_back(south);
      south->candidate = true;
    }

    CELL* southeast = south->eastNeighbor();
    if (southeast && !southeast->candidate) {
      _candidates.push_back(southeast);
      southeast->candidate = true;
    }

    CELL* southwest = south->westNeighbor();
    if (southwest && !southwest->candidate) {
      _candidates.push_back(southwest);
      southwest->candidate = true;
    }
  }

  CELL* west = cell->westNeighbor();
  if (west && !west->candidate) {
    _candidates.push_back(west);
    west->candidate = true;
  }
}

//////////////////////////////////////////////////////////////////////
// add particle to the aggregate
//////////////////////////////////////////////////////////////////////
bool QUAD_DBM_2D::addParticle()
{
  static float invSqrtTwo = 1.0f / sqrt(2.0f);
  static int totalParticles = 0;
  static int skipSolve = 0;
 
  // compute the potential
  int iterations = 0;
  if (!skipSolve)
    iterations = _quadPoisson->solve();
  skipSolve++;
  if (skipSolve == _skips) skipSolve = 0;

  // construct probability distribution
  vector<float> probabilities;
  float totalPotential = 0.0f;
  for (int x = 0; x < _candidates.size(); x++)
  {
    if (_candidates[x]->candidate)
    {
      probabilities.push_back(_candidates[x]->potential);
      totalPotential += _candidates[x]->potential;
    }
    else
      probabilities.push_back(0.0f);
  }
  
  // get all the candidates
  // if none are left, stop
  if (_candidates.size() == 0) {
    return false;
  }
 
  // if there is not enough potential, go Brownian
  int toAddIndex = 0;
  if (totalPotential < 1e-8)
    toAddIndex = _candidates.size() * _twister.getDoubleLR();
  // else follow DBM algorithm
  else
  {
    // add a neighbor
    float random = _twister.getDoubleLR();
    float invTotalPotential = 1.0f / totalPotential;
    float potentialSeen = probabilities[0] * invTotalPotential;
    while ((potentialSeen < random) && (toAddIndex < _candidates.size()))
    {
      toAddIndex++;
      potentialSeen += probabilities[toAddIndex] * invTotalPotential;
    }
  }
  _candidates[toAddIndex]->boundary = true;
  _candidates[toAddIndex]->potential = 0.0f;
  _candidates[toAddIndex]->state = NEGATIVE;

  CELL* neighbor = NULL;
  CELL* added = _candidates[toAddIndex];
  CELL* north = added->northNeighbor();
  if (north)
  {
    if (north->state == NEGATIVE)
      neighbor = north;
    CELL* northeast = north->eastNeighbor();
    if (northeast && northeast->state == NEGATIVE)
      neighbor = northeast;
    CELL* northwest = north->westNeighbor();
    if (northwest && northwest->state == NEGATIVE)
      neighbor = northwest;
  }
  CELL* east = added->eastNeighbor();
  if (east && east->state == NEGATIVE)
    neighbor = east;
  
  CELL* south = added->southNeighbor();
  if (south)
  {
    if (south->state == NEGATIVE)
      neighbor = south;
    CELL* southeast = south->eastNeighbor();
    if (southeast && southeast->state == NEGATIVE)
      neighbor = southeast;
    CELL* southwest = south->westNeighbor();
    if (southwest && southwest->state == NEGATIVE)
      neighbor = southwest;
  }
  CELL* west = added->westNeighbor();
  if (west && west->state == NEGATIVE)
    neighbor = west;
  
  // insert it as a node for bookkeeping
  _quadPoisson->insert(added->center[0], added->center[1]);
  checkForCandidates(added);

  // insert into the dag
  int newIndex = (int)(added->center[0] * _xRes) +
                 (int)(added->center[1] * _yRes) * _xRes;
  int neighborIndex = (int)(neighbor->center[0] * _xRes) +
                      (int)(neighbor->center[1] * _yRes) * _xRes;
  _dag->addSegment(newIndex, neighborIndex);

  totalParticles++;
  if (!(totalParticles % 200))
    cout << " " << totalParticles;
 
  hitGround(added);
  
  return true;
}

//////////////////////////////////////////////////////////////////////
// hit ground yet?
//////////////////////////////////////////////////////////////////////
bool QUAD_DBM_2D::hitGround(CELL* cell)
{
  if (_bottomHit)
    return true;

  if (!cell)
    return false;

  bool hit = false;
  if (cell->northNeighbor())
  {
    CELL* north = cell->northNeighbor();
    if (north->state == POSITIVE)
      hit = true;
    if (north->eastNeighbor()->state == POSITIVE)
      hit = true;
    if (north->westNeighbor()->state == POSITIVE)
      hit = true;
  }
  if (cell->eastNeighbor())
    if (cell->eastNeighbor()->state == POSITIVE)
      hit = true;
  if (cell->southNeighbor())
  {
    CELL* south = cell->southNeighbor();
    if (south->state == POSITIVE)
      hit = true;
    if (south->eastNeighbor()->state == POSITIVE)
      hit = true;
    if (south->westNeighbor()->state == POSITIVE)
      hit = true;
  }
  if (cell->westNeighbor())
    if (cell->westNeighbor()->state == POSITIVE)
      hit = true;
  
  if (hit)
  {
    _bottomHit = (int)(cell->center[0] * _xRes) +
                 (int)(cell->center[1] * _yRes) * _xRes;
    _dag->buildLeader(_bottomHit);
    return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////
// drawing functions
////////////////////////////////////////////////////////////////////
void QUAD_DBM_2D::draw() {
  glPushMatrix();
  glTranslatef(-0.5, -0.5, 0);
  
  list<CELL*> leaves;
  list<CELL*>::iterator cellIterator = leaves.begin();
  _quadPoisson->getAllLeaves(leaves);
  for (cellIterator = leaves.begin(); cellIterator != leaves.end(); cellIterator++)
  {
    float color = (*cellIterator)->potential;

    if ((*cellIterator)->boundary) {
      if (color <= 0.0f)
        _quadPoisson->drawCell(*cellIterator, 0,0,0);
      else
        _quadPoisson->drawCell(*cellIterator, 0,0,color);
    }
    else
      _quadPoisson->drawCell(*cellIterator, color, 0,0);
  }
  _quadPoisson->draw(NULL);
  
  glPopMatrix();
}

////////////////////////////////////////////////////////////////////
// read in attractors from an image
////////////////////////////////////////////////////////////////////
bool QUAD_DBM_2D::readImage(unsigned char* initial,
                            unsigned char* attractors, 
                            unsigned char* repulsors,
                            unsigned char* terminators, 
                            int xRes, int yRes)
{
  _dag->inputWidth() = xRes;
  _dag->inputHeight() = yRes;
 
  bool initialFound = false;
  bool terminateFound = false;
  
  int index = 0;
  for (int y = 0; y < yRes; y++)
    for (int x = 0; x < xRes; x++, index++)
    {
      // insert initial condition
      if (initial[index])
      {
        // insert something
        CELL* negative = _quadPoisson->insert(x, y);
        negative->boundary  = true;
        negative->potential = 0.0f;
        negative->state = NEGATIVE;
        negative->candidate = true;

        checkForCandidates(negative);

        initialFound = true;
      }
      
      // insert attractors
      if (attractors[index])
      {
        // insert something
        CELL* positive = _quadPoisson->insert(x, y);
        positive->boundary  = true;
        positive->potential = 1.0f;
        positive->state = ATTRACTOR;
        positive->candidate = true;
      }
      
      // insert repulsors
      if (repulsors[index])
        {
          // only insert the repulsor if it is the edge of a repulsor
          bool edge = false;

          if (x != 0)
          {
            if (!repulsors[index - 1]) edge = true;
            if (y != 0         && !repulsors[index - xRes - 1]) edge = true;
            if (y != yRes - 1 && !repulsors[index + xRes - 1]) edge = true;
          }
          if (x != _xRes - 1)
          {
            if (!repulsors[index + 1]) edge = true;
            if (y != 0         && !repulsors[index - xRes + 1]) edge = true;
            if (y != yRes - 1 && !repulsors[index + xRes + 1]) edge = true;
          }
          if (y != 0         && !repulsors[index - xRes]) edge = true;
          if (y != yRes - 1 && !repulsors[index + xRes]) edge = true;
          
          if (edge)
          {
            // insert something
            CELL* negative = _quadPoisson->insert(x, y);
            negative->boundary  = true;
            negative->potential = 0.0f;
            negative->state = REPULSOR;
            negative->candidate = true;
          }
        }

      // insert terminators
      if (terminators[index])
      {
        // insert something
        CELL* positive = _quadPoisson->insert(x, y);
        positive->boundary  = true;
        positive->potential = 1.0f;
        positive->state = POSITIVE;
        positive->candidate = true;

        terminateFound = true;
      }
    }
  
  if (!initialFound) {
    cout << " The lightning does not start anywhere! " << endl;
    return false;
  }
  if (!terminateFound) {
    cout << " The lightning does not end anywhere! " << endl;
    return false;
  }
 
  return true;
}
