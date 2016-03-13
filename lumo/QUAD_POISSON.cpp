///////////////////////////////////////////////////////////////////////////////////
// File : QUAD_POISSON.cpp
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

#include "QUAD_POISSON.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

QUAD_POISSON::QUAD_POISSON(int xRes, int yRes, int iterations) :
  _root(new CELL(1.0f, 1.0f, 0.0f, 0.0f)),
  _noise()
{
  _root->refine();
 
  // figure out the max depth needed
  float xMax = log((float)xRes) / log(2.0f);
  float yMax = log((float)yRes) / log(2.0f);
 
  float max = (xMax > yMax) ? xMax : yMax;
  if (max - floor(max) > 1e-7)
    max = max + 1;
  max = floor(max);
  
  _maxRes = pow(2.0f, (float)max);
  _maxDepth = max;

  // create the blue noise
  _noiseFunc = new BLUE_NOISE(5.0f / (float)_maxRes);
  _noise = new bool[_maxRes * _maxRes];
	_noiseFunc->complete();
	_noiseFunc->maximize();
  _noiseFunc->writeToBool(_noise, _maxRes);

  _solver = new CG_SOLVER(_maxDepth, iterations);
}

QUAD_POISSON::~QUAD_POISSON()
{
  deleteGhosts();
  delete _root;
  delete _solver;
  delete _noiseFunc;
  delete[] _noise;
}

//////////////////////////////////////////////////////////////////////
// draw boundaries to OGL
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::draw(CELL* cell)
{
  // see if it's the root
  if (cell == NULL) {
    draw(_root);
    return;
  }

  // draw the current cell
  glColor4f(1,1,1,0.1);

  glBegin(GL_LINE_STRIP);
    glVertex2f(cell->bounds[1], 1.0f - cell->bounds[0]);
    glVertex2f(cell->bounds[1], 1.0f - cell->bounds[2]);
    glVertex2f(cell->bounds[3], 1.0f - cell->bounds[2]);
    glVertex2f(cell->bounds[3], 1.0f - cell->bounds[0]);
    glVertex2f(cell->bounds[1], 1.0f - cell->bounds[0]);
  glEnd();
  
  // draw the children
  for (int x = 0; x < 4; x++)
    if (cell->children[x] != NULL)
      draw(cell->children[x]);
}

//////////////////////////////////////////////////////////////////////
// draw the given cell
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::drawCell(CELL* cell, float r, float g, float b)
{
  // draw the current cell
  glColor4f(r,g,b,1.0f);
  glBegin(GL_QUADS);
    glVertex2f(cell->bounds[1], 1.0f - cell->bounds[0]);
    glVertex2f(cell->bounds[1], 1.0f - cell->bounds[2]);
    glVertex2f(cell->bounds[3], 1.0f - cell->bounds[2]);
    glVertex2f(cell->bounds[3], 1.0f - cell->bounds[0]);
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// subdivide quadtree to max level for (xPos, yPos)
// returns a pointer to the cell that was created
//////////////////////////////////////////////////////////////////////
CELL* QUAD_POISSON::insert(float xPos, float yPos)
{
  int currentDepth = 0;
  CELL* currentCell = _root;
  bool existed = true;
 
  while (currentDepth < _maxDepth) {
    // find quadrant of current point
    float diff[2];
    diff[0] = xPos - currentCell->center[0];
    diff[1] = yPos - currentCell->center[1];
    int quadrant = 1;
    if (diff[0] > 0.0f) {
      if (diff[1] < 0.0f)
        quadrant = 2;
    }
    else if (diff[1] < 0.0f)
      quadrant = 3;
    else
      quadrant = 0;
    
    // check if it exists
    if (currentCell->children[quadrant] == NULL) {
      existed = false;
      currentCell->refine();
    }
    
    // recurse to next level
    currentCell = currentCell->children[quadrant];

    // increment depth
    currentDepth++;
  }
  // if we had to subdivide to get the cell, add them to the list
  if (!existed)
    for (int i = 0; i < 4; i++)
    {
      _smallestLeaves.push_back(currentCell->parent->children[i]);
      setNoise(currentCell->parent->children[i]);
    }
  
  ///////////////////////////////////////////////////////////////////
  // force orthogonal neighbors to be same depth
  // I have commented the first block, the rest follow the same flow
  ///////////////////////////////////////////////////////////////////

  // see if neighbor exists
  CELL* north = currentCell->northNeighbor();
  if (north && north->depth != _maxDepth) {
    // while the neighbor needs to be refined
    while (north->depth != _maxDepth) {

      // refine it 
      north->refine();
      
      // set to the newly refined neighbor
      north = currentCell->northNeighbor();
    }
    // add newly created nodes to the list
    for (int i = 0; i < 4; i++)
    {
      _smallestLeaves.push_back(north->parent->children[i]);
      setNoise(north->parent->children[i]);
    }
  }
  CELL* south = currentCell->southNeighbor();
  if (south && south->depth != _maxDepth) {
    while (south->depth != _maxDepth) {
      south->refine();
      south = currentCell->southNeighbor();
    }
    for (int i = 0; i < 4; i++)
    {
      _smallestLeaves.push_back(south->parent->children[i]);
      setNoise(south->parent->children[i]);
    }
  }
  CELL* west = currentCell->westNeighbor();
  if (west && west->depth != _maxDepth) {
    while (west->depth != _maxDepth) {
      west->refine();
      west = currentCell->westNeighbor();
    }
    for (int i = 0; i < 4; i++)
    {
      _smallestLeaves.push_back(west->parent->children[i]);
      setNoise(west->parent->children[i]);
    }
  }
  CELL* east = currentCell->eastNeighbor();
  if (east && east->depth != _maxDepth) {
    while (east->depth != _maxDepth) {
      east->refine();
      east = currentCell->eastNeighbor();
    }
    for (int i = 0; i < 4; i++)
    {
      _smallestLeaves.push_back(east->parent->children[i]);
      setNoise(east->parent->children[i]);
    }
  }

  ///////////////////////////////////////////////////////////////////
  // force diagonal neighbors to be same depth
  // The same flow follows as above, except that it makes sure that
  // the 'north' and 'south' neighbors already exist
  ///////////////////////////////////////////////////////////////////
  
  if (north) {
    CELL* northwest = north->westNeighbor();
    if (northwest && northwest->depth != _maxDepth) {
      while (northwest->depth != _maxDepth) {
        northwest->refine();
        northwest = northwest->children[2];
      }
      for (int i = 0; i < 4; i++)
      {
        _smallestLeaves.push_back(northwest->parent->children[i]);
        setNoise(northwest->parent->children[i]);
      }
    }
    CELL* northeast = north->eastNeighbor();
    if (northeast && northeast->depth != _maxDepth) {
      while (northeast->depth != _maxDepth) {
        northeast->refine();
        northeast= northeast->children[3];
      }
      for (int i = 0; i < 4; i++)
      {
        _smallestLeaves.push_back(northeast->parent->children[i]);
        setNoise(northeast->parent->children[i]);
      }
    }
  }
  if (south) {
    CELL* southwest = south->westNeighbor();
    if (southwest && southwest->depth != _maxDepth) {
      while (southwest->depth != _maxDepth) {
        southwest->refine();
        southwest = southwest->children[1];
      }
      for (int i = 0; i < 4; i++)
      {
        _smallestLeaves.push_back(southwest->parent->children[i]);
        setNoise(southwest->parent->children[i]);
      }
    }
    CELL* southeast = south->eastNeighbor();
    if (southeast && southeast->depth != _maxDepth) {
      while (southeast->depth != _maxDepth) {
        southeast->refine();
        southeast= southeast->children[0];
      }
      for (int i = 0; i < 4; i++)
      {
        _smallestLeaves.push_back(southeast->parent->children[i]);
        setNoise(southeast->parent->children[i]);
      }
    }
  }
  
  return currentCell;
}

//////////////////////////////////////////////////////////////////////
// check if a cell hits a noise node
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::setNoise(CELL* cell)
{
  if (!(cell->state == EMPTY))
    return;
  
  int x = cell->center[0] * _maxRes;
  int y = cell->center[1] * _maxRes;

  if (_noise[x + y * _maxRes])
  {
    cell->boundary = true;
    cell->state = ATTRACTOR;
    cell->potential = 0.5f;
    cell->candidate = true;
  }
}

//////////////////////////////////////////////////////////////////////
// insert all leaves into a list
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::getAllLeaves(list<CELL*>& leaves, CELL* currentCell)
{
  // if we're at the root
  if (currentCell == NULL)
  {
    getAllLeaves(leaves, _root);
    return;
  }
  
  // if we're at a leaf, add it to the list
  if (currentCell->children[0] == NULL)
  {
    leaves.push_back(currentCell);
    return;
  }
 
  // if children exist, call recursively
  for (int x = 0; x < 4; x++)
    getAllLeaves(leaves, currentCell->children[x]);
}

//////////////////////////////////////////////////////////////////////
// insert all leaves not on the boundary into a list
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::getEmptyLeaves(list<CELL*>& leaves, CELL* currentCell)
{
  // if we're at the root
  if (currentCell == NULL) {
    getEmptyLeaves(leaves, _root);
    return;
  }
  
  // if we're at a leaf, check if it's a boundary and then
  // add it to the list
  if (currentCell->children[0] == NULL) {
    if (!(currentCell->boundary))
      leaves.push_back(currentCell);
    return;
  }
  
  // if children exist, call recursively
  for (int x = 0; x < 4; x++)
    getEmptyLeaves(leaves, currentCell->children[x]);
}

//////////////////////////////////////////////////////////////////////
// balance the current tree
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::balance()
{
  // collect all the leaf nodes
  list<CELL*> leaves;
  getAllLeaves(leaves);

  // while the list is not empty
  list<CELL*>::iterator cellIterator = leaves.begin();
  for (cellIterator = leaves.begin(); cellIterator != leaves.end(); cellIterator++) {
    CELL* currentCell = *cellIterator;

    // if a north neighbor exists
    CELL* north = currentCell->northNeighbor();
    if (north != NULL)
      // while the neighbor is not balanced
      while (north->depth < currentCell->depth - 1) {
        // refine it
        north->refine();

        // add the newly refined nodes to the list of
        // those to be checked
        for (int x = 0; x < 4; x++)
          leaves.push_back(north->children[x]);

        // set the cell to the newly created one
        north = currentCell->northNeighbor();
      }

    // the rest of the blocks flow the same as above
    CELL* south = currentCell->southNeighbor();
    if (south!= NULL)
      while (south->depth < currentCell->depth - 1) {
        south->refine();
        for (int x = 0; x < 4; x++)
          leaves.push_back(south->children[x]);
        south = currentCell->southNeighbor();
      }
  
    CELL* west = currentCell->westNeighbor();
    if (west != NULL)
      while (west->depth < currentCell->depth - 1) {
        west->refine();
        for (int x = 0; x < 4; x++)
          leaves.push_back(west->children[x]);
        west = currentCell->westNeighbor();
      }

    CELL* east = currentCell->eastNeighbor();
    if (east != NULL)
      while (east->depth < currentCell->depth - 1) {
        east->refine();
        for (int x = 0; x < 4; x++)
          leaves.push_back(east->children[x]);
        east = currentCell->eastNeighbor();
      }
  }
}

//////////////////////////////////////////////////////////////////////
// build the neighbor lists of the current quadtree
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::buildNeighbors()
{
  balance();

  // collect all the leaf nodes
  list<CELL*> leaves;
  getAllLeaves(leaves);

  list<CELL*>::iterator cellIterator = leaves.begin();
  for (cellIterator = leaves.begin(); cellIterator != leaves.end(); cellIterator++)
  {
    CELL* currentCell = *cellIterator;
    
    // build north neighbors
    CELL* north = currentCell->northNeighbor();
    if (north != NULL) {
      if (north->children[0] == NULL) {
        currentCell->neighbors[0] = north;
        currentCell->neighbors[1] = NULL;
      }
      else {
        currentCell->neighbors[0] = north->children[3];
        currentCell->neighbors[1] = north->children[2];
      }
    }
    // else build a ghost cell
    else 
      currentCell->neighbors[0] = new CELL(currentCell->depth);

    // build east neighbors
    CELL* east = currentCell->eastNeighbor();
    if (east != NULL) {
      if (east->children[0] == NULL) {
        currentCell->neighbors[2] = east;
        currentCell->neighbors[3] = NULL;
      }
      else {
        currentCell->neighbors[2] = east->children[0];
        currentCell->neighbors[3] = east->children[3];
      }
    }
    // else build a ghost cell
    else 
      currentCell->neighbors[2] = new CELL(currentCell->depth);

    // build south neighbors
    CELL* south = currentCell->southNeighbor();
    if (south != NULL) {
      if (south->children[0] == NULL) {
        currentCell->neighbors[4] = south;
        currentCell->neighbors[5] = NULL;
      }
      else {
        currentCell->neighbors[4] = south->children[1];
        currentCell->neighbors[5] = south->children[0];
      }
    }
    // else build a ghost cell
    else 
      currentCell->neighbors[4] = new CELL(currentCell->depth);

    // build west neighbors
    CELL* west = currentCell->westNeighbor();
    if (west != NULL) {
      if (west->children[0] == NULL) {
        currentCell->neighbors[6] = west;
        currentCell->neighbors[7] = NULL;
      }
      else {
        currentCell->neighbors[6] = west->children[2];
        currentCell->neighbors[7] = west->children[1];
      }
    }
    // else build a ghost cell
    else 
      currentCell->neighbors[6] = new CELL(currentCell->depth);
  }
}

//////////////////////////////////////////////////////////////////////
// delete ghost cells
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::deleteGhosts(CELL* currentCell)
{
  // if at the base, call on the root
  if (currentCell == NULL)
  {
    deleteGhosts(_root);
    return;
  }

  // if there are children, delete those too
  if (currentCell->children[0]) {
    // call recursively
    for (int x = 0; x < 4; x++)
      deleteGhosts(currentCell->children[x]);
    return;
  }

  // check the neighbors for stuff to delete
  for (int x = 0; x < 8; x++) {
    // if the neighbor exists
    if (currentCell->neighbors[x])
      // and if it is a ghost cell, delete it
      if (currentCell->neighbors[x]->parent == NULL)
        delete currentCell->neighbors[x];
  }
}

//////////////////////////////////////////////////////////////////////
// solve the Poisson problem
//////////////////////////////////////////////////////////////////////
int QUAD_POISSON::solve() {
  // maintain the quadtree
  balance();
  buildNeighbors();

  // retrieve leaves at the lowest level
  _emptyLeaves.clear();
  getEmptyLeaves(_emptyLeaves);
 
  static bool firstSolve = true;
  static int iterations = _solver->iterations();

  // do a full precision solve the first time
  if (firstSolve)
  {
    iterations = _solver->iterations();
    _solver->iterations() = 10000;
    firstSolve = false;
  }
  else
    _solver->iterations() = iterations;
 
  // return the number of iterations
  return _solver->solve(_emptyLeaves);
};


//////////////////////////////////////////////////////////////////////
// get the leafnode that corresponds to the coordinate
//////////////////////////////////////////////////////////////////////
CELL* QUAD_POISSON::getLeaf(float xPos, float yPos)
{
  CELL* currentCell = _root;
  
  while (currentCell->children[0] != NULL)
  {
    // find quadrant of current point
    float diff[2];
    diff[0] = xPos - currentCell->center[0];
    diff[1] = yPos - currentCell->center[1];
    int quadrant = 1;
    if (diff[0] > 0.0f)
    {
      if (diff[1] < 0.0f)
        quadrant = 2;
    }
    else if (diff[1] < 0.0f)
      quadrant = 3;
    else
      quadrant = 0;
    
    // check if it exists
    if (currentCell->children[quadrant] != NULL)
      currentCell = currentCell->children[quadrant];
  }
  return currentCell;
}
