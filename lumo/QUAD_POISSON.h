///////////////////////////////////////////////////////////////////////////////////
// File : QUAD_POISSON.h
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

#ifndef QUAD_POISSON_H
#define QUAD_POISSON_H

#include <gl/glut.h>
#include <cstdlib>
#include "CELL.h"
#include <list>
#include "CG_SOLVER.h"
#include "CG_SOLVER_SSE.h"
#include "BlueNoise/BLUE_NOISE.h"

#include <iostream>

using namespace std;

//////////////////////////////////////////////////////////////////////
/// \brief Quadtree Poisson solver
//////////////////////////////////////////////////////////////////////
class QUAD_POISSON  
{
public:
  /// \brief quadtree constructor 
  ///
  /// \param xRes         maximum x resolution
  /// \param yRes         maximum y resolution
  /// \param iterations   maximum conjugate gradient iterations
	QUAD_POISSON(int xRes, 
               int yRes,
               int iterations = 10);
  
  //! destructor
	virtual ~QUAD_POISSON();
 
  /// \brief OpenGL drawing function
  /// 
  /// \param cell         internally used param, should always be NULL externally
  void draw(CELL* cell = NULL);

  /// \brief Draw a single OpenGL cell  
  ///
  /// \param cell         quadtree cell to draw
  /// \param r            red intensity to draw
  /// \param g            green intensity to draw
  /// \param b            blue intensity to draw
  void drawCell(CELL* cell, 
                float r = 1.0f, 
                float g = 0.0f, 
                float b = 0.0f);

  //! Solve the Poisson problem
  int solve();  
 
  /// \brief insert point at maximum subdivision level
  ///
  /// \param xPos         x position to insert at
  /// \param yPos         y position to insert at
  CELL* insert(float xPos, float yPos);

  /// \brief insert point at maximum subdivision level
  /// 
  /// \param xPos         grid cell x index to insert
  /// \param yPos         grid cell y index to insert
  CELL* insert(int xPos, int yPos) { 
    return insert((float)xPos / _maxRes, (float)yPos / _maxRes);
  };
  
  /// \brief get all the leaf nodes
  /// \return Leaves are returned in the 'leaves' param
  void getAllLeaves(list<CELL*>& leaves, CELL* currentCell = NULL);
  
  /// \brief  get all the leaf nodes at finest subdivision level
  /// \return Leaves are  returned in the 'leaves' param
  list<CELL*>& getSmallestLeaves() { return _smallestLeaves; };

  //! maximum resolution accessor
  int& maxRes() { return _maxRes; }

  //! maximum depth accessor
  int& maxDepth() { return _maxDepth; };
  
  //! get leaf at coordinate (x,y)
  CELL* getLeaf(float xPos, float yPos);
  
private:
  //! root of the quadtree
  CELL* _root;

  //! maximum resolution of quadtree
  int _maxRes;

  //! maxmimum depth of quadtree
  int _maxDepth;
  
  //! dependant leaves
  list<CELL*> _emptyLeaves;

  //! smallest leaves
  list<CELL*> _smallestLeaves;
  
  //! current Poisson solver
  CG_SOLVER* _solver;
  
  //! balance quadtree
  void balance();

  //! get the leaf nodes not on the boundary
  void getEmptyLeaves(list<CELL*>& leaves, CELL* currentCell = NULL);
  
  //! build the neighbor lists of the cells
  void buildNeighbors();

  //! delete ghost cells
  void deleteGhosts(CELL* currentCell = NULL);

  //! Blue noise function
  BLUE_NOISE* _noiseFunc;

  //! Blue noise sample locations
  bool* _noise;

  //! check if a cell hits a noise node
  void setNoise(CELL* cell);
};

#endif
