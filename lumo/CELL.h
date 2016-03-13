///////////////////////////////////////////////////////////////////////////////////
// File : CELL.h
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

#ifndef CELL_H
#define CELL_H
#include <cstdlib>

//////////////////////////////////////////////////////////////////////
/// \enum Possible states of the cell in the DBM simulation
//////////////////////////////////////////////////////////////////////
enum CELL_STATE {EMPTY, NEGATIVE, POSITIVE, REPULSOR, ATTRACTOR};

//////////////////////////////////////////////////////////////////////
/// \brief Basic cell data structure of the quadtree
//////////////////////////////////////////////////////////////////////
class CELL  
{
public:
  //! normal cell constructor  
  CELL(float north, 
       float east, 
       float south, 
       float west, 
       CELL* parent = NULL, 
       int depth = 0);
  
  //! ghost cell constructor  
  CELL(int depth = 0);

  //! destructor
  ~CELL();
  
  //! The children of the node in the quadtree
  /*! 
      Winding order of children is:

      \verbatim
        _________
        |   |   |
        | 0 | 1 |  
        |___|___|
        |   |   |
        | 3 | 2 |   
        |___|___|
      \endverbatim */
  CELL* children[4]; 
  
  //! The physical bounds of the current grid cell
  /*! 
      Winding order of bounds is:

      \verbatim
        0 - north
        1 - east
        2 - south
        3 - west
      \endverbatim */
  float bounds[4];

  //! The neighbors in the balanced quadtree
  /*!
    winding order of the neighbors is:
  
    \verbatim
         | 0  | 1  |
     ____|____|____|_____
         |         |
       7 |         |  2
     ____|         |_____
         |         |
       6 |         |  3
     ____|_________|_____
         |    |    |
         | 5  |  4 |
    \endverbatim
 
    Neighbors 0,2,4,6 should always exist. Depending on
    if the neighbor is on a lower refinement level,
    neighbors 1,3,5,7 may or may not exist. If they are not
    present, the pointer value should ne NULL.  */
  CELL* neighbors[8];

  //! Poisson stencil coefficients
  /*!
    winding order of the stencil coefficients:
  
    \verbatim
         | 0  | 1  |
     ____|____|____|_____
         |         |
      7  |         |  2
     ____|    8    |_____
         |         |
      6  |         |  3
     ____|_________|_____
         |    |    |
         | 5  | 4  |
    \endverbatim
    Stencils 0,2,4,6 should always exist. Depending on
    if the neighbor is on a lower refinement level,
    stencils 1,3,5,7 may or may not exist. If they are not
    present, the pointer value should ne NULL.    */
  float stencil[9];
  
  float center[2];    ///< center of the cell
  int depth;          ///< current tree depth
  bool candidate;     ///< already a member of candidate list?

  CELL* parent;       ///< parent node in the quadtree
  CELL_STATE state;   ///< DBM state of the cell

  void refine();      ///< subdivide the cell

  ////////////////////////////////////////////////////////////////
  // solver-related variables
  ////////////////////////////////////////////////////////////////
  bool boundary;      ///< boundary node to include in the solver?
  float potential;    ///< current electric potential
  float b;            ///< rhs of the linear system
  float residual;     ///< residual in the linear solver
  int index;          ///< lexicographic index for the solver

  ////////////////////////////////////////////////////////////////
  // neighbor lookups
  ////////////////////////////////////////////////////////////////
  CELL* northNeighbor();  ///< lookup northern neighbor
  CELL* southNeighbor();  ///< lookup southern neighbor
  CELL* westNeighbor();   ///< lookup western neighbor
  CELL* eastNeighbor();   ///< lookup eastern neighbor
};

#endif
