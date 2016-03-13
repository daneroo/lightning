///////////////////////////////////////////////////////////////////////////////////
// File : CELL.cpp
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

#include "CELL.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

// normal cell constructor
CELL::CELL(float north, float east, float south, float west, CELL* parent, int depth) : 
  parent(parent), depth(depth), index(-1), candidate(false),
  boundary(false), potential(0.0f), state(EMPTY)
{
  for (int x = 0; x < 4; x++)
    children[x] = NULL;
  for (int x = 0; x < 8; x++)
    neighbors[x] = NULL;

  bounds[0] = north; bounds[1] = east; bounds[2] = south; bounds[3] = west;

  center[0] = (bounds[1] + bounds[3]) * 0.5f;
  center[1] = (bounds[0] + bounds[2]) * 0.5f;
}

// ghost cell constructor
CELL::CELL(int depth) : parent(NULL), depth(depth), index(-1), candidate(false),
  boundary(true), potential(0.0f), state(EMPTY)
{
  for (int x = 0; x < 4; x++)
    children[x] = NULL;
  for (int x = 0; x < 8; x++)
    neighbors[x] = NULL;
  
  bounds[0] = 0.0f; bounds[1] = 0.0f; bounds[2] = 0.0f; bounds[3] = 0.0f; 
  center[0] = 0.0f; center[1] = 0.0f;
}

CELL::~CELL() {
  int x;

  for (x = 0; x < 4; x++)
    if (children[x] != NULL)
    {
      delete children[x];
      children[x] = NULL;
    }
}

//////////////////////////////////////////////////////////////////////
// refine current cell
//////////////////////////////////////////////////////////////////////
void CELL::refine() {
  if (children[0] != NULL) return;
  float center[] = {(bounds[0] + bounds[2]) * 0.5f, (bounds[1] + bounds[3]) * 0.5f};
 
  children[0] = new CELL(bounds[0], center[1], center[0], bounds[3], this, depth + 1);
  children[1] = new CELL(bounds[0], bounds[1], center[0], center[1], this, depth + 1);
  children[2] = new CELL(center[0], bounds[1], bounds[2], center[1], this, depth + 1);
  children[3] = new CELL(center[0], center[1], bounds[2], bounds[3], this, depth + 1);

  children[0]->potential = potential;
  children[1]->potential = potential;
  children[2]->potential = potential;
  children[3]->potential = potential;
}

//////////////////////////////////////////////////////////////////////
// return north neighbor to current cell
//////////////////////////////////////////////////////////////////////
CELL* CELL::northNeighbor()
{
  // if it is the root
  if (this->parent == NULL) return NULL;

  // if it is the southern child of the parent
  if (parent->children[3] == this) return parent->children[0];
  if (parent->children[2] == this) return parent->children[1];

  // else look up higher
  CELL* mu = parent->northNeighbor();

  // if there are no more children to look at,
  // this is the answer
  if (mu == NULL || mu->children[0] == NULL) return mu;
  // if it is the NW child of the parent
  else if (parent->children[0] == this) return mu->children[3];
  // if it is the NE child of the parent
  else return mu->children[2];
}

//////////////////////////////////////////////////////////////////////
// return north neighbor to current cell
//////////////////////////////////////////////////////////////////////
CELL* CELL::southNeighbor()
{
  // if it is the root
  if (this->parent == NULL) return NULL;

  // if it is the northern child of the parent
  if (parent->children[0] == this) return parent->children[3];
  if (parent->children[1] == this) return parent->children[2];

  // else look up higher
  CELL* mu = parent->southNeighbor();

  // if there are no more children to look at,
  // this is the answer
  if (mu == NULL || mu->children[0] == NULL) return mu;
  // if it is the SW child of the parent
  else if (parent->children[3] == this) return mu->children[0];
  // if it is the SE child of the parent
  else return mu->children[1];
}

//////////////////////////////////////////////////////////////////////
// return north neighbor to current cell
//////////////////////////////////////////////////////////////////////
CELL* CELL::westNeighbor()
{
  // if it is the root
  if (this->parent == NULL) return NULL;

  // if it is the eastern child of the parent
  if (parent->children[1] == this) return parent->children[0];
  if (parent->children[2] == this) return parent->children[3];

  // else look up higher
  CELL* mu = parent->westNeighbor();

  // if there are no more children to look at,
  // this is the answer
  if (mu == NULL || mu->children[0] == NULL) return mu;
  // if it is the NW child of the parent
  else if (parent->children[0] == this) return mu->children[1];
  // if it is the SW child of the parent
  else return mu->children[2];
}

//////////////////////////////////////////////////////////////////////
// return north neighbor to current cell
//////////////////////////////////////////////////////////////////////
CELL* CELL::eastNeighbor()
{
  // if it is the root
  if (this->parent == NULL) return NULL;

  // if it is the western child of the parent
  if (parent->children[0] == this) return parent->children[1];
  if (parent->children[3] == this) return parent->children[2];

  // else look up higher
  CELL* mu = parent->eastNeighbor();

  // if there are no more children to look at,
  // this is the answer
  if (mu == NULL || mu->children[0] == NULL) return mu;
  // if it is the NE child of the parent
  else if (parent->children[1] == this) return mu->children[0];
  // if it is the SE child of the parent
  else return mu->children[3];
}

