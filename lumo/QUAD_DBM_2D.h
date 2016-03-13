///////////////////////////////////////////////////////////////////////////////////
// File : QUAD_DBM_2D.h
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

#ifndef QUAD_DBM_2D_H 
#define QUAD_DBM_2D_H

#include <vector>
#include <gl/glut.h>
#include "DAG.h"
#include "QUAD_POISSON.h"

////////////////////////////////////////////////////////////////////
/// \brief Quadtree DBM solver. This is the highest level class.
////////////////////////////////////////////////////////////////////
class QUAD_DBM_2D  
{
public:
  /// \brief DBM constructor 
  ///
  /// \param xRes         maximum x resolution
  /// \param yRes         maximum y resolution
  /// \param iterations   maximum conjugate gradient iterations
	QUAD_DBM_2D(int xRes = 128, int yRes = 128, int iterations = 10);

  //! destructor
	virtual ~QUAD_DBM_2D();

  //! add to aggregate
  bool addParticle();
  
  /// \brief Hit ground yet?
  /// \return returns true if a terminator as already been hit
  bool hitGround(CELL* cell = NULL);
 
  //! draw the quadtree cells to OpenGL
  void draw();

  //! draw the DAG to OpenGL
  void drawSegments()  {
    glLineWidth(1.0f);
    glPushMatrix();
    glTranslatef(-0.5f, -0.5f, 0.0f);
    _dag->draw();
    glPopMatrix();
  };
  
  ////////////////////////////////////////////////////////////////
  // file IO
  ////////////////////////////////////////////////////////////////

  //! write everything to a file
  void writeFields(const char* filename);

  //! read everything from a file
  void readFields(const char* filename);

  /// \brief read in control parameters from an input file
  ///
  /// \param initial        initial pixels of lightning
  /// \param attractors     pixels that attract the lightning
  /// \param repulsors      pixels that repulse the lightning
  /// \param terminators    pixels that halt the simulation if hit
  /// \param xRes           x resolution of the image
  /// \param yRes           y resolution of the image
  ///
  /// \return Returns false if it finds something wrong with the images
  bool readImage(unsigned char* initial, 
                 unsigned char* attractors,
                 unsigned char* repulsors,
                 unsigned char* terminators,
                 int xRes, int yRes);

  //! read in a new DAG
  void readDAG(const char* filename)     { _dag->read(filename); };

  //! write out the current DAG
  void writeDAG(const char* filename)    { _dag->write(filename); };
  
  /// \brief render to a software-only buffer
  ///
  /// \param scale      a (scale * xRes) x (scale * yRes) image is rendered
  float*& renderOffscreen(int scale = 1) { return _dag->drawOffscreen(scale); };
  
  //! access the DBM x resolution 
  int xRes() { return _xRes; };
  //! access the DBM y resolution
  int yRes() { return _yRes; };
  //! access the DAG x resolution
  int xDagRes() { return _dag->xRes(); };
  //! access the DAG y resolution
  int yDagRes() { return _dag->yRes(); };
  //! access the x resolution of the input image
  int inputWidth() { return _dag->inputWidth(); };
  //! access the y resolution of the input image
  int inputHeight() { return _dag->inputHeight(); };

private:
  void allocate();
  void deallocate();
  
  ////////////////////////////////////////////////////////////////////
  // dielectric breakdown model components
  ////////////////////////////////////////////////////////////////////
  
  // field dimensions
  int _xRes;
  int _yRes;
  int _maxRes;
  float _dx;
  float _dy;
  int _iterations;
  
  // which cell did it hit bottom with?
  int _bottomHit;

  DAG* _dag;

  QUAD_POISSON* _quadPoisson;

  // current candidate list
  vector<CELL*> _candidates;

  // check if any of the neighbors of cell should be added to the
  // candidate list
  void checkForCandidates(CELL* cell);

  // number of particles to add before doing another Poisson solve
  int _skips;

  // Mersenne Twister
  RNG _twister;
};

#endif
