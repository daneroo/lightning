///////////////////////////////////////////////////////////////////////////////////
// File : DAG.h
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

#ifndef DAG_H
#define DAG_H

#include <map>
#include <vector>
#include <cmath>
#include <gl/glut.h>
#include <iostream>

using namespace std;

////////////////////////////////////////////////////////////////////
/// \brief Directed acyclic graph that renders the final lightning
////////////////////////////////////////////////////////////////////
class DAG  
{
public:
  //! constructor
	DAG(int xRes, int yRes);
  //! destructor
	virtual ~DAG();

  //! build the stepped ladder
  void buildLeader(int bottomHit);

  //! add DAG segment
  bool addSegment(int index, int neighbor);

  //! draw to OpenGL
  void draw() { drawNode(_root); };

  //! draw to an offscreen buffer
  float*& drawOffscreen(int scale = 1);

  //! read in a new DAG
  void read(const char* filename);
  //! write out the current DAG
  void write(const char* filename);

  //! quadtree x resolution accessor
  int xRes() { return _xRes; };
  //! quadtree y resolution accessor
  int yRes() { return _yRes; };
 
  //! input image x resolution accessor
  int& inputWidth() { return _inputWidth; };
  //! input image y resolution accessor
  int& inputHeight() { return _inputHeight; };
  
private:
  //! x resolution of quadtree
  int _xRes;
  //! y resolution of quadtree
  int _yRes;
  //! physical length of one grid cell
  float _dx;
  //! physical length of one grid cell
  float _dy;
  
  ////////////////////////////////////////////////////////////////////
  /// \brief node for line segment tree
  ////////////////////////////////////////////////////////////////////
  struct NODE {
    int index;
    vector<NODE*> neighbors;
    NODE* parent;
    bool leader;
    bool secondary;
    int depth;
    NODE* maxDepthNode;
    float intensity;

    NODE(int indexIn) { 
      index = indexIn; 
      parent = NULL;
      leader = false;
      depth = 0;
      maxDepthNode = NULL;
    };
  };
  //! recursive destructor
  void deleteNode(NODE* root);

  //! root of DAG
  NODE* _root;

  //! hash table of DAG nodes
  map<int, NODE*> _hash;

  //! build side branch
  void buildBranch(NODE* node, int depth);
  //! draw a node to OpenGL
  void drawNode(NODE* root);
  //! find the deepest node in a given subtree
  void findDeepest(NODE* root, NODE*& deepest);

  //! read in a DAG node 
  void readNode(FILE* file);
  //! write out a DAG node 
  void writeNode(NODE* root, FILE* file);

  //! total number of nodes in scene
  int _totalNodes;

  //! node that finally hit bottom
  int _bottomHit;

  //! set the line segment intensities
  void buildIntensity(NODE* root);
  
  //! brightness of secondary branch
  float _secondaryIntensity;
  //! brightness of primary branch
  float _leaderIntensity;

  ////////////////////////////////////////////////////////////////
  // offscreen buffer variables
  ////////////////////////////////////////////////////////////////
  //! offscreen buffer
  float* _offscreenBuffer;
  
  //! width of offscreen buffer
  int _width;
  
  //! height of offscreen buffer
  int _height;
  
  //! scale of offscreen buffer compared to original image
  int _scale;
  
  //! draw a given node in the DAG
  void drawOffscreenNode(NODE* root);

  //! rasterize a single line to the offscreen buffer
  void drawLine(int begin[], int end[], float intensity);

  //! input image x resolution
  int _inputWidth;
  //! input image y resolution
  int _inputHeight;
};

#endif
