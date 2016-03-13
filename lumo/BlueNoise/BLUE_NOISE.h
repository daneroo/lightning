///////////////////////////////////////////////////////////////////////////////////
// File : BLUE_NOISE.h
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
///////////////////////////////////////////////////////////////////////////////////
//
// This class is a very thin wrapper to Daniel Dunbar's blue noise generator.
// With the exception of BLUE_NOISE.h and BLUE_NOISE.cpp, the other files
// in this directory are unmodified copies of his code.
//
// For the original, untainted code, see: 
//   http://www.cs.virginia.edu/~gfx/pubs/antimony/
//
///////////////////////////////////////////////////////////////////////////////////

// $Id: PDSampling.h,v 1.6 2006/07/06 23:13:18 zr Exp $

#include "RNG.h"
#include <cmath>
#include <vector>

#define kMaxPointsPerCell 9

class RangeList;
class ScallopedRegion;

class Vec2 {
public:
	Vec2() {};
	Vec2(float _x, float _y) : x(_x), y(_y) {};

	float x,y;

	float length() { return sqrt(x*x + y*y); }

	bool operator ==(const Vec2 &b) const { return x==b.x && y==b.y; }
	Vec2 operator +(Vec2 b) { return Vec2(x+b.x, y+b.y); }
	Vec2 operator -(Vec2 b) { return Vec2(x-b.x, y-b.y); }
	Vec2 operator *(Vec2 b) { return Vec2(x*b.x, y*b.y); }
	Vec2 operator /(Vec2 b) { return Vec2(x/b.x, y*b.y); }

	Vec2 operator +(float n) { return Vec2(x+n, y+n); }
	Vec2 operator -(float n) { return Vec2(x-n, y-n); }
	Vec2 operator *(float n) { return Vec2(x*n, y*n); }
	Vec2 operator /(float n) { return Vec2(x/n, y*n); }

	Vec2 &operator +=(Vec2 b) { x+=b.x; y+=b.y; return *this; }
	Vec2 &operator -=(Vec2 b) { x-=b.x; y-=b.y; return *this; }
	Vec2 &operator *=(Vec2 b) { x*=b.x; y*=b.y; return *this; }
	Vec2 &operator /=(Vec2 b) { x/=b.x; y/=b.y; return *this; }

	Vec2 &operator +=(float n) { x+=n; y+=n; return *this; }
	Vec2 &operator -=(float n) { x-=n; y-=n; return *this; }
	Vec2 &operator *=(float n) { x*=n; y*=n; return *this; }
	Vec2 &operator /=(float n) { x/=n; y/=n; return *this; }
};

/// \brief Daniel Dunbar's blue noise generator
/// 
/// The original code has been modified so that the 'boundary sampling'
/// method is the only one available.
class BLUE_NOISE {
protected:
	RNG m_rng;
	std::vector<int> m_neighbors;
	
	int (*m_grid)[kMaxPointsPerCell];
	int m_gridSize;
	float m_gridCellSize;
	
public:
	std::vector<Vec2> points;
	float radius;
	bool isTiled;

public:
	BLUE_NOISE(float radius, bool isTiled=true, bool usesGrid=true);
	virtual ~BLUE_NOISE() { };

	//

	bool pointInDomain(Vec2 &a);

		// return shortest distance between _a_ 
		// and _b_ (accounting for tiling)
	float getDistanceSquared(Vec2 &a, Vec2 &b) { Vec2 v = getTiled(b-a); return v.x*v.x + v.y*v.y; }
	float getDistance(Vec2 &a, Vec2 &b) { return sqrt(getDistanceSquared(a, b)); }

		// generate a random point in square
	Vec2 randomPoint();

		// return tiled coordinates of _v_
	Vec2 getTiled(Vec2 v);

		// return grid x,y for point
	void getGridXY(Vec2 &v, int *gx_out, int *gy_out);

		// add _pt_ to point list and grid
	void addPoint(Vec2 pt);

		// populate m_neighbors with list of
		// all points within _radius_ of _pt_
		// and return number of such points
	int findNeighbors(Vec2 &pt, float radius);

		// return distance to closest neighbor within _radius_
	float findClosestNeighbor(Vec2 &pt, float radius);

		// find available angle ranges on boundary for candidate 
		// by subtracting occluded neighbor ranges from _rl_
	void findNeighborRanges(int index, RangeList &rl);

		// extend point set by boundary sampling until domain is 
		// full
	void maximize();

		// apply one step of Lloyd relaxation
	void relax();

	void complete();

  void writeToBool(bool* noise, int size);
};
