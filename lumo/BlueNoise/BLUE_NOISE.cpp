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

// $Id: PDSampling.cpp,v 1.12 2006/07/11 16:45:22 zr Exp $

#define _USE_MATH_DEFINES
#include <cmath>

#include <map>

#include "BLUE_NOISE.h"
#include "RangeList.h"
#include "ScallopedSector.h"
#include "WeightedDiscretePDF.h"

typedef std::vector<int> IntVector;

///

BLUE_NOISE::BLUE_NOISE(float _radius, bool _isTiled, bool usesGrid) :
	m_rng(123456),
	radius(_radius),
	isTiled(_isTiled)
{
	if (usesGrid) {
			// grid size is chosen so that 4*radius search only 
			// requires searching adjacent cells, this also
			// determines max points per cell
		m_gridSize = (int) ceil(2./(4.*_radius));
		if (m_gridSize<2) m_gridSize = 2;

		m_gridCellSize = 2.0f/m_gridSize;
		m_grid = new int[m_gridSize*m_gridSize][kMaxPointsPerCell];

		for (int y=0; y<m_gridSize; y++) {
			for (int x=0; x<m_gridSize; x++) {
				for (int k=0; k<kMaxPointsPerCell; k++) {
					m_grid[y*m_gridSize + x][k] = -1;
				}
			}
		}
	} else {
		m_gridSize = 0;
		m_gridCellSize = 0;
		m_grid = 0;
	}
}
	
bool BLUE_NOISE::pointInDomain(Vec2 &a)
{
	return -1<=a.x && -1<=a.y && 1>=a.x && 1>=a.y;
}

Vec2 BLUE_NOISE::randomPoint()
{
	return Vec2(2*m_rng.getFloatL()-1, 2*m_rng.getFloatL()-1);
}

Vec2 BLUE_NOISE::getTiled(Vec2 v)
{
	float x = v.x, y = v.y;

	if (isTiled) {
		if (x<-1) x += 2;
		else if (x>1) x -= 2;

		if (y<-1) y += 2;
		else if (y>1) y -= 2;
	}

	return Vec2(x,y);
}

void BLUE_NOISE::getGridXY(Vec2 &v, int *gx_out, int *gy_out)
{
	int gx = *gx_out = (int) floor(.5*(v.x + 1)*m_gridSize);
	int gy = *gy_out = (int) floor(.5*(v.y + 1)*m_gridSize);
	if (gx<0 || gx>=m_gridSize || gy<0 || gy>=m_gridSize) {
		printf("Internal error, point outside grid was generated, ignoring.\n");
	}
}

void BLUE_NOISE::addPoint(Vec2 pt)
{
	int i, gx, gy, *cell;

	points.push_back(pt);

	if (m_grid) {
		getGridXY(pt, &gx, &gy);
		cell = m_grid[gy*m_gridSize + gx];
		for (i=0; i<kMaxPointsPerCell; i++) {
			if (cell[i]==-1) {
				cell[i] = (int) points.size()-1;
				break;
			}
		}
		if (i==kMaxPointsPerCell) {
			printf("Internal error, overflowed max points per grid cell. Exiting.\n");
			exit(1);
		}
	}
}

int BLUE_NOISE::findNeighbors(Vec2 &pt, float distance)
{
	if (!m_grid) {
		printf("Internal error, sampler cannot search without grid.\n");
		exit(1);
	}

	float distanceSqrd = distance*distance;
	int i, j, k, gx, gy, N = (int) ceil(distance/m_gridCellSize);
	if (N>(m_gridSize>>1)) N = m_gridSize>>1;
	
	m_neighbors.clear();
	getGridXY(pt, &gx, &gy);
	for (j=-N; j<=N; j++) {
		for (i=-N; i<=N; i++) {
			int cx = (gx+i+m_gridSize)%m_gridSize;
			int cy = (gy+j+m_gridSize)%m_gridSize;
			int *cell = m_grid[cy*m_gridSize + cx];

			for (k=0; k<kMaxPointsPerCell; k++) {
				if (cell[k]==-1) {
					break;
				} else {
					if (getDistanceSquared(pt, points[cell[k]])<distanceSqrd)
						m_neighbors.push_back(cell[k]);
				}
			}
		}
	}

	return (int) m_neighbors.size();
}

float BLUE_NOISE::findClosestNeighbor(Vec2 &pt, float distance)
{
	if (!m_grid) {
		printf("Internal error, sampler cannot search without grid.\n");
		exit(1);
	}

	float closestSqrd = distance*distance;
	int i, j, k, gx, gy, N = (int) ceil(distance/m_gridCellSize);
	if (N>(m_gridSize>>1)) N = m_gridSize>>1;
	
	getGridXY(pt, &gx, &gy);
	for (j=-N; j<=N; j++) {
		for (i=-N; i<=N; i++) {
			int cx = (gx+i+m_gridSize)%m_gridSize;
			int cy = (gy+j+m_gridSize)%m_gridSize;
			int *cell = m_grid[cy*m_gridSize + cx];

			for (k=0; k<kMaxPointsPerCell; k++) {
				if (cell[k]==-1) {
					break;
				} else {
					float d = getDistanceSquared(pt, points[cell[k]]);

					if (d<closestSqrd)
						closestSqrd = d;
				}
			}
		}
	}

	return sqrt(closestSqrd);
}

void BLUE_NOISE::findNeighborRanges(int index, RangeList &rl)
{
	if (!m_grid) {
		printf("Internal error, sampler cannot search without grid.\n");
		exit(1);
	}

	Vec2 &candidate = points[index];
	float rangeSqrd = 4*4*radius*radius;
	int i, j, k, gx, gy, N = (int) ceil(4*radius/m_gridCellSize);
	if (N>(m_gridSize>>1)) N = m_gridSize>>1;
	
	getGridXY(candidate, &gx, &gy);

	int xSide = (candidate.x - (-1 + gx*m_gridCellSize))>m_gridCellSize*.5;
	int ySide = (candidate.y - (-1 + gy*m_gridCellSize))>m_gridCellSize*.5;
	int iy = 1;
	for (j=-N; j<=N; j++) {
		int ix = 1;

		if (j==0) iy = ySide;
		else if (j==1) iy = 0;

		for (i=-N; i<=N; i++) {
			if (i==0) ix = xSide;
			else if (i==1) ix = 0;

				// offset to closest cell point
			float dx = candidate.x - (-1 + (gx+i+ix)*m_gridCellSize);
			float dy = candidate.y - (-1 + (gy+j+iy)*m_gridCellSize);

			if (dx*dx+dy*dy<rangeSqrd) {
				int cx = (gx+i+m_gridSize)%m_gridSize;
				int cy = (gy+j+m_gridSize)%m_gridSize;
				int *cell = m_grid[cy*m_gridSize + cx];

				for (k=0; k<kMaxPointsPerCell; k++) {
					if (cell[k]==-1) {
						break;
					} else if (cell[k]!=index) {
						Vec2 &pt = points[cell[k]];
						Vec2 v = getTiled(pt-candidate);
						float distSqrd = v.x*v.x + v.y*v.y;

						if (distSqrd<rangeSqrd) {
							float dist = sqrt(distSqrd);
							float angle = atan2(v.y,v.x);
							float theta = acos(.25f*dist/radius);

							rl.subtract(angle-theta, angle+theta);
						}
					}
				}
			}
		}
	}
}

void BLUE_NOISE::maximize()
{
	RangeList rl(0,0);
	int i, N = (int) points.size();

	for (i=0; i<N; i++) {
		Vec2 &candidate = points[i];

		rl.reset(0, (float) M_PI*2);
		findNeighborRanges(i, rl);
		while (rl.numRanges) {
			RangeEntry &re = rl.ranges[m_rng.getInt31()%rl.numRanges];
			float angle = re.min + (re.max-re.min)*m_rng.getFloatL();
			Vec2 pt = getTiled(Vec2(candidate.x + cos(angle)*2*radius,
									candidate.y + sin(angle)*2*radius));

			addPoint(pt);
			rl.subtract(angle - (float) M_PI/3, angle + (float) M_PI/3);
		}
	}
}

void BLUE_NOISE::relax()
{
	FILE *tmp = fopen("relaxTmpIn.txt","w");
	int dim, numVerts, numFaces;
	Vec2 *verts = 0;
	int numPoints = (int) points.size();

		// will overwrite later
	fprintf(tmp, "2                  \n");
	for (int i=0; i<(int) points.size(); i++) {
		Vec2 &pt = points[i];
		fprintf(tmp, "%f %f\n", pt.x, pt.y);
	}
	for (int y=-1; y<=1; y++) {
		for (int x=-1; x<=1; x++) {
			if (x || y) {
				for (int i=0; i<(int) points.size(); i++) {
					Vec2 &pt = points[i];
					if (fabs(pt.x+x*2)-1<radius*4 || fabs(pt.y+y*2)-1<radius*4) {
						fprintf(tmp, "%f %f\n", pt.x+x*2, pt.y+y*2);
						numPoints++;
					}
				}
			}
		}
	}
	fseek(tmp, 0, 0);
	fprintf(tmp, "2 %d", numPoints);
	fclose(tmp);

	tmp = fopen("relaxTmpOut.txt", "w");
	fclose(tmp);
	system("qvoronoi p FN < relaxTmpIn.txt > relaxTmpOut.txt");

	tmp = fopen("relaxTmpOut.txt", "r");
	fscanf(tmp, "%d\n%d\n", &dim, &numVerts);

	if (dim!=2) {
		printf("Error calling out to qvoronoi, skipping relaxation.\n");
		goto exit;
	}

	verts = new Vec2[numVerts];
	for (int i=0; i<numVerts; i++) {
		fscanf(tmp, "%f %f\n", &verts[i].x, &verts[i].y);
	}

	fscanf(tmp, "%d\n", &numFaces);

	for (int i=0; i<(int) points.size(); i++) {
		Vec2 center(0,0);
		int N, skip=0;

		fscanf(tmp, "%d", &N);
		for (int j=0; j<N; j++) {
			int index;

			fscanf(tmp, "%d", &index);
			if (index<0) {
				skip = 1;
			} else {
				center += verts[index];
			}
		}

		if (!skip) {
			center *= (1.0f/N);
			points[i] = getTiled(center);
		}
	}

exit:
	if (verts) delete verts;
}

///

typedef std::map<int, ScallopedRegion*> RegionMap;

void BLUE_NOISE::complete()
{
	RangeList rl(0,0);
	IntVector candidates;

	addPoint(randomPoint());
	candidates.push_back((int) points.size()-1);

	while (candidates.size()) {
		int c = m_rng.getInt32()%candidates.size();
		int index = candidates[c];
		Vec2 candidate = points[index];
		candidates[c] = candidates[candidates.size()-1];
		candidates.pop_back();

		rl.reset(0, (float) M_PI*2);
		findNeighborRanges(index, rl);
		while (rl.numRanges) {
			RangeEntry &re = rl.ranges[m_rng.getInt32()%rl.numRanges];
			float angle = re.min + (re.max-re.min)*m_rng.getFloatL();
			Vec2 pt = getTiled(Vec2(candidate.x + cos(angle)*2*radius,
									candidate.y + sin(angle)*2*radius));

			addPoint(pt);
			candidates.push_back((int) points.size()-1);

			rl.subtract(angle - (float) M_PI/3, angle + (float) M_PI/3);
		}
	}
}

void BLUE_NOISE::writeToBool(bool* noise, int size)
{
  // wipe
  int index = 0;
  for (index = 0; index < size * size; index++)
    noise[index] = false;

  for (int x = 0; x < points.size(); x++)
  {
    int i = (points[x].x + 1.0f) * 0.5f * size;
    int j = (points[x].y + 1.0f) * 0.5f * size;
    index = i + j * size;
    noise[index] = true;
  }
}
