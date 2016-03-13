///////////////////////////////////////////////////////////////////////////////////
// File : main.cpp
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
///
/// \mainpage Fast Animation of Lightning Using An Adaptive Mesh
/// \section Introduction
///
/// This project is an implementation of the paper
/// <b><em>Fast Animation of Lightning Using An Adaptive Mesh</em></b>. It
/// includes both the simulation and rendering components described in that paper.
///
/// Several pieces of software are used in this project that the respective
/// authors were kind enough to make freely available:
///
/// <UL>
/// <LI> <A HREF = "http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html">Mersenne twister</A>
///      (Thanks to Makoto Matsumoto)
/// <LI> <A HREF = "http://www.fftw.org/">FFTW</A>
///      (Thanks to Matteo Frigo and Steven G. Johnson)
/// <LI> <A HREF = "http://www.openexr.com/">OpenEXR</A>
///      (Thanks to ILM)
/// <LI> <A HREF = "http://www.cs.unc.edu/~walk/software/glvu/">GLVU</A>
///      (Thanks to Walkthru)
/// <LI> <A HREF = "http://www.cs.virginia.edu/~gfx/pubs/antimony/">Antimony</A>
///      (Thanks to Daniel Dunbar and Greg Humphreys)
/// </UL>
///
/// <em>Theodore Kim, kim@cs.unc.edu, October 2006</em>
///
///////////////////////////////////////////////////////////////////////////////////

#define COMMAND_LINE_VERSION 1

#include <iostream>
#include "ppm\ppm.hpp"
#include "APSF.h"
#include "FFT.h"
#include "QUAD_DBM_2D.h"
#include "EXR.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////
// globals
////////////////////////////////////////////////////////////////////////////
int iterations = 10;
static QUAD_DBM_2D* potential = new QUAD_DBM_2D(256, 256, iterations);
APSF apsf(512);

// input image info
int inputWidth = -1;
int inputHeight = -1;

// input params
string inputFile;
string outputFile;

// image scale
float scale = 5;

// pause the simulation?
bool pause = false;

////////////////////////////////////////////////////////////////////////////
// render the glow
////////////////////////////////////////////////////////////////////////////
void renderGlow(string filename, int scale = 1)
{
  int w = potential->xDagRes() * scale;
  int h = potential->yDagRes() * scale;

  // draw the DAG
  float*& source = potential->renderOffscreen(scale);
  
  // if there is no input dimensions specified, else there were input
  // image dimensions, so crop it
  if (inputWidth == -1)
  {
    inputWidth  = potential->inputWidth();
    inputHeight = potential->inputHeight();
  }

  // copy out the cropped version
  int wCropped = inputWidth * scale;
  int hCropped = inputHeight * scale;
  float* cropped = new float[wCropped * hCropped];
  cout << endl << " Generating EXR image width: " << wCropped << " height: " << hCropped << endl;
  for (int y = 0; y < hCropped; y++)
    for (int x = 0; x < wCropped; x++)
    {
      int uncroppedIndex = x + y * w;
      int croppedIndex = x + y * wCropped;
      cropped[croppedIndex] = source[uncroppedIndex];
    }

  // create the filter
  apsf.generateKernelFast();
 
  // convolve with FFT
  bool success = FFT::convolve(cropped, apsf.kernel(), wCropped, hCropped, apsf.res(), apsf.res());
   
  if (success) {
    EXR::writeEXR(filename.c_str(), cropped, wCropped, hCropped);
    cout << " " << filename << " written." << endl;
  }
  else
    cout << " Final image generation failed." << endl;
    
  delete[] cropped;
}

////////////////////////////////////////////////////////////////////////////
// load image file into the DBM simulation
////////////////////////////////////////////////////////////////////////////
bool loadImages(string inputFile)
{
  // load the files
  unsigned char* input = NULL;
  LoadPPM(inputFile.c_str(), input, inputWidth, inputHeight);

  unsigned char* start       = new unsigned char[inputWidth * inputHeight];
  unsigned char* repulsor    = new unsigned char[inputWidth * inputHeight];
  unsigned char* attractor   = new unsigned char[inputWidth * inputHeight];
  unsigned char* terminators = new unsigned char[inputWidth * inputHeight];
  
  // composite RGB channels into one
  for (int x = 0; x < inputWidth * inputHeight; x++)
  {
    start[x]     = (input[3 * x] == 255)     ? 255 : 0;
    repulsor[x]  = (input[3 * x + 1] == 255) ? 255 : 0;
    attractor[x] = (input[3 * x + 2] == 255) ? 255 : 0;
    terminators[x] = 0;
    
    if (input[3 * x] + input[3 * x + 1] + input[3 * x + 2] == 255 * 3)
    {
      terminators[x] = 255;
      start[x] = repulsor[x] = attractor[x] = 0;
    }
  }

  if (potential) delete potential;
  potential = new QUAD_DBM_2D(inputWidth, inputHeight, iterations);
  bool success = potential->readImage(start, attractor, repulsor, terminators, inputWidth, inputHeight);
  
  // delete the memory
  delete[] input;
  delete[] start;
  delete[] repulsor;
  delete[] attractor;
  delete[] terminators;

  return success;
}

int width  = 600;
int height = 600;
bool animate = false;
float camera[2];
float translate[2];

////////////////////////////////////////////////////////////////////////////
// window Reshape function 
////////////////////////////////////////////////////////////////////////////
void Reshape(int w, int h)
{
  if (h == 0) h = 1;
  
  glViewport(0, 0, w, h);
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-camera[0] - translate[0], camera[0] + translate[0], -camera[1] - translate[1], camera[1] + translate[1]);
  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

////////////////////////////////////////////////////////////////////////////
// GLUT Display callback
////////////////////////////////////////////////////////////////////////////
void Display()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-camera[0] + translate[0], camera[0] + translate[0], -camera[1] + translate[1], camera[1] + translate[1]);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  
  potential->draw();
  potential->drawSegments();
  
  glutSwapBuffers();
}

////////////////////////////////////////////////////////////////////////////
// GLUT Keyboard callback
////////////////////////////////////////////////////////////////////////////
void Keyboard(unsigned char key, int x, int y)
{
  switch(key) 
  {
    case 'p':
      pause = !pause;
      break;
    
    case 'q':
      cout << " You terminated the simulation prematurely." << endl;
      exit(0);
      break;
  }
  
  glutPostRedisplay();
}

////////////////////////////////////////////////////////////////////////////
// window Reshape function 
////////////////////////////////////////////////////////////////////////////
void Idle()
{
  if (!pause)
    for (int x = 0; x < 100; x++)
    {
      bool success = potential->addParticle();

      if (!success)
      {
        cout << " No nodes left to add! Is your terminator reachable?" << endl;
        //exit(1);
        return;
      }
      
      if (potential->hitGround())
      {
        glutPostRedisplay();
        cout << endl << endl;

        // write out the DAG file
        string lightningFile = inputFile.substr(0, inputFile.size() - 3) + string("lightning");
        cout << " Intermediate file " << lightningFile << " written." << endl;
        potential->writeDAG(lightningFile.c_str());
        
        // render the final EXR file
        renderGlow(outputFile, scale);
        delete potential;
        exit(0);
      }
    }
  glutPostRedisplay();
}

////////////////////////////////////////////////////////////////////////////
// GLUT Main 
////////////////////////////////////////////////////////////////////////////
int glutMain()
{ 
  float smaller = 1.0f;
  camera[0] = smaller * 0.5f;
  camera[1] = smaller * 0.5f;
  translate[0] = 0.0f;
  translate[1] = 0.0f;
 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(width, height);
  glutCreateWindow("Lumos: A Lightning Generator v0.1");

  glutDisplayFunc(Display);
  glutKeyboardFunc(Keyboard);
  glutIdleFunc(Idle);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE);
  
  Reshape(width, height); 

  glClearColor(0.0, 0.0, 0.0, 1.0);
  glShadeModel(GL_SMOOTH);
  
  // Go!
  glutMainLoop();

  return 0;
}

////////////////////////////////////////////////////////////////////////////
// Main 
////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  if (argc < 3)
  {
    cout << endl;
    cout << "   LumosQuad <input file> <output file> <scale (optional)>" << endl;
    cout << "   =========================================================" << endl;
    cout << "      <input file>  - *.ppm file with input colors" << endl;
    cout << "                      --OR--" << endl;
    cout << "                      *.lightning file from a previous run" << endl;
    cout << "      <output file> - The OpenEXR file to output" << endl;
    cout << "      <scale>       - Scaling constant for final image." << endl;
    cout << "   Press 'q' to terminate the simulation prematurely." << endl;
    cout << "   Send questions and comments to kim@cs.unc.edu" << endl;
    return 1;
  }

  cout << endl << "Lumos: A lightning generator v0.1" << endl;
  cout << "------------------------------------------------------" << endl;

  // store the input params
  inputFile = string(argv[1]);
  outputFile = string(argv[2]);
  if (argc > 3) scale = atoi(argv[3]);
 
  // see if the input is a *.lightning file
  if (inputFile.size() > 10)
  {
    string postfix = inputFile.substr(inputFile.size() - 9, inputFile.size());

    cout << " Using intermediate file " << inputFile << endl;
    if (postfix == string("lightning"))
    {
      potential->readDAG(inputFile.c_str());
      renderGlow(outputFile, scale);
      delete potential;
      return 0;
    }
  }
  
  // read in the *.ppm input file
  if (!loadImages(inputFile))
  {
    cout << " ERROR: " << inputFile.c_str() << " is not a valid PPM file." << endl;
    return 1;
  }
  cout << " " << inputFile << " read." << endl << endl;


  // loop simulation until it hits a terminator
  cout << " Total particles added: ";
  glutMain();

  return 0;
}
