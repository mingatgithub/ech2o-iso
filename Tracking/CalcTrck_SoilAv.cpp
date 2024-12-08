/*******************************************************************************
 * Ech2o, a spatially-distributed, ecohydrologic simulator
 * Copyright (c) 2016 Marco Maneta <marco.maneta@umontana.edu>
 *
 *     This file is part of ech2o, a hydrologic model developed at the 
 *     University of Montana.
 *
 *     Ech2o is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     Ech2o is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with Ech2o.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors:
 *    Marco Maneta, Sylvain Kuppel
 *******************************************************************************/
/*
 * CalcTrck_L1L2.cpp
 *
 *  Created on: Jun 27, 2017
 *      Author: Sylvain Kuppel
 */

#include "Tracking.h"
#include "Grid.h"

// Calculates isotopes, chloride or age weighted average over the top two soil layers
int Tracking::CalcTrcksoil_12(Basin &bsn, int ic){
  
  double d1, d2;
  double theta1, theta2;
  int r, c;
#pragma omp parallel default(shared) private(d1, d2, theta1, theta2, r,c)
  {
#pragma omp for nowait
    for (unsigned int j = 0; j < bsn.getSortedGrid().cells.size(); j++) {
      r = bsn.getSortedGrid().cells[j].row;
      c = bsn.getSortedGrid().cells[j].col;
      d1 = bsn.getSoilDepth1()->matrix[r][c];
      d2 = bsn.getSoilDepth2()->matrix[r][c];
      theta1 = bsn.getSoilMoist1()->matrix[r][c];
      theta2 = bsn.getSoilMoist2()->matrix[r][c];

      if(ic == 1)
	_d2Hsoil_12->matrix[r][c] = (_d2Hsoil1->matrix[r][c] * d1 * theta1
				     + _d2Hsoil2->matrix[r][c] * d2* theta2) /
	  (d1*theta1+d2*theta2);
      if(ic == 2)
	_d18Osoil_12->matrix[r][c] = (_d18Osoil1->matrix[r][c] * d1 * theta1
				      + _d18Osoil2->matrix[r][c] * d2 * theta2) /
	  (d1*theta1+d2*theta2);
      if(ic == 3)
	_cClsoil_12->matrix[r][c] = (_cClsoil1->matrix[r][c] * d1 * theta1
				   + _cClsoil2->matrix[r][c] * d2 * theta2) /
	  (d1*theta1+d2*theta2);
      if(ic == 4)
	_Agesoil_12->matrix[r][c] = (_Agesoil1->matrix[r][c] * d1 * theta1
				   + _Agesoil2->matrix[r][c] * d2  *theta2) /
	  (d1*theta1+d2*theta2);
    }
  }
  return EXIT_SUCCESS;
}

// Calculates isotopes, chloride or age weighted average over the soil layers
int Tracking::CalcTrcksoil_Av(Basin &bsn, int ic){
  
  double depth, d1, d2, d3;
  //double fc; 
  double theta1, theta2, theta3;
  int r, c;
#pragma omp parallel default(shared) private(r,c,depth, d1, d2, d3, \
					     theta1, theta2, theta3)
  {
#pragma omp for nowait
    for (unsigned int j = 0; j < bsn.getSortedGrid().cells.size(); j++) {
      r = bsn.getSortedGrid().cells[j].row;
      c = bsn.getSortedGrid().cells[j].col;
      depth = bsn.getSoilDepth()->matrix[r][c];
      d1 = bsn.getSoilDepth1()->matrix[r][c];
      d2 = bsn.getSoilDepth2()->matrix[r][c];
      d3 = depth - d1 - d2;
      //fc = bsn.getFieldCapacity()->matrix[r][c];
      theta1 = bsn.getSoilMoist1()->matrix[r][c];
      theta2 = bsn.getSoilMoist2()->matrix[r][c];
      theta3 = bsn.getSoilMoist3()->matrix[r][c];


      if(ic == 1)
	_d2HsoilAv->matrix[r][c] = (_d2Hsoil1->matrix[r][c] * d1 * theta1
				    + _d2Hsoil2->matrix[r][c] * d2 * theta2
				    + _d2Hsoil3->matrix[r][c] * d3 * theta3 )/ 
	  (d1*theta1+d2*theta2+d3*theta3);
      //std::min<double>(fc,theta3)
      //+ _d2Hgroundwater->matrix[r][c] * d3 * std::max<double>(0.0,theta3-fc))
      if(ic == 2)
	_d18OsoilAv->matrix[r][c] = (_d18Osoil1->matrix[r][c] * d1 * theta1
				    + _d18Osoil2->matrix[r][c] * d2 * theta2
				    + _d18Osoil3->matrix[r][c] * d3 * theta3 )/ 
	  (d1*theta1+d2*theta2+d3*theta3);
      if(ic == 3)
	_cClsoilAv->matrix[r][c] = (_cClsoil1->matrix[r][c] * d1 * theta1
				    + _cClsoil2->matrix[r][c] * d2 * theta2
				    + _cClsoil3->matrix[r][c] * d3 * theta3 )/ 
	  (d1*theta1+d2*theta2+d3*theta3);
      if(ic == 4)
	_AgesoilAv->matrix[r][c] = (_Agesoil1->matrix[r][c] * d1 * theta1
				    + _Agesoil2->matrix[r][c] * d2 * theta2
				    + _Agesoil3->matrix[r][c] * d3 * theta3 )/ 
	  (d1*theta1+d2*theta2+d3*theta3);
      
    }
  }
  return EXIT_SUCCESS;
}

// Calculates isotopes, chloride or age weighted averaged over the saturated layers
// If there is a perched aquifer, the singature is only averaged over it
// (it is considered to be disconnected from lower saturated pools)
int Tracking::CalcTrcksoil_GW(Basin &bsn, int ic){
  
  
  double depth, d1, d2, d3;
  double fc1, fc2, fc3; 
  double theta1, theta2, theta3;
  double w1, w2, w3; // weights for the weighted average
  int r, c;
  
#pragma omp parallel default(shared) private(r,c,depth, d1, d2, d3,	\
					     fc1, fc2, fc3, theta1, theta2, theta3, \
					     w1, w2, w3)
  {
#pragma omp for nowait
    for (unsigned int j = 0; j < bsn.getSortedGrid().cells.size(); j++) {
      r = bsn.getSortedGrid().cells[j].row;
      c = bsn.getSortedGrid().cells[j].col;
      depth = bsn.getSoilDepth()->matrix[r][c];
      d1 = bsn.getSoilDepth1()->matrix[r][c];
      d2 = bsn.getSoilDepth2()->matrix[r][c];
      d3 = depth - d1 - d2;
      fc1 = bsn.getFieldCapacityL1()->matrix[r][c];
      fc2 = bsn.getFieldCapacityL2()->matrix[r][c];
      fc3 = bsn.getFieldCapacityL3()->matrix[r][c];
      theta1 = bsn.getSoilMoist1()->matrix[r][c];
      theta2 = bsn.getSoilMoist2()->matrix[r][c];
      theta3 = bsn.getSoilMoist3()->matrix[r][c];
      w1 = 0;
      w2 = 0;
      w3 = 0 ;
      
      if(theta1 > fc1){
	if(theta2 < fc2) // Perched aquifer in L1 only
	  w1 = 1 ;
	else{
	  if(theta3 < fc3) {// Perched aquifer in L1+L2
	    w1 = d1 * (theta1-fc1);
	    w2 = d2 * (theta2-fc2);
	  }
	  else{ // Saturated zone over L1, L2 and L3
	     w1 = d1 * (theta1-fc1);
	     w2 = d2 * (theta2-fc2);
	     w3 = d3 * (theta3-fc3);
	  }
	}
      } else {
	if(theta2 < fc2) // Aquifer (if present) in L3 only
	  w3 = 1 ;
	else{
	  if(theta3 < fc3) // Perched aquifer in L2 only
	    w2 = 1;
	  else { // Saturated zone in L2+L3
	    w2 = d2 * (theta2-fc2);
	    w3 = d3 * (theta3-fc3);
	  }
	}
      }

      
      if(ic == 1)
	_d2Hgroundwater->matrix[r][c] =	(_d2Hsoil1->matrix[r][c] * w1 + 
					 _d2Hsoil2->matrix[r][c] * w2 +
					 _d2Hsoil3->matrix[r][c] * w3)/ 
	  (w1 + w2 + w3);	          
      
      if(ic == 2)
	_d18Ogroundwater->matrix[r][c] = (_d18Osoil1->matrix[r][c] * w1 +
					  _d18Osoil2->matrix[r][c] * w2 +
					  _d18Osoil3->matrix[r][c] * w3 )/
	  (w1 + w2 + w3);
      if(ic == 3)
	_cClgroundwater->matrix[r][c] = (_cClsoil1->matrix[r][c] * w1 +
					 _cClsoil2->matrix[r][c] * w2 +
					 _cClsoil3->matrix[r][c] * w3 ) /
	  (w1 + w2 + w3);
      if(ic == 4)
	_Agegroundwater->matrix[r][c] = (_Agesoil1->matrix[r][c] * w1 +
					 _Agesoil2->matrix[r][c] * w2 +
					 _Agesoil3->matrix[r][c] * w3 )/ 
	  (w1 + w2 + w3);
      
    }
  }
  return EXIT_SUCCESS;
}

