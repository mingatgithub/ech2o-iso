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
 * CalculateGroundwater.cpp
 *
 *  Created on: Apr 9, 2020
 *      Author: Sylvain Kuppel
 */

#include"Basin.h"

int Basin::CalculateGrndWaterVol() {

  // Summed volume (height) of saturated zones (may include perched aquifers)
  
  int r, c;
  double fc1, fc2, fc3, d1, d2, d3 ;
  double theta1, theta2, theta3 ;
  
#pragma omp parallel for						\
  default(none) private(r, c, fc1, fc2, fc3, d1, d2, d3,		\
			theta1, theta2, theta3)			\
  shared(std::cout)
  
  for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
    
    r = _vSortedGrid.cells[j].row;
    c = _vSortedGrid.cells[j].col;
    fc1 = _fieldcapL1->matrix[r][c];
    fc2 = _fieldcapL2->matrix[r][c];
    fc3 = _fieldcapL3->matrix[r][c];
    d1 = _depth_layer1->matrix[r][c] ;
    d2 = _depth_layer2->matrix[r][c];
    d3 = _soildepth->matrix[r][c] - d1 - d2 ;
    theta1 = _soilmoist1->matrix[r][c];
    theta2 = _soilmoist2->matrix[r][c];
    theta3 = _soilmoist3->matrix[r][c];

    _GrndWater->matrix[r][c] = \
      std::max<double>(0, theta1 - fc1)* d1 +
      std::max<double>(0, theta2 - fc2)* d2 +
      std::max<double>(0, theta3 - fc3)* d3 ;

    // Volume of "continuous" groundwater starting from the surface
    // (different from GrndWater if there's is a discontnuity in L2)
    if(theta1 > fc1 and theta2 < fc2) // Perched aquifer in L1 only
      _GrndWater2->matrix[r][c] = (theta1 - fc1)* d1 ;
    else // Other cases -> equal to GrndWater
      _GrndWater2->matrix[r][c] = _GrndWater->matrix[r][c];
    
  } // end loop on pixels
  
  return EXIT_SUCCESS;
  
} // end pragma


int Basin::CalculateWaterTableDepth(Control &ctrl) {

  // Works from top to bottom, so that WTD may reflect perched, transient
  // water table
  
  int r, c, i;
  double theta, eta1, eta2, eta3, fc1, fc2, fc3 ;
  double d, d1, d2;
  double nstep, step, H, H_old, hGW, fact, phi_z, thr_z;
  double th_r, phi0, kphi;
#pragma omp parallel for						\
  default(none) private(r, c, i, fc1, fc2, fc3,				\
			theta, eta1, eta2, eta3, d, d1, d2,		\
			nstep,step, H, H_old, hGW, fact, phi_z, thr_z ,	\
			th_r, phi0, kphi)				\
  shared(ctrl, std::cout)
  
  for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
    
    r = _vSortedGrid.cells[j].row;
    c = _vSortedGrid.cells[j].col;
    fc1 = _fieldcapL1->matrix[r][c];
    fc2 = _fieldcapL2->matrix[r][c];
    fc3 = _fieldcapL3->matrix[r][c];
    eta1 = _porosityL1->matrix[r][c];
    eta2 = _porosityL2->matrix[r][c];
    eta3 = _porosityL3->matrix[r][c];
    d = _soildepth->matrix[r][c];
    d1 = _depth_layer1->matrix[r][c];
    d2 = _depth_layer2->matrix[r][c];
  
    // If the theta1 is above porosity, then the water table intersects the surface
    if(eta1 - _soilmoist1->matrix[r][c] < RNDOFFERR) {
      _WaterTableDepth->matrix[r][c] = 0.0 ;
      //cout << "no WTD..." ;
    }
    else{   

      if(ctrl.toggle_Poros!=1) {
	// Constant or layer-prescribed porosity
	// use linear assumptions (poros constant across the layer)
	if(_soilmoist1->matrix[r][c] - fc1 > RNDOFFERR) {
	  // (Perched) water table within L1
	  _WaterTableDepth->matrix[r][c] =
	    d1*std::min<double>(1.0,(eta1-_soilmoist1->matrix[r][c])/(eta1-fc1));
	}
	else if(_soilmoist2->matrix[r][c] - fc2 > RNDOFFERR) {
	  // (Perched) water table within L2
	  _WaterTableDepth->matrix[r][c] =
	    d1 + d2*std::min<double>(1.0,
				     (eta2 - _soilmoist2->matrix[r][c])/(eta2 -fc2));
	} else if(_soilmoist3->matrix[r][c] - fc3 > RNDOFFERR) {
	  _WaterTableDepth->matrix[r][c] =
	    d1+d2+ (d-d1-d2)*std::min<double>(1.0,
					      (eta3-_soilmoist3->matrix[r][c])/
					      (eta3 - fc3));
	} else
	  _WaterTableDepth->matrix[r][c] = d;
      }
      
      // Exponential porosity profile case:
      // Calculate taking into account the variable porosity along depth
      // iteratively goes from bottom (n bins) until reaching GW volume
      else{
	
	fact = 1 - powl(_psi_ae->matrix[r][c]/3.36,1/_BClambda->matrix[r][c]);
	phi0 = _porosity0->matrix[r][c];
	kphi = _kporos->matrix[r][c];
	H_old = 0.0;  //  iterative WTD
	H = 0.0;  // iterative WTD
	i = 0;      
	  
	if(_soilmoist1->matrix[r][c] - fc1 > RNDOFFERR) {
	  // Water table within first layer
	  th_r = _theta_rL1->matrix[r][c];
	  theta = _soilmoist1->matrix[r][c];
	  hGW = d1 * std::max<double>(0.0,theta - fc1) ; // saturated "volume"
	  // step min: 5mm, max step numbers: 1000
	  nstep = std::min<double>(1000, std::ceil(d1/0.005)); 
	  step = d1 / nstep;  // bins height
	  do {
	    H_old = H ;
	    phi_z = phi0*expl(-(d1-(i+0.5)*step)/kphi);
	    thr_z = std::min<double>(th_r , phi_z / 4);
	    H += step * fact*(phi_z-thr_z);
	    i += 1;
	  } while ( H <= hGW and i <= 1000) ;
	  // Interpolate between current and previous bins' depths
	  if(i == 1) // close to bottom
	    _WaterTableDepth->matrix[r][c] =
	      std::min<double>(d1, d1 - step*0.5*hGW/H);
	  else	    
	    _WaterTableDepth->matrix[r][c] =
	      std::min<double>(d1, d1 - step*(i-0.5-(H - hGW)/(H - H_old)));
	}
	else if(_soilmoist2->matrix[r][c] - fc2 > RNDOFFERR) {
	  // Wate table within L2
	  th_r = _theta_rL2->matrix[r][c];
	  theta = _soilmoist2->matrix[r][c];
	  hGW = d2 * std::max<double>(0.0,theta - fc2) ; // saturated "volume"
	  // step min: 5mm, max step numbers: 1000
	  nstep = std::min<double>(1000, std::ceil(d2/0.005)); 
	  step = d2 / nstep;  // bins height
	  do {
	    H_old = H ;
	    phi_z = phi0*expl(-(d1+d2-(i+0.5)*step)/kphi);
	    thr_z = std::min<double>(th_r , phi_z / 4);
	    H += step * fact*(phi_z-thr_z);
	    i += 1;
	  } while ( H <= hGW and i <= nstep) ;
	  // Interpolate between current and previous bins' depths
	  if(i == 1) // close to bottom
	    _WaterTableDepth->matrix[r][c] =
	      std::min<double>(d1 + d2, d1 + d2 - step*0.5*hGW/H);
	  else	    
	    _WaterTableDepth->matrix[r][c] =
	      std::min<double>(d1 + d2, d1 + d2 - step*(i-0.5-(H - hGW)/(H - H_old))); 	  
	}
	else if (_soilmoist3->matrix[r][c] - fc3 > RNDOFFERR) {
	  // Water table within L3
	  th_r = _theta_rL3->matrix[r][c];
	  theta = _soilmoist3->matrix[r][c];
	  hGW = (d -d1 - d2)*(theta - fc3) ; // saturated "volume"
	  //cout << "hGW : " << hGW << " theta3 : " << theta << endl ;
	  // step min: 5mm, max step numbers: 1000
	  nstep = std::min<double>(1000, std::ceil((d-d1-d2)/0.005)); 
	  step = (d - d1 - d2) / nstep;  // bins depth
	  //cout << "r " << r << " c " << c << " fact : " << fact
	  // << " hGW : " << hGW << "theta_r : " << th_r << endl;
	  do { 
	    H_old = H ;
	    phi_z = phi0*expl(-(d-(i+0.5)*step)/kphi);
	    thr_z = std::min<double>(th_r , phi_z / 4);
	    H += step * fact*(phi_z-thr_z);
	    i += 1;
	    //if(r==10 and c==22 and H - H_old < 0 and i==1)
	    //  cout << " H_old : " << H_old << " H : " << H <<
	    // " phi_z : " << phi_z << " thr_z : " << thr_z << endl ;
	  } while ( H <= hGW and i <= nstep) ;
	  // Interpolate between current and previous bins' depths
	  if(i == 1) // close to bottom
	    _WaterTableDepth->matrix[r][c] =
	      std::min<double>(d, d - step*0.5*hGW/H);
	  else	    
	    _WaterTableDepth->matrix[r][c] =
	      std::min<double>(d ,d - step*(i-0.5-(H - hGW)/(H - H_old)));
	  //if(r==25 and c==35){
	  //  cout << "phi3 : " << eta3 << " theta3 : " << theta << endl ;	       
	  //  cout << "hGW (m) " << hGW << " H_old : " << H_old << " H : " << H << endl ;
	  //  cout << "depth: " << d - step*(i-0.5-(H - hGW)/(H - H_old)) << " i" << i << endl ;
	  //}
	}
	else
	  _WaterTableDepth->matrix[r][c] = d;

      } // End porosity cases      

    } // end cases where there is non-zero

  } // end loop on pixels
  
  return EXIT_SUCCESS;
  
} // end pragma
