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
 * CheckMaps.cpp
 *
 *  Created on: Feb 10, 2011
 *      Author: Marco.Maneta
 */

#include "Basin.h"

void Basin::CheckMaps(Control &ctrl) {

  UINT4 r, c;
  UINT4 j = 0;
  UINT4 excep_thrown = 0; //  poor man way  to rethrow string exception outside omp pragma
  UINT4 length = _vSortedGrid.cells.size();

#pragma omp parallel for						\
  default(shared) private(r,c,j) //shared(length, cout, excep_thrown)
  for (j = 0; j < length; j++) {
    r = _vSortedGrid.cells[j].row;
    c = _vSortedGrid.cells[j].col;
    try {

      if (_vSortedGrid.cells[j].dir < 1
	  || _vSortedGrid.cells[j].dir > 9)
	throw string("The drain direction map is not a valid ldd map...\n");

      if (_slope->matrix[r][c] == _slope->nodata)
	throw string("Slope map contains no data values inside the valid domain...\n");

      if (_slope->matrix[r][c] == 0)
	_slope->matrix[r][c] = MIN_SLOPE;
      if (_slope->matrix[r][c] <= 0)
	throw string("Slope map contains negative or zero values inside the valid domain...\n");

      if(ctrl.toggle_Ksat!=1){
	if (_KsatL1->matrix[r][c] == _KsatL1->nodata)
	  throw string("Layer 1 hydraulic conductivity map contains no data values inside the valid domain...\n");
	if (_KsatL1->matrix[r][c] <= 0)
	  throw string("Layer 1 hydraulic conductivity map contains negative or zero values inside the valid domain...\n");
      }

      if(ctrl.toggle_Ksat==1){ // Exponential Kh profile
	if (_Ksat0->matrix[r][c] == _Ksat0->nodata)
	  throw string("Top-of-profile Ksat map contains no data values inside the valid domain...\n");
	if (_Ksat0->matrix[r][c] <= 0)
	  throw string("Top-of-profile Ksat map contains negative or zero values inside the valid domain...\n");
	if (_kKsat->matrix[r][c] == _kKsat->nodata)
	  throw string("Ksat profile map contains no data values inside the valid domain...\n");
	if (_kKsat->matrix[r][c] <= 0)
	  throw string("Ksat profile map contains negative or zero values inside the valid domain...\n");
      }
      if(ctrl.toggle_Ksat==2){ // Layer-defined Khs
	if (_KsatL2->matrix[r][c] == _KsatL2->nodata)
	  throw string("Layer 2 Ksat map contains no data values inside the valid domain...\n");
	if (_KsatL3->matrix[r][c] == _KsatL3->nodata)
	  throw string("Layer 3 Ksat map contains no data values inside the valid domain...\n");
	if (_KsatL2->matrix[r][c] <= 0)
	  throw string("Layer 2 Ksat map contains negative or zero values inside the valid domain...\n");
	if (_KsatL3->matrix[r][c] <= 0)
	  throw string("Layer 3 Ksat map contains negative or zero values inside the valid domain...\n");
      }

      if(ctrl.toggle_Poros!=1){
	if (_porosityL1->matrix[r][c] == _porosityL1->nodata)
	  throw string("Layer 1 porosity map contains no data values inside the valid domain...\n");
	if (_porosityL1->matrix[r][c] <= 0)
	  throw string("Layer 1 porosity map contains negative or zero values inside the valid domain...\n");
      }
      if(ctrl.toggle_Poros==1){ // Exponential porosity profile
	if (_porosity0->matrix[r][c] == _porosity0->nodata)
	  throw string("Top-of-profile porosity map contains no data values inside the valid domain...\n");
	if (_porosity0->matrix[r][c] <= 0)
	  throw string("Top-of-profile porosity map contains negative or zero values inside the valid domain...\n");
	if (_kporos->matrix[r][c] == _kporos->nodata)
	  throw string("Porosity profile map contains no data values inside the valid domain...\n");
	if (_kporos->matrix[r][c] <= 0)
	  throw string("Porosity profile map contains negative or zero values inside the valid domain...\n");
      }
      if(ctrl.toggle_Poros==2){ // Layer-defined phis
	if (_porosityL2->matrix[r][c] == _porosityL2->nodata)
	  throw string("Layer 2 porosity map contains no data values inside the valid domain...\n");
	if (_porosityL3->matrix[r][c] == _porosityL3->nodata)
	  throw string("Layer 3 porosity map contains no data values inside the valid domain...\n");
	if (_porosityL2->matrix[r][c] <= 0)
	  throw string("Layer 2 porosity map contains negative or zero values inside the valid domain...\n");
	if (_porosityL3->matrix[r][c] <= 0)
	  throw string("Layer 3 porosity map contains negative or zero values inside the valid domain...\n");
      }

      if (_psi_ae->matrix[r][c] == _psi_ae->nodata)
	throw string("Air entry pressure map contains no data values inside the valid domain...\n");
      if (_psi_ae->matrix[r][c] <= 0)
	throw string("Air entry pressure map contains negative or zero values inside the valid domain...\n");

      if (_BClambda->matrix[r][c] == _BClambda->nodata)
	throw string("Brooks and Cory lambda map contains no data values inside the valid domain...\n");
      if (_BClambda->matrix[r][c] <= 2) {
	_BClambda->matrix[r][c] = 2;
	throw string("WARNING: Brooks and Corey lambda map is too small, switching to minimum value of 2...\n");
      }

      // Initial check of theta_r using L1 (it is the same for all 3 layers,
      // at least before comparing to _porosityL*) 
      if (_theta_rL1->matrix[r][c] == _theta_rL1->nodata)
	throw string("residual moisture map (L1) contains no data values inside the valid domain...\n");
      if (_theta_rL1->matrix[r][c] <= 0)
	throw string("Residual moisture map (L1) contains negative or zero values inside the valid domain...\n");
      //if (_theta_rL1->matrix[r][c] > _porosityL1->matrix[r][c])
      //throw string("Residual soil moisture (L1) map is larger than porosity inside the valid domain...\n");

      // Compare theta_rL* and porosityL* (and nudge if needed, to avoid inconsistencies when
      // exponential profile or layer-specific toggle are on)
      if (_theta_rL1->matrix[r][c] > 0.5 * _porosityL1->matrix[r][c]) {
	string e("WARNING: Topsoil residual soil moisture is > 0.5*topsoil porosity, let's tone it down...\n");
	//cout << e;
	_theta_rL1->matrix[r][c] = 0.5 * _porosityL1->matrix[r][c];
      }
      if (_theta_rL2->matrix[r][c] > 0.5 * _porosityL2->matrix[r][c]) {
	string e("WARNING: Residual soil moisture in layer 2 is > 0.5*(L2's) porosity, let's tone it down...\n");
	//cout << e;
	_theta_rL2->matrix[r][c] = 0.5 * _porosityL2->matrix[r][c];
	// 
      }
      if (_theta_rL3->matrix[r][c] > 0.5 * _porosityL3->matrix[r][c]) {
	string e("WARNING: Residual soil moisture in layer 3 is > 0.5*(L3's) porosity, let's tone it down...\n");
	//cout << e;
	_theta_rL3->matrix[r][c] = 0.5 * _porosityL3->matrix[r][c];
	// 
      }

      if (_soildepth->matrix[r][c] == _soildepth->nodata)
	throw string("soil depth map contains no data values inside the valid domain...\n");

      if (_soildepth->matrix[r][c] < 0)
	throw string("soil depth map contains negative values inside the valid domain...\n");

      if (_paramWc->matrix[r][c] == _paramWc->nodata)
	throw string("Soil parameter Wc map contains no data values inside the valid domain...\n");

      if (_paramWc->matrix[r][c] <= 0)
	throw string("Soil parameter Wc map contains negative or zero values inside the valid domain...\n");

      if (_paramWp->matrix[r][c] == _paramWp->nodata)
	throw string("Soil parameter Wp map contains no data values inside the valid domain...\n");

      if (_paramWp->matrix[r][c] <= 0)
	throw string("Soil parameter Wp map contains negative or zero values inside the valid domain...\n");

      if (_meltCoeff->matrix[r][c] == _meltCoeff->nodata)
	throw string("Snowmelt coefficient map contains no data values inside the valid domain...\n");

      if (_meltCoeff->matrix[r][c] <= 0)
	throw string("Snowmelt coefficient map contains negative or zero values inside the valid domain...\n");

      if (_snow->matrix[r][c] == _snow->nodata)
	throw string("Initial SWE map contains no data values inside the valid domain...\n");

      if (_snow->matrix[r][c] < 0)
	throw string("Initial SWE map contains negative values inside the valid domain...\n");

      if (_albedo->matrix[r][c] == _albedo->nodata)
	throw string("Albedo map contains no data values inside the valid domain...\n");

      if (_albedo->matrix[r][c] <= 0)
	throw string("Albedo map contains negative or zero values inside the valid domain...\n");

      if (_emiss_surf->matrix[r][c] == _emiss_surf->nodata)
	throw string("Surface emissivity map contains no data values inside the valid domain...\n");

      if (_emiss_surf->matrix[r][c] <= 0)
	throw string("Surface emissivity map contains negative or zero values inside the valid domain...\n");

      if (_soil_dry_heatcap->matrix[r][c] == _soil_dry_heatcap->nodata)
	throw string("Soil dry heat capacity map contains no data values inside the valid domain...\n");

      if (_soil_dry_heatcap->matrix[r][c] <= 0)
	throw string("Soil dry heat capacity map contains negative or zero values inside the valid domain...\n");

      if (_soil_dry_thermcond->matrix[r][c] == _soil_dry_thermcond->nodata)
	throw string("Dry soil thermal conductivity map contains no data values inside the valid domain...\n");

      if (_soil_dry_thermcond->matrix[r][c] <= 0)
	throw string("Dry soil thermal conductivity map contains negative or zero values inside the valid domain...\n");

      if (_dampdepth->matrix[r][c] == _dampdepth->nodata)
	throw string("Thermal soil damping depth map contains no data values inside the valid domain...\n");

      if (_dampdepth->matrix[r][c] <= 0)
	throw string("Thermal soil damping depth map contains negative or zero values inside the valid domain...\n");

      if (_Temp_d->matrix[r][c] == _Temp_d->nodata)
	throw string("Soil temp at damp depth map contains no data values inside the valid domain...\n");

      // Soil moisture L1
      if (_soilmoist1->matrix[r][c] == _soilmoist1->nodata)
	throw string("Initial soil moisture map contains no data values inside the valid domain...\n");
      if (_soilmoist1->matrix[r][c] <= 0)
	throw string("Topsoil Initial soil moisture map contains negative or zero values inside the valid domain...\n");
      if (_soilmoist1->matrix[r][c] > _porosityL1->matrix[r][c]) {
	string e("WARNING: Topsoil Initial soil moisture map is larger than porosity, let's bring it down...\n");
	//cout << e;
	_soilmoist1->matrix[r][c] = 0.75*_porosityL1->matrix[r][c];
      }
      if (_soilmoist1->matrix[r][c] < _theta_rL1->matrix[r][c]) {
	string e("WARNING: Topsoil Initial soil moisture map is lower than residual moisture, let's jump it up...\n");
	//cout << e;
	_soilmoist1->matrix[r][c] = 1.5*_theta_rL1->matrix[r][c];
      }

      // Soil moisture L2
      if (_soilmoist2->matrix[r][c] == _soilmoist2->nodata)
	throw string("Initial soil moisture map in layer 2 contains no data values inside the valid domain...\n");
      if (_soilmoist2->matrix[r][c] <= 0)
	throw string("Initial soil moisture map in layer 2 contains negative or zero values inside the valid domain...\n");
      if (_soilmoist2->matrix[r][c] <= _theta_rL2->matrix[r][c])
	throw string("Initial soil moisture map in soil layer 2 is lower or equal than residual moisture inside the valid domain...\n");
      if (_soilmoist2->matrix[r][c] > _porosityL2->matrix[r][c]) {
	string e("WARNING: Initial soil moisture in layer 2 is larger than porosity, let's bring it down...\n");
	//cout << e;
	_soilmoist2->matrix[r][c] = 0.75*_porosityL2->matrix[r][c];
      }
      if (_soilmoist2->matrix[r][c] < _theta_rL2->matrix[r][c]) {
	string e("WARNING: Initial soil moisture map in layer 2 is lower than residual moisture, let's jump it up...\n");
	//cout << e;
	_soilmoist2->matrix[r][c] = 1.5*_theta_rL2->matrix[r][c];
      }

      // Soil moisture L3
      if (_soilmoist3->matrix[r][c] == _soilmoist3->nodata)
	throw string("Initial soil moisture map in layer 3 contains no data values inside the valid domain...\n");
      if (_soilmoist3->matrix[r][c] <= 0)
	throw string("Initial soil moisture map in layer 3 contains negative or zero values inside the valid domain...\n");
      if (_soilmoist3->matrix[r][c] > _porosityL3->matrix[r][c]) {
	string e("WARNING: Initial soil moisture in layer 3 is larger than porosity, let's bring it down...\n");
	//cout << e;
	_soilmoist3->matrix[r][c] = 0.75*_porosityL3->matrix[r][c];
	//
      }
      if (_soilmoist3->matrix[r][c] < _theta_rL3->matrix[r][c]) {
	string e("WARNING: Initial soil moisture map in layer 3 is lower than residual moisture, let's jump it up...\n");
	// cout << e;
	_soilmoist3->matrix[r][c] = 1.5*_theta_rL3->matrix[r][c];
	//
      }
      
      if (_channelwidth->matrix[r][c] < 0)
	throw string("The channel width map contains negative values\n");

    } catch (string &e) {
      std::cout << e << std::endl;
      std::cout << "In row " << r << " col " << c << std::endl;
#pragma omp atomic
      ++excep_thrown;
    }

  }//for

  if (excep_thrown != 0)
    throw;
}
