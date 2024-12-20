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
 * MassBalanceError.cpp
 *
 *  Created on: Mar 18, 2010
 *      Author: Marco Maneta
 */

#include "Budget.h"

void Budget::MassBalanceError(const Control *ctrl)
{
  double inputs = 0.0;
  double outputs = 0.0;
  double ds = 0.0;

  inputs = precipitation + initsnowpack + initponding + initL1 + initL2 + initL3 ; //+ initGW;
  outputs = evaporationS + evaporationI + transpiration + ovlndflow + gwtrflow + leakage;
  ds = canopy + snowpack + ponding + soilL1 + soilL2 + soilL3 ; //+ grndwater;

  if(inputs>0)
    MBErr = 100/inputs*(inputs-outputs - ds);
  else
    MBErr = 0;

  // Tracking -------------
  
  //Deuterium balance
  if(ctrl->sw_2H){
    
    double inputs_d2H = 0.0;
    double outputs_d2H = 0.0;
    double ds_d2H = 0.0;
    
    inputs_d2H = precipitation_d2H + initsnowpack_d2H + initponding_d2H + 
      initL1_d2H + initL2_d2H + initL3_d2H ; //+ initGW_d2H;
    
    outputs_d2H = evaporationS_d2H + evaporationI_d2H + transpiration_d2H +
      ovlndflow_d2H + gwtrflow_d2H + leakage_d2H;
    
    ds_d2H = canopy_d2H + snowpack_d2H + ponding_d2H + soilL1_d2H + soilL2_d2H + soilL3_d2H ; //+ 
    //grndwater_d2H;
    
    if(inputs>0) 
      MBErr_d2H = 100/inputs_d2H*(inputs_d2H - outputs_d2H - ds_d2H);
    else 
      MBErr_d2H = 0;
    
  }
	
  // Oxygen 18 balance
  if(ctrl->sw_18O){

    double inputs_d18O = 0.0;
    double outputs_d18O = 0.0;
    double ds_d18O = 0.0;
  
    inputs_d18O = precipitation_d18O + initsnowpack_d18O + initponding_d18O + 
      initL1_d18O + initL2_d18O + initL3_d18O ; //+ initGW_d18O;
    
    outputs_d18O = evaporationS_d18O + evaporationI_d18O + transpiration_d18O +
      ovlndflow_d18O + gwtrflow_d18O + leakage_d18O;
	  
    ds_d18O = canopy_d18O + snowpack_d18O + ponding_d18O + soilL1_d18O + soilL2_d18O + soilL3_d18O ;
    //+ grndwater_d18O;	  
	  
    if(inputs>0) 
      MBErr_d18O = 100/inputs_d18O*(inputs_d18O - outputs_d18O - ds_d18O);
    else 
      MBErr_d18O = 0;
  }

  // Chloride balance
  if(ctrl->sw_Cl){

    double inputs_cCl = 0.0;
    double outputs_cCl = 0.0;
    double ds_cCl = 0.0;
  
    inputs_cCl = precipitation_cCl + initsnowpack_cCl + initponding_cCl + 
      initL1_cCl + initL2_cCl + initL3_cCl ; //+ initGW_cCl;
	  
    outputs_cCl = ovlndflow_cCl + gwtrflow_cCl + leakage_cCl;
	  
    ds_cCl = canopy_cCl + snowpack_cCl + ponding_cCl + soilL1_cCl + soilL2_cCl + soilL3_cCl ;
    //+ grndwater_cCl;	  
	  
    if(inputs > 0) 
      MBErr_cCl = 100/inputs_cCl*(inputs_cCl - outputs_cCl - ds_cCl);
    else 
      MBErr_cCl = 0;
  }

  // Age mass balance
  if(ctrl->sw_Age){
	  
    double inputs_Age = 0.0;
    double outputs_Age = 0.0;
    double ds_Age = 0.0;
    
    // Ageing initial storage
    initsnowpack_Age += initsnowpack * ctrl->dt / 86400 ;
    initponding_Age += initponding * ctrl->dt / 86400;
    initL1_Age += initL1 * ctrl->dt / 86400 ;
    initL2_Age += initL2 * ctrl->dt / 86400 ;
    initL3_Age += initL3 * ctrl->dt / 86400 ;
    //initGW_Age += initGW * ctrl->dt / 86400;	  

    inputs_Age = precipitation_Age + // gradually aging precip (incremented in TotalPrecipitation.cpp)
      initsnowpack_Age + initponding_Age + 
      initL1_Age + initL2_Age + initL3_Age ; //+ initGW_Age;

    outputs_Age = evaporationS_Age + evaporationI_Age + transpiration_Age +
      ovlndflow_Age + gwtrflow_Age + leakage_Age;
	  
    ds_Age = canopy_Age + snowpack_Age + ponding_Age + 
      soilL1_Age + soilL2_Age + soilL3_Age ; //+ grndwater_Age;
	  
    if(inputs_Age > RNDOFFERR) 
      MBErr_Age = 100/inputs_Age*(inputs_Age - outputs_Age - ds_Age);
    else 
      MBErr_Age = 0;
  }

}
