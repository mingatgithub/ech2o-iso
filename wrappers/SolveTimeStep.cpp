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
 * SolveTimeStep.cpp
 *
 *  Created on: Aug 2, 2010
 *      Author: Marco.Maneta, Sylvain Kuppel
 */
#include <iostream>
#include "Sativa.h"


int SolveTimeStep(){

  //oBasin->UpdateSnowPack(*oAtmosphere, *oControl);
  oBasin->SolveCanopyFluxes(*oAtmosphere, *oControl, *oTracking);

  oBasin->SolveSurfaceFluxes(*oAtmosphere, *oControl, *oTracking);

  // if(oControl->sw_veg_dyn)
  // commented out the "if", otherwise GPP and NPP are not calculated in "static"
  // (toggle_veg_dyn=0) and "LAI-forced" (toggle_veg_dyn=2) modes??
    oBasin->CalculateGrowForest(*oAtmosphere, *oControl);

  //oBasin->DailySurfaceRouting(*oAtmosphere, *oControl);
  //if(oControl->toggle_soil_water_profile < 2)
  oBasin->DailyGWRouting(*oAtmosphere, *oControl, *oTracking);

  // Groundwater volume and table depth (if necessary)
  oBasin->CalculateGrndWaterVol();
  if(oControl->Rep_WaterTableDepth || oControl->RepTs_WaterTableDepth)
    oBasin->CalculateWaterTableDepth(*oControl);

  // Saturated area
  oBasin->CalculateSatArea(*oControl);
  
 // Tracking	  
  if(oControl->sw_trck){
    // If Two-pore...
    if(oControl->sw_TPD){
      // Calculate the soil-averaged tracer values
      oTracking->CalcTPDtoLayers(*oBasin, *oControl);
      // Calculate the relative fraction of tightly-bound domain
      oBasin->CalcFracMobileWater();
    }

    if(oControl->sw_2H){ 
      if(oControl->Rep_d2HsoilUp || oControl->RepTs_d2HsoilUp)
	oTracking->CalcTrcksoil_12(*oBasin, 1);
      if(oControl->Rep_d2HsoilAv || oControl->RepTs_d2HsoilAv)
	oTracking->CalcTrcksoil_Av(*oBasin, 1);
      // In any case, get the groundwater signature
      oTracking->CalcTrcksoil_GW(*oBasin, 1);
    }
    
    if(oControl->sw_18O){
      if(oControl->Rep_d18OsoilUp || oControl->RepTs_d18OsoilUp)
	oTracking->CalcTrcksoil_12(*oBasin, 2);
      if(oControl->Rep_d18OsoilAv || oControl->RepTs_d18OsoilAv)
	oTracking->CalcTrcksoil_Av(*oBasin, 2);
      // In any case, get the groundwater signature
      oTracking->CalcTrcksoil_GW(*oBasin, 2);
    }

    if(oControl->sw_Cl){
      // Reported quantities
      if(oControl->Rep_cClsoilUp || oControl->RepTs_cClsoilUp)
	oTracking->CalcTrcksoil_12(*oBasin, 3);
      if(oControl->Rep_cClsoilAv || oControl->RepTs_cClsoilAv)
	oTracking->CalcTrcksoil_Av(*oBasin, 3);
      // In any case, get the groundwater signature
      oTracking->CalcTrcksoil_GW(*oBasin, 3);
    }
    
    if(oControl->sw_Age){
      // Increment age by one time step duration
      oTracking->IncrementAge(*oBasin, *oControl);
      // Reported quantities
      if(oControl->Rep_AgesoilUp || oControl->RepTs_AgesoilUp)
	oTracking->CalcTrcksoil_12(*oBasin, 4);
      if(oControl->Rep_AgesoilAv || oControl->RepTs_AgesoilAv)
	oTracking->CalcTrcksoil_Av(*oBasin, 4);
      // In any case, get the groundwater signature
      oTracking->CalcTrcksoil_GW(*oBasin, 4);
    }
  }
  
  return EXIT_SUCCESS;
}
