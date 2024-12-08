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
 * MixingV_snow.cpp
 *
 *  Created on: Jun 21, 2017
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"

void Tracking::MixingV_snow(Atmosphere &atm, Basin &bsn, Control &ctrl,
			    double &h, double &dh, int r, int c) //time step
{
//DEBUG  if(r==8 and c==1)
//    cout << endl << "----" ;

  double h_eff; // Effective snowpack (before snowfall) SWE used for mixing
  double dh_eff; // Effective snowfall-to-snowpack used for mixing
  double snow_in = bsn.getFluxCnptoSnow()->matrix[r][c] * ctrl.dt ;

  h_eff = h - snow_in + dh ;

  // - in snowpack (snowfall in + snowmelt out),
  // considering that snowmelt "flushes" the most recent snowfall first, without mixing
  /*if(h_eff < 0)
    cout << "ERROR in snowpack mixing calculation!!" << endl ;*/

  //  No initial snow - no need to mix
  if(h_eff < RNDOFFERR){

//  if(r==1 and c==1)
  //    cout << "No initial snow (h_eff=0): h_t=" << h_eff << ", h_t+1=" << h <<", snowfall=" << snow_in << ", melt=" << dh << endl ;

    if(ctrl.sw_2H){
      _d2Hsnowmelt->matrix[r][c] = dh > RNDOFFERR ? atm.getd2Hprecip()->matrix[r][c] : sqrt(-1);
      if(snow_in > RNDOFFERR)
				_d2Hsnowpack->matrix[r][c] = atm.getd2Hprecip()->matrix[r][c] ;
    }

    if(ctrl.sw_18O){
      _d18Osnowmelt->matrix[r][c] =  dh > RNDOFFERR ? atm.getd18Oprecip()->matrix[r][c] : sqrt(-1);
      if(snow_in > RNDOFFERR)
	_d18Osnowpack->matrix[r][c] = atm.getd18Oprecip()->matrix[r][c];
    }

    if(ctrl.sw_Cl){
      _cClsnowmelt->matrix[r][c] =  dh > RNDOFFERR ? atm.getcClprecip()->matrix[r][c] : sqrt(-1);
      if(snow_in > RNDOFFERR)
	_cClsnowpack->matrix[r][c] = atm.getcClprecip()->matrix[r][c];
    }

    if(ctrl.sw_Age){
      _Agesnowmelt->matrix[r][c] = dh > RNDOFFERR ? 0.0 : sqrt(-1);
      if(snow_in > RNDOFFERR)
	_Agesnowpack->matrix[r][c] = 0.0;
    }

  } else if (h > RNDOFFERR and snow_in > dh){

    // Case where there was initially a snowpack
    // if there is more snowfall than snowmelt:
    // --> snowpack mixed, snowmelt has snowfall signature

    dh_eff = snow_in - dh;

// DEBUG   if(r==8 and c==1)
  //    cout << "Remaining snowpack, snowfall > snowmelt: h_t=" << h_eff << ", h_t+1=" << h <<", snowfall=" << snow_in << ", melt=" << dh << endl ;


    if(ctrl.sw_2H){
      _d2Hsnowmelt->matrix[r][c] = dh > RNDOFFERR ? atm.getd2Hprecip()->matrix[r][c] : sqrt(-1);
      _d2Hsnowpack->matrix[r][c] = InputMix(h_eff, _d2Hsnowpack->matrix[r][c],
					    dh_eff, atm.getd2Hprecip()->matrix[r][c]);
    }

    if(ctrl.sw_18O){
      _d18Osnowmelt->matrix[r][c] = dh > RNDOFFERR ? atm.getd18Oprecip()->matrix[r][c] : sqrt(-1);
      _d18Osnowpack->matrix[r][c] = InputMix(h_eff, _d18Osnowpack->matrix[r][c],
					    dh_eff, atm.getd18Oprecip()->matrix[r][c]);
    }

    if(ctrl.sw_Cl){
      _cClsnowmelt->matrix[r][c] = dh > RNDOFFERR ? atm.getcClprecip()->matrix[r][c] : sqrt(-1);
      _cClsnowpack->matrix[r][c] = InputMix(h_eff, _cClsnowpack->matrix[r][c],
					     dh_eff, atm.getcClprecip()->matrix[r][c]);
    }

    if(ctrl.sw_Age){
      _Agesnowmelt->matrix[r][c] = dh > RNDOFFERR ? 0.0 : sqrt(-1);
      _Agesnowpack->matrix[r][c] = InputMix(h_eff, _Agesnowpack->matrix[r][c], dh_eff, 0.0);
    }

  } else { // snow_in < dh

    // Case where there is more snowmelt than snowfall:
    // --> no mixing in snowpack, snowmelt has mixed signature
    dh_eff = dh - snow_in;

// DEBUG    if(r==8 and c==1)
//   cout << "Remaining snowpack, snowfall < snowmelt: h_t=" << h_eff << ", h_t+1=" << h <<", snowfall=" << snow_in << ", melt=" << dh << endl ;

    // Snowpack: no change (all recent snow has melted)
    // Snowmelt: mix of snowpack and throughfall
    if(ctrl.sw_2H)
      _d2Hsnowmelt->matrix[r][c] = InputMix(dh_eff, _d2Hsnowpack->matrix[r][c],
					    snow_in, atm.getd2Hprecip()->matrix[r][c]);
    if(ctrl.sw_18O)
      _d18Osnowmelt->matrix[r][c] =InputMix(dh_eff, _d18Osnowpack->matrix[r][c],
					    snow_in, atm.getd18Oprecip()->matrix[r][c]);
    if(ctrl.sw_Cl)
      _cClsnowmelt->matrix[r][c] = InputMix(dh_eff, _cClsnowpack->matrix[r][c],
					    snow_in, atm.getcClprecip()->matrix[r][c]);
    if(ctrl.sw_Age)
      _Agesnowmelt->matrix[r][c] = InputMix(dh_eff, _Agesnowpack->matrix[r][c], snow_in, 0.0);

  }

  // If no snowpack, nan values for tracers
  if(abs(h) < RNDOFFERR){
    if(ctrl.sw_2H)
      _d2Hsnowpack->matrix[r][c] = sqrt(-1);
    if(ctrl.sw_18O)
      _d18Osnowpack->matrix[r][c] = sqrt(-1);
    if(ctrl.sw_Cl)
      _cClsnowpack->matrix[r][c] = sqrt(-1);
    if(ctrl.sw_Age)
      _Agesnowpack->matrix[r][c] = sqrt(-1);
  }

/*//DEBUG if(r==8 and c==1){
    if(ctrl.sw_2H)
      cout << "d2H_snowpack=" << _d2Hsnowpack->matrix[r][c] << 
      ", d2H_snowmelt=" << _d2Hsnowmelt->matrix[r][c] << 
      ", d2H_precip=" << atm.getd2Hprecip()->matrix[r][c] << endl ;
    if(ctrl.sw_Age)
        cout << "Age_snowpack=" << _Agesnowpack->matrix[r][c] << 
      ", Age_snowmelt=" << _Agesnowmelt->matrix[r][c] << endl ;*/
    

}
