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
 * MixingV_down.cpp
 *
 *  Created on: Jun 21, 2017
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"

void Tracking::MixingV_down(Basin &bsn, Control &ctrl,
			    double &d1, double &d2, double &d3, double &fc,
			    int r, int c,
			    int step) // 0 for SolveSurface, 1 for GWrouting
{

  double theta1_old = bsn.getSoilMoist1()->matrix[r][c];
  double theta2_old = bsn.getSoilMoist2()->matrix[r][c];
  double theta3_old = bsn.getSoilMoist3()->matrix[r][c];

  double LattoSrf = 0;
  double SnowtoSrf = 0;
  double FinSrf = 0;
  double pond_old = 0;
  double d2Hin = 0 ;
  double d18Oin = 0;
  double cClin = 0;
  double Agein = 0;

  double SrftoL1 = bsn.getFluxSrftoL1()->matrix[r][c] * ctrl.dt ;
  double L1toL2 = bsn.getFluxL1toL2()->matrix[r][c] * ctrl.dt ;
  double L2toL3 = bsn.getFluxL2toL3()->matrix[r][c] * ctrl.dt ;
  double Leak = bsn.getFluxLeakage()->matrix[r][c] * ctrl.dt ;

  int mixmod = ctrl.toggle_mix;

  double theta_MW1 = 0;
  double theta_MW2 = 0;
  //double theta_r = 0;
  //double porosity = 0;
  double SrftoTB1 = 0;
  double SrftoMW1 = 0;
  double MW1toTB2 = 0;
  double MW1toMW2 = 0;

  double d_old = 0;

  // Initialize
  if(step == 0){
	  if(ctrl.sw_2H){
	  	_d2Hsnowmelt->matrix[r][c] = sqrt(-1) ;
      _d2HGWtoChn->matrix[r][c] = sqrt(-1) ;
      _d2HSrftoChn->matrix[r][c] = sqrt(-1) ;
      _d2HRecharge->matrix[r][c] = sqrt(-1) ;
    }
  	if(ctrl.sw_18O){
			_d18Osnowmelt->matrix[r][c] = sqrt(-1) ;
    	_d18OGWtoChn->matrix[r][c] = sqrt(-1) ;
    	_d18OSrftoChn->matrix[r][c] = sqrt(-1) ;
    	_d18ORecharge->matrix[r][c] = sqrt(-1) ;
  	}
 		if(ctrl.sw_Cl){
			_cClsnowmelt->matrix[r][c] = sqrt(-1) ;
    	_cClGWtoChn->matrix[r][c] = sqrt(-1) ;
    	_cClSrftoChn->matrix[r][c] = sqrt(-1) ;
    	_cClRecharge->matrix[r][c] = sqrt(-1) ;
  	}
  	if(ctrl.sw_Age){
			_Agesnowmelt->matrix[r][c] = sqrt(-1) ;
    	_AgeGWtoChn->matrix[r][c] = sqrt(-1) ;
    	_AgeSrftoChn->matrix[r][c] = sqrt(-1) ;
    	_AgeRecharge->matrix[r][c] = sqrt(-1) ;
  	}
	
  }
  
  if(ctrl.sw_TPD){
    //theta_r = bsn.getSoilMoistR()->matrix[r][c];
    theta_MW1 = bsn.getMoistureMW1()->matrix[r][c];
    theta_MW2 = bsn.getMoistureMW2()->matrix[r][c];
    //porosity = bsn.getPorosity()->matrix[r][c];
    // Routing to tightly-bound domain: fraction from relative volume of "inactive"-pore domain,
    // limited by the "available space" in tightly-bound domain max(0,d1*(theta_MW-theta1_old))
    SrftoTB1 = std::min<double>(SrftoL1,//*(theta_MW-theta_r)/(porosity-theta_r),
				std::max<double>(0,d1*(theta_MW1-theta1_old)));
    // Remainder directly goes to mobile water
    SrftoMW1 = std::max<double>(0,SrftoL1-SrftoTB1);
    // L2: Percolation only from mobile water in L1
    MW1toTB2 = std::min<double>(L1toL2,//*(theta_MW-theta_r)/(porosity-theta_r),
				std::max<double>(0,d2*(theta_MW2-theta2_old)));
    MW1toMW2 = std::max<double>(0,L1toL2-MW1toTB2);
    // L3: Percolation only from mobile water water in L2 (MW2toL3 = L2toL3).
  }

  // Surface (if reinfiltrating) ----------------------------------------------------
  if(step == 1){

    // Fluxes in
    SnowtoSrf = bsn.getFluxSnowmelt()->matrix[r][c] * ctrl.dt;
    LattoSrf = bsn.getFluxLattoSrf()->matrix[r][c] * ctrl.dt ;
    FinSrf = SnowtoSrf + LattoSrf;

    if(FinSrf > RNDOFFERR) {

      pond_old = bsn.getPondingWater()->matrix[r][c] - FinSrf;

      if(ctrl.sw_2H){
	d_old = _d2Hsurface->matrix[r][c];
	d2Hin = 0 ; // add component only if flux > 0 (to avoid nan propagation)
	if(LattoSrf > RNDOFFERR)
	  d2Hin += _Fd2HLattoSrf->matrix[r][c] / FinSrf ;
	if(SnowtoSrf > RNDOFFERR)
	  d2Hin += SnowtoSrf*_d2Hsnowmelt->matrix[r][c] / FinSrf ;

	_d2Hsurface->matrix[r][c] = pond_old > RNDOFFERR ?
	  InOutMix(pond_old, _d2Hsurface->matrix[r][c], FinSrf, d2Hin, SrftoL1, mixmod) : d2Hin ;

	  //cout << "d2Hsurf_t=" << d_old << ", d2Hsurf_t+1=" << _d2Hsurface->matrix[r][c] << endl ;

	/*if(abs(_d2Hsurface->matrix[r][c])>100)
	  cout << r << " " << c << "| d2H-MixingRouting1 " <<
	  "| dSrf_new:" << _d2Hsurface->matrix[r][c] << "|dSrf_old:" << d_old <<
	  "|d2Hin:" << d2Hin << endl;*/
      }

      if(ctrl.sw_18O){
	// Update surface
	d_old = _d18Osurface->matrix[r][c];
	d18Oin = 0 ; // add component only if flux > 0 (to avoid nan propagation)
	if(LattoSrf > RNDOFFERR)
	  d18Oin += _Fd18OLattoSrf->matrix[r][c] / FinSrf ;
	if(SnowtoSrf > RNDOFFERR)
	  d18Oin += SnowtoSrf*_d18Osnowmelt->matrix[r][c] / FinSrf ;

	_d18Osurface->matrix[r][c] = pond_old > RNDOFFERR ?
	  InOutMix(pond_old, _d18Osurface->matrix[r][c], FinSrf, d18Oin, SrftoL1, mixmod) : d18Oin;
	/*if(abs(_d18Osurface->matrix[r][c]>100))
	  cout << r << " " << c << "| d18O-MixingRouting1 " <<
	  "| dSrf_new:" << _d18Osurface->matrix[r][c] << "|dSrf_old:" << d_old <<
	  "|d18Oin:" << d18Oin << endl;*/
      }

      if(ctrl.sw_Cl){
	// Update surface
	cClin = 0 ; // add component only if flux > 0 (to avoid nan propagation)
	if(LattoSrf > RNDOFFERR)
	  cClin += _FcClLattoSrf->matrix[r][c] / FinSrf ;
	if(SnowtoSrf > RNDOFFERR)
	  cClin += SnowtoSrf*_cClsnowmelt->matrix[r][c] / FinSrf ;

	_cClsurface->matrix[r][c] = pond_old > RNDOFFERR ?
	  InOutMix(pond_old, _cClsurface->matrix[r][c], FinSrf, cClin, SrftoL1, mixmod) : cClin;
      }

      if(ctrl.sw_Age){
	Agein = 0 ; // add component only if flux > 0 (to avoid nan propagation)
		d_old = _Agesurface->matrix[r][c];
	if(LattoSrf > RNDOFFERR)
	  Agein += _FAgeLattoSrf->matrix[r][c] / FinSrf ;
	if(SnowtoSrf > RNDOFFERR)
	  Agein += SnowtoSrf*_Agesnowmelt->matrix[r][c] / FinSrf ;

	_Agesurface->matrix[r][c] = pond_old > RNDOFFERR ?
	  InOutMix(pond_old, _Agesurface->matrix[r][c], FinSrf, Agein, SrftoL1, mixmod) : Agein ;

	  //cout << "Agesurf_t=" << d_old << ", Agesurf_t+1=" << _Agesurface->matrix[r][c] << endl ;
      }
    }
  }

  // === Only if initial infitlration or reinfilt + non-channel
  if(step == 0 or (step==1 and ctrl.sw_reinfilt and
		   !(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0))){

    // Layer 1 ------------------------------------------------------------------------

    // If two-pore domain activated, and different groundwater origin (MW2 instead of soil2)
    if(ctrl.sw_TPD and SrftoL1>RNDOFFERR){
      if(ctrl.sw_2H){
	d_old = _d2H_TB1->matrix[r][c];
	_d2H_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1,
					  _d2H_TB1->matrix[r][c],
					  SrftoTB1, _d2Hsurface->matrix[r][c]);

	if(abs(_d2H_TB1->matrix[r][c])>300)
	  cout << r << " " << c << "| d2H " <<
	    "| dTB1_new:" << _d2H_TB1->matrix[r][c] << "| dTB1_old:" << d_old <<
	    "| dSrf:" << _d2Hsurface->matrix[r][c] << "| SrftoTB1:" << SrftoTB1 << endl;

	d_old = _d2H_MW1->matrix[r][c];
	_d2H_MW1->matrix[r][c] = std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR ?
	  InOutMix(std::max<double>(0,theta1_old-theta_MW1)*d1, _d2H_MW1->matrix[r][c],
		   SrftoMW1, _d2Hsurface->matrix[r][c], L1toL2, mixmod) :
	  _d2Hsurface->matrix[r][c] ;

	if(abs(_d2H_MW1->matrix[r][c])>300) // or (r==40 and c==92))
	  cout << r << " " << c << "| d2H " <<
	    "| dMW1_new:" << _d2H_MW1->matrix[r][c] << "| dMW1_old:" << d_old <<
	    "| dSrf:" << _d2Hsurface->matrix[r][c] <<
	    "| SrftoMW1:" << SrftoMW1 << "| MW1toL2:" << L1toL2 << endl;
      }

      if(ctrl.sw_18O){
	d_old = _d18O_TB1->matrix[r][c];
	_d18O_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1,
					   _d18O_TB1->matrix[r][c],
					   SrftoTB1, _d18Osurface->matrix[r][c]);

	if(abs(_d18O_TB1->matrix[r][c])>100)
	  cout << r << " " << c << "| d18O " <<
	    "| dTB1_new:" << _d18O_TB1->matrix[r][c] << "| dTB1_old:" << d_old <<
	    "| dSrf:" << _d18Osurface->matrix[r][c] << "| SrftoTB1:" << SrftoTB1 << endl;

	d_old = _d18O_MW1->matrix[r][c];

	_d18O_MW1->matrix[r][c] = std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR ?
	  InOutMix(std::max<double>(0,theta1_old-theta_MW1)*d1, _d18O_MW1->matrix[r][c],
		   SrftoMW1, _d18Osurface->matrix[r][c], L1toL2, mixmod) :
	  _d18Osurface->matrix[r][c] ;

	if(abs(_d18O_MW1->matrix[r][c])>100) // or (r==40 and c==92))
	  cout << r << " " << c << "| d18O " <<
	    "| dMW1_new:" << _d18O_MW1->matrix[r][c] << "| dMW1_old:" << d_old <<
	    "| dSrf:" << _d18Osurface->matrix[r][c] <<
	    "| SrftoMW1:" << SrftoMW1 << "| MW1toL2:" << L1toL2 << endl;

      }

      if(ctrl.sw_Cl){
	_cCl_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1,
					  _cCl_TB1->matrix[r][c],
					  SrftoTB1, _cClsurface->matrix[r][c]);
      	_cCl_MW1->matrix[r][c] = std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR ?
	  InOutMix(std::max<double>(0,theta1_old-theta_MW1)*d1, _cCl_MW1->matrix[r][c],
		   SrftoMW1, _cClsurface->matrix[r][c], L1toL2, mixmod) :
	  _cClsurface->matrix[r][c] ;

      }

      if(ctrl.sw_Age){
	_Age_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1,
					  _Age_TB1->matrix[r][c],
					  SrftoTB1, _Agesurface->matrix[r][c]);
	_Age_MW1->matrix[r][c] = std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR ?
	  InOutMix(std::max<double>(0,theta1_old-theta_MW1)*d1, _Age_MW1->matrix[r][c],
		   SrftoMW1, _Agesurface->matrix[r][c], L1toL2, mixmod) :
	  _Agesurface->matrix[r][c] ;

      }
    } else if (SrftoL1>RNDOFFERR){ // Soil-averaged
      if(ctrl.sw_2H){
      	d_old = _d2Hsoil1->matrix[r][c];
	_d2Hsoil1->matrix[r][c] = InOutMix(theta1_old*d1, _d2Hsoil1->matrix[r][c],
					   SrftoL1, _d2Hsurface->matrix[r][c], L1toL2, mixmod);

	//cout << "d2Hsoil1_t=" << d_old << ", d2Hsoil_t+1=" << _d2Hsoil1->matrix[r][c] <<
	//			", d2Hsurface=" << _d2Hsurface->matrix[r][c] << ", SrftoL1=" << SrftoL1 << endl ;
}

      if(ctrl.sw_18O)
	_d18Osoil1->matrix[r][c] = InOutMix(theta1_old*d1, _d18Osoil1->matrix[r][c],
					    SrftoL1, _d18Osurface->matrix[r][c], L1toL2, mixmod);

      if(ctrl.sw_Cl)
	_cClsoil1->matrix[r][c] = InOutMix(theta1_old*d1, _cClsoil1->matrix[r][c],
					   SrftoL1, _cClsurface->matrix[r][c], L1toL2, mixmod);

      if(ctrl.sw_Age){
      	d_old = _Agesoil1->matrix[r][c];
	_Agesoil1->matrix[r][c] = InOutMix(theta1_old*d1, _Agesoil1->matrix[r][c],
					   SrftoL1, _Agesurface->matrix[r][c], L1toL2, mixmod);

	//cout << "Agesoil1_t=" << d_old << ", Agesoil_t+1=" << _Agesoil1->matrix[r][c] <<
	//			", Agesurface=" << _Agesurface->matrix[r][c] << ", SrftoL1=" << SrftoL1 << endl ;
}

    }

    // Layer 2 ------------------------------------------------------------------------

    // If two-pore domain activated: only MW1 percolates
    if(ctrl.sw_TPD and L1toL2>RNDOFFERR){
      if(ctrl.sw_2H){
	_d2H_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2,
					  _d2H_TB2->matrix[r][c],
					  MW1toTB2, _d2H_MW1->matrix[r][c]);
	_d2H_MW2->matrix[r][c] = std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR ?
	  InOutMix(std::max<double>(0,theta2_old-theta_MW2)*d2, _d2H_MW2->matrix[r][c],
		   MW1toMW2, _d2H_MW1->matrix[r][c], L2toL3, mixmod) : _d2H_MW1->matrix[r][c] ;

      }
      if(ctrl.sw_18O){
	_d18O_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2,
					   _d18O_TB2->matrix[r][c],
					   MW1toTB2, _d18O_MW1->matrix[r][c]);
	_d18O_MW2->matrix[r][c] = std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR ?
	  InOutMix(std::max<double>(0,theta2_old-theta_MW2)*d2, _d18O_MW2->matrix[r][c],
		   MW1toMW2, _d18O_MW1->matrix[r][c], L2toL3, mixmod) : _d18O_MW1->matrix[r][c] ;
      }

      if(ctrl.sw_Cl){
	_cCl_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2,
					  _cCl_TB2->matrix[r][c],
					  MW1toTB2, _cCl_MW1->matrix[r][c]);
	_cCl_MW2->matrix[r][c] = std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR ?
	  InOutMix(std::max<double>(0,theta2_old-theta_MW2)*d2, _cCl_MW2->matrix[r][c],
		   MW1toMW2, _cCl_MW1->matrix[r][c], L2toL3, mixmod) : _cCl_MW1->matrix[r][c] ;
      }

      if(ctrl.sw_Age){
	_Age_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2,
					  _Age_TB2->matrix[r][c],
					  MW1toTB2, _Age_MW1->matrix[r][c]);
	_Age_MW2->matrix[r][c] = std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR ?
	  InOutMix(std::max<double>(0,theta2_old-theta_MW2)*d2, _Age_MW2->matrix[r][c],
		   MW1toMW2, _Age_MW1->matrix[r][c], L2toL3, mixmod) : _Age_MW1->matrix[r][c] ;
      }
    } else if (L1toL2 > RNDOFFERR){
      if(ctrl.sw_2H)
	_d2Hsoil2->matrix[r][c] = InOutMix(theta2_old*d2, _d2Hsoil2->matrix[r][c],
					   L1toL2, _d2Hsoil1->matrix[r][c], L2toL3, mixmod);

      if(ctrl.sw_18O)
	_d18Osoil2->matrix[r][c] = InOutMix(theta2_old*d2, _d18Osoil2->matrix[r][c],
					    L1toL2, _d18Osoil1->matrix[r][c], L2toL3, mixmod);

      if(ctrl.sw_Cl)
	_cClsoil2->matrix[r][c] = InOutMix(theta2_old*d2, _cClsoil2->matrix[r][c],
					   L1toL2, _cClsoil1->matrix[r][c], L2toL3, mixmod);

      if(ctrl.sw_Age)
	_Agesoil2->matrix[r][c] = InOutMix(theta2_old*d2, _Agesoil2->matrix[r][c],
					   L1toL2, _Agesoil1->matrix[r][c], L2toL3, mixmod);

    }

    // Layer 3 ------------------------------------------------------------------------

    // If two-pore domain activated: only MW2 percolates
    if(ctrl.sw_TPD and L2toL3>RNDOFFERR){
      if(ctrl.sw_2H){
	_d2Hsoil3->matrix[r][c] = InOutMix(theta3_old*d3, _d2Hsoil3->matrix[r][c],
					   L2toL3, _d2H_MW2->matrix[r][c], Leak, mixmod);
	_d2HRecharge->matrix[r][c] = step == 0 ? _d2H_MW2->matrix[r][c] :
	  InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c]*ctrl.dt - L2toL3),
		   _d2HRecharge->matrix[r][c], L2toL3, _d2H_MW2->matrix[r][c]);
      }

      if(ctrl.sw_18O){
	_d18Osoil3->matrix[r][c] = InOutMix(theta3_old*d3, _d18Osoil3->matrix[r][c],
					    L2toL3, _d18O_MW2->matrix[r][c], Leak, mixmod);
	_d18ORecharge->matrix[r][c] = step == 0 ? _d18O_MW2->matrix[r][c] :
	  InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c]*ctrl.dt - L2toL3),
		   _d18ORecharge->matrix[r][c], L2toL3, _d18O_MW2->matrix[r][c]);
      }

      if(ctrl.sw_Cl){
	_cClsoil3->matrix[r][c] = InOutMix(theta3_old*d3, _cClsoil3->matrix[r][c],
					   L2toL3, _cCl_MW2->matrix[r][c], Leak, mixmod);
	_cClRecharge->matrix[r][c] = step == 0 ? _cCl_MW2->matrix[r][c] :
	  InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c]*ctrl.dt - L2toL3),
		   _cClRecharge->matrix[r][c], L2toL3, _cCl_MW2->matrix[r][c]);
      }

      if(ctrl.sw_Age){
	_Agesoil3->matrix[r][c] = InOutMix(theta3_old*d3, _Agesoil3->matrix[r][c],
					   L2toL3, _Age_MW2->matrix[r][c], Leak, mixmod);
	_AgeRecharge->matrix[r][c] = step == 0 ? _Age_MW2->matrix[r][c] :
	  InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c]*ctrl.dt - L2toL3),
		   _AgeRecharge->matrix[r][c], L2toL3, _Age_MW2->matrix[r][c]);
      }

    } else if (L2toL3>RNDOFFERR){

      if(ctrl.sw_2H){
	_d2Hsoil3->matrix[r][c] = InOutMix(theta3_old*d3, _d2Hsoil3->matrix[r][c],
					   L2toL3, _d2Hsoil2->matrix[r][c], Leak, mixmod);

	_d2HRecharge->matrix[r][c] = step == 0 ? _d2Hsoil2->matrix[r][c] :
	  InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c]*ctrl.dt - L2toL3),
		   _d2HRecharge->matrix[r][c], L2toL3, _d2Hsoil2->matrix[r][c]);
      }
      if(ctrl.sw_18O){
	_d18Osoil3->matrix[r][c] = InOutMix(theta3_old*d3, _d18Osoil3->matrix[r][c],
					    L2toL3, _d18Osoil2->matrix[r][c], Leak, mixmod);
	_d18ORecharge->matrix[r][c] = step == 0 ? _d18Osoil2->matrix[r][c] :
	  InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c]*ctrl.dt - L2toL3),
		   _d18ORecharge->matrix[r][c], L2toL3, _d18Osoil2->matrix[r][c]);
      }

      if(ctrl.sw_Cl){
	_cClsoil3->matrix[r][c] = InOutMix(theta3_old*d3, _cClsoil3->matrix[r][c],
					   L2toL3, _cClsoil2->matrix[r][c], Leak, mixmod);
	_cClRecharge->matrix[r][c] = step == 0 ? _cClsoil2->matrix[r][c] :
	  InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c]*ctrl.dt - L2toL3),
		   _cClRecharge->matrix[r][c], L2toL3, _cClsoil2->matrix[r][c]);
      }

      if(ctrl.sw_Age){
	_Agesoil3->matrix[r][c] = InOutMix(theta3_old*d3, _Agesoil3->matrix[r][c],
					   L2toL3, _Agesoil2->matrix[r][c], Leak, mixmod);
	_AgeRecharge->matrix[r][c] = step == 0 ? _Agesoil2->matrix[r][c] :
	  InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c]*ctrl.dt - L2toL3),
		   _AgeRecharge->matrix[r][c], L2toL3, _Agesoil2->matrix[r][c]);
      }
    }

    // Groundwater ------------------------------------------------------------------------
    /*
      if(ctrl.sw_2H){
      d_old = _d2Hgroundwater->matrix[r][c] ;
      _d2Hgroundwater->matrix[r][c] = _d2Hsoil3->matrix[r][c];

      if(abs(_d2Hgroundwater->matrix[r][c])>300)
      cout << r << " " << c << "d2H : | L3_old:" << theta3_old*d3 << //"| L2toGW:" << L2toGW <<
      "| d_GWnew:" << _d2Hgroundwater->matrix[r][c] << "| dGW_old:" << d_old <<
      "| d_soil2:" << _d2Hsoil2->matrix[r][c] << endl;
      }

      if(ctrl.sw_18O){
      d_old = _d18Ogroundwater->matrix[r][c] ;
      _d18Ogroundwater->matrix[r][c] = _d18Osoil3->matrix[r][c];

      if(abs(_d18Ogroundwater->matrix[r][c])>100)
      cout << r << " " << c << "d18O : | L3_old:" << theta3_old*d3 << //"| L2toGW:" << L2toGW <<
      "| d_GWnew:" << _d18Ogroundwater->matrix[r][c] << "| dGW_old:" << d_old <<
      "| d_soil2:" << _d18Osoil2->matrix[r][c] << endl;
      }

      if(ctrl.sw_Cl)
      _cClgroundwater->matrix[r][c] = _cClsoil3->matrix[r][c];

      if(ctrl.sw_Age)
      _Agegroundwater->matrix[r][c] = _Agesoil3->matrix[r][c];
    */
  }

  // -- Leakage : only if first infiltration (in SolveSurfaceFluxes)
  if(step == 0){
    /*
      if(ctrl.sw_2H)
      //_d2Hleakage->matrix[r][c] = _d2Hgroundwater->matrix[r][c];
      _d2Hleakage->matrix[r][c] = _d2Hsoil3->matrix[r][c];
      if(ctrl.sw_18O)
      _d18Oleakage->matrix[r][c] = _d18Ogroundwater->matrix[r][c];
      if(ctrl.sw_Cl)
      _cClleakage->matrix[r][c] = _cClgroundwater->matrix[r][c];
      if(ctrl.sw_Age)
      _Ageleakage->matrix[r][c] = _Agegroundwater->matrix[r][c];
    */
    if(ctrl.sw_2H)
      _d2Hleakage->matrix[r][c] = _d2Hsoil3->matrix[r][c];
    if(ctrl.sw_18O)
      _d18Oleakage->matrix[r][c] = _d18Osoil3->matrix[r][c];
    if(ctrl.sw_Cl)
      _cClleakage->matrix[r][c] = _cClsoil3->matrix[r][c];
    if(ctrl.sw_Age)
      _Ageleakage->matrix[r][c] = _Agesoil3->matrix[r][c];
  }
}
