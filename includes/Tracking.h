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
 * Tracking.h
 *
 *  Created on: Nov 14, 2016
 *      Author: Sylvain Kuppel
 *
 */

#ifndef TRACKING_H_
#define TRACKING_H_

#include "Grid.h"
#include "SortGrid.h"
#include "Basin.h"
#include "Atmosphere.h"
#include "InitConf.h"

#include <omp.h>
#include <errno.h>

using namespace std;

class Basin;
class Tracking {

  //Spatial distribution of tracked signatures (d2H, d18O, age) in the compartments
  grid *_d2Hprecip, *_d18Oprecip, *_cClprecip;
  // Canopy (avg over veg fractions)
  grid *_d2Hcanopy_sum, *_d18Ocanopy_sum, *_cClcanopy_sum, *_Agecanopy_sum;
  // Snowpack
  grid *_d2Hsnowpack, *_d18Osnowpack, *_cClsnowpack, *_Agesnowpack;
  // Snowmelt
  grid *_d2Hsnowmelt, *_d18Osnowmelt, *_cClsnowmelt, *_Agesnowmelt;
  // Ponding
  grid *_d2Hsurface, *_d18Osurface, *_cClsurface, *_Agesurface;
  // Vadose layer 1
  grid *_d2Hsoil1, *_d18Osoil1, *_cClsoil1, *_Agesoil1;
  // Vadose layer 2
  grid *_d2Hsoil2, *_d18Osoil2, *_cClsoil2, *_Agesoil2;
  // Weighted average L1+L2
  grid *_d2Hsoil_12, *_d18Osoil_12, *_cClsoil_12, *_Agesoil_12;
  // Vadose layer 3
  grid *_d2Hsoil3, *_d18Osoil3, *_cClsoil3, *_Agesoil3;
  // Vadose weighted average
  grid *_d2HsoilAv, *_d18OsoilAv, *_cClsoilAv, *_AgesoilAv;
  // Mobile water water 1
  grid *_d2H_MW1, *_d18O_MW1, *_cCl_MW1, *_Age_MW1;
  // Mobile water water 2
  grid *_d2H_MW2, *_d18O_MW2, *_cCl_MW2, *_Age_MW2;
  // Tightly-bound water 1
  grid *_d2H_TB1, *_d18O_TB1, *_cCl_TB1, *_Age_TB1;
  // Tightly-bound water 2
  grid *_d2H_TB2, *_d18O_TB2, *_cCl_TB2, *_Age_TB2; 
  // Tightly-bound water 1+2
  grid *_Age_MW12, *_Age_TB12;
  // Groundwater
  grid *_d2Hgroundwater, *_d18Ogroundwater, *_cClgroundwater, *_Agegroundwater; 
  // Signature of evaporating water (not for Cl)
  grid *_d2HevapS_sum, *_d18OevapS_sum, *_AgeevapS_sum;
  grid *_d2HevapI_sum, *_d18OevapI_sum, *_AgeevapI_sum;
  grid *_d2HevapT_sum, *_d18OevapT_sum, *_AgeevapT_sum;
  // Signature of outgoing water
  grid *_Fd2HLattoSrf, *_Fd18OLattoSrf, *_FcClLattoSrf, *_FAgeLattoSrf;
  grid *_Fd2HLattoChn, *_Fd18OLattoChn, *_FcClLattoChn, *_FAgeLattoChn;
  grid *_Fd2HLattoGW, *_Fd18OLattoGW, *_FcClLattoGW, *_FAgeLattoGW;
  grid *_d2Hleakage, *_d18Oleakage, *_cClleakage, *_Ageleakage;
  // Internal age contributions
  grid *_d2HGWtoChn, *_d2HSrftoChn, *_d2HRecharge;
  grid *_d18OGWtoChn, *_d18OSrftoChn, *_d18ORecharge;
  grid *_cClGWtoChn, *_cClSrftoChn, *_cClRecharge;
  grid *_AgeGWtoChn, *_AgeSrftoChn, *_AgeRecharge;

  //vectors containing signature of water output for each cell with no drainage (ldd value of 5). 
  // The vectCell structure contains the row and col  
  //of that cell with no output and the area draining to that cell
  // surface output
  vectCells _d2HOvlndOutput, _d18OOvlndOutput, _cClOvlndOutput, _AgeOvlndOutput;
  // groundwater output
  vectCells _d2HGwtrOutput, _d18OGwtrOutput, _cClGwtrOutput, _AgeGwtrOutput; 

  // Control volume (area) from transit/residence times tracking: only where _Age_domain > 0
  // are water ages incremented. By default, it is the whole simulation domain
  grid *_AgeDomain ;
  
  //check maps mainly to make sure no nodata values are in the domain.
  void CheckMapsTrck(Control &ctrl, Basin &bsn);

 public:
  
  //Constructors
  Tracking();
  Tracking(Control &ctrl, Basin &bsn);
  
  //Destructor
  ~Tracking();
  
  void MixingV_down(Basin &bsn, Control &ctrl, 
		    double &d1, double &d2, double &d3, double &fc,
		    int r, int c, int step);
  
  void MixingV_latup(Basin &bsn, Control &ctrl, 
		     double &d1, double &d2, double &d3, double &fc,
		     double &Qk1, double &dtdx, double &dx, int r, int c);
  
  void MixingV_evapS(Atmosphere &atm, Basin &bsn, Control &ctrl, 
		     double &d1, double &theta_new,
		     double Ts, double &etp, double &beta,
		     double &d2Hevap, double &d18Oevap, double &Agevap,
		     int r, int c);
  
  void MixingV_snow(Atmosphere &atm, Basin &bsn, Control &ctrl, double &h, double &dh, int r, int c);
  
  void MixingV_through(Atmosphere &atm, Basin &bsn, Control &ctrl, double &rain, double &p, int r, int c);
  
  void FCdownstream(Basin &bsn, Control &ctrl, 
		    double &Qk1, double &dtdx, double &dx, int r, int c, int rr, int cc);

  void MixingTPD_postET(Basin &bsn, Control &ctrl,
			double &dtheta, double &dtheta2,
			double &kTB_L1, double &kTB_L2,
			double &kMW_L1, double &kMW_L2,
			int r, int c);
  
  // Outlet lateral fluxes' signatures
  void OutletVals(Control &ctrl, int mode, int r, int c);

  // Generic function for both isotopes, using the 'iso' toggle: 0=deuterium, 1=oxygen18
  int Frac_Esoil(Atmosphere &atm, Basin &bsn, Control &ctrl,
		 REAL8 V_old, REAL8 V_new, REAL8 &theta,
		 REAL8 &di_old, REAL8 &di_new, REAL8 &di_evap,
		 REAL8 &Ts, int r, int c, int iso);
  
  // Age increment at end of time step
  int IncrementAge(Basin &bsn, Control &ctrl);
  
  // Soil averaged quantities
  int CalcTrcksoil_12(Basin &bsn, int mode);
  int CalcTrcksoil_Av(Basin &bsn, int mode);
  int CalcTrcksoil_GW(Basin &bsn, int mode);
 
  // Convert grid from isotopic ratios to isotopic deltas
  void Ratio2DeltaGrid(const Basin &bsn, const grid &m, grid &mO, int iso);
  void Ratio2DeltaGrid_L1L2(const Basin &bsn, const grid &mL1, const grid &mL2, grid &mO, int iso);
  void Ratio2DeltaGrid_SoilAv(const Basin &bsn, const grid &mL1, const grid &mL2,
			      const grid &mL3, const grid &mGW, grid &mO, int iso);
  //int ReadConfigTrck(Control &ctrl, string confilename = "configTrck.ini");
  
  // Conversion from soil to two-pore and vice-versa
  int CalcInitTPD(Basin &bsn, Control &ctrl);
  int CalcTPDtoLayers(Basin &bsn, Control &ctrl);

  //Getters 2H (ratios)
  //grid *getd2Hprecip() const {
  //  return _d2Hprecip;
  //}
  grid *getd2Hcanopy_sum() const {
    return _d2Hcanopy_sum;
  }
  grid *getd2Hsnowpack() const {
    return _d2Hsnowpack;
  }
  grid *getd2Hsnowmelt() const {
    return _d2Hsnowmelt;
  }
  grid *getd2Hsurface() const {
    return _d2Hsurface;
  }
  grid *getd2Hsoil1() const {
    return _d2Hsoil1;
  }
  grid *getd2Hsoil2() const {
    return _d2Hsoil2;
  }
  grid *getd2Hsoil_12() const {
    return _d2Hsoil_12;
  }
  grid *getd2Hsoil3() const {
    return _d2Hsoil3;
  }
  grid *getd2Hsoil_Av() {
    return _d2HsoilAv;
  }
  grid *getd2Hgroundwater() const {
    return _d2Hgroundwater;
  }
  grid *getd2HevapS_sum() const {
    return _d2HevapS_sum;
  }
  grid *getd2HevapI_sum() const {
    return _d2HevapI_sum;
  }
  grid *getd2HevapT_sum() const {
    return _d2HevapT_sum;
  }
  grid *getd2Hleakage() const {
    return _d2Hleakage;
  }
  grid *getd2HGWtoChn() const {
    return _d2HGWtoChn;
  }
  grid *getd2HSrftoChn() const {
    return _d2HSrftoChn;
  }
  grid *getd2HRecharge() const {
    return _d2HRecharge;
  }
  const vectCells *getd2HOvlndOutput() const {
    return &_d2HOvlndOutput;
  }  
  const vectCells *getd2HGwtrOutput() const {
    return &_d2HGwtrOutput;
  }
  grid *getd2H_MW1() const {
    return _d2H_MW1;
  }
  grid *getd2H_MW2() const {
    return _d2H_MW2;
  }
  grid *getd2H_TB1() const {
    return _d2H_TB1;
  }
  grid *getd2H_TB2() const {
    return _d2H_TB2;
  }
  // -------
  
  
  // 18O
  //grid *getd18Oprecip() const {
  //  return _d18Oprecip;
  //}
  grid *getd18Ocanopy_sum() const {
    return _d18Ocanopy_sum;
  }
  grid *getd18Osnowpack() const {
    return _d18Osnowpack;
  }
  grid *getd18Osnowmelt() const {
    return _d18Osnowmelt;
  }
  grid *getd18Osurface() const {
    return _d18Osurface;
  }
  grid *getd18Osoil1() const {
    return _d18Osoil1;
  }
  grid *getd18Osoil2() const {
    return _d18Osoil2;
  }
  grid *getd18Osoil_12() const {
    return _d18Osoil_12;
  }
  grid *getd18Osoil3() const {
    return _d18Osoil3;
  }
  grid *getd18Osoil_Av() const {
    return _d18OsoilAv;
  }
  grid *getd18Ogroundwater() const {
    return _d18Ogroundwater;
  }
  grid *getd18OevapS_sum() const {
    return _d18OevapS_sum;
  }
  grid *getd18OevapI_sum() const {
    return _d18OevapI_sum;
  }
  grid *getd18OevapT_sum() const {
    return _d18OevapT_sum;
  }
  grid *getd18Oleakage() const {
    return _d18Oleakage;
  }
  grid *getd18OGWtoChn() const {
    return _d18OGWtoChn;
  }
  grid *getd18OSrftoChn() const {
    return _d18OSrftoChn;
  }
  grid *getd18ORecharge() const {
    return _d18ORecharge;
  }
  const vectCells *getd18OOvlndOutput() const {
    return &_d18OOvlndOutput;
  }  
  const vectCells *getd18OGwtrOutput() const {
    return &_d18OGwtrOutput;
  }
  grid *getd18O_MW1() const {
    return _d18O_MW1;
  }
  grid *getd18O_MW2() const {
    return _d18O_MW2;
  }
  grid *getd18O_TB1() const {
    return _d18O_TB1;
  }
  grid *getd18O_TB2() const {
    return _d18O_TB2;
  }
  // -------

  // Chloride
  //grid *getcClprecip() const {
  //  return _cClprecip;
  //}
  grid *getcClcanopy_sum() const {
    return _cClcanopy_sum;
  }
  grid *getcClsnowpack() const {
    return _cClsnowpack;
  }
  grid *getcClsnowmelt() const {
    return _cClsnowmelt;
  }
  grid *getcClsurface() const {
    return _cClsurface;
  }
  grid *getcClsoil1() const {
    return _cClsoil1;
  }
  grid *getcClsoil2() const {
    return _cClsoil2;
  }
  grid *getcClsoil_12() const {
    return _cClsoil_12;
  }
  grid *getcClsoil3() const {
    return _cClsoil3;
  }
  grid *getcClsoil_Av() const {
    return _cClsoilAv;
  }
  grid *getcClgroundwater() const {
    return _cClgroundwater;
  }
  grid *getcClleakage() const {
    return _cClleakage;
  }
  grid *getcClGWtoChn() const {
    return _cClGWtoChn;
  }
  grid *getcClSrftoChn() const {
    return _cClSrftoChn;
  }
  grid *getcClRecharge() const {
    return _cClRecharge;
  }
  const vectCells *getcClOvlndOutput() const {
    return &_cClOvlndOutput;
  }  
  const vectCells *getcClGwtrOutput() const {
    return &_cClGwtrOutput;
  }
  grid *getcCl_MW1() const {
    return _cCl_MW1;
  }
  grid *getcCl_MW2() const {
    return _cCl_MW2;
  }
  grid *getcCl_TB1() const {
    return _cCl_TB1;
  }
  grid *getcCl_TB2() const {
    return _cCl_TB2;
  }
  // -------

  // Age
  grid *getAgecanopy_sum() const {
    return _Agecanopy_sum;
  }
  grid *getAgesnowpack() const {
    return _Agesnowpack;
  }
  grid *getAgesnowmelt() const {
    return _Agesnowmelt;
  }
  grid *getAgesurface() const {
    return _Agesurface;
  }
  grid *getAgesoil1() const {
    return _Agesoil1;
  }
  grid *getAgesoil2() const {
    return _Agesoil2;
  }
  grid *getAgesoil_12() const {
    return _Agesoil_12;
  }
  grid *getAgesoil3() const {
    return _Agesoil3;
  }
  grid *getAgesoil_Av() const {
    return _AgesoilAv;
  }
  grid *getAgegroundwater() const {
    return _Agegroundwater;
  }
  grid *getAgeevapS_sum() const {
    return _AgeevapS_sum;
  }
  grid *getAgeevapI_sum() const {
    return _AgeevapI_sum;
  }
  grid *getAgeevapT_sum() const {
    return _AgeevapT_sum;
  }
  grid *getAgeleakage() const {
    return _Ageleakage;
  }
  grid *getAgeGWtoChn() const {
    return _AgeGWtoChn;
  }
  grid *getAgeSrftoChn() const {
    return _AgeSrftoChn;
  }
  grid *getAgeRecharge() const {
    return _AgeRecharge;
  }

  // Age domain
  grid *getAgeDomain() const {
    return _AgeDomain;
  }

  const vectCells *getAgeOvlndOutput() const {
    return &_AgeOvlndOutput;
  }  
  const vectCells *getAgeGwtrOutput() const {
    return &_AgeGwtrOutput;
  }
  grid *getAge_MW1() const {
    return _Age_MW1;
  }
  grid *getAge_MW2() const {
    return _Age_MW2;
  }
  grid *getAge_MW12() const {
    return _Age_MW12;
  }
  grid *getAge_TB1() const {
    return _Age_TB1;
  }
  grid *getAge_TB2() const {
    return _Age_TB2;
  }
  grid *getAge_TB12() const {
    return _Age_TB12;
  }
  // ---
  
  // --- Setters
  // - 2H
  void setd2Hcanopy_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hcanopy_sum->matrix[row][col] = value;
  }
  void setd2Hsnowpack(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hsnowpack->matrix[row][col] = value;
  }
  void setd2Hsurface(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hsurface->matrix[row][col] = value;
  }
  void setd2Hsoil1(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hsoil1->matrix[row][col] = value;
  }
  void setd2Hsoil2(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hsoil2->matrix[row][col] = value;
  }
  void setd2Hsoil3(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hsoil3->matrix[row][col] = value;
  }
  void setd2Hgroundwater(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hgroundwater->matrix[row][col] = value;
  }
  void setd2HevapS_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d2HevapS_sum->matrix[row][col] = value;
  }
  void setd2HevapI_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d2HevapI_sum->matrix[row][col] = value;
  }
  void setd2HevapT_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d2HevapT_sum->matrix[row][col] = value;
  }
  void resetFd2HLat() { // reset summed lateral contributions
    _Fd2HLattoGW->reset();
    _Fd2HLattoChn->reset();
    _Fd2HLattoSrf->reset();
  }

  // - 18O
  void setd18Ocanopy_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d18Ocanopy_sum->matrix[row][col] = value;
  }  
  void setd18Osnowpack(UINT4 row, UINT4 col, REAL8 value) {
    _d18Osnowpack->matrix[row][col] = value;
  }
  void setd18Osurface(UINT4 row, UINT4 col, REAL8 value) {
    _d18Osurface->matrix[row][col] = value;
  }
  void setd18Osoil1(UINT4 row, UINT4 col, REAL8 value) {
    _d18Osoil1->matrix[row][col] = value;
  }
  void setd18Osoil2(UINT4 row, UINT4 col, REAL8 value) {
    _d18Osoil2->matrix[row][col] = value;
  }
  void setd18Osoil3(UINT4 row, UINT4 col, REAL8 value) {
    _d18Osoil3->matrix[row][col] = value;
  }
  void setd18Ogroundwater(UINT4 row, UINT4 col, REAL8 value) {
    _d18Ogroundwater->matrix[row][col] = value;
  }
  void setd18OevapS_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d18OevapS_sum->matrix[row][col] = value;
  }
  void setd18OevapI_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d18OevapI_sum->matrix[row][col] = value;
  }
  void setd18OevapT_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d18OevapT_sum->matrix[row][col] = value;
  }
  void resetFd18OLat() { // reset summed lateral contributions
    _Fd18OLattoGW->reset();
    _Fd18OLattoChn->reset();
    _Fd18OLattoSrf->reset();
  }

  // - Chloride
  void setcClcanopy_sum(UINT4 row, UINT4 col, REAL8 value) {
    _cClcanopy_sum->matrix[row][col] = value;
  }  
  void setcClsnowpack(UINT4 row, UINT4 col, REAL8 value) {
    _cClsnowpack->matrix[row][col] = value;
  }
  void setcClsurface(UINT4 row, UINT4 col, REAL8 value) {
    _cClsurface->matrix[row][col] = value;
  }
  void setcClsoil1(UINT4 row, UINT4 col, REAL8 value) {
    _cClsoil1->matrix[row][col] = value;
  }
  void setcClsoil2(UINT4 row, UINT4 col, REAL8 value) {
    _cClsoil2->matrix[row][col] = value;
  }
  void setcClsoil3(UINT4 row, UINT4 col, REAL8 value) {
    _cClsoil3->matrix[row][col] = value;
  }
  void setcClgroundwater(UINT4 row, UINT4 col, REAL8 value) {
    _cClgroundwater->matrix[row][col] = value;
  }
  void resetFcClLat() { // reset summed lateral contributions
    _FcClLattoGW->reset();
    _FcClLattoChn->reset();
    _FcClLattoSrf->reset();
  }

  // - Age
  void setAgecanopy_sum(UINT4 row, UINT4 col, REAL8 value) {
    _Agecanopy_sum->matrix[row][col] = value;
  }    
  void setAgesnowpack(UINT4 row, UINT4 col, REAL8 value) {
    _Agesnowpack->matrix[row][col] = value;
  }
  void setAgesurface(UINT4 row, UINT4 col, REAL8 value) {
    _Agesurface->matrix[row][col] = value;
  }
  void setAgesoil1(UINT4 row, UINT4 col, REAL8 value) {
    _Agesoil1->matrix[row][col] = value;
  }
  void setAgesoil2(UINT4 row, UINT4 col, REAL8 value) {
    _Agesoil2->matrix[row][col] = value;
  }
  void setAgesoil3(UINT4 row, UINT4 col, REAL8 value) {
    _Agesoil3->matrix[row][col] = value;
  }
  void setAgegroundwater(UINT4 row, UINT4 col, REAL8 value) {
    _Agegroundwater->matrix[row][col] = value;
  }
  void setAgeevapS_sum(UINT4 row, UINT4 col, REAL8 value) {
    _AgeevapS_sum->matrix[row][col] = value;
  }
  void setAgeevapI_sum(UINT4 row, UINT4 col, REAL8 value) {
    _AgeevapI_sum->matrix[row][col] = value;
  }
  void setAgeevapT_sum(UINT4 row, UINT4 col, REAL8 value) {
    _AgeevapT_sum->matrix[row][col] = value;
  }
  void resetFAgeLat() { // reset summed lateral contributions
    _FAgeLattoGW->reset();
    _FAgeLattoChn->reset();
    _FAgeLattoSrf->reset();
  }

  // ---- Mixing equations ------------------------------------------------------------------
  // When there's only input
  double InputMix(double hold, double iold, double qin, double iin){

    double inew = 0;
    // Assign fluxes
    if(hold + qin > RNDOFFERR)
      inew = (iold*hold + iin*qin)/ (hold + qin);
    else
      inew = iold;
    return inew;
  }
  
  // Mixing when one input one output: two modes to apporximate what aggregated h and i
  // values the DE approximation uses.
  // If mode = 0, h=hold and i=inew. Mixing without change of volume due to output.
  // If mode = 1, h=0.5*(hold+hnew) and i=0.5*(inew+iold). A somewhat linear model of mixing.
  // A more complex, and stable weighted-linear model of mixing.
  double InOutMix(double hold, double iold, double qin, double iin, double qout, int mode){

    double inew = 0;
    double hsum = 0;

    if(mode==0)
      inew = hold + qin > RNDOFFERR ?
	(iold*hold + iin*qin)/ (hold + qin) : iold;
    /* else if (mode==1) */
    /*   inew = hold + qin - 0.5*qout > RNDOFFERR ?  */
    /* 	((iold+1000)*(hold-0.5*qout) + (iin+1000)*qin)/ (hold - 0.5*qout + qin) -1000 : iold; */
    /* else if (mode==2){ */
    /*   hsum = 2*hold+qin-qout; */
    /*   inew = powl(hsum,2)/2 + hold*qin > RNDOFFERR ?  */
    /* 	((iold+1000)*(powl(hsum,2)/2-(hsum-hold)*qin) + (iin+1000)*qin*hsum) / */
    /* 	(powl(hsum,2)/2 + hold*qin) -1000 : iold; */
    /*   //(iold*(powl(hsum,2)/2-(hsum-hold)*qin) + iin*qin*hsum) / */
    /*   //(powl(hsum,2)/2 + hold*qin) : iold; */
    //}
    else if(mode ==1) {
      hsum = 0.5*(hold+qin+std::max<double>(0,hold-qout));
      inew = 0.5*hsum+qin > RNDOFFERR ? 
	//((iold+1000)*(hsum-0.5*qin) + (iin+1000)*qin)/ (hsum + 0.5*qin) -1000 : iold;
      	(iold*(hsum-0.5*qin) + iin*qin)/ (hsum + 0.5*qin) : iold;
    }

    return inew;
  }
  // ------------------------------------------------------------------------------------------
};

#endif /* TRACKING_H_ */
