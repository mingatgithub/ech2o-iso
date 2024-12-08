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
 * Budget.h
 *
 *  Created on: Mar 8, 2010
 *      Author: Marco Maneta
 */

#ifndef BUDGET_H_
#define BUDGET_H_

#include "Basin.h"
#include "Atmosphere.h"
#include "InitConf.h"

struct Budget {


  double dt;
  double MBErr;

  //init storages
  double initsnowpack; //m3
  double initponding; //m3
  //double initvadose; //m3
  double initL1, initL2, initL3; //m3
  //double initgravwater; //m3
  //double initgrndwater; //m3
  double initGW;
  //storages
  REAL8 canopy; //m3
  double snowpack; //m3
  double ponding; //m3
  double vadose; //m3
  double soilL1, soilL2, soilL3; //m3
  double rootzone; //m3
  //double gravwater; //m3
  double grndwater; //m3

  //water inputs
  double precipitation;
  double snowfall;


  //water outputs
  double evaporation; //m3s-1
  double evaporationS; //m3s-1
  double evaporationI; //m3s-1
  double transpiration; //m3s-1
  double leakage; //m3s-1
  double ovlndflow; //m3s-1
  double gwtrflow; //m3s-1

  // internal fluxes
  double gwtochn; //m3
  double srftochn; //m3
  double recharge; //m3

  //saturation extent
  double satarea; //%

  double AccountFluxes(const grid *map, const Basin *b);
  double AccountFluxes(const grid *map, const Atmosphere *b);
  double AccountFluxes(const vectCells *timeseries, const Basin *b);
  double AccountStorages(const grid *map, const Basin *b);
  double AccountStorages(const grid *map1, const grid *map2, const Basin *b);
  double AccountRelArea(const grid *map, const Basin *b);

  void TotalPrecipitation(const grid* map, const Atmosphere *b);
  void TotalEvaporation(const grid* map, const Basin *b);
  void TotalEvaporationS(const grid* map, const Basin *b);
  void TotalEvaporationI(const grid* map, const Basin *b);
  void TotalTranspiration(const grid* map, const Basin *b);
  void TotalBedrockLeakage(const grid* map, const Basin*b);
  void TotalOvlndFlow(const vectCells *timeseries, const Basin *b);
  void TotalGrndFlow(const vectCells *timeseries, const Basin*b);
  void TotalStorage(	const grid *Canopy,
			const grid *Snow,
			const grid *Ponding,
			//const grid *Vadose,
			const grid *SoilL1,
			const grid *SoilL2,
			const grid *SoilL3,
			//const grid *GravWater,
			const grid *GrndWater,
			const grid *ProotL1, const grid *ProotL2, const grid *ProotL3,
			const Basin *b);
  void TotalSaturationArea(const grid* map, const Basin*b);
  void TotalGWtoChn(const grid* map, const Basin*b);
  void TotalSrftoChn(const grid* map, const Basin*b);
  void TotalRecharge(const grid* map, const Basin*b);

  void MassBalanceError(const Control *ctrl);

  // Tracking --------------------------------------------------------------------------
  double MBErr_d2H, MBErr_d18O, MBErr_cCl, MBErr_Age;

  double AccountTrckFluxes(const grid *map1, const grid *map2, const Basin *b);
  double AccountTrckFluxes2(const grid *map1, const grid *map2, const Basin *b);
  double AccountTrckFluxes(const grid *map1, const grid *map2, const Atmosphere *b);
  double AccountTrckFluxes(const grid *map, const Atmosphere *b);//, const Control *ctrl);
  double AccountTrckFluxes2(const grid *map1, const grid *map2, const Atmosphere *b);
  double AccountTrckFluxes(const vectCells *timeseries1, const vectCells *timeseries2);
  double AccountTrckFluxes2(const vectCells *timeseries1, const vectCells *timeseries2);
  double AccountTrckStorages(const grid *map1, const grid *map2, const Basin *b);
  double AccountTrckStorages2(const grid *map1, const grid *map2, const Basin *b);
  double AccountTrckVadose(const grid *mapL1, const grid *mapCL1,
			   const grid *mapL2, const grid *mapCL2,
			   const grid *mapL3, const grid *mapCL3,
			   //const grid *mapGW, const grid *mapCGW,
			   const Basin *b);
  double AccountTrckRootZone(const grid *mapL1, const grid *mapCL1, const grid *mappL1,
			     const grid *mapL2, const grid *mapCL2, const grid *mappL2,
			     const grid *mapL3, const grid *mapCL3, const grid *mappL3,
			     //const grid *mapGW, const grid *mapCGW,
			     const Basin *b);
  double AccountTrckDomain(const grid *mapCanopy, const grid *mapCCanopy,
			   const grid *mapSnow, const grid *mapCSnow,
			   const grid *mapSurface, const grid *mapCSurface,
			   const grid *mapL1, const grid *mapCL1,
			   const grid *mapL2, const grid *mapCL2,
			   const grid *mapL3, const grid *mapCL3,
			   //const grid *mapGW, const grid *mapCGW,
			   const Basin *b);
  double AccountTrckET(const grid* evapS, const grid* CevapS,
		       const grid* evapI, const grid* CevapI,
		       const grid* evapT, const grid* CevapT,
		       const Basin *b);
  double AccountTrckOut(const grid* evapS, const grid* CevapS,
			const grid* evapI, const grid* CevapI,
			const grid* evapT, const grid* CevapT,
			const grid* leakage, const grid* Cleakage,
			const vectCells *OvlndOut, const vectCells *COvlndOut,
			const vectCells *GWOut, const vectCells *CGWOut,
			const Basin *b);
  double AccountTrckOut(const grid* leakage, const grid* Cleakage,
			const vectCells *OvlndOut, const vectCells *COvlndOut,
			const vectCells *GWOut, const vectCells *CGWOut,
			const Basin *b);

  void TotalPrecipitation_d2H(const grid* map1, const grid* map2, const Atmosphere *b);
  void TotalEvaporationS_d2H(const grid* map1, const grid* map2, const Basin *b);
  void TotalEvaporationI_d2H(const grid* map1, const grid* map2, const Basin *b);
  void TotalTranspiration_d2H(const grid* map1, const grid* map2, const Basin *b);
  void TotalBedrockLeakage_d2H(const grid* map1, const grid* map2, const Basin*b);
  void TotalOvlndFlow_d2H(const vectCells *timeseries1, const vectCells *timeseries2);
  void TotalGrndFlow_d2H(const vectCells *timeseries1, const vectCells *timeseries2);
  void TotalStorage_d2H( const grid *Canopy, const grid *Canopy_d2H,
			const grid *Snow, const grid *Snow_d2H,
			const grid *Ponding, const grid *Ponding_d2H,
			const grid *SoilL1_Water, const grid *SoilL1_d2H,
			 const grid *ProotzoneL1,
			 const grid *SoilL2, const grid *SoilL2_d2H,
			 const grid *ProotzoneL2,
			 const grid *SoilL3, const grid *SoilL3_d2H,
			 const grid *ProotzoneL3,
			 const grid *GWater, const grid *GWater_d2H,
			const Basin *b);//, const Control *ctrl);
  void InstPrecipitation_d2H(const grid* precip, const grid* Cprecip,
			     const Basin *b);
  void InstPrecipitation_d2H(const grid* precip, const grid* Cprecip,
         		const Atmosphere *b);
  void InstEvaporation_d2H(const grid* evapS, const grid* CevapS,
			    const grid* evapI, const grid* CevapI,
			    const grid* evapT, const grid* CevapT,
			    const Basin *b);
  void InstEvaporationS_d2H(const grid* map1, const grid* map2, const Basin *b);
  void InstEvaporationI_d2H(const grid* map1, const grid* map2, const Basin *b);
  void InstTranspiration_d2H(const grid* map1, const grid* map2, const Basin *b);
  void InstBedrockLeakage_d2H(const grid* map1, const grid* map2, const Basin*b);
  void InstOvlndFlow_d2H(const vectCells *timeseries1, const vectCells *timeseries2);
  void InstGrndFlow_d2H(const vectCells *timeseries1, const vectCells *timeseries2);
  void InstOut_d2H(const grid* evapS, const grid* CevapS,
			   const grid* evapI, const grid* CevapI,
			   const grid* evapT, const grid* CevapT,
			   const grid* leakage, const grid* Cleakage,
			   const vectCells *OvlndOut, const vectCells *COvlndOut,
			   const vectCells *GWOut, const vectCells *CGWOut,
			   const Basin *b);
  void InstSrftoChn_d2H(const grid* map1, const grid* map2, const Basin*b);
  void InstGWtoChn_d2H(const grid* map1, const grid* map2, const Basin*b);
  void InstRecharge_d2H(const grid* map1, const grid* map2, const Basin*b);

  void TotalPrecipitation_d18O(const grid* map1, const grid* map2, const Atmosphere *b);
  void TotalEvaporationS_d18O(const grid* map1, const grid* map2, const Basin *b);
  void TotalEvaporationI_d18O(const grid* map1, const grid* map2, const Basin *b);
  void TotalTranspiration_d18O(const grid* map1, const grid* map2, const Basin *b);
  void TotalBedrockLeakage_d18O(const grid* map1, const grid* map2, const Basin*b);
  void TotalOvlndFlow_d18O(const vectCells *timeseries1, const vectCells *timeseries2);
  void TotalGrndFlow_d18O(const vectCells *timeseries1, const vectCells *timeseries2);
  void TotalStorage_d18O( const grid *Canopy, const grid *Canopy_d18O,
			  const grid *Snow, const grid *Snow_d18O,
			  const grid *Ponding, const grid *Ponding_d18O,
			  const grid *SoilL1_Water, const grid *SoilL1_d18O,
			  const grid *ProotzoneL1,
			  const grid *SoilL2, const grid *SoilL2_d18O,
			  const grid *ProotzoneL2,
			  const grid *SoilL3, const grid *SoilL3_d18O,
			  const grid *ProotzoneL3,
 			  const grid *GWater, const grid *GWater_d18O,
			  const Basin *b);//, const Control *ctrl);
  void InstPrecipitation_d18O(const grid* precip, const grid* Cprecip,
			      const Basin *b);
  void InstPrecipitation_d18O(const grid* precip, const grid* Cprecip,
         		const Atmosphere *b);
  void InstEvaporation_d18O(const grid* evapS, const grid* CevapS,
			    const grid* evapI, const grid* CevapI,
			    const grid* evapT, const grid* CevapT,
			    const Basin *b);
  void InstEvaporationS_d18O(const grid* map1, const grid* map2, const Basin *b);
  void InstEvaporationI_d18O(const grid* map1, const grid* map2, const Basin *b);
  void InstTranspiration_d18O(const grid* map1, const grid* map2, const Basin *b);
  void InstBedrockLeakage_d18O(const grid* map1, const grid* map2, const Basin*b);
  void InstOvlndFlow_d18O(const vectCells *timeseries1, const vectCells *timeseries2);
  void InstGrndFlow_d18O(const vectCells *timeseries1, const vectCells *timeseries2);
  void InstOut_d18O(const grid* evapS, const grid* CevapS,
			   const grid* evapI, const grid* CevapI,
			   const grid* evapT, const grid* CevapT,
			   const grid* leakage, const grid* Cleakage,
			   const vectCells *OvlndOut, const vectCells *COvlndOut,
			   const vectCells *GWOut, const vectCells *CGWOut,
			   const Basin *b);
  void InstSrftoChn_d18O(const grid* map1, const grid* map2, const Basin*b);
  void InstGWtoChn_d18O(const grid* map1, const grid* map2, const Basin*b);
  void InstRecharge_d18O(const grid* map1, const grid* map2, const Basin*b);

  void TotalPrecipitation_cCl(const grid* map1, const grid* map2, const Atmosphere *b);
  void TotalBedrockLeakage_cCl(const grid* map1, const grid* map2, const Basin*b);
  void TotalOvlndFlow_cCl(const vectCells *timeseries1, const vectCells *timeseries2);
  void TotalGrndFlow_cCl(const vectCells *timeseries1, const vectCells *timeseries2);
  void TotalStorage_cCl( const grid *Canopy, const grid *Canopy_cCl,
			  const grid *Snow, const grid *Snow_cCl,
			  const grid *Ponding, const grid *Ponding_cCl,
			  const grid *SoilL1_Water, const grid *SoilL1_cCl,
			  const grid *ProotzoneL1,
			  const grid *SoilL2, const grid *SoilL2_cCl,
			  const grid *ProotzoneL2,
			  const grid *SoilL3, const grid *SoilL3_cCl,
			  const grid *ProotzoneL3,
 			  const grid *GWater, const grid *GWater_cCl,
			  const Basin *b);//, const Control *ctrl);
  void InstPrecipitation_cCl(const grid* precip, const grid* Cprecip,
			      const Basin *b);
  void InstPrecipitation_cCl(const grid* precip, const grid* Cprecip,
           const Atmosphere *b);
  void InstBedrockLeakage_cCl(const grid* map1, const grid* map2, const Basin*b);
  void InstOvlndFlow_cCl(const vectCells *timeseries1, const vectCells *timeseries2);
  void InstGrndFlow_cCl(const vectCells *timeseries1, const vectCells *timeseries2);
  void InstOut_cCl(const grid* leakage, const grid* Cleakage,
		   const vectCells *OvlndOut, const vectCells *COvlndOut,
		   const vectCells *GWOut, const vectCells *CGWOut,
		   const Basin *b);
  void InstSrftoChn_cCl(const grid* map1, const grid* map2, const Basin*b);
  void InstGWtoChn_cCl(const grid* map1, const grid* map2, const Basin*b);
  void InstRecharge_cCl(const grid* map1, const grid* map2, const Basin*b);

  void TotalPrecipitation_Age();//const grid* map, const Atmosphere *b);//, const Control *ctrl);
  void TotalEvaporationS_Age(const grid* map1, const grid* map2, const Basin *b);
  void TotalEvaporationI_Age(const grid* map1, const grid* map2, const Basin *b);
  void TotalTranspiration_Age(const grid* map1, const grid* map2, const Basin *b);
  void TotalBedrockLeakage_Age(const grid* map1, const grid* map2, const Basin*b);
  void TotalOvlndFlow_Age(const vectCells *timeseries1, const vectCells *timeseries2);
  void TotalGrndFlow_Age(const vectCells *timeseries1, const vectCells *timeseries2);
  void TotalStorage_Age( const grid *Canopy, const grid *Canopy_Age,
			 const grid *Snow, const grid *Snow_Age,
			 const grid *Ponding, const grid *Ponding_Age,
			 const grid *SoilL1_Water, const grid *SoilL1_Age,
			 const grid *ProotzoneL1,
			 const grid *SoilL2, const grid *SoilL2_Age,
			 const grid *ProotzoneL2,
			 const grid *SoilL3, const grid *SoilL3_Age,
			 const grid *ProotzoneL3,
			 const grid *GWater, const grid *GWater_Age,
			 const Basin *b);//, const Control *ctrl);
  void InstEvaporation_Age(const grid* evapS, const grid* CevapS,
			    const grid* evapI, const grid* CevapI,
			    const grid* evapT, const grid* CevapT,
			    const Basin *b);
  void InstEvaporationS_Age(const grid* map1, const grid* map2, const Basin *b);
  void InstEvaporationI_Age(const grid* map1, const grid* map2, const Basin *b);
  void InstTranspiration_Age(const grid* map1, const grid* map2, const Basin *b);
  void InstBedrockLeakage_Age(const grid* map1, const grid* map2, const Basin*b);
  void InstOvlndFlow_Age(const vectCells *timeseries1, const vectCells *timeseries2);
  void InstGrndFlow_Age(const vectCells *timeseries1, const vectCells *timeseries2);
  void InstOut_Age(const grid* evapS, const grid* CevapS,
			   const grid* evapI, const grid* CevapI,
			   const grid* evapT, const grid* CevapT,
			   const grid* leakage, const grid* Cleakage,
			   const vectCells *OvlndOut, const vectCells *COvlndOut,
			   const vectCells *GWOut, const vectCells *CGWOut,
			   const Basin *b);
  void InstSrftoChn_Age(const grid* map1, const grid* map2, const Basin*b);
  void InstGWtoChn_Age(const grid* map1, const grid* map2, const Basin*b);
  void InstRecharge_Age(const grid* map1, const grid* map2, const Basin*b);


  double initcanopy_d2H, initcanopy_d18O, initcanopy_cCl, initcanopy_Age; //m3.'tracer-unit'
  double initsnowpack_d2H, initsnowpack_d18O, initsnowpack_cCl, initsnowpack_Age; //m3.'tracer-unit'
  double initponding_d2H, initponding_d18O, initponding_cCl, initponding_Age; //m3.'tracer-unit'
  double initL1_d2H, initL1_d18O, initL1_cCl, initL1_Age; //m3.'tracer-unit'
  double initL2_d2H, initL2_d18O, initL2_cCl, initL2_Age; //m3.'tracer-unit'
  double initL3_d2H, initL3_d18O, initL3_cCl, initL3_Age; //m3.'tracer-unit'
  double initGW_d2H, initGW_d18O, initGW_cCl, initGW_Age; //m3.'tracer-unit'
  //storages
  double canopy_d2H, canopy_d18O, canopy_cCl, canopy_Age; //m3.'tracer-unit'
  double snowpack_d2H, snowpack_d18O, snowpack_cCl, snowpack_Age; //m3.'tracer-unit'
  double ponding_d2H, ponding_d18O, ponding_cCl, ponding_Age; //m3.'tracer-unit'
  double soilL1_d2H, soilL1_d18O, soilL1_cCl, soilL1_Age; //m3.'tracer-unit'
  double soilL2_d2H, soilL2_d18O, soilL2_cCl, soilL2_Age; //m3.'tracer-unit'
  double soilL3_d2H, soilL3_d18O, soilL3_cCl, soilL3_Age; //m3.'tracer-unit'
  double grndwater_d2H, grndwater_d18O, grndwater_cCl, grndwater_Age; //m3.'tracer-unit'
  //water inputs
  double precipitation_d2H, precipitation_d18O, precipitation_cCl, precipitation_Age;
  //water outputs
  double evaporationS_d2H, evaporationS_d18O, evaporationS_Age; //m3
  double evaporationI_d2H, evaporationI_d18O, evaporationI_Age; //m3
  double transpiration_d2H, transpiration_d18O, transpiration_Age; //m3
  double leakage_d2H, leakage_d18O, leakage_cCl, leakage_Age; //m3
  double ovlndflow_d2H, ovlndflow_d18O, ovlndflow_cCl, ovlndflow_Age; //m3
  double gwtrflow_d2H, gwtrflow_d18O, gwtrflow_cCl, gwtrflow_Age; //m3

  // Values specificially for Basind2HSummary.txt
  double d2Hprecip, d2HTot, d2Hvadose, d2Hrootzone;
  double d2Hcanopy, d2Hsnowpack, d2Hponding;
  double d2HsoilL1, d2HsoilL2, d2HsoilL3, d2Hgrndwater;
  // "instantaneous" output ages
  double d2HET, d2HevapS, d2HevapI, d2HevapT;
  double d2Hgwtochn, d2Hsrftochn, d2Hrecharge;
  double d2Hleakage, d2HOvlndOut, d2HGWOut, d2HOut;

  // Values specificially for Basind18OSummary.txt
  double d18Oprecip, d18OTot, d18Ovadose, d18Orootzone;
  double d18Ocanopy, d18Osnowpack, d18Oponding;
  double d18OsoilL1, d18OsoilL2, d18OsoilL3, d18Ogrndwater;
  // "instantaneous" output ages
  double d18OET, d18OevapS, d18OevapI, d18OevapT;
  double d18Ogwtochn, d18Osrftochn, d18Orecharge;
  double d18Oleakage, d18OOvlndOut, d18OGWOut, d18OOut;

  // Values specificially for BasincClSummary.txt
  double cClprecip, cClTot, cClvadose, cClrootzone;
  double cClcanopy, cClsnowpack, cClponding;
  double cClsoilL1, cClsoilL2, cClsoilL3, cClgrndwater;
  // "instantaneous" output ages
  double cClgwtochn, cClsrftochn, cClrecharge;
  double cClleakage, cClOvlndOut, cClGWOut, cClOut;

  // Values specificially for BasinAgeSummary.txt
  double AgeTot, Agevadose, Agerootzone;
  double Agecanopy, Agesnowpack, Ageponding;
  double AgesoilL1, AgesoilL2, AgesoilL3, Agegrndwater;
  // "instantaneous" output ages
  double AgeET, AgeevapS, AgeevapI, AgeevapT;
  double Agegwtochn, Agesrftochn, Agerecharge;
  double Ageleakage, AgeOvlndOut, AgeGWOut, AgeOut;


  void TrckBalanceError(const Control *ctrl);
  // -----------------------------------------------------------------------------------------

  //constructor inline
  Budget(const Basin *b, const Control *ctrl, const Tracking *trck)
  {

    dt = ctrl->dt;
    canopy = 0;
    snowpack = 0;
    ponding = 0;
    soilL1 = 0;
    soilL2 = 0;
    soilL3 = 0;
    rootzone = 0;
    //vadose = 0;
    //gravwater = 0;
    grndwater = 0;
    precipitation = 0;
    snowfall = 0;

    evaporation = 0;
    evaporationS = 0;
    evaporationI = 0;
    transpiration = 0;
    leakage = 0;
    ovlndflow = 0;
    gwtrflow = 0;

    recharge = 0;
    srftochn = 0;
    gwtochn = 0;
    satarea = 0;

    //calculate initial storages
    initsnowpack = AccountStorages(b->getSnowWaterEquiv(), b);
    initponding = AccountStorages(b->getPondingWater(), b);
    //initvadose = AccountStorages(b->getSoilWaterDepth(), b);
    initL1 = AccountStorages(b->getSoilWaterDepthL1(), b);
    initL2 = AccountStorages(b->getSoilWaterDepthL2(), b);
    initL3 = AccountStorages(b->getSoilWaterDepthL3(), b);
    initGW = AccountStorages(b->getInitGroundwater(), b);
    //initgravwater = AccountStorages(b->getGravityWater(), b);
    //initgrndwater = AccountStorages(b->getGrndWater(),b);

    // Tracking ---------------------------------------------------------------------
    if(ctrl->sw_trck && ctrl->sw_2H){
      canopy_d2H = 0;
      snowpack_d2H = 0;
      ponding_d2H = 0;
      soilL1_d2H = 0;
      soilL2_d2H = 0;
      soilL3_d2H = 0;
      grndwater_d2H = 0;
      precipitation_d2H = 0;
      evaporationS_d2H = 0;
      evaporationI_d2H = 0;
      transpiration_d2H = 0;
      leakage_d2H = 0;
      ovlndflow_d2H = 0;
      gwtrflow_d2H = 0;
      // Values for Basind2HSummary.txt
      d2Hprecip = 0;
      d2Hcanopy = 0 ;
      d2Hsnowpack = 0;
      d2Hponding = 0;
      d2HsoilL1 = 0;
      d2HsoilL2 = 0;
      d2HsoilL3 = 0;
      d2Hgrndwater = 0;
      d2HTot = 0;
      d2Hvadose = 0;
      d2Hrootzone = 0;
      d2HET = 0;
      d2HevapS = 0;
      d2HevapI = 0;
      d2HevapT = 0;
      d2Hleakage = 0;
      d2HOvlndOut = 0;
      d2HGWOut = 0;
      d2Hgwtochn = 0;
      d2Hsrftochn = 0;
      d2Hrecharge = 0;
      d2HOut = 0;
      //calculate initial storages
      initcanopy_d2H = 0;
      initsnowpack_d2H = AccountTrckStorages(b->getSnowWaterEquiv(), trck->getd2Hsnowpack(), b);
      initponding_d2H = AccountTrckStorages(b->getPondingWater(), trck->getd2Hsurface(), b);
      initL1_d2H = AccountTrckStorages(b->getSoilWaterDepthL1(), trck->getd2Hsoil1(), b);
      initL2_d2H = AccountTrckStorages(b->getSoilWaterDepthL2(), trck->getd2Hsoil2(), b);
      initL3_d2H = AccountTrckStorages(b->getSoilWaterDepthL3(), trck->getd2Hsoil3(), b);
      initGW_d2H = AccountTrckStorages(b->getInitGroundwater(), trck->getd2Hgroundwater(), b);
    }

    if(ctrl->sw_trck && ctrl->sw_18O){
      canopy_d18O = 0;
      snowpack_d18O = 0;
      ponding_d18O = 0;
      soilL1_d18O = 0;
      soilL2_d18O = 0;
      soilL3_d18O = 0;
      grndwater_d18O = 0;
      precipitation_d18O = 0;
      evaporationS_d18O = 0;
      evaporationI_d18O = 0;
      transpiration_d18O = 0;
      leakage_d18O = 0;
      ovlndflow_d18O = 0;
      gwtrflow_d18O = 0;
      // Values for Basind18OSummary.txt
      d18Oprecip = 0;
      d18Ocanopy = 0 ;
      d18Osnowpack = 0;
      d18Oponding = 0;
      d18OsoilL1 = 0;
      d18OsoilL2 = 0;
      d18OsoilL3 = 0;
      d18Ogrndwater = 0;
      d18OTot = 0;
      d18Ovadose = 0;
      d18Orootzone = 0;
      d18OET = 0;
      d18OevapS = 0;
      d18OevapI = 0;
      d18OevapT = 0;
      d18Oleakage = 0;
      d18OOvlndOut = 0;
      d18OGWOut = 0;
      d18Ogwtochn = 0;
      d18Osrftochn = 0;
      d18Orecharge = 0;
      d18OOut = 0;
      //calculate initial storages
      initcanopy_d18O = 0;
      initsnowpack_d18O = AccountTrckStorages(b->getSnowWaterEquiv(), trck->getd18Osnowpack(), b);
      initponding_d18O = AccountTrckStorages(b->getPondingWater(), trck->getd18Osurface(), b);
      initL1_d18O = AccountTrckStorages(b->getSoilWaterDepthL1(), trck->getd18Osoil1(), b);
      initL2_d18O = AccountTrckStorages(b->getSoilWaterDepthL2(), trck->getd18Osoil2(), b);
      initL3_d18O = AccountTrckStorages(b->getSoilWaterDepthL3(), trck->getd18Osoil3(), b);
      initGW_d18O = AccountTrckStorages(b->getInitGroundwater(), trck->getd18Ogroundwater(), b);
    }

    if(ctrl->sw_trck && ctrl->sw_Cl){
      canopy_cCl = 0;
      snowpack_cCl = 0;
      ponding_cCl = 0;
      soilL1_cCl = 0;
      soilL2_cCl = 0;
      soilL3_cCl = 0;
      grndwater_cCl = 0;
      precipitation_cCl = 0;
      leakage_cCl = 0;
      ovlndflow_cCl = 0;
      gwtrflow_cCl = 0;
      // Values for BasincClSummary.txt
      cClprecip = 0;
      cClcanopy = 0 ;
      cClsnowpack = 0;
      cClponding = 0;
      cClsoilL1 = 0;
      cClsoilL2 = 0;
      cClsoilL3 = 0;
      cClgrndwater = 0;
      cClTot = 0;
      cClvadose = 0;
      cClrootzone = 0;
      cClleakage = 0;
      cClOvlndOut = 0;
      cClGWOut = 0;
      cClgwtochn = 0;
      cClsrftochn = 0;
      cClrecharge = 0;
      cClOut = 0;
      //calculate initial storages
      initcanopy_cCl = 0;
      initsnowpack_cCl = AccountTrckStorages(b->getSnowWaterEquiv(), trck->getcClsnowpack(), b);
      initponding_cCl = AccountTrckStorages(b->getPondingWater(), trck->getcClsurface(), b);
      initL1_cCl = AccountTrckStorages(b->getSoilWaterDepthL1(), trck->getcClsoil1(), b);
      initL2_cCl = AccountTrckStorages(b->getSoilWaterDepthL2(), trck->getcClsoil2(), b);
      initL3_cCl = AccountTrckStorages(b->getSoilWaterDepthL3(), trck->getcClsoil3(), b);
      initGW_cCl = AccountTrckStorages(b->getInitGroundwater(), trck->getcClgroundwater(), b);
    }

    if(ctrl->sw_trck && ctrl->sw_Age){
      precipitation_Age = 0;
      canopy_Age = 0 ;
      snowpack_Age = 0;
      ponding_Age = 0;
      soilL1_Age = 0;
      soilL2_Age = 0;
      soilL3_Age = 0;
      grndwater_Age = 0;
      evaporationS_Age = 0;
      evaporationI_Age = 0;
      transpiration_Age = 0;
      leakage_Age = 0;
      ovlndflow_Age = 0;
      gwtrflow_Age = 0;
      // Values for BasinAgeSummary.txt
      Agecanopy = 0 ;
      Agesnowpack = 0;
      Ageponding = 0;
      AgesoilL1 = 0;
      AgesoilL2 = 0;
      AgesoilL3 = 0;
      Agegrndwater = 0;
      AgeTot = 0;
      Agevadose = 0;
      Agerootzone = 0;
      AgeET = 0;
      AgeevapS = 0;
      AgeevapI = 0;
      AgeevapT = 0;
      Ageleakage = 0;
      AgeOvlndOut = 0;
      AgeGWOut = 0;
      Agegwtochn = 0;
      Agesrftochn = 0;
      Agerecharge = 0;
      AgeOut = 0;
      //calculate initial storages
      initcanopy_Age = 0;
      initsnowpack_Age = AccountTrckStorages(b->getSnowWaterEquiv(), trck->getAgesnowpack(), b);
      initponding_Age = AccountTrckStorages(b->getPondingWater(), trck->getAgesurface(), b);
      initL1_Age = AccountTrckStorages(b->getSoilWaterDepthL1(), trck->getAgesoil1(), b);
      initL2_Age = AccountTrckStorages(b->getSoilWaterDepthL2(), trck->getAgesoil2(), b);
      initL3_Age = AccountTrckStorages(b->getSoilWaterDepthL3(), trck->getAgesoil3(), b);
      initGW_Age = AccountTrckStorages(b->getInitGroundwater(), trck->getAgegroundwater(), b);
    }
    // ---------------------------------------------------------------------------------------

  }

};

#endif /* BUDGET_H_ */
