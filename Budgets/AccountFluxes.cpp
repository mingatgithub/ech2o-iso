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
 * AccountInputFluxes.cpp
 *
 *  Created on: Mar 8, 2010
 *      Author: Marco Maneta
 */
#include "Budget.h"

double Budget::AccountFluxes(const grid *map, const Basin *b) {

  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();

#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)

  for (UINT4 i = 0; i < length; i++) {

    r = b->getSortedGrid().cells[i].row;
    c = b->getSortedGrid().cells[i].col;

    result += (map->matrix[r][c] * dx * dx * dt);
  }

  return result;
}

double Budget::AccountFluxes(const grid *map, const Atmosphere *b) {

  UINT4 zones = b->getSortedGrid().size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();

#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)

  for (UINT4 i = 0; i < zones; i++)
    for (UINT4 j = 0; j < b->getSortedGrid()[i].cells.size(); j++) {

      r = b->getSortedGrid()[i].cells[j].row;
      c = b->getSortedGrid()[i].cells[j].col;

      result += (map->matrix[r][c] * dx * dx * dt);
    }

  return result;
}

double Budget::AccountFluxes(const vectCells *timeseries, const Basin *b) {

  UINT4 length = timeseries->cells.size(); //b->getSortedGrid().cells.size();

  REAL8 result = 0;

#pragma omp parallel for			\
  reduction (+:result)

  for (UINT4 i = 0; i < length; i++) {

    result += timeseries->cells[i].val * dt;
  }

  return result;
}

// --- Tracking : uses two maps (or vectors) to multiply -----------------------------------

double Budget::AccountTrckFluxes(const grid *map1, const grid *map2, const Basin *b) {

  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();

#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)

  for (UINT4 i = 0; i < length; i++) {

    r = b->getSortedGrid().cells[i].row;
    c = b->getSortedGrid().cells[i].col;

    if(map1->matrix[r][c] > RNDOFFERR)
    result +=  (map1->matrix[r][c] * map2->matrix[r][c] *dx*dx*dt);

  }

  return result;
}

// Precip isotopes
double Budget::AccountTrckFluxes(const grid *map1, const grid *map2, const Atmosphere *a) {

  UINT4 zones = a->getSortedGrid().size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = a->getCellSize();

#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)

  for (UINT4 i = 0; i < zones; i++)
    for (UINT4 j = 0; j < a->getSortedGrid()[i].cells.size(); j++) {

      r = a->getSortedGrid()[i].cells[j].row;
      c = a->getSortedGrid()[i].cells[j].col;

      if(map1->matrix[r][c] > RNDOFFERR)
        result += (map1->matrix[r][c] * map2->matrix[r][c] * dx * dx * dt);
    }

  return result;
}

// Precip age: it's 0 at entry, but the budgets must account for aging previously-input precip!
double Budget::AccountTrckFluxes(const grid *map, const Atmosphere *a){

  UINT4 zones = a->getSortedGrid().size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = a->getCellSize();

#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)

  for (UINT4 i = 0; i < zones; i++)
    for (UINT4 j = 0; j < a->getSortedGrid()[i].cells.size(); j++) {

      r = a->getSortedGrid()[i].cells[j].row;
      c = a->getSortedGrid()[i].cells[j].col;

      if(map->matrix[r][c] == map->matrix[r][c])
      result += map->matrix[r][c] * dx * dx * dt * dt / 86400;
    }

  return result;
}

double Budget::AccountTrckFluxes(const vectCells *timeseries1, const vectCells *timeseries2) {

  UINT4 length = timeseries1->cells.size(); //b->getSortedGrid().cells.size();
  REAL8 result = 0;

#pragma omp parallel for default(shared)	\
  reduction (+:result)

  for (UINT4 i = 0; i < length; i++){
    if(timeseries1->cells[i].val > RNDOFFERR)
    result +=  (timeseries1->cells[i].val * timeseries2->cells[i].val * dt);
  }

  return result;
}

// == AgeReporting stuff ------------------------------------------------------------------------------

double Budget::AccountTrckFluxes2(const grid *map1, const grid *map2, const Basin *b) {

  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 numer = 0;
  REAL8 denom = 0;

#pragma omp parallel default(shared) private(r,c)	\

  {
#pragma omp for reduction (+:numer, denom)

    for (UINT4 i = 0; i < length; i++) {

      r = b->getSortedGrid().cells[i].row;
      c = b->getSortedGrid().cells[i].col;

      if(map1->matrix[r][c] > RNDOFFERR){
        numer += map1->matrix[r][c] * map2->matrix[r][c];
        denom += map1->matrix[r][c];
      }

      // cout << r << " " << c << " " <<
      //    map1->matrix[r][c] << " " << map2->matrix[r][c] << " " <<
      //  numer << " " << denom << " " << numer/ denom << endl ;

    }
  }

  result = denom > RNDOFFERR ? numer / denom : sqrt(-1) ;

  return result;
}

// Precip isotopes
double Budget::AccountTrckFluxes2(const grid *map1, const grid *map2, const Atmosphere *a) {

  UINT4 zones = a->getSortedGrid().size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 numer = 0;
  REAL8 denom = 0;

#pragma omp parallel default(shared) private(r,c) \

  {
#pragma omp for reduction (+:numer, denom)

  for (UINT4 i = 0; i < zones; i++)
    for (UINT4 j = 0; j < a->getSortedGrid()[i].cells.size(); j++) {

      r = a->getSortedGrid()[i].cells[j].row;
      c = a->getSortedGrid()[i].cells[j].col;

      if(map1->matrix[r][c] > RNDOFFERR){
        numer += map1->matrix[r][c] * map2->matrix[r][c];
        denom += map1->matrix[r][c];
      }

    }
  }
  result = denom > RNDOFFERR ? numer / denom : sqrt(-1) ;

  return result ;
}

double Budget::AccountTrckFluxes2(const vectCells *timeseries1, const vectCells *timeseries2) {

  UINT4 length = timeseries1->cells.size(); //b->getSortedGrid().cells.size();
  REAL8 result = 0;
  REAL8 numer = 0;
  REAL8 denom = 0;

  //cout << length << endl;

#pragma omp parallel default(shared)
  {

#pragma omp for reduction (+:numer, denom)

    for (UINT4 i = 0; i < length; i++){
      if(timeseries1->cells[i].val > RNDOFFERR){
      numer +=  timeseries1->cells[i].val * timeseries2->cells[i].val;
      //cout << timeseries1->cells[i].val << endl;
      //cout << timeseries2->cells[i].val << endl;
      denom +=  timeseries1->cells[i].val;
    }
  }
  }
  result = denom > RNDOFFERR ? numer / denom : sqrt(-1) ;

  return result;
}

// -- Flux-weighted average across evaporative losses
double Budget::AccountTrckET(const grid* evapS, const grid* CevapS,
			     const grid* evapI, const grid* CevapI,
			     const grid* evapT, const grid* CevapT,
			     const Basin *b)
{

  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 numer = 0;
  REAL8 denom = 0;

#pragma omp parallel default(shared) private(r,c)

  {
#pragma omp for reduction (+:numer, denom)
    for (UINT4 i = 0; i< length; i++){

      r = b->getSortedGrid().cells[i].row;
      c = b->getSortedGrid().cells[i].col;

      if(evapS->matrix[r][c] > RNDOFFERR){
        numer += evapS->matrix[r][c]* CevapS->matrix[r][c] ;
        denom += evapS->matrix[r][c] ;
      }
      if(evapI->matrix[r][c] > RNDOFFERR){
        numer += evapI->matrix[r][c]* CevapI->matrix[r][c] ;
        denom += evapI->matrix[r][c] ;
      }
      if(evapT->matrix[r][c] > RNDOFFERR){
        numer += evapT->matrix[r][c]* CevapT->matrix[r][c] ;
        denom += evapT->matrix[r][c] ;
      }
    }
  }

  result = denom > RNDOFFERR ? numer / denom : sqrt(-1) ;
  return result;
}

// -- Flux-weighted average across outputs
double Budget::AccountTrckOut(const grid* evapS, const grid* CevapS,
			      const grid* evapI, const grid* CevapI,
			      const grid* evapT, const grid* CevapT,
			      const grid* leakage, const grid* Cleakage,
			      const vectCells *OvlndOut, const vectCells *COvlndOut,
			      const vectCells *GWOut, const vectCells *CGWOut,
			      const Basin *b)
{

  UINT4 length1 = b->getSortedGrid().cells.size();
  UINT4 length2 = OvlndOut->cells.size(); //b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 numer1 = 0;
  REAL8 numer2 = 0;
  REAL8 denom1 = 0;
  REAL8 denom2 = 0;
  REAL8 dx = b->getCellSize();

#pragma omp parallel default(shared) private(r,c)

  {
#pragma omp for reduction (+:numer1, denom1)

    for (UINT4 i = 0; i< length1; i++){

      r = b->getSortedGrid().cells[i].row;
      c = b->getSortedGrid().cells[i].col;

      if(evapS->matrix[r][c] > RNDOFFERR){
        numer1 += evapS->matrix[r][c]* CevapS->matrix[r][c] *dx*dx;
        denom1 += evapS->matrix[r][c] *dx*dx;
      }
      if(evapI->matrix[r][c] > RNDOFFERR){
        numer1 += evapI->matrix[r][c]* CevapI->matrix[r][c] *dx*dx;
        denom1 += evapI->matrix[r][c]*dx*dx ;
      }
      if(evapT->matrix[r][c] > RNDOFFERR){
        numer1 += evapT->matrix[r][c]* CevapT->matrix[r][c] *dx*dx;
        denom1 += evapT->matrix[r][c] *dx*dx;
      }
      if(leakage->matrix[r][c] > RNDOFFERR){
        numer1 += leakage->matrix[r][c]* Cleakage->matrix[r][c] *dx*dx;
        denom1 += leakage->matrix[r][c] *dx*dx;
      }

    }

    //cout << numer1 << " " << denom1 << endl;
    //cout << numer2 << " " << denom2 << endl;

#pragma omp for reduction (+:numer2, denom2)

    for (UINT4 j = 0; j < length2; j++) {
      if(OvlndOut->cells[j].val > RNDOFFERR){
        numer2 +=  OvlndOut->cells[j].val * COvlndOut->cells[j].val ;
        denom2 +=  OvlndOut->cells[j].val ;
      }
      if(GWOut->cells[j].val > RNDOFFERR){
        numer2 +=  GWOut->cells[j].val * CGWOut->cells[j].val;
        denom2 +=  GWOut->cells[j].val ;
      }

    }
  }
  //cout << numer1 << " " << denom1 << endl;
  //cout << numer2 << " " << denom2 << endl;

  result = (denom1 + denom2) > RNDOFFERR ?
    (numer1 + numer2) / (denom1 + denom2) : sqrt(-1);

  //cout << result << endl;

 return result;
}

// -- Flux-weighted average across outputs
double Budget::AccountTrckOut(const grid* leakage, const grid* Cleakage,
			      const vectCells *OvlndOut, const vectCells *COvlndOut,
			      const vectCells *GWOut, const vectCells *CGWOut,
			      const Basin *b)
{

  UINT4 length1 = b->getSortedGrid().cells.size();
  UINT4 length2 = OvlndOut->cells.size(); //b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 numer1 = 0;
  REAL8 numer2 = 0;
  REAL8 denom1 = 0;
  REAL8 denom2 = 0;
  REAL8 dx = b->getCellSize();

#pragma omp parallel default(shared) private(r,c)

  {
#pragma omp for reduction (+:numer1, denom1)

    for (UINT4 i = 0; i< length1; i++){

      r = b->getSortedGrid().cells[i].row;
      c = b->getSortedGrid().cells[i].col;

      if(leakage->matrix[r][c] > RNDOFFERR){
      numer1 += leakage->matrix[r][c]* Cleakage->matrix[r][c] *dx*dx ;
      denom1 += leakage->matrix[r][c] * dx * dx ;
    }
    }

#pragma omp for reduction (+:numer2, denom2)

    for (UINT4 j = 0; j < length2; j++) {
      if(OvlndOut->cells[j].val > RNDOFFERR){
        numer2 +=  OvlndOut->cells[j].val * COvlndOut->cells[j].val ;
        denom2 +=  OvlndOut->cells[j].val ;
      }
      if(GWOut->cells[j].val > RNDOFFERR){
        numer2 +=  GWOut->cells[j].val * CGWOut->cells[j].val;
        denom2 +=  GWOut->cells[j].val ;
      }

    }
  }

  result = (denom1 + denom2) > RNDOFFERR ?
    (numer1 + numer2) / (denom1 + denom2) : sqrt(-1);

  return result;
}
