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
 * AdvanceClimateMaps.cpp
 *
 *  Created on: Oct 18, 2009
 *      Author: Marco Maneta
 */

#include "Atmosphere.h"



int Atmosphere::AdvanceClimateMaps(Control &ctrl){

// #ifdef _OPENMP
// UINT4 number_threads;
//  number_threads = omp_get_num_threads();
// #endif
size_t errCount = 0;

// #pragma omp parallel num_threads(8) if (number_threads > 1)
// {
try {
    if (UpdateClimateMap(ifLdown, *_Ldown) != _vSsortedGridTotalCellNumber) {
        throw string("error advancing long wave time step");
    }
    if (UpdateClimateMap(ifSdown, *_Sdown) != _vSsortedGridTotalCellNumber) {
        throw string("error advancing short wave time step");
    }
    if (UpdateClimateMap(ifTp, *_Tp) != _vSsortedGridTotalCellNumber) {
        throw string("error advancing av air temp time step");
    }
    if (UpdateClimateMap(ifMaxTp, *_MaxTp) != _vSsortedGridTotalCellNumber) {
        throw string("error advancing max air temp time step");
    }
    if (UpdateClimateMap(ifMinTp, *_MinTp) != _vSsortedGridTotalCellNumber) {
        throw string("error advancing min air temp time step");
    }

    if (ctrl.sw_trck && ctrl.sw_2H) {
        if (UpdateClimateMap(ifd2Hprecip, *_d2Hprecip) != _vSsortedGridTotalCellNumber) {
            throw string("error advancing input d2H time step");
        }
    }

    if (UpdateClimateMap(ifPrecip, *_Precip) != _vSsortedGridTotalCellNumber) {
        throw string("error advancing precipitation time step");
    }
    AdjustPrecip();

    if (UpdateClimateMap(ifRelHumid, *_Rel_humid) != _vSsortedGridTotalCellNumber) {
        throw string("error advancing RH time step");
    }

    if (ctrl.sw_trck && ctrl.sw_18O) {
        if (UpdateClimateMap(ifd18Oprecip, *_d18Oprecip) != _vSsortedGridTotalCellNumber) {
            throw string("error advancing input d18O time step");
        }
    }

    if (UpdateClimateMap(ifWindSpeed, *_Wind_speed) != _vSsortedGridTotalCellNumber) {
        throw string("error advancing wind speed time step");
    }

    if (ctrl.sw_trck && ctrl.sw_Cl) {
        if (UpdateClimateMap(ifcClprecip, *_cClprecip) != _vSsortedGridTotalCellNumber) {
            throw string("error advancing input chloride time step");
        }
    }
} catch (string &s) {
    cout << s;
    ++errCount;
}//catch
// } //parallel region

 if (errCount != 0)
   throw;

 return EXIT_SUCCESS;

}
