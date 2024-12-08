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
 *    Marco Maneta
 *******************************************************************************/
/*
 * Exfiltration.cpp
 *
 *  Created on: Jul 13, 2010
 *      Author: Marco.Maneta
 */

#include "Basin.h"

REAL8 Basin::ExfiltrationCapacity(	REAL8 theta,//soil moisture
									REAL8 dt,
									UINT4 r,
									UINT4 c)

{

	REAL8 S; //relative soil moisture
	REAL8 thetar = _theta_r->matrix[r][c];
	REAL8 poros = _porosity->matrix[r][c];
	REAL8 Ksat = _Ksat->matrix[r][c];
	REAL8 psiae = _psi_ae->matrix[r][c];
	REAL8 poresizeind = _BClambda->matrix[r][c];


	S = (theta - thetar)/(poros - thetar);
	if (S<=0)
		S = 0.01;

	REAL8 num = 8 * poros * Ksat * 86400 * psiae; //m2d-1
	REAL8 denom = 3 * (1 + 3*poresizeind) * (1 + 4*poresizeind); //dimensionless


	return (powl(S, 1/(2*poresizeind) + 2) * sqrt(num/denom) * sqrt(1/86400)) / sqrt(1/dt); // returns water limited exfiltration volume in ms-1


}
