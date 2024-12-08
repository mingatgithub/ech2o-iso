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
 * CalculateInitialStreamStorage.cpp
 *
 *  Created on: Aug 8, 2015
 *      Author: marco
 */


#include "Basin.h"

int Basin::CalcInitialStreamStorage(Control & ctrl){

	UINT4 r, c;
	REAL8 w, n, a, sqrtS, Q;

	UINT4 length = _vSortedGrid.cells.size();

#pragma omp parallel default(none)\
  private(  r,c, w, n, a, sqrtS,Q), shared(length, ctrl)
	for (UINT4 j = 0; j < length ; j++)
	{
					r = _vSortedGrid.cells[j].row;
					c = _vSortedGrid.cells[j].col;

      _ponding->matrix[r][c] = 0.0;

      if(ctrl.sw_channel){
					sqrtS = powl(_slope->matrix[r][c], 0.5);
					w = _channelwidth->matrix[r][c];
			        n = _Manningn->matrix[r][c];
					a = powl(powl(w,0.67)*n/sqrtS, 0.6); //wetted perimeter is approximated with channel width
					Q = _Disch_old->matrix[r][c];

	_Schannel->matrix[r][c] = (w > 0) && (Q > 0) ? a * powl(Q, 0.6)/_dx : 0.0;
      } else
	_Schannel->matrix[r][c] = 0.0;

	}

	return EXIT_SUCCESS;
}
