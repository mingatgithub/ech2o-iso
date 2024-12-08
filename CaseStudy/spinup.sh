#!/bin/bash

COUNT=1
NITER=11

while [ $COUNT -lt $NITER ]
do
	echo "Running iteration ${COUNT}"
	xterm -e ech2o_iso config.ini
	
	echo "finished run, copying files"
	cp -f ./Outputs/root0_00.365 ./Spatial/root_0.map
	cp -f ./Outputs/p0_00000.365 ./Spatial/p_0.map
	cp -f ./Outputs/ntr0_000.365 ./Spatial/ntr_0.map
	cp -f ./Outputs/lai0_000.365 ./Spatial/lai_0.map
	cp -f ./Outputs/hgt0_000.365 ./Spatial/hgt_0.map
	cp -f ./Outputs/bas0_000.365 ./Spatial/bas_0.map
	cp -f ./Outputs/age0_000.365 ./Spatial/age_0.map

	cp -f ./Outputs/SWE00000.365 ./Spatial/Init_swe.map
	cp -f ./Outputs/SWC1_000.365 ./Spatial/Init_soilmoisture.L1.map
	cp -f ./Outputs/SWC2_000.365 ./Spatial/Init_soilmoisture.L2.map
	cp -f ./Outputs/SWC3_000.365 ./Spatial/Init_soilmoisture.L3.map
	cp -f ./Outputs/Ts000000.365 ./Spatial/Init_soiltemp.map
	cp -f ./Outputs/Q0000000.365 ./Spatial/Init_streamflow.map

	cat ./Outputs/lai_0.tab >> ./Outputs/laiaccum.txt
	cat ./Outputs/NPP_0.tab >> ./Outputs/NPPaccum.txt
	cat ./Outputs/SoilMoistureAv.tab >> ./Outputs/SWCaccum.txt


	let "COUNT++"
done
	echo "Finished simulation\n"


