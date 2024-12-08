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
 * GenerateConfigTemplate.cpp
 *
 *  Created on: Jan 12, 2016
 *      Author: marcomaneta
 */

#include  <fstream>
#include <unistd.h>

#include "Sativa.h"

void GenerateConfigTemplate(const char *fn){

  ofstream ofOut;

  try{

    if (access(fn, F_OK) != -1) {

      cout << "File exists. Do you want to overwrite? (y, n):  " << endl;
      char c;
      cin.get(c);
      switch (c) {
      case 'y':
	break;
      case 'n':
	exit(EXIT_SUCCESS);
	break;
      default:
	cout << "Not a valid option. Bye" << endl;
	exit(EXIT_SUCCESS);

      }
    }

    ofOut.open(fn);
    if(!ofOut)
      throw std::ios::failure("Error opening file ");


    ofOut << "#ECH2O configuration file v2" << std::endl << std::endl;
    ofOut << "# Please, check Appendix A of Documentation" << std::endl;
    ofOut << "# for units of parameters and variables  " << std::endl;
    ofOut << "# (http://ech2o-iso.readthedocs.io/en/latest/Keywords.html)" << std::endl << std::endl;

    ofOut << "#" << endl << "#Folder section" << endl << "#" << endl << endl;

    ofOut << "Maps_Folder = ./Spatial" << endl;
    ofOut << "Clim_Maps_Folder = ./Climate" << endl;
    ofOut << "Output_Folder = ./Outputs" << endl << endl;

    ofOut << "#" << endl << "#Water tracking (isotopes and/or ages)" << endl;
    ofOut << "Tracking = 1" << endl ;
    ofOut << "TrackingConfig = ./configTrck.ini" << endl << endl; 

    ofOut << "#" << endl << "# Options section" << endl << "#" << endl << endl;
    
    ofOut << "MapTypes = csf" << endl;
    ofOut << "Species_State_Variable_Input_Method = tables # maps or tables" << endl << endl;

    ofOut << "# -- Boolean switches" << endl;
    ofOut << "Reinfiltration = 1" << endl;
    ofOut << "Channel = 1" << endl ;
    ofOut << "Channel_Infiltration = 0" << endl << endl ;

    ofOut << "# -- Toggle switches" << endl << endl ;
    ofOut << "# Hydraulic conductivity and porosity profiles --> 3 options each: " << endl ;
    ofOut << "# 0=constant, 1=exponential profile, 2=defined for each layer" << endl ;
    ofOut << "Hydraulic_Conductivity_profile = 0" << endl ;
    ofOut << "Porosity_profile = 0 " << endl << endl;
    
    ofOut << "# Vegetation dynamics (via allocation): 3 modes" << endl;
    ofOut << "# 0 -> desactivated, no dynamic allocation and constant LAI to initial value" <<endl;
    ofOut << "# 1 -> fully activated" << endl;
    ofOut << "# 2 -> partially activated, except that LAI is prescribed via an input file" << endl;
    ofOut << "Vegetation_dynamics = 1" << endl;
    ofOut << "# Used only if Vegetation_dynamics = 2. Files names for each species is" << endl;
    ofOut << "# name below + '_'+ species number (starting at 0) + '.bin'" << endl;
    ofOut<< "TimeSeries_LAI = LAI" << endl << endl;
    
    ofOut << "# Aerodynamic resistance choices: " << endl;
    ofOut << "# 0 = Penman Monteith option " << endl;
    ofOut << "# 1 = Thom and Oliver 1977 " << endl;
    ofOut << "Aerodyn_resist_opt = 0 " << endl << endl;

    ofOut << "# Soil resistance to vapor diffusion choices: " << endl;
    ofOut << "# 0 = No resistance" << endl;
    ofOut << "# 1 = Passerat de Silans et al. 1989" << endl;
    ofOut << "# 2 = Sellers et al. 1992" << endl;
    ofOut << "# 3 = Sakaguchi and Zeng 2009 (CLM 3.5)" << endl;
    ofOut << "Soil_resistance_opt = 3 " << endl << endl;

    ofOut << "#" << endl << "# Time variables section" << endl << "#" << endl;
    ofOut << "Simul_start = 0 #" << endl;
    ofOut << "Simul_end = 31536000 # seconds (365 days)" << endl;
    ofOut << "Simul_tstep = 86400 # seconds (daily)" << endl;
    ofOut << "Clim_input_tstep = 86400 # seconds (daily)" << endl;
    ofOut << "Report_interval = 86400 # seconds (daily)" << endl ;
    ofOut << "ReportMap_interval = 86400 # seconds (daily)" << endl ;
    ofOut << "ReportMap_starttime = 86400 # seconds (from first time step)" << endl << endl;

    ofOut << "# -------------------------------------------------" << endl;
    ofOut << "# == Climate input information" << endl;
    ofOut << "# (files must be folder pointed by Clim_Maps_Folder)" << endl << endl ;

    ofOut << "Precipitation = Precip.bin " << " # Precip rate in meters/second"<< endl;
    ofOut << "AirTemperature = Tavg.bin " << " # Average air temperature in degC" << endl;
    ofOut << "MaxAirTemp = Tmax.bin " << " # Maximum air temperature in degC" << endl;
    ofOut << "MinAirTemp = Tmin.bin " << " # Minimum air temperature in degC"<< endl;
    ofOut << "RelativeHumidity = RH.bin " << " # air relative humidity in kPa/kPa"<< endl;
    ofOut << "WindSpeed = windspeed.bin " << " # Wind speed in meters/second" << endl;
    ofOut << "IncomingLongWave = Ldown.bin " << " # Downwelling longwave radiation in W/sq.meter" << endl;
    ofOut << "IncomingShortWave = Sdown.bin " << " # Solar radiation in W/sq.meter" << endl << endl;
    ofOut << "ClimateZones = ClimZones.map  # Climate zones applied to all above forcing fields" << endl;
    //ofOut << "Snow_rain_temp_threshold = 2  # Snow to rain temperatures threshold in degC" << endl;
    ofOut << "Isohyet_map = isohyet.map  # Precipitation multiplier map"<< endl;

    
    ofOut << "# -------------------------------------------------" << endl;
    ofOut << "# == Spatial and vegetation input information" << endl;
    ofOut << "# (files must be folder pointed by Maps_Folder)" << endl << endl ;

    ofOut << "# -- Site geometry" << endl << endl ;

    ofOut << "DEM = DEM.map" << endl;
    ofOut << "Slope = slope.map " << endl;
    ofOut << "# Drainage network" << endl ;
    ofOut << "local_drain_direc = ldd.map" << endl;

    ofOut << "# -- Soil parameters" << endl << endl ;

    ofOut << "# Kh: full profile / top-of-profile / L1 value (resp.) " << endl ;
    ofOut << "# if Hydraulic_Conductivity_profile = 0 / 1 / 2 (resp.)" << endl;
    ofOut << "Horiz_Hydraulic_Conductivity = Khsat.map " << endl;
    ofOut << "# if Hydraulic_Conductivity_profile = 1" << endl;
    ofOut << "Horiz_Hydraulic_Conductivity_Profile_Coeff = kKhsat.map " << endl;
    ofOut << "# if Hydraulic_Conductivity_profile = 2" << endl;
    ofOut << "Horiz_Hydraulic_Conductivity_Layer2 = Khsat.L2.map # if Ksat profile = 2" << endl;
    ofOut << "Horiz_Hydraulic_Conductivity_Layer3 = Khsat.L3.map # if Ksat profile = 2" << endl << endl ;
    
    ofOut << "Vert_Horz_Anis_ratio = KvKh.map " << endl;
    ofOut << "Terrain_Random_Roughness = randrough.map " << endl << endl ;

    ofOut << "# porosity: full profile / top-of-profile / L1 value (resp.) " << endl ;
    ofOut << "# if Porosity_profile = 0 / 1 / 2 (resp.)" << endl;
    ofOut << "Porosity = poros.map " << endl;
    ofOut << "# if Porosity_profile = 1" << endl;
    ofOut << "Porosity_Profile_Coeff = kporos.map # if poros profile = 1" << endl;
    ofOut << "# if Porosity_profile = 2" << endl;
    ofOut << "Porosity_Layer2 = poros.L2.map # if poros profile = 2" << endl;
    ofOut << "Porosity_Layer3 = poros.L3.map # if poros profile = 2" << endl << endl ;
    
    ofOut << "Air_entry_pressure = psi_ae.map " << endl;
    ofOut << "Brooks_Corey_lambda = BClambda.map " << endl;
    ofOut << "Residual_soil_moisture = theta_r.map " << endl;
    ofOut << "Soil_depth = soildepth.map " << endl;
    ofOut << "Depth_soil_layer_1 = soildepth.L1.map " << endl;
    ofOut << "Depth_soil_layer_2 = soildepth.L2.map " << endl;
    ofOut << "Veget_water_use_param1 = Wc.map " << endl;
    ofOut << "Veget_water_use_param2 = Wp.map " << endl;
    //ofOut << "Fraction_roots_soil_layer_1 = rootfrac1.map " << endl;
    //ofOut << "Fraction_roots_soil_layer_2 = rootfrac2.map " << endl;
    //ofOut << "Root_profile_coeff = Kroot.map " << endl;
    ofOut << "Soil_bedrock_leakance = leakance.map " << endl << endl;
    ofOut << "Albedo = albedo.map" << endl;
    ofOut << "Surface_emissivity = emissivity.map" << endl;
    ofOut << "Dry_Soil_Heat_Capacity = soilheatcap.map" << endl;
    ofOut << "Dry_Soil_Therm_Cond = soilthermalK.map" << endl;
    ofOut << "Damping_depth = dampdepth.map" << endl;
    ofOut << "Temp_at_damp_depth = temp_damp.map" << endl;
    ofOut << "Snow_rain_temp_threshold = SnowRainTemp.map " << endl;
    ofOut << "Snow_Melt_Coeff = snowmeltCoeff.map" << endl << endl;

    ofOut << "# -- Channel parameters" << endl << endl ;
    
    ofOut << "channel_width = chanwidth.map" << endl;
    ofOut << "channel_gw_transfer_param = chanparam.map" << endl;
    ofOut << "mannings_n = chanmanningn.map" << endl;

    ofOut << "# -- Initial hydrological/temperature conditions  " << endl << endl ;
    
    ofOut << "Streamflow = Init_streamflow.map " << endl;
    ofOut << "snow_water_equivalent = Init_SWE.map " << endl;
    ofOut << "Soil_moisture_1 = Init_soilmoisture.L1.map " << endl;
    ofOut << "Soil_moisture_2 = Init_soilmoisture.L2.map " << endl;
    ofOut << "Soil_moisture_3 = Init_soilmoisture.L3.map " << endl;
    ofOut << "Soil_temperature = Init_soiltemp.map " << endl << endl;

    ofOut << "# -- Vegetation distribution and initial states " << endl << endl ;

    ofOut << "Species_Parameters = SpeciesParams.tab " << endl << endl;
    ofOut << "Number_of_Species = 1  " << endl;
    ofOut << "# the above number can be smaller than the rows found in Species_Parameters" << endl;
    ofOut << "ForestPatches = patches.map  " << endl;
    ofOut << "# needed in all cases, in 'maps' mode ForestPatches can be a uniform map" << endl;
    
    ofOut << "# Initial states: 2 options" << endl;
    ofOut << "# 1. if Species_State_Variable_Input_Method = tables, the tables below are needed" << endl;
    ofOut << "Species_Proportion_Table = SpecsProp.tab " << endl;
    ofOut << "Species_StemDensity_Table = SpecsStemDens.tab " << endl;
    ofOut << "Species_LAI_Table = SpecsLAI.tab  # (not needed of Vegetation_dynamics = 2)" << endl;
    ofOut << "Species_AGE_Table = SpecsAge.tab " << endl;
    ofOut << "Species_BasalArea_Table = SpeciesBasalArea.tab " << endl;
    ofOut << "Species_Height_table = SpeciesHeight.tab " << endl;
    ofOut << "Species_RootMass_table = SpecsRootDensity.tab " << endl ;
    ofOut << "# 2. if Species_State_Variable_Input_Method = maps, then initial conditions" << endl;
    ofOut << "#    are read from maps with pre-established names: " << endl;
    ofOut << "#    - p_i.map for species fraction cover" << endl;
    ofOut << "#    - ntr_i.map for stem density (tree m-2)" << endl;
    ofOut << "#    - lai_i.map for LAI (not needed of Vegetation_dynamics = 2)" << endl ;
    ofOut << "#    - age_i.map for vegetation age (years)" << endl;
    ofOut << "#    - bas_i.map for stem basal area (m2)" << endl;
    ofOut << "#    - hgt_i.map for vegetation height (m)" << endl;
    ofOut << "#    - root_i.map for root density (g m-2)" << endl << endl << endl ; 


    ofOut << "# --------------------------------------------------------------" << endl;
    ofOut << "# == Report flags for output files" << endl;
    ofOut << "# (files will be written in be folder pointed by Output_Folder)" << endl << endl ;

    ofOut << "# ------------------ " << endl;
    ofOut << "# Report map section " << endl;
    ofOut << "#   " << endl << endl ;
    
    ofOut << "Report_Long_Rad_Down = 0 " << endl;
    ofOut << "Report_Short_Rad_Down = 0 " << endl;
    ofOut << "Report_Precip = 0 " << endl;
    ofOut << "Report_Rel_Humidity = 0 " << endl;
    ofOut << "Report_Wind_Speed = 0 " << endl;
    ofOut << "Report_AvgAir_Temperature = 0 " << endl;
    ofOut << "Report_MinAir_Temperature = 0 " << endl;
    ofOut << "Report_MaxAir_Temperature = 0 " << endl << endl;

    ofOut << "Report_SWE = 0 " << endl;
    ofOut << "Report_Snowmelt = 0 " << endl;
    ofOut << "Report_Infilt_Cap = 0 " << endl;
    ofOut << "Report_Streamflow = 0 " << endl;
    ofOut << "Report_Saturation_Area = 0 " << endl;
    ofOut << "Report_Ponding = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_Average = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_Up = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_L1 = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_L2 = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_L3 = 0 " << endl;
    ofOut << "Report_WaterTableDepth = 0 " << endl;

    ofOut << "# Time-constant variables (reported only once)" << endl;
    ofOut << "Report_RootZone_in_L1 = 0 " << endl;
    ofOut << "Report_RootZone_in_L2 = 0 " << endl;
    ofOut << "Report_RootZone_in_L3 = 0 " << endl;
    ofOut << "Report_Field_Capacity_L1 = 0 " << endl;
    ofOut << "Report_Field_Capacity_L2 = 0 " << endl;
    ofOut << "Report_Field_Capacity_L3 = 0 " << endl << endl ;
    ofOut << "# -------------------------------------------" << endl;
    ofOut << "Report_Soil_Sat_Deficit = 0 " << endl;
    ofOut << "Report_Ground_Water = 0 " << endl;
    ofOut << "Report_Surface_Net_Rad = 0 " << endl;
    ofOut << "Report_Vegetation_Net_Rad = 0 " << endl;
    ofOut << "Report_Total_Net_Rad = 0 " << endl;
    ofOut << "Report_Surface_Latent_Heat = 0 " << endl;
    ofOut << "Report_Vegetation_Latent_Heat = 0 " << endl;
    ofOut << "Report_Total_Latent_Heat = 0 " << endl;
    ofOut << "Report_Surface_Sens_Heat = 0 " << endl;
    ofOut << "Report_Vegetation_Sens_Heat = 0 " << endl;
    ofOut << "Report_Total_Sens_Heat = 0 " << endl;
    ofOut << "Report_Grnd_Heat = 0 " << endl;
    ofOut << "Report_Snow_Heat = 0 " << endl;
    ofOut << "Report_Soil_Temperature = 0 " << endl;
    ofOut << "Report_Skin_Temperature = 0 " << endl << endl;

    ofOut << "Report_Total_ET = 0 " << endl;
    ofOut << "Report_Transpiration_sum = 0 " << endl;
    ofOut << "Report_Transpiration_Layer1 = 0 " << endl;
    ofOut << "Report_Transpiration_Layer2 = 0 " << endl;
    ofOut << "Report_Transpiration_Layer3 = 0 " << endl;
    ofOut << "Report_Einterception_sum = 0 " << endl;
    ofOut << "Report_Esoil_sum = 0 " << endl;
    ofOut << "Report_Canopy_Water_Stor_sum = 0 " << endl << endl;

    ofOut << "Report_Veget_frac = 0 " << endl;
    ofOut << "Report_Stem_Density = 0 " << endl;
    ofOut << "Report_RootFracL1_species = 0 " << endl;
    ofOut << "Report_RootFracL2_species = 0 " << endl;
    ofOut << "Report_Leaf_Area_Index = 0 " << endl;
    ofOut << "Report_Stand_Age = 0 " << endl;
    ofOut << "Report_Canopy_Conductance = 0 " << endl;
    ofOut << "Report_GPP = 0 " << endl;
    ofOut << "Report_NPP = 0 " << endl;
    ofOut << "Report_Basal_Area = 0 " << endl;
    ofOut << "Report_Tree_Height = 0 " << endl;
    ofOut << "Report_Root_Mass = 0 " << endl;
    ofOut << "Report_Canopy_Temp = 0 " << endl;
    ofOut << "Report_Canopy_NetR = 0 " << endl;
    ofOut << "Report_Canopy_LE_E = 0 " << endl;
    ofOut << "Report_Canopy_LE_T = 0 " << endl;
    ofOut << "Report_Canopy_Sens_Heat = 0 " << endl;
    ofOut << "Report_Canopy_Water_Stor = 0 " << endl;
    ofOut << "Report_species_ET = 0 " << endl;
    ofOut << "Report_Transpiration = 0 " << endl;
    ofOut << "Report_Einterception = 0 " << endl;
    ofOut << "Report_Esoil = 0 " << endl << endl;

    ofOut << "Report_GW_to_Channel = 0 " << endl;
    ofOut << "Report_Surface_to_Channel = 0 " << endl;
    ofOut << "Report_Infiltration = 0" << endl ;
    ofOut << "Report_Return_Flow_Surface = 0" << endl ;
    ofOut << "Report_Percolation_to_Layer2 = 0" << endl ;
    ofOut << "Report_Return_Flow_to_Layer1 = 0" << endl ;
    ofOut << "Report_Percolation_to_Layer3 = 0" << endl ;
    ofOut << "Report_Groundwater_Recharge = 0" << endl ;
    ofOut << "Report_Leakage_Out_of_System = 0" << endl ;
    ofOut << "Report_Return_Flow_to_Layer2 = 0" << endl ;
	
    ofOut << "Report_Overland_Inflow = 0" << endl ;
    ofOut << "Report_Stream_Inflow = 0" << endl;
    ofOut << "Report_Groundwater_Inflow = 0 " << endl ;
    ofOut << "Report_Overland_Outflow = 0" << endl ;
    ofOut << "Report_Stream_Outflow = 0" << endl;
    ofOut << "Report_Groundwater_Outflow = 0" << endl ;

    ofOut << "Report_Infiltration_acc = 0" << endl ;
    ofOut << "Report_Return_Flow_Surface_acc = 0" << endl ;
    ofOut << "Report_Return_Flow_Surface_acc = 0" << endl ;
    ofOut << "Report_Percolation_to_Layer2_acc = 0" << endl ;
    ofOut << "Report_Return_Flow_to_Layer1_acc = 0" << endl ;
    ofOut << "Report_Percolation_to_Layer3_acc = 0" << endl ;
    ofOut << "Report_Groundwater_Recharge_acc = 0" << endl ;
    ofOut << "Report_Leakage_Out_of_System_acc = 0" << endl ;
    ofOut << "Report_Return_Flow_to_Layer2_acc = 0" << endl ;
    ofOut << "Report_Soil_Evaporation_acc = 0" << endl ;
    ofOut << "Report_Transpiration_Layer1_acc = 0" << endl ;
    ofOut << "Report_Transpiration_Layer2_acc = 0" << endl ;
    ofOut << "Report_Transpiration_Layer3_acc = 0" << endl ;
    ofOut << "Report_Overland_Inflow_acc = 0" << endl ;
    ofOut << "Report_Stream_Inflow_acc = 0" << endl ;
    ofOut << "Report_Groundwater_Inflow_acc = 0" << endl ;
    ofOut << "Report_Overland_Outflow_acc = 0" << endl ;
    ofOut << "Report_Stream_Outflow_acc = 0" << endl ;
    ofOut << "Report_Groundwater_Outflow_acc = 0" << endl << endl;
    ofOut << "Report_GW_to_Channel_acc = 0 " << endl;
    ofOut << "Report_Surface_to_Channel_acc = 0 " << endl;
    
    ofOut << "# ---------------------------------------------------- " << endl;
    ofOut << "# Report time section (at locations set by TS_mask map) " << endl;
    ofOut << "#   " << endl << endl ;
    ofOut << "TS_mask = Tsmask.map " << endl <<"#" << endl;
    ofOut << "Ts_OutletDischarge = 1 " << endl;
    ofOut << "Ts_Long_Rad_Down = 0 " << endl;
    ofOut << "Ts_Short_Rad_Down = 0 " << endl;
    ofOut << "Ts_Precip = 0 " << endl;
    ofOut << "Ts_Rel_Humidity = 0 " << endl;
    ofOut << "Ts_Wind_Speed = 0 " << endl;
    ofOut << "Ts_AvgAir_Temperature = 0 " << endl;
    ofOut << "Ts_MinAir_Temperature = 0 " << endl;
    ofOut << "Ts_MaxAir_Temperature = 0 " << endl;
    ofOut << "Ts_SWE = 1 " << endl;
    ofOut << "Ts_Snowmelt = 1 " << endl;
    ofOut << "Ts_Infilt_Cap = 0 " << endl;
    ofOut << "Ts_Streamflow = 0 " << endl;
    ofOut << "Ts_Ponding = 0 " << endl;
    ofOut << "Ts_Soil_Water_Content_Average = 1 " << endl;
    ofOut << "Ts_Soil_Water_Content_Up = 0 " << endl;
    ofOut << "Ts_Soil_Water_Content_L1 = 1 " << endl;
    ofOut << "Ts_Soil_Water_Content_L2 = 1 " << endl;
    ofOut << "Ts_Soil_Water_Content_L3 = 1 " << endl;
    ofOut << "Ts_WaterTableDepth = 0 " << endl;
    ofOut << "Ts_Field_Capacity_L1 = 0 " << endl;
    ofOut << "Ts_Field_Capacity_L2 = 0 " << endl;
    ofOut << "Ts_Field_Capacity_L3 = 0 " << endl;
    ofOut << "Ts_Soil_Sat_Deficit = 0 " << endl;
    ofOut << "Ts_Ground_Water = 0 " << endl;
    ofOut << "Ts_Surface_Net_Rad = 0 " << endl;
    ofOut << "Ts_Vegetation_Net_Rad = 0 " << endl;
    ofOut << "Ts_Total_Net_Rad = 0 " << endl;
    ofOut << "Ts_Surface_Latent_Heat = 0 " << endl;
    ofOut << "Ts_Vegetation_Latent_Heat = 0 " << endl;
    ofOut << "Ts_Total_Latent_Heat = 0 " << endl;
    ofOut << "Ts_Surface_Sens_Heat = 0 " << endl;
    ofOut << "Ts_Vegetation_Sens_Heat = 0 " << endl;
    ofOut << "Ts_Total_Sens_Heat = 0 " << endl;
    ofOut << "Ts_Grnd_Heat = 0 " << endl;
    ofOut << "Ts_Snow_Heat = 0 " << endl;
    ofOut << "Ts_Soil_Temperature = 0 " << endl;
    ofOut << "Ts_Skin_Temperature = 0 " << endl << endl ;

    ofOut << "Ts_Total_ET = 1 " << endl;
    ofOut << "Ts_Transpiration_sum = 1 " << endl;
    ofOut << "Ts_Transpiration_Layer1 = 0 " << endl;
    ofOut << "Ts_Transpiration_Layer2 = 0 " << endl;
    ofOut << "Ts_Transpiration_Layer3 = 0 " << endl;
    ofOut << "Ts_Einterception_sum = 1 " << endl;
    ofOut << "Ts_Esoil_sum = 1 " << endl;
    ofOut << "Ts_Canopy_Water_Stor_sum = 0 " << endl << endl;

    ofOut << "Ts_Veget_frac = 0 " << endl;
    ofOut << "Ts_Stem_Density = 0 " << endl;
    ofOut << "Ts_RootFracL1_species = 0 " << endl;
    ofOut << "Ts_RootFracL2_species = 0 " << endl;
    ofOut << "Ts_Leaf_Area_Index = 1 " << endl;
    ofOut << "Ts_Stand_Age = 0 " << endl;
    ofOut << "Ts_Canopy_Conductance = 1 " << endl;
    ofOut << "Ts_GPP = 0 " << endl;
    ofOut << "Ts_NPP = 1 " << endl;
    ofOut << "Ts_Basal_Area = 0 " << endl;
    ofOut << "Ts_Tree_Height = 0 " << endl;
    ofOut << "Ts_Root_Mass = 0 " << endl;
    ofOut << "Ts_Canopy_Temp = 0 " << endl;
    ofOut << "Ts_Canopy_NetR = 0 " << endl;
    ofOut << "Ts_Canopy_LE_E = 0 " << endl;
    ofOut << "Ts_Canopy_LE_T = 0 " << endl;
    ofOut << "Ts_Canopy_Sens_Heat = 0 " << endl;
    ofOut << "Ts_Canopy_Water_Stor = 0 " << endl;
    ofOut << "Ts_species_ET = 0 " << endl;
    ofOut << "Ts_Transpiration = 0 " << endl;
    ofOut << "Ts_Einterception = 0 " << endl;
    ofOut << "Ts_Esoil = 0 " << endl << endl;

    ofOut << "Ts_GW_to_Channel = 0 " << endl;
    ofOut << "Ts_Surface_to_Channel = 0 " << endl;
    ofOut << "Ts_Infiltration = 0" << endl ;
    ofOut << "Ts_Return_Flow_Surface = 0" << endl ;
    ofOut << "Ts_Percolation_to_Layer2 = 0" << endl ;
    ofOut << "Ts_Return_Flow_to_Layer1 = 0" << endl ;
    ofOut << "Ts_Percolation_to_Layer3 = 0" << endl ;
    ofOut << "Ts_Groundwater_Recharge = 0" << endl ;
    ofOut << "Ts_Leakage_Out_of_System = 0" << endl ;
    ofOut << "Ts_Return_Flow_to_Layer2 = 0" << endl ;
    ofOut << "Ts_Overland_Inflow = 0" << endl ;
    ofOut << "Ts_Stream_Inflow = 0" << endl;
    ofOut << "Ts_Groundwater_Inflow = 0 " << endl ;
    ofOut << "Ts_Overland_Outflow = 0" << endl ;
    ofOut << "Ts_Stream_Outflow = 0" << endl;
    ofOut << "Ts_Groundwater_Outflow = 0" << endl ;

    if (ofOut)
      ofOut.close();
  }
  catch(const std::exception &e){
    cout << "Failure writing configuration template file with  " << e.what() << endl;
    exit(EXIT_FAILURE);
  }
}

