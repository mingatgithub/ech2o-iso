################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Hydro/AdvanceLAIMaps.cpp \
../Hydro/CalcCatchArea.cpp \
../Hydro/CalcFracMobileWater.cpp \
../Hydro/CalcInitialStreamStorage.cpp \
../Hydro/CalcSoilResist.cpp \
../Hydro/CalcPropLayers.cpp \
../Hydro/CalcRootDistrib.cpp \
../Hydro/CalcTPDMoisture.cpp \
../Hydro/CalculateForestGrowth.cpp \
../Hydro/CalculateGroundwater.cpp \
../Hydro/CalculateSatArea.cpp \
../Hydro/CanopyInterception.cpp \
../Hydro/GWrouting.cpp \
../Hydro/Infilt_GreenAmpt.cpp \
../Hydro/KinematicWave.cpp \
../Hydro/SnowOutputPhase.cpp \
../Hydro/SoilEvapotranspiration.cpp \
../Hydro/SoilWaterRedistribution.cpp \
../Hydro/SolveCanopyFluxes.cpp \
../Hydro/SolveSurfaceEnergyBalance.cpp \
../Hydro/UpdateSnowPack.cpp 

OBJS += \
./Hydro/AdvanceLAIMaps.o \
./Hydro/CalcCatchArea.o \
./Hydro/CalcFracMobileWater.o \
./Hydro/CalcInitialStreamStorage.o \
./Hydro/CalcPropLayers.o \
./Hydro/CalcRootDistrib.o \
./Hydro/CalcSoilResist.o \
./Hydro/CalcTPDMoisture.o \
./Hydro/CalculateForestGrowth.o \
./Hydro/CalculateGroundwater.o \
./Hydro/CalculateSatArea.o \
./Hydro/CanopyInterception.o \
./Hydro/GWrouting.o \
./Hydro/Infilt_GreenAmpt.o \
./Hydro/KinematicWave.o \
./Hydro/SnowOutputPhase.o \
./Hydro/SoilEvapotranspiration.o \
./Hydro/SoilWaterRedistribution.o \
./Hydro/SolveCanopyFluxes.o \
./Hydro/SolveSurfaceEnergyBalance.o \
./Hydro/UpdateSnowPack.o 

CPP_DEPS += \
./Hydro/AdvanceLAIMaps.d \
./Hydro/CalcCatchArea.d \
./Hydro/CalcFracMobileWater.d \
./Hydro/CalcInitialStreamStorage.d \
./Hydro/CalcSoilResist.d \
./Hydro/CalcRootDistrib.d \
./Hydro/CalcPropLayers.d \
./Hydro/CalcTPDMoisture.d \
./Hydro/CalculateForestGrowth.d \
./Hydro/CalculateSatArea.d \
./Hydro/CalculateGroundwater.d \
./Hydro/CanopyInterception.d \
./Hydro/GWrouting.d \
./Hydro/Infilt_GreenAmpt.d \
./Hydro/KinematicWave.d \
./Hydro/SnowOutputPhase.d \
./Hydro/SoilEvapotranspiration.d \
./Hydro/SoilWaterRedistribution.d \
./Hydro/SolveCanopyFluxes.d \
./Hydro/SolveSurfaceEnergyBalance.d \
./Hydro/UpdateSnowPack.d 


# Each subdirectory must supply rules for building sources it contributes
Hydro/%.o: ../Hydro/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) -DCPU_LITTLE_ENDIAN -I"../includes" -I"/usr/local/Cellar/boost/include" -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


