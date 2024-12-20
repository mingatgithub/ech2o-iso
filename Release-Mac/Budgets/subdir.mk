################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Budgets/AccountFluxes.cpp \
../Budgets/AccountRelArea.cpp \
../Budgets/AccountStorages.cpp \
../Budgets/MassBalanceError.cpp \
../Budgets/TotalOutputs.cpp \
../Budgets/TotalEvaporationS.cpp \
../Budgets/TotalEvaporationI.cpp \
../Budgets/TotalGWtoChn.cpp \
../Budgets/TotalTranspiration.cpp \
../Budgets/TotalLeakage.cpp \
../Budgets/TotalOvlndFlow.cpp \
../Budgets/TotalPrecipitation.cpp \
../Budgets/TotalRecharge.cpp \
../Budgets/TotalSaturationArea.cpp \
../Budgets/TotalSrftoChn.cpp \
../Budgets/TotalStorage.cpp \
../Budgets/totalGrndFlow.cpp 

OBJS += \
./Budgets/AccountFluxes.o \
./Budgets/AccountRelArea.o \
./Budgets/AccountStorages.o \
./Budgets/MassBalanceError.o \
./Budgets/TotalOutputs.o \
./Budgets/TotalEvaporationS.o \
./Budgets/TotalEvaporationI.o \
./Budgets/TotalGWtoChn.o \
./Budgets/TotalTranspiration.o \
./Budgets/TotalLeakage.o \
./Budgets/TotalOvlndFlow.o \
./Budgets/TotalPrecipitation.o \
./Budgets/TotalRecharge.o \
./Budgets/TotalSaturationArea.o \
./Budgets/TotalSrftoChn.o \
./Budgets/TotalStorage.o \
./Budgets/totalGrndFlow.o 

CPP_DEPS += \
./Budgets/AccountFluxes.d \
./Budgets/AccountRelArea.d \
./Budgets/AccountStorages.d \
./Budgets/MassBalanceError.d \
./Budgets/TotalOutputs.d \
./Budgets/TotalEvaporationS.d \
./Budgets/TotalEvaporationI.d \
./Budgets/TotalGWtoChn.d \
./Budgets/TotalTranspiration.d \
./Budgets/TotalLeakage.d \
./Budgets/TotalOvlndFlow.d \
./Budgets/TotalPrecipitation.d \
./Budgets/TotalRecharge.d \
./Budgets/TotalSaturationArea.d \
./Budgets/TotalSrftoChn.d \
./Budgets/TotalStorage.d \
./Budgets/totalGrndFlow.d 


# Each subdirectory must supply rules for building sources it contributes
Budgets/%.o: ../Budgets/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) -DCPU_LITTLE_ENDIAN -I"../includes" -I"/usr/local/Cellar/boost/include" -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


