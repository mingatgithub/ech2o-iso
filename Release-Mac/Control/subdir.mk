################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Control/GenerateConfigTemplate.cpp \
../Control/GenerateConfigTrckTemplate.cpp \
../Control/ReadConfigFile.cpp \
../Control/ReadConfigTrck.cpp 

OBJS += \
./Control/GenerateConfigTemplate.o \
./Control/GenerateConfigTrckTemplate.o \
./Control/ReadConfigFile.o \
./Control/ReadConfigTrck.o 

CPP_DEPS += \
./Control/GenerateConfigTemplate.d \
./Control/GenerateConfigTrckTemplate.d \
./Control/ReadConfigFile.d \
./Control/ReadConfigTrck.d


# Each subdirectory must supply rules for building sources it contributes
Control/%.o: ../Control/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) -DCPU_LITTLE_ENDIAN -I"../includes" -I"/usr/local/Cellar/boost/include" -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


