################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../ech2o.cpp 

OBJS += \
./ech2o.o 

CPP_DEPS += \
./ech2o.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking2: Clang C++ Compiler'
	$(CXX) -DCPU_LITTLE_ENDIAN -I"../includes" -I"/usr/local/Cellar/boost/include" -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<" 
	@echo 'Finished building: $<'
	@echo ' '


