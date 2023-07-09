################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CoordTransform.cpp \
../src/Mesh.cpp \
../src/MeshDeformationTool.cpp \
../src/MeshQuality.cpp \
../src/ReadConfigFile.cpp \
../src/SPDS.cpp \
../src/WriteResults.cpp \
../src/getNodeType.cpp \
../src/greedy.cpp \
../src/rbfGenFunc.cpp \
../src/rbfds.cpp \
../src/rbfps.cpp \
../src/rbfstd.cpp 

CPP_DEPS += \
./src/CoordTransform.d \
./src/Mesh.d \
./src/MeshDeformationTool.d \
./src/MeshQuality.d \
./src/ReadConfigFile.d \
./src/SPDS.d \
./src/WriteResults.d \
./src/getNodeType.d \
./src/greedy.d \
./src/rbfGenFunc.d \
./src/rbfds.d \
./src/rbfps.d \
./src/rbfstd.d 

OBJS += \
./src/CoordTransform.o \
./src/Mesh.o \
./src/MeshDeformationTool.o \
./src/MeshQuality.o \
./src/ReadConfigFile.o \
./src/SPDS.o \
./src/WriteResults.o \
./src/getNodeType.o \
./src/greedy.o \
./src/rbfGenFunc.o \
./src/rbfds.o \
./src/rbfps.o \
./src/rbfstd.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"C:\Users\floyd\OneDrive\Documenten\Thesis" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src

clean-src:
	-$(RM) ./src/CoordTransform.d ./src/CoordTransform.o ./src/Mesh.d ./src/Mesh.o ./src/MeshDeformationTool.d ./src/MeshDeformationTool.o ./src/MeshQuality.d ./src/MeshQuality.o ./src/ReadConfigFile.d ./src/ReadConfigFile.o ./src/SPDS.d ./src/SPDS.o ./src/WriteResults.d ./src/WriteResults.o ./src/getNodeType.d ./src/getNodeType.o ./src/greedy.d ./src/greedy.o ./src/rbfGenFunc.d ./src/rbfGenFunc.o ./src/rbfds.d ./src/rbfds.o ./src/rbfps.d ./src/rbfps.o ./src/rbfstd.d ./src/rbfstd.o

.PHONY: clean-src

