################################################################################
# Automatically-generated file. Do not edit!
################################################################################
CXX := gcc
NVCC := nvcc
CFLAGS := -static -std=c++11 -Wall -Wextra
NVCCFLAGS := 
IFLAGS = -I/afs/crc.nd.edu/x86_64_linux/c/cuda/10.0/include/

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../IncrementExternalForceOnDevice.cu \
../CalculateEquilibrium.cu \
../TorsionSolveOnDevice.cu \
../WLCSolveOnDevice.cu \
../DPDParticle.cu \
../AdvancePositionOnDevice.cu \
../LinkNodesOnDevice.cu \
../BucketSchemeOnDevice.cu \
../NodeSystemDevice.cu \
../NodeSystemBuilder.cpp \
../ForceDiagramStorage.cu \
../main.cpp 


# this is a variable
OBJS += \
./IncrementExternalForceOnDevice.o \
./CalculateEquilibrium.o \
./TorsionSolveOnDevice.o \
./WLCSolveOnDevice.o \
./DPDParticle.o \
./AdvancePositionOnDevice.o \
./LinkNodesOnDevice.o \
./BucketSchemeOnDevice.o \
./NodeSystemDevice.o \
./NodeSystemBuilder.o \
./ForceDiagramStorage.o \
./main.o 

 
CPP_DEPS += \
./IncrementExternalForceOnDevice.d \
./CalculateEquilibrium.d \
./TorsionSolveOnDevice.d \
./WLCSolveOnDevice.d \
./DPDParticle.d \
./AdvancePositionOnDevice.d \
./LinkNodesOnDevice.d \
./BucketSchemeOnDevice.d \
./NodeSystemDevice.d \
./NodeSystemBuilder.d \
./ForceDiagramStorage.d \
./main.d 


#cpp files
%.o : ./%.cpp 
	$(CXX) $(CFLAGS) $(IFLAGS) -o  $@ -c $^

	
#cuda files
%.o : ./%.cu 
	$(NVCC) $(NVCCFLAGS) -dc -o $@ $^ 



