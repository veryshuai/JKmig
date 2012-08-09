################################################################################
# MAKEFILE FOR DYN 
#
################################################################################

INC_PATH := -I$$HOME/libs/eigen-eigen-ca142d0540d3 
#LIB_PATH := -L$$HOME/work/lib -L$$HOME/work/src/casadi/build/lib

# All Target
all: DavidsTarget

nothing:
	@echo ' ' 
	@echo $(notdir $(CPP_SRCS:.cpp=.d))
	@echo ' ' 

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS = \
../dyn.cpp \
../mig.cpp

CPP_DEPS = $(notdir $(CPP_SRCS:.cpp=.d))
OBJS = $(notdir $(CPP_SRCS:.cpp=.o))

# Each subdirectory must supply rules for building sources it contributes
%.o:    ../%.cpp
	@echo 'Building file: $<'
	g++ -O0 -g3 -Wall -Wno-sign-compare -c -fmessage-length=0 $(INC_PATH) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo ' '

#LIBS := -lknitro -ldl -lxerces-c -lcasadi -larprec

# Tool invocations
DavidsTarget: $(OBJS)
	@echo 'Building target: $@'
	g++  -o"dyn" $(LIB_PATH) $(OBJS) $(LIBS)
	@echo ' '

.PHONY: all 

