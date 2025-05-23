HOME    :=/work/d448/d448/yichongchen/Proposed_Knudsen/Rarefied/L100
SRC     :=${HOME}/postProcessing
BIN     :=${HOME}/bin
DataDir :=${HOME}/DataDir

# CXX      := g++
# CXXFLAGS := -O2 -Wall $(INCL)

# .SUFFIXES: .cpp
# # It informs make that you will be using this special suffixes to make your own rules.
# # Implicit rule for the construction of .o (object) files out of .cpp (source files)
# %.o : ${SRC}/%.cpp
# $(CXX) $(CXXFLAGS) -c -o $@ $<
# # The -c flag says to generate the objective file. 
# # The -o $@ says to put the output of the compilation in the file named on the left side of the :	

TARGETS := pp_Knudsen_ACs \
           pp_Knudsen_Bins \
		   pp_Knudsen_doubleBins \
		   pp_Knudsen_GasOverTime \
		   pp_Knudsen_VDFPr \
		   pp_Knudsen_xBins \
		   pp_reservoir_Bins \
		   pp_vaccum_Bins

# Pattern rule for compiling the codes
$(TARGETS): % : ${SRC}/%.cpp
	${CXX} -o $@ $^
	./$@
	mv $@ ${BIN}

# Pattern rule for building the objec files
#${SRC}/%.o: ${SRC}/%.cpp
#	${CXX} $(CXXFLAGS) -c -o $@ $<		

# ------------------------------------------------------------
# Phony target, that it does not build anything
.PHONY :PoiseData Penetration AccommodationCoefficients ScatteringData collect clean q grep help

## PoiseData	: Compile the source codes for Poiseuille Flow simulation.
PoiseData: pp_meas_Bins pp_meas_GasOverTime
Penetration: pp_meas_penACs pp_meas_penBeam pp_meas_penOverall
AccommodationCoefficients: pp_meas_ACs pp_meas_Collisions pp_meas_penACs
VDF: pp_meas_VDFF pp_meas_VDFPi pp_meas_VDFPr
ScatteringData: $(TARGETS)

## collect	: Collect the auto-generated files
collect :
	mv -f *.txt ${DataDir}
	mv -f *.dat ${DataDir}
	mv ${DataDir}/Specification.dat ./
	mv ${DataDir}/data.dat ./

## clean	: Remove the auto-generated files
clean :
	$(RM) *.txt *.dat *.o

## q		: Monitoring jobs
q : 
	squeue -u $(USER)

# Target to extract the timesteps in the dumpfile
grep :
	grep -o 'TIMESTEP' dump_meas_gas.lammpstrj | wc -l

# New target to update Specification.dat
update :
	@echo "nTimeStep = $$(grep -o 'TIMESTEP' dump_meas_gas.lammpstrj | wc -l)" && \
	nTimeSteps=$$(grep -o 'TIMESTEP' dump_meas_gas.lammpstrj | wc -l) && \
	sed -i "1s/.*/$$nTimeSteps                      # nTimeSteps/" Specification.dat


# sed stands for 'stream editor'. sed reads in some text, does some filtering, and writes out the filtered text.
help : Makefile
	@sed -n 's/^##//p' $<
