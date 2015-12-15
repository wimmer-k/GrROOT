.EXPORT_ALL_VARIABLES:

.PHONY: clean all docs

ifndef GRROOT_DIR
GRROOT_DIR := $(PWD)/..
endif

ifdef GRROOT_LIBDIR
LIB_DIR := $(GRROOT_LIBDIR)
else
LIB_DIR := $(GRROOT_DIR)/lib
endif

ifdef GRROOT_BINDIR
BIN_DIR := $(GRROOT_BINDIR)
else
BIN_DIR := $(GRROOT_DIR)/bin
endif

DEBUG = -g
SWITCH = -DS800_DETAILEDTREE -USPECTCL_MODE

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --glibs) -lSpectrum
INCLUDES      = -I./

CPP             = g++
CFLAGS		= -pedantic -Wall -Wno-long-long -O3 -fPIC -D_FILE_OFFSET_BITS=64 -MMD \
                   $(DEBUG) $(ROOTCFLAGS) $(INCLUDES) $(SWITCH)

CFLAGS += -Wl,--no-as-needed

LIBS  =  $(ROOTLIBS) -L$(LIB_DIR) -lCommandLineInterface -lGretina -lS800 -lScaler -lRunInfo -lFit -lSettings -lPeaks

O_FILES = UnpackedEvent.o \
	Calibration.o \
	Simulation.o \
	RawHistograms.o \
	CalHistograms.o \
	SimHistograms.o \
	Tracking.o

TO_FILES = Calibration.o Tracking.o

all: libs HFC AnalysisProgs SimulationProgs CalibrationProgs MiscProgs Track

libs: libSettings.so libGretina.so libS800.so libScaler.so libRunInfo.so libCommandLineInterface.so libPeaks.so libFit.so

AnalysisProgs: GrROOT Calculate Histos Cal_histos ScalerAnalysis TreeSplitter MakeMode2

SimulationProgs:  SimCalculate Sim_histos 

CalibrationProgs: GammaCal ICCal PadCal

MiscProgs: RawEventLoop CalEventLoop

define EXE_COMMANDS
@echo "Compiling $@"
@$(CPP) $(CFLAGS) $(LIBS) $< $(filter %.o,$^) -o $@ && cp $@ $(BIN_DIR)
endef

GrROOT: GrROOT.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

HFC:
	@cd hfc; make
	@echo "GEB_HFC compiled"

Calculate: Calculate.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

SimCalculate: SimCalculate.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

TreeSplitter: TreeSplitter.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

Histos: Histos.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

Cal_histos: Cal_histos.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

Sim_histos: Sim_histos.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

Track: Track.cc $(TO_FILES) | libs
	$(EXE_COMMANDS)

ScalerAnalysis: ScalerAnalysis.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

MakeMode2: MakeMode2.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

RawEventLoop: RawEventLoop.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

CalEventLoop: CalEventLoop.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

GammaCal: GammaCal.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

ICCal: ICCal.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

PadCal: PadCal.cc $(O_FILES) | libs
	$(EXE_COMMANDS)

define LIB_COMMANDS
@echo "Making $@"
@$(CPP) $(CFLAGS) -fPIC -shared -Wl,-soname,$@ -o $@ $^ -lc
@cp $@ $(LIB_DIR)
endef

libGretina.so: Gretina.o GretinaDictionary.o GretinaCalc.o GretinaCalcDictionary.o GretinaTrack.o GretinaTrackDictionary.o GretinaSim.o GretinaSimDictionary.o Mode3Calc.o Mode3CalcDictionary.o Trace.o TraceDictionary.o
	$(LIB_COMMANDS)

libS800.so: S800.o S800Dictionary.o S800Calc.o S800CalcDictionary.o
	$(LIB_COMMANDS)

libScaler.so: Scaler.o ScalerDictionary.o
	$(LIB_COMMANDS)

libSettings.so: Settings.o SettingsDictionary.o SimSettings.o SimSettingsDictionary.o  TrackSettings.o TrackSettingsDictionary.o 
	$(LIB_COMMANDS)

libRunInfo.so: RunInfo.o RunInfoDictionary.o
	$(LIB_COMMANDS)

lib%.so: %.o
	$(LIB_COMMANDS)

docs: *.cc *.hh
	doxygen doxy-config

-include $(wildcard *.d)

%.o: %.cc %.hh
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) -c $< -o $@

CommandLineInterface.o: CommandLineInterface.cc CommandLineInterface.hh | libFit.so
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) -c $< -o $@

%Dictionary.o: %Dictionary.cc %Dictionary.h
	@echo "Compiling $@ as a Dictionary"
	@$(CPP) -p -fPIC $(CFLAGS) -c $<

%Dictionary.cc %Dictionary.h: %.hh %LinkDef.h
	@echo Building $@
	@rootcint -f $*Dictionary.cc -c $(SWITCH) $*.hh $*LinkDef.h

clean:
	@echo Cleaning up
	@rm -f $(BIN_DIR)/* $(LIB_DIR)/* *.o *Dictionary.cc *Dictionary.h lib*.so *.d
	@rm -f Calculate CalEventLoop Cal_histos GammaCal GrROOT Histos ICCal PadCal RawEventLoop ScalerAnalysis SimCalculate Sim_histos TreeSplitter MakeMode2 Track
	@cd hfc; make clean
