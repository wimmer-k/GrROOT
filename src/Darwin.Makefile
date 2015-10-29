.EXPORT_ALL_VARIABLES:

.PHONY: clean all

#GRROOT_DIR := $(shell if [ -z $(GRROOT_DIR) ]; then echo ".."; else echo $(GRROOT_DIR); fi)
GRROOT_DIR := $(shell if [ -z $(GRROOT_DIR) ]; then echo $(PWD)"/.."; else echo $(GRROOT_DIR); fi)

LIB_DIR = $(shell if [ -z $(GRROOT_LIBDIR) ]; then echo "$(GRROOT_DIR)/lib/lib32"; else echo $(GRROOT_LIBDIR); fi)
BIN_DIR = $(shell if [ -z $(GRROOT_BINDIR) ]; then echo "$(GRROOT_DIR)/bin/bin32"; else echo $(GRROOT_BINDIR); fi)

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTINC      := -I$(shell root-config --incdir)

SPECTRUM     := $(shell root-config --libdir)/libSpectrum.so

CPP             = g++
CFLAGS		= -pedantic -Wall -Wno-long-long -g -O3 $(ROOTCFLAGS) -fPIC
DFLAGS		= -Wall -Wno-long-long -g -O3 $(ROOTCFLAGS) -fPIC
#CFLAGS		= -pedantic -Wall -Wno-long-long -g -O3 $(ROOTCFLAGS) -fPIC -m32
#DFLAGS		= -Wall -Wno-long-long -g -O3 $(ROOTCFLAGS) -fPIC -m32

BASELIBS = -lm $(ROOTLIBS) $(ROOTGLIBS) -L$(LIB_DIR) -lSpectrum
ALLIBS  =  $(BASELIBS) -lCommandLineInterface -lGretina -lS800 -lScaler -lRunInfo -lFit -lSettings -lPeaks
LIBS 		= $(ALLIBS)

INCLUDES        = -I./
LFLAGS		= -g -fPIC
#LFLAGS		= -g -fPIC -m32
SOFLAGS	        = -dynamiclib -single_module -undefined dynamic_lookup

SWITCH = -DS800_DETAILEDTREE -USPECTCL_MODE

CFLAGS += $(SWITCH)
LFLAGS += $(SWITCH)
DFLAGS += $(SWITCH)
CFLAGS += -Wl,--no-as-needed
LFLAGS += -Wl,--no-as-needed
DFLAGS += -Wl,--no-as-needed



O_FILES = UnpackedEvent.o \
	Calibration.o \
	Settings.o \
	Gretina.o \
	RunInfo.o \
	RawHistograms.o \
	CalHistograms.o

RO_FILES = S800.o \
	S800Calc.o \
	Gretina.o \
	GretinaCalc.o \
	RunInfo.o \
	Trace.o \
	Mode3Calc.o \
	Scaler.o

TO_FILES = GretinaCalc.o \
	GretinaTrack.o \
	Calibration.o \
	TrackSettings.o \
	Tracking.o

all: libs HFC Analysis Calibration Misc
	echo Done

libs: Settings Gretina S800 Scaler RunInfo libCommandLineInterface.so libPeaks.so
	echo Done

Analysis: GrROOT Calculate Histos Cal_histos ScalerAnalysis
	echo Done

Calibration: GammaCal ICCal PadCal
	echo Done

Misc: RawEventLoop CalEventLoop
	echo Done

GrROOT: GrROOT.cc $(O_FILES) $(RO_FILES) libs
	$(CPP) $(CFLAGS) $(INCLUDES) GrROOT.cc $(LIBS) $(O_FILES) -o $(BIN_DIR)/$@

HFC: 
	cd hfc; make

Calculate: Calculate.cc $(O_FILES) $(RO_FILES) Gretinadefs.h
	$(CPP) $(CFLAGS) $(INCLUDES) Calculate.cc $(LIBS) $(O_FILES) -o $(BIN_DIR)/$@

TreeSplitter: TreeSplitter.cc $(RO_FILES)
	$(CPP) $(CFLAGS) $(INCLUDES) TreeSplitter.cc $(LIBS) $(RO_FILES) -o $(BIN_DIR)/$@

Histos: Histos.cc $(O_FILES) $(RO_FILES)
	$(CPP) $(CFLAGS) $(INCLUDES) Histos.cc $(LIBS) $(O_FILES) -o $(BIN_DIR)/$@

Cal_histos: Cal_histos.cc $(O_FILES) $(RO_FILES)
	$(CPP) $(CFLAGS) $(INCLUDES) Cal_histos.cc $(LIBS) $(O_FILES) -o $(BIN_DIR)/$@

Track: Track.cc $(TO_FILES)
	$(CPP) $(CFLAGS) $(INCLUDES) Track.cc $(LIBS) $(TO_FILES) -o $(BIN_DIR)/$@

ScalerAnalysis: ScalerAnalysis.cc
	$(CPP) $(CFLAGS) $(INCLUDES) $^ $(LIBS) $(O_FILES) $(SPECTRUM) -o $(BIN_DIR)/$@

RawEventLoop: RawEventLoop.cc $(O_FILES) $(RO_FILES)
	$(CPP) $(CFLAGS) $(INCLUDES) RawEventLoop.cc $(LIBS) $(O_FILES) -o $(BIN_DIR)/$@

CalEventLoop: CalEventLoop.cc $(O_FILES) $(RO_FILES)
	$(CPP) $(CFLAGS) $(INCLUDES) CalEventLoop.cc $(LIBS) $(O_FILES) -o $(BIN_DIR)/$@

GammaCal: GammaCal.cc libCommandLineInterface.so libPeaks.so
	$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(O_FILES) $(SPECTRUM) -o $(BIN_DIR)/$@

ICCal: ICCal.cc libCommandLineInterface.so libPeaks.so
	$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(O_FILES) $(SPECTRUM) -o $(BIN_DIR)/$@

PadCal: PadCal.cc libCommandLineInterface.so libPeaks.so
	$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(O_FILES) $(SPECTRUM) -o $(BIN_DIR)/$@


Gretina: Gretina.o GretinaDictionary.o GretinaCalc.o GretinaCalcDictionary.o GretinaTrack.o GretinaTrackDictionary.o Mode3Calc.o Mode3CalcDictionary.o Trace.o TraceDictionary.o
	$(CPP) $(LFLAGS) $(SOFLAGS) -install_name lib$@.so -o $(LIB_DIR)/lib$@.so $^ -lc

S800: S800.o S800Dictionary.o S800Calc.o S800CalcDictionary.o
	$(CPP) $(LFLAGS) $(SOFLAGS) -install_name lib$@.so -o $(LIB_DIR)/lib$@.so $^ -lc

Scaler: Scaler.o ScalerDictionary.o
	$(CPP) $(LFLAGS) $(SOFLAGS) -install_name lib$@.so -o $(LIB_DIR)/lib$@.so $^ -lc

Settings: Settings.o SettingsDictionary.o
	$(CPP) $(LFLAGS) $(SOFLAGS) -install_name lib$@.so -o $(LIB_DIR)/lib$@.so $^ -lc

TrackSettings.o: TrackSettings.cc TrackSettings.hh Settings.hh
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

RunInfo: RunInfo.o RunInfoDictionary.o
	$(CPP) $(LFLAGS) $(SOFLAGS) -install_name lib$@.so -o $(LIB_DIR)/lib$@.so $^ -lc


docs: *.cc *.hh
	doxygen doxy-config

%.o: %.cc %.hh
	@echo Default .o rule
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

lib%.so: %.o
	$(CPP) $(LFLAGS) $(SOFLAGS) -install_name $@ -o $(LIB_DIR)/$@ $^ $(BASELIBS) -lc

GretinaCalc.o: GretinaCalc.cc GretinaCalc.hh\
		Gretina.hh Gretinadefs.h
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

Fit.o: Fit.cc Fit.hh
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

Peaks.o: Peaks.cc Peaks.hh
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

CommandLineInterface.o: CommandLineInterface.cc CommandLineInterface.hh libFit.so
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

RawHistograms.o: RawHistograms.cc RawHistograms.hh\
	          S800.hh Gretina.hh Trace.hh Settings.hh
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

CalHistograms.o: CalHistograms.cc CalHistograms.hh\
                  CommandLineInterface.hh S800Calc.hh Gretina.hh Mode3Calc.hh Scaler.hh Settings.hh
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

UnpackedEvent.o: UnpackedEvent.cc UnpackedEvent.hh\
		  Calibration.hh Gretina.hh Gretinadefs.h Trace.hh\
	  	  S800.hh S800defs.h Scalerdefs.h Scaler.hh S800Calc.hh Mode3Calc.hh\
	          Settings.hh RawHistograms.hh CalHistograms.hh
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

Calibration.o: Calibration.cc Calibration.hh Gretina.hh S800.hh Trace.hh Gretinadefs.h Settings.hh S800Calc.hh Mode3Calc.hh Fit.hh
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

GretinaDictionary.o: GretinaDictionary.cc GretinaDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

GretinaDictionary.cc: Gretina.hh GretinaLinkDef.h 
	 rm -f GretinaDictionary.cc GretinaDictionary.h; rootcint -f $@ -c $(SWITCH) Gretina.hh GretinaLinkDef.h 

TraceDictionary.o: TraceDictionary.cc TraceDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

TraceDictionary.cc: Trace.hh TraceLinkDef.h 
	 rm -f TraceDictionary.cc TraceDictionary.h; rootcint -f $@ -c $(SWITCH) Trace.hh TraceLinkDef.h 

S800Dictionary.o: S800Dictionary.cc S800Dictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

S800Dictionary.cc: S800.hh S800LinkDef.h 
	 rm -f S800Dictionary.cc S800Dictionary.h; rootcint -f $@ -c $(SWITCH) S800.hh S800LinkDef.h 

ScalerDictionary.o: ScalerDictionary.cc ScalerDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

ScalerDictionary.cc: Scaler.hh ScalerLinkDef.h 
	 rm -f ScalerDictionary.cc ScalerDictionary.h; rootcint -f $@ -c $(SWITCH) Scaler.hh ScalerLinkDef.h 

SettingsDictionary.o: SettingsDictionary.cc SettingsDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

SettingsDictionary.cc: Settings.hh SettingsLinkDef.h 
	 rm -f SettingsDictionary.cc SettingsDictionary.h; rootcint -f $@ -c $(SWITCH) Settings.hh SettingsLinkDef.h 

TrackSettingsDictionary.o: TrackSettingsDictionary.cc TrackSettingsDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

TrackSettingsDictionary.cc: TrackSettings.hh TrackSettingsLinkDef.h 
	 rm -f TrackSettingsDictionary.cc TrackSettingsDictionary.h; rootcint -f $@ -c $(SWITCH) TrackSettings.hh TrackSettingsLinkDef.h 

RunInfoDictionary.o: RunInfoDictionary.cc RunInfoDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

RunInfoDictionary.cc: RunInfo.hh RunInfoLinkDef.h 
	 rm -f RunInfoDictionary.cc RunInfoDictionary.h; rootcint -f $@ -c $(SWITCH) RunInfo.hh RunInfoLinkDef.h 

S800CalcDictionary.o: S800CalcDictionary.cc S800CalcDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

S800CalcDictionary.cc: S800Calc.hh S800CalcLinkDef.h 
	 rm -f S800CalcDictionary.cc S800CalcDictionary.h; rootcint -f $@ -c $(SWITCH) S800Calc.hh S800CalcLinkDef.h 

FitDictionary.o: FitDictionary.cc FitDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

FitDictionary.cc: Fit.hh FitLinkDef.h 
	 rm -f FitDictionary.cc FitDictionary.h; rootcint -f $@ -c $(SWITCH) Fit.hh FitLinkDef.h 

PeaksDictionary.o: PeaksDictionary.cc PeaksDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

PeaksDictionary.cc: Peaks.hh PeaksLinkDef.h 
	 rm -f PeaksDictionary.cc PeaksDictionary.h; rootcint -f $@ -c $(SWITCH) Peaks.hh PeaksLinkDef.h 

Mode3CalcDictionary.o: Mode3CalcDictionary.cc Mode3CalcDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

Mode3CalcDictionary.cc: Mode3Calc.hh Mode3CalcLinkDef.h 
	 rm -f Mode3CalcDictionary.cc Mode3CalcDictionary.h; rootcint -f $@ -c $(SWITCH) Mode3Calc.hh Mode3CalcLinkDef.h 

GretinaCalcDictionary.o: GretinaCalcDictionary.cc GretinaCalcDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

GretinaCalcDictionary.cc: GretinaCalc.hh GretinaCalcLinkDef.h 
	 rm -f GretinaCalcDictionary.cc GretinaCalcDictionary.h; rootcint -f $@ -c $(SWITCH) GretinaCalc.hh GretinaCalcLinkDef.h 

GretinaTrackDictionary.o: GretinaTrackDictionary.cc GretinaTrackDictionary.h
	 $(CPP) -p -fPIC $(DFLAGS) -c $<

GretinaTrackDictionary.cc: GretinaTrack.hh GretinaTrackLinkDef.h 
	 rm -f GretinaTrackDictionary.cc GretinaTrackDictionary.h; rootcint -f $@ -c $(SWITCH) GretinaTrack.hh GretinaTrackLinkDef.h 

clean:
	rm  -f $(BIN_DIR)/* $(LIB_DIR)/* *.o *Dictionary.cc *Dictionary.h
	cd hfc; make clean

tar:
	@echo "creating zipped tar-ball ... "
	tar -chvzf GUnpack.tar.gz Makefile *LinkDef.h \
	*.hh *.cc *defs.h crmat.LINUX
