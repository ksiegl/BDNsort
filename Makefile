ifeq ($(ROOTSYS),)
$(error ROOTSYS is not define, do ssetup root)
endif
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
#INCDIRS := -I/opt/scarlet-3.x/include -I/opt/ktool-2.0/include -I/home/scaldwell/code/include 
CXX = g++
CXXFLAGS = -c -g -O2 $(ROOTCFLAGS) -I/opt/scarlet-3.x/include -I/opt/ktool-2.0/include -I/home/scaldwell/code/include #$(INCDIRS)
LIBS = -Wl,-rpath,/opt/scarlet-3.x/lib:/opt/ktool-2.0/lib \
 -L/opt/scarlet-3.x/lib -lscarlet -L/opt/ktool-2.0/lib -lktool

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $< -o $@

#%.d: %.cxx
#	$(CXX) -MM $(INCDIRS) $< -o $@

cxxsrcs = $(wildcard *.cxx)

.PHONY: all clean

targets = tof_cuts gate_on_low_tof_noise tof_from_E cooling no_spikes_sb135 draw_no_spikes_loop write_metadata no_spikes_diagnostic betas_vs_cycle_time betas_vs_cycle_time_i137 tof_official beta_gamma mcp_cal mcp_cal_i137 rf_phase gammas_vs_cycle_time beta_gamma_0 beta_gamma_1 bdn_sort_20130903 bdn_sort_20130923 bdn_sort_20130924 bdn_sort_20130925 bdn_sort_Ge_only bdn_sort_20131029 bdn_sort_empty bdn_sort_20131112 bdn_sort_ADC1_only bdn_sort_ADC1_TDC1_only bdn_sort_20131119 bdn_sort_20131120 bdn_sort_20131120_noLiveTime bdn_sort_20131125 bdn_sort_20131203 bdn_Sort_09272012_for_2013_run_grtrthan_1681 bdn_Sort_09272012_for_2013_run_lessthan_1682 bdn_sort_20131210 bdn_sort_20140104 bdn_Sort_09272012 bdn_Sort_09272012_for_137i02_run00002 BFit Metadata bdn_sort_20140308 mcp_cal_pedSubtract bdn_sort_20140417 DeadtimeCorrection bdn_sort_20140515 ExampleProgram bdn_sort_20140527 bdn_sort_20140613

all: $(targets)

tof_cuts: tof_cuts.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
gate_on_low_tof_noise: gate_on_low_tof_noise.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
tof_from_E: tof_from_E.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
cooling: cooling.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
no_spikes_sb135: no_spikes_sb135.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
draw_no_spikes_loop: draw_no_spikes_loop.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
write_metadata: write_metadata.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
no_spikes_diagnostic: no_spikes_diagnostic.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
betas_vs_cycle_time: betas_vs_cycle_time.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
betas_vs_cycle_time_i137: betas_vs_cycle_time_i137.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
tof_official: tof_official.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
beta_gamma: beta_gamma.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
mcp_cal: mcp_cal.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
mcp_cal_i137: mcp_cal_i137.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
rf_phase: rf_phase.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
gammas_vs_cycle_time: gammas_vs_cycle_time.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
beta_gamma_0: beta_gamma_0.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
beta_gamma_1: beta_gamma_1.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20130903: bdn_sort_20130903.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20130923: bdn_sort_20130923.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20130924: bdn_sort_20130924.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20130925: bdn_sort_20130925.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_Ge_only: bdn_sort_Ge_only.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20131029: bdn_sort_20131029.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_empty: bdn_sort_empty.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20131112: bdn_sort_20131112.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_ADC1_only: bdn_sort_ADC1_only.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_ADC1_TDC1_only: bdn_sort_ADC1_TDC1_only.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20131119: bdn_sort_20131119.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20131120: bdn_sort_20131120.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20131120_noLiveTime: bdn_sort_20131120_noLiveTime.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20131125: bdn_sort_20131125.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20131203: bdn_sort_20131203.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_Sort_09272012_for_2013_run_grtrthan_1681: bdn_Sort_09272012_for_2013_run_grtrthan_1681.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_Sort_09272012_for_2013_run_lessthan_1682: bdn_Sort_09272012_for_2013_run_lessthan_1682.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20131210: bdn_sort_20131210.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
#bdn_sort_20131211: bdn_sort_20131211.o
#	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
#En: En.o
#	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20140104: bdn_sort_20140104.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_Sort_09272012: bdn_Sort_09272012.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_Sort_09272012_for_137i02_run00002: bdn_Sort_09272012_for_137i02_run00002.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
#B_fit: B_fit.o bdn_cases.o B_fit_cases.o
#	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
BFit: BFit.o CSVtoStruct.o BFitModel.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
Metadata: Metadata.o CSVtoStruct.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20140308: bdn_sort_20140308.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
mcp_cal_pedSubtract: mcp_cal_pedSubtract.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20140417: bdn_sort_20140417.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
DeadtimeCorrection: DeadtimeCorrection.o CSVtoStruct.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20140515: bdn_sort_20140515.o bdn_histograms.o bdn_trees.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
ExampleProgram: ExampleProgram.o CSVtoStruct.o mcpGridCorrection.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20140527: bdn_sort_20140527.o bdn_histograms.o bdn_trees.o CSVtoStruct.o mcpGridCorrection.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
bdn_sort_20140613: bdn_sort_20140613.o bdn_histograms.o bdn_trees_20140613.o CSVtoStruct.o mcpGridCorrection.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
-include $(cxxsrcs:.cxx=.d)

clean:
	rm -f $(targets) *.o
	
