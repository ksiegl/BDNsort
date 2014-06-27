/*
// 2012-10-06 Shane Caldwell
// My version of sort code for BDN experiment.
// Made from combining bdn_Sort_09272012.cxx and bdn_online_10032012.cxx
// 2012-10-11 Made T_mcpX, etc., doubles in bdn struct and tree
// 2012-10-15 Added static_cast<double> to assignments to T_mcpX, etc.
// 2012-10-15 Added event_good flag to confirm that each event has all markers where expected
// 2012-10-23:
//	Adding TEventLists list_LT, list_LR, etc., for plastic-mcp coincidences with a minimal set of cuts
//  Also forming trees tree_LT, tree_LR, etc., in case I don't get the lists working soon....
// 2012-10-25:
//	Adding metadata tree.  Now reading file #.  n_hits variables record scaler counts of actual detector hits,
//	not just the ones the trigger the system.
// 2012-10-28:
//	Adding En spectra.
// 2012-10-31:
//	Fixed a problem with the definition of the trees
// 2012-11-03:
//	Added the scaler counting
//	Fixed the determination of s_capt_state to be correct for all 135Sb runs (via function "fix_capt_state")
//	Added flag "bkgd_good" to metadata_Tree to show when we were accidentally adding a NIS shot to the bkgd msmt (via function "get_bkgd_good_flag")
//	Both functions defined in sb135.h
//	Added month to metadata via header functions get_start_month and get_stop_month
// 2012-11-06:
//	Rearranged which scaler channels correspond to which detectors after finding the wiring completely jumbled
// 2012-11-07
//	Added s_SiX4_ts, sCyc_SiX4_ts, and sCyc_SiX4_ts to trig sync function
// 2013-01-23:
//	Taken from /disks/3/bpttrap/138-I/i138_sort_2012-11-18.cxx
//	Changed header i137.h (yes 137) to sb135.h.  The header files still need updating.
//	Now saving bdn.root to ./ROOT/files2
//	There should be a comment reflecting that the 2012-11-18 version completely omits bad events
//	On second thought I have made this go to ./ROOT/files as well, so it should be identical to sb135_sort_2012-11-18.cxx.
// 2013-09-03: 
//	Adding ADC2 readout and ha_X_ge_highE and ha_X_mcpE histos
// 2013-09-23:
//	Adding ha_X_dEsum plots to show dEa+dEb data
//	Also uncommented the tree assignment of Ge high-E and MCP-E variables... meant to do this on 09-03
//	Changed Ge Top and Right from ADC1 channel 9 and 1 to ADC2 ch. 9 and 10
// 2013-09-24:
//	Added filling of h_[T/R]_mcp[X/Y] histos. (They were already declared (as ha_...).)
// 2013-09-25
//	Adding a lot of diagnostic histograms to be filled at the end of the TRIGGERED event loop
//	New tree variable: a_X_dEsum	= a_X_dEa + a_X_dEb
//	New tree variable: t_X_dE		= 0.5*(t_X_dEa + t_X_dEb)
//	Now use cut on MCP sum (instead of each post) for the dE-MCP TOF histos
//	Lots of new coincident timing (TOF) histos
//	Made t_X_dE (the mean of dEa and dEb) a double instead of an int
// 2013-09-29: Adding zero_time tree for events in the zero-time peak
// 2013-10-29
//	Some new histos
//	Right Ge now on ADC2 channel 8 instead of channel 10
// 2013-11-12
//	New MCP map defintions (cut on sum instead of each post)
//  Moved MCP E channels from ADC2 to ADC1
// 2013-11-13
//	New MCP maps that show and slow recoils only
// 2013-11-19
//  Added live time scaler variables: tot_runTime_us (metadata), tot_liveTime_us (metadata), deadTime_us (bdn), all_trigs (metadata)
// 2013-11-20
//  Subtract zero-time values from TOF spectra
//	Add beta_recoil_tree
//	Changed placeholder values to -900060001 to avoid t_dE's averaging into physical range
// 2013-11-22
//  Make h_betas_vs_cycle_time
// 2013-11-25
//  Make h_all_vs_cycle_time, h_slow_vs_cycle_time, h_fast_vs_cycle_time, h_zero_vs_cycle_time, h_lowTOF_vs_cycle_time
//	Define event_deadtime and use it together with h_all_vs_cycle_time to correct other _vs_cycle_time plots.
// 2013-11-26 Rename the _vs_cycle_time histos to _vs_cycle_time_observed, since they are to be corrected for deadtime
// 2013-11-27 Adding calibrated Ge spectra in both singles (he_X_ge) and beta-gamma (he_beta_ge)
// 2013-12-02
//	Switching Top Ge ADC channel based on a flag on n_run
//	Switching scaler readout based on n_run as well
// 2013-12-03
//	Switched the vs_cycle_time binning from (302000,-999.5,301000.5) to (302000,-1000.5,300999.5)
//	Adding many new Ge singles and Beta-Gamma histos and reorganizing how they are filled
// 2013-12-10
//	Before the change to v. _20131210
//	...chaged the vs_cycle_time binning from (302000,-1000.5,300999.5) to (302000,-1000.0,301000.0)
//	After the change to v. 20131210
//	...subtracting E-MCP zerotimes from E-MCP TOF spectra
//	Created bdn_histograms.h and bdn_histograms.cxx to declare and construct histograms -- not working yet
//  *** Adding tof_LT, tof_LR, tof_BT, and tof_BR variables to beta_recoils_tree ***
//	Moved the filling of the beta-recoil _vs_cycle_time histos (zero, lowTOF, fast, slow, oops) out of the capt_state cut.
// 2014-01-04
//	Adding beta_gamma_tree for beta-gamma coincidences
// 2014-03-08
//	- Corrections for events where ADC range was exceeded. (These are the same as the events
//	  where there was no hit on the hit register.) This correction depends on the shape of the
//	  pulse height distribution, so alpha-source data must be corrected differently than ion data.
//	  So take argv[2] from the command line to control which correction is made.
//	- Adding Rndm() to numerator and denominator in calculation of MCP x and y coords.
//	- Adding count of missing MCP hits (eg. na_R_mcpA_missing) to metadata_tree
// 2014-04-17
//	* New non-tree variables for corrected MCP-post ADC values; these are important for calculating the X and Y coords on the MCP
//		- Example defn: a_T_mcpA_corr = a_T_mcpA - pedestal_T_mcpA - 0.5 + rand(0,1)
//		- This is done in the ADC readout
//		- These data are shown in a new set of histograms, eg. ha_T_mcpA_corr
//		- MCP coords are calculated using these corrected values
//		- 
//	* Imposing fiducial areas on MCP maps, -23.0mm to +23.0mm in each coord
//		- This cut must be applied to all data where we are counting MCP hits, esp TOF spectra
//		- ...
// 2014-04-24
//	Adding the reconstruction of a missing fourth post, based on code from current version of mcp_cal.cxx
//	Eight new tree vars (eg. miss_R_mcpA) to tell you which posts, if any, were missing.
//	Changing what the second cmd line argument does. Now "posts" turns on the one-post reconstruction. Missing posts are still counted.
// 2014-04-26
//	Removing the MCP map requirement that "all four posts get hit" (>0) from the maps
//	Renaming the map with no sum cut from 'h_R_mcpMap' to 'h_R_mcpMap_all'
//	Now 'h_R_mcpMap' is the map with the 'official cut' (a_mcp_lo found in bdn.h); fill this just below the calculated tree values
// 2014-05-05
//	Attempting to correct the map for one-missing-post events; the corrected points are in h_R_mcpMap_3post
//	New tree variables, eg. bdn.R_mcpPhysX, bdn.R_mcpPhysY, to hold the "small" MCP coords; from now on, bdn.R_mcpX and Y hold estimates of the physical coords
// 2014-05-12
//	Adding rf_phase to trees. This is a function of t_rf and new var rf_freq. In future rf_freq should be read in from an external source.
// 2014-05-15 New version
//	(Finally) moving histogram and tree declarations and definitions into external code (header and source).
// 2014-05-22
//	Added bdn.fid_area_hit_X_mcp to indicate events when the MCP fiducial areas were hit. Makes it easier to include this in cuts.
// 2014-05-27
//	Using CSVtoStruct to pass metadata to the sort code via the BDNCase_t struct.
//	Computing ion speed v, 1/v, and En in the Beta-Recoil coincidence section now.
//	New function for computing tof-to-grid from BDN case parameters (grid accel correction) in mcpGridCorrection.cxx.
//	To make the calculations easy to read I'm using new vars t1, t2, z1, x2, y2, s2, etc. This follows my notation from my thesis appendix.
//	Assigning placeholder values to beta-recoil variables tof, v, vInv, and En. These should really live only in the beta_recoil_tree but I don't want to spend the time to change this.
//	Deleted the commented-out block of tree and histo definitions.
//	Your BDN Case Code goes in the third command line argument, argv[3]
// 2014-06-10
//	Removed MCP Fiducial Area cut from selection of beta_recoil_tree events. It is still imposed on the counting of MCP data within the beta-recoil section.
// 2014-06-13
//	TOF cuts now use values from stBDNCase instead of the generic ones found in bdn.h. (Look for tof_T_fast_lo, tof_R_fast_hi, etc.; defs under comment "TOF bounds for this case")
//	Resurrecting the fast ion and slow ion counts for the metadata tree.
// 2014-06-17
//	The TOF counting appears to be working correctly now. I have it being done in two ways:
//	  1) "Counting" -- In the beta-recoil section of the event loop, I increment the metadata_Tree variables nFastCount[LT], nSlowCount[BR], etc., whenever the appropriate cuts are passed.
//	  2) "Integral" -- After the event loop, integrate the h_tof_LT/LR/BT/BR spectra over the relevant ranges.
//	The histos and the counters are incremented under the same cuts, so their exact agreement is a confirmation that things are working correctly.
//  Sort code now reports net recoils which are intended to be correct/final!
//
//////////////////////////////////////////////////////////////////////////////////////////////
//
//	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	GUIDE TO THE NAMES OF VARIABLES
//	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  "n_"	...	integer # of something, or a running total
//	"t_"	...	TDC data, or something in a linear relation to TDC data
//	"a_"	... ADC data, or something in a linear relation to ADC data
//	"e_"	... calibrated energy value
//	"s_"	... Scaler data
//	"r_"	... a rate (n_ or s_ divided by units of time
//	"h_"	... histogram
//	"hx_"	... histogram of x data, eg. ha_, ht_
//	"he_"	... histogram of calibrated energy values eg. he_ge
//	"BLRT"	... detector location (bottom, left, right, top)
//	"mcp, E, dEa, dEb, ge" ... detector type
//	"bg"	... beta-gamma coincidence data
//	"_observed" ... histo that needs to be corrected for deadtime but is not yet corrected (see DeadtimeCorrection.cxx)
//
///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Description of program:
//	The main function takes the name of a Scarlet data file (eg. run00241) from the cmd line.
//	It then goes into the data file and loops over the event headers in the file. For each event
//	type (ie. ACQUIRE, TRIGGERED, SYNC, and STOP) it executes code appropriate to that event.
//	For example for ACQUIRE and STOP events it reads the time stamp at the beginning and end of
//	the file. For TRIGGERED events it reads out ADC, TDC, and scalers data into a TTree and
//	fills histograms. For SYNC events it reads out the trig sync scaler and calculates some
//	rates. The trees and histos are saved in a ROOT file called bdn.root.
//
//	To execute:
//	  ./bdn_sort_<version> <run12345> <mcp_corr>
//	<version> is part of the name of this file
//	<run12345> is the runfile
//	<mcp_corr> is an optional second argument that fills in the missing above-range MCP pulse heights
//		mcp_corr = alpha -- gives the correction for alphas
//		mcp_corr = ion   -- gives the correction for recoil ions
//
/////////////////////////////////////////////////////////////////////////////////////////////// 
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "libgen.h"
#include "ScarletEvntSrc.h"
#include "ScarletEvnt.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TEventList.h"
#include "TRandom3.h"
#include "TMath.h"
#include "bdn.h"
#include "bdn_trees_20140613.h"
#include "bdn_histograms.h"
//#include "sb135.h"
#include "CSVtoStruct.h"

// Declare functions:
int countbit(int x);
int time_in_seconds(int, int, int, int);
Double_t tofToMCPGrid (BDNCase_t, char, Double_t);
int tofCounting();
int ReadADC1(int**, int, int*, int*);
int ReadADC2(int**, int, int*, int*);
int ReadTDC1(int**, int, int*, int*);
int ReadTDC2(int**, int, int*, int*);
int ReadScalers(int**, int, int, int*, int*, int*)

using namespace std;
using namespace TMath;
//using std::vector;
TRandom3 *randgen = new TRandom3(1);

int main(int argc, char *argv[]) {
	
	//printf("LT Zerotime = %f",LT_zeroTime[0]);
	
	// Get run # from command line and set ROOT filename
	char *argdup;
	argdup = strdup(argv[1]);
	char *filename;
	filename = basename(argdup);
	free(argdup);
	int n_run = atoi(&filename[3]);
	cout << endl << "Sorting " << filename;
	char *rootFileName;
	rootFileName ="bdn.root";
	
//	// To allow building ./sorted/runxxxxx.root directly instead of building bdn.root and renaming in the shell
//	// As is, TFile complains when you give it the rootFileName this produces
//	// So it doesn't work yet
//	std::stringstream ss;
//	ss << "./sorted/" << string(filename) << ".root";
//	std::string rFN = ss.str();
//	
//	//rootFileName = new char [rFN.size()+1];
//	char *rootFileName[rFN.size()+1];
//	strcpy (rootFileName, rFN.c_str());
//	printf(rootFileName);	
	
// Metadata structure
	BDNCase_t	stBDNCases[FILE_ROWS_BDN];
	int			iBDNCaseIndex;
	int			iNumStructs_BDN;
	char *csvBDNCases;
	csvBDNCases = "BDNCases.csv_transposed";
	cout << endl << "Importing metadata from CSV files..." << endl;
	iNumStructs_BDN  = CSVtoStruct_BDN  (csvBDNCases, stBDNCases);
	cout << "Imported " << iNumStructs_BDN << " BDN cases" << endl;
	iBDNCaseIndex		 = FindStructIndex ( stBDNCases,  sizeof(BDNCase_t),  iNumStructs_BDN,  argv[3] );
	BDNCase_t  stBDNCase = stBDNCases[iBDNCaseIndex];
	// Optional error catching
	if ( iBDNCaseIndex == -1 )
	{ // One of the read-ins failed and already printed a message about it
		cout << "How to run this program:" << endl;
		cout << "'./ExampleProgram <BDN case code>'" << endl;
		cout << "where valid case codes are listed in the CSV files." << endl << endl;
		return -1; // error return
	}
// tofCounting()	
// TOF bounds for this case	
	Double_t	tof_R_fast_lo		= 1000.0 * stBDNCase.dRightMCPMinFastIonTOF;
	Double_t	tof_R_fast_hi		= 1000.0 * stBDNCase.dRightMCPMaxFastIonTOF;
	Double_t	tof_T_fast_lo		= 1000.0 * stBDNCase.dTopMCPMinFastIonTOF;
	Double_t	tof_T_fast_hi		= 1000.0 * stBDNCase.dTopMCPMaxFastIonTOF;
	Double_t	tof_R_slow_lo		= tof_R_fast_hi;
	Double_t	tof_R_slow_hi		= tof_oops_lo; // from bdn.h
	Double_t	tof_T_slow_lo		= tof_T_fast_hi;
	Double_t	tof_T_slow_hi		= tof_oops_lo; // from bdn.h
	Float_t		tofBinsPerNs		= TOFBins/(TOFMax-TOFMin); // from bdn_histograms.h
	Int_t		tofBin_R_fast_lo	= Nint( tofBinsPerNs * (tof_R_fast_lo - TOFMin) + 1 );
	Int_t		tofBin_R_fast_hi	= Nint( tofBinsPerNs * (tof_R_fast_hi - TOFMin) + 1 );
	Int_t		tofBin_T_fast_lo	= Nint( tofBinsPerNs * (tof_T_fast_lo - TOFMin) + 1 );
	Int_t		tofBin_T_fast_hi	= Nint( tofBinsPerNs * (tof_T_fast_hi - TOFMin) + 1 );
	Int_t		tofBin_R_slow_lo	= Nint( tofBinsPerNs * (tof_R_slow_lo - TOFMin) + 1 );
	Int_t		tofBin_R_slow_hi	= Nint( tofBinsPerNs * (tof_R_slow_hi - TOFMin) + 1 );
	Int_t		tofBin_T_slow_lo	= Nint( tofBinsPerNs * (tof_T_slow_lo - TOFMin) + 1 );
	Int_t		tofBin_T_slow_hi	= Nint( tofBinsPerNs * (tof_T_slow_hi - TOFMin) + 1 );
	Int_t		tofBin_oops_lo		= Nint( tofBinsPerNs * (tof_oops_lo   - TOFMin) + 1 );
	Int_t		tofBin_oops_hi		= Nint( tofBinsPerNs * (tof_oops_hi   - TOFMin) + 1 );
	Int_t		tofBin_zero_lo		= Nint( tofBinsPerNs * (tof_zero_lo   - TOFMin) + 1 );
	Int_t		tofBin_zero_hi		= Nint( tofBinsPerNs * (tof_zero_hi   - TOFMin) + 1 );
	Int_t		tofBin_lowTOF_lo	= Nint( tofBinsPerNs * (tof_lowTOF_lo - TOFMin) + 1 );
	Int_t		tofBin_lowTOF_hi	= Nint( tofBinsPerNs * (tof_lowTOF_hi - TOFMin) + 1 );
	printf("\nCounting   Top MCP recoils using TOF ranges:\n\tFast\t\t= (%8.2f,%8.2f) ns\t= bins (%6d,%6d);\n\tSlow\t\t= (%8.2f,%8.2f) ns\t= bins (%6d,%6d);\n\tAccidentals\t= (%8.2f,%8.2f) ns\t= bins (%6d,%6d).", tof_T_fast_lo, tof_T_fast_hi, tofBin_T_fast_lo, tofBin_T_fast_hi, tof_T_slow_lo, tof_T_slow_hi, tofBin_T_slow_lo, tofBin_T_slow_hi, tof_oops_lo, tof_oops_hi, tofBin_oops_lo, tofBin_oops_hi);
	printf("\nCounting Right MCP recoils using TOF ranges:\n\tFast\t\t= (%8.2f,%8.2f) ns\t= bins (%6d,%6d);\n\tSlow\t\t= (%8.2f,%8.2f) ns\t= bins (%6d,%6d);\n\tAccidentals\t= (%8.2f,%8.2f) ns\t= bins (%6d,%6d).", tof_R_fast_lo, tof_R_fast_hi, tofBin_R_fast_lo, tofBin_R_fast_hi, tof_R_slow_lo, tof_R_slow_hi, tofBin_R_slow_lo, tofBin_R_slow_hi, tof_oops_lo, tof_oops_hi, tofBin_oops_lo, tofBin_oops_hi);
//	printf("\nIntegrating   Top MCP TOF spectra over ranges :\n\tFast =\t\tbin(%d,%d);\n\tSlow =\t\tbin(%d,%d);\n\tAccidentals =\tbin(%d,%d).", , , );
//	printf("\nCounting    Right MCP recoils using TOF ranges:\n\tFast\t\t= (%g,%g) us;\n\tSlow\t\t= (%g,%g) us;\n\tAccidentals\t= (%g,%g) us.", tof_R_fast_lo, tof_R_fast_hi, tof_R_slow_lo, tof_R_slow_hi, tof_oops_lo, tof_oops_hi);
//	printf("\nIntegrating Right MCP TOF spectra over ranges :\n\tFast\t\tbin(%d,%d);\n\tSlow =\t\tbin(%d,%d);\n\tAccidentals =\tbin(%d,%d).\n", tofBin_R_fast_lo, tofBin_R_fast_hi, tofBin_R_slow_lo, tofBin_R_slow_hi, tofBin_oops_lo, tofBin_oops_hi);
/*	
//	Float_t nZeroTOFCount[4]		= {0.0, 0.0, 0.0, 0.0};
//	Float_t nLowTOFCount[4]			= {0.0, 0.0, 0.0, 0.0};
//	Float_t nOopsCount[4]			= {0.0, 0.0, 0.0, 0.0};
//	Float_t nFastCount[4]			= {0.0, 0.0, 0.0, 0.0};
//	Float_t nSlowCount[4]			= {0.0, 0.0, 0.0, 0.0};
//	Float_t nZeroTOFBkgdCount[4]	= {0.0, 0.0, 0.0, 0.0};
//	Float_t nLowTOFBkgdCount[4]		= {0.0, 0.0, 0.0, 0.0};
//	Float_t nOopsBkgdCount[4]		= {0.0, 0.0, 0.0, 0.0};
//	Float_t nFastBkgdCount[4]		= {0.0, 0.0, 0.0, 0.0};
//	Float_t nSlowBkgdCount[4]		= {0.0, 0.0, 0.0, 0.0};
//	Float_t nZeroTOFIntegral[4]		= {0.0, 0.0, 0.0, 0.0};
//	Float_t nLowTOFIntegral[4]		= {0.0, 0.0, 0.0, 0.0};
//	Float_t nOopsIntegral[4]		= {0.0, 0.0, 0.0, 0.0};
//	Float_t nFastIntegral[4]		= {0.0, 0.0, 0.0, 0.0};
//	Float_t nSlowIntegral[4]		= {0.0, 0.0, 0.0, 0.0};
//	Float_t nZeroTOFBkgdIntegral[4]	= {0.0, 0.0, 0.0, 0.0};
//	Float_t nLowTOFBkgdIntegral[4]	= {0.0, 0.0, 0.0, 0.0};
//	Float_t nOopsBkgdIntegral[4]	= {0.0, 0.0, 0.0, 0.0};
//	Float_t nFastBkgdIntegral[4]	= {0.0, 0.0, 0.0, 0.0};
//	Float_t nSlowBkgdIntegral[4]	= {0.0, 0.0, 0.0, 0.0};
////	std::array<float,4> nFastCount	= {0.0, 0.0, 0.0, 0.0};
*/	
// Scarlet variables
	ScarletEvntSrc *esrc;
	ScarletEvntHdr *h;
	ScarletEvnt e0, e1;
	int *p;
	
	// Data structs:
	// struct Globals{int t_eject, t_capt, s_capt_state, s_capt, rf_phase, clock, clock_tot, run, event;};
	// struct MCP{int t, A, B, C, D, E, x, y;};
	// struct Plastic{int t_dEa, t_dEb, t_E, a_dEa, a_dEb, a_E;};
	// struct Ge{int t, a;};
	
	// ADC values:
	int a_placeholder = -900060001;
	int a_B_dEa, a_B_dEb, a_B_E;
	int a_L_dEa, a_L_dEb, a_L_E;
	int a_R_mcpA, a_R_mcpB, a_R_mcpC, a_R_mcpD, a_R_mcpE, a_R_ge, a_R_ge_highE;
	int a_T_mcpA, a_T_mcpB, a_T_mcpC, a_T_mcpD, a_T_mcpE, a_T_ge, a_T_ge_highE;
	double a_R_mcpA_corr, a_R_mcpB_corr, a_R_mcpC_corr, a_R_mcpD_corr, a_R_mcpE_corr, a_R_mcpSum_corr;
	double a_T_mcpA_corr, a_T_mcpB_corr, a_T_mcpC_corr, a_T_mcpD_corr, a_T_mcpE_corr, a_T_mcpSum_corr;
		// "_ge" 		= HPGe 0-3.6 MeV range
		// "_ge_highE"	= HPGe 0-9.2 MeV range
	
	// TDC values:
	int t_placeholder = -900060001;
	int t_B_dEa, t_B_dEb, t_B_E;
	int t_L_dEa, t_L_dEb, t_L_E;
	int t_R_mcp, t_R_ge;
	int t_T_mcp, t_T_ge;
	int t_rf;
	
	// Scaler values:
	int s_B_dEa, s_B_dEb, s_B_E;
	int s_L_dEa, s_L_dEb, s_L_E;
	int s_R_mcp, s_R_ge;
	int s_T_mcp, s_T_ge;
	int s_ms_since_capt;			// time since last capt pulse in ms
	int s_ms_since_eject	= 0;	// time since last eject pulse in ms
	int s_capt_state;			// 0 = trap full, 1 = trap empty
	int s_capt;				// # of capturesin current trap cycle
	int s_SiX4		= 0;
	int s_liveTime_us = 0;
	int deadTime_us = 0;
	int all_trigs = 0;
	long tot_trigs = 0;
	
	// Energy values (new 2013-12-01)
	
	Double_t y, e_B_E, e_L_E, e_R_ge, e_R_ge_highE, e_T_ge, e_T_ge_highE;
	
	// Event counts for ADC (na_), TDC (nt_):
	int na_B_dEa		= 0;
	int na_B_dEb		= 0;
	int na_B_E			= 0;
	int na_L_dEa		= 0;
	int na_L_dEb		= 0;
	int na_L_E			= 0;
	int na_R_mcpA		= 0; int na_R_mcpA_missing = 0;
	int na_R_mcpB		= 0; int na_R_mcpB_missing = 0;
	int na_R_mcpC		= 0; int na_R_mcpC_missing = 0;
	int na_R_mcpD		= 0; int na_R_mcpD_missing = 0;
	int na_R_mcpE		= 0; int na_R_mcpE_missing = 0;
	int na_R_ge			= 0;
	int na_R_ge_highE	= 0;
	int na_T_mcpA		= 0; int na_T_mcpA_missing = 0;
	int na_T_mcpB		= 0; int na_T_mcpB_missing = 0;
	int na_T_mcpC		= 0; int na_T_mcpC_missing = 0;
	int na_T_mcpD		= 0; int na_T_mcpD_missing = 0;
	int na_T_mcpE		= 0; int na_T_mcpE_missing = 0;
	int na_T_ge			= 0;
	int na_T_ge_highE	= 0;
	
	int nt_B_dEa	= 0;
	int nt_B_dEb	= 0;
	int nt_B_E		= 0;
	int nt_L_dEa	= 0;
	int nt_L_dEb	= 0;
	int nt_L_E		= 0;
	int nt_R_mcp	= 0;
	int nt_R_ge		= 0;
	int nt_T_mcp	= 0;
	int nt_T_ge		= 0;
	int nt_all;
	
	// Variables for ACQUIRE events:
	int start_day, start_hour, start_min, start_sec;
	
	// Variables for TRIGGERED events:
	int n_trig		= 0; // # of TRIGGERED events
	int trig_bin	= 0; // bin # for histos on TRIGGERED data
	int s_ms_since_eject_last;
	int s_SiX4_last;
	int sTot_SiX4	= 0;
	int event_good	= 1; // 0 means a data marker was not found where expected
	int n_bad_events= 0;
	//int trig_day, trig_hour, trig_min, trig_sec;	// not used
	
	// Variables for SYNC events
	int n_sync		= 0;	// # of SYNC events
	int n_sync_update = 1;	// # of syncs between updates of sync data
	int sync_bin	= 0;	// bin # for histos on SYNC data
	int sCyc_B_dEa	= 0;	int sTot_B_dEa	= 0;
	int sCyc_B_dEb	= 0;	int sTot_B_dEb	= 0;
	int sCyc_B_E	= 0;	int sTot_B_E	= 0;
	int sCyc_L_dEa	= 0;	int sTot_L_dEa	= 0;
	int sCyc_L_dEb	= 0;	int sTot_L_dEb	= 0;
	int sCyc_L_E	= 0;	int sTot_L_E	= 0;
	int sCyc_R_mcp	= 0;	int sTot_R_mcp	= 0;
	int sCyc_R_ge	= 0;	int sTot_R_ge	= 0;
	int sCyc_T_mcp	= 0;	int sTot_T_mcp	= 0;
	int sCyc_T_ge	= 0;	int sTot_T_ge	= 0;
	int sCyc_SiX4_ts= 0;	int sTot_SiX4_ts= 0; // ts for trig sync
	int s_SiX4_ts	= 0;
	int sTot_all;
	int sync_day_last, sync_day, sync_hour, sync_min, sync_sec;
	
	// Variables for STOP events:
	int stop_flag = 0;
	int stop_day, stop_hour, stop_min, stop_sec, stop_time_sec;
	long tot_liveTime_us = 0;
	long tot_runTime_us = 0;
	long s_runTime = 1;
		// Other variables:
	int fake_day; // to stand in for next day when the day goes back to 0, eg. Oct 32
	int bkgd_good	= 1;//get_bkgd_good_flag(n_run); // 0 = not good = we were accidentally sending in 1 shot from the NIS
//	int start_month	= get_start_month(n_run);
//	int stop_month	= get_stop_month(n_run);
	
	Double_t t_TrapFull	= stBDNCase.dCycleTime - stBDNCase.dBackgroundTime;
	Double_t t_cycle	= stBDNCase.dCycleTime;
	int now_day, now_hour, now_min, now_sec, now_temp;
	int start_time_sec	= 0;
	int trig_time_sec	= 0;
	int sync_time_sec	= 0;
	int now_time_sec;
	int run_time_min, run_time_sec, run_remainder_sec;
	int clock = 0;
	int lastClock = 0;
	int event = 0;
	int first_cycle_flag = 0;	// Indicates that first (partial) trapping cycle is complete
								// Set to 1 upon first ejection pulse
//	double mcp_x, mcp_y, mcp_r, tof, vInv, En;
	double n_fast_LT = 0.0;
	double n_fast_LR = 0.0;
	double n_fast_BT = 0.0;
	double n_fast_BR = 0.0;
	double n_slow_LT = 0.0;
	double n_slow_LR = 0.0;
	double n_slow_BT = 0.0;
	double n_slow_BR = 0.0;
	
	Double_t xx = 0.0;
	Double_t yy = 0.0;
	Double_t xx1 = 0.0;
	Double_t yy1 = 0.0;
	Double_t th = 0.0;
	Double_t x0 = 0.0;
	Double_t y0 = 0.0;
	Double_t a0 = 0.0;
	Double_t b0 = 0.0;
	Double_t a2 = 0.0;
	Double_t b2 = 0.0;
	
	Double_t t1, z1, t2, x2, y2, s2, vs, vz;
	const Double_t c	= 299792.46; // speed of light in mm/us
	
	int x;	  //data
	int adc_ch;
	int tdc_ch;
	int i=1;
	int j;
	int wordc;	  //word counts in a sequential readout
	
	// ROOT and Scarlet variables:
	TFile *f = new TFile(rootFileName, "recreate");
	book_trees();
	book_histograms();
//	extern bdn_struct bdn;
//	extern metadata_struct metadata;
	
// Procedure:
	cout << endl;
	if (argc <= 1) return 0;
	try {esrc = new ScarletFileSrc(argv[1]);}
	catch (std::exception &x){
		std::cerr << "exception opening " << argv[1] << std::endl;
		return 0;
	}
	
//	if (!strcmp(argv[2],"alpha")) {
//		printf("Correcting MCP pulse heights for missing data, assuming alphas.\n");
//	}
//	if (!strcmp(argv[2],"ion")) {
//		printf("Correcting MCP pulse heights for missing data, assuming ions.\n");
//	}
	if (!strcmp(argv[2],"posts")) {
		printf("Reconstructing missing MCP post values when only one is missing.\n");
	}
	
	while ((h = esrc->getevent()) != 0) {
		
		switch (h->type) {
			
			case SE_TYPE_TRIGGERED:
				
				n_trig++;
				//if (n_trig%1000==0) printf("event %d",n_trig);
				
			// Store latest values:
				s_ms_since_eject_last	= s_ms_since_eject;
				s_SiX4_last 			= s_SiX4;
				
			// ADC placeholder values:
		        a_B_dEa			= a_placeholder;
		        a_B_dEb			= a_placeholder;
		        a_B_E			= a_placeholder;
				a_L_dEa			= a_placeholder;
				a_L_dEb			= a_placeholder;
				a_L_E			= a_placeholder;
				a_R_mcpA		= a_placeholder;
				a_R_mcpB		= a_placeholder;
				a_R_mcpC		= a_placeholder;
				a_R_mcpD		= a_placeholder;
				a_R_mcpE		= a_placeholder;
				a_R_ge			= a_placeholder;
				a_R_ge_highE	= a_placeholder;
				a_T_mcpA		= a_placeholder;
				a_T_mcpB		= a_placeholder;
				a_T_mcpC		= a_placeholder;
				a_T_mcpD		= a_placeholder;
				a_T_mcpE		= a_placeholder;
				a_T_ge			= a_placeholder;
				a_T_ge_highE	= a_placeholder;
				e_R_ge			= a_placeholder;
				e_T_ge			= a_placeholder;
				e_R_ge_highE	= a_placeholder;
				e_T_ge_highE	= a_placeholder;
				e_B_E			= a_placeholder;
				e_L_E			= a_placeholder;
				a_R_mcpA_corr	= a_placeholder;
				a_R_mcpB_corr	= a_placeholder;
				a_R_mcpC_corr	= a_placeholder;
				a_R_mcpD_corr	= a_placeholder;
				a_R_mcpE_corr	= a_placeholder;
				a_T_mcpA_corr	= a_placeholder;
				a_T_mcpB_corr	= a_placeholder;
				a_T_mcpC_corr	= a_placeholder;
				a_T_mcpD_corr	= a_placeholder;
				a_T_mcpE_corr	= a_placeholder;
				
			// TDC placeholder values:
				t_B_dEa			= t_placeholder;
				t_B_dEb			= t_placeholder;
				t_B_E			= t_placeholder;
				t_L_dEa			= t_placeholder;
				t_L_dEb			= t_placeholder;
				t_L_E			= t_placeholder;
				t_R_mcp			= t_placeholder;
				t_R_ge			= t_placeholder;
				t_T_mcp			= t_placeholder;
				t_T_ge			= t_placeholder;
				t_rf			= t_placeholder;
			
			// Coincidence values
				bdn.tof_LT		= t_placeholder;
				bdn.tof_LR		= t_placeholder;
				bdn.tof_BT		= t_placeholder;
				bdn.tof_BR		= t_placeholder;
				bdn.v_LT		= t_placeholder;
				bdn.v_LR		= t_placeholder;
				bdn.v_BT		= t_placeholder;
				bdn.v_BR		= t_placeholder;
				bdn.En_LT		= t_placeholder;
				bdn.En_LR		= t_placeholder;
				bdn.En_BT		= t_placeholder;
				bdn.En_BR		= t_placeholder;
				
			// Missing post flags
				bdn.miss_R_mcpA = 0;
				bdn.miss_R_mcpB = 0;
				bdn.miss_R_mcpC = 0;
				bdn.miss_R_mcpD = 0;
				bdn.miss_T_mcpA = 0;
				bdn.miss_T_mcpB = 0;
				bdn.miss_T_mcpC = 0;
				bdn.miss_T_mcpD = 0;
				bdn.fid_area_hit_R_mcp = 0;
				bdn.fid_area_hit_T_mcp = 0;
				
				deadTime_us = 0;
				s_liveTime_us = 0;
				s_runTime = 1; // this is because clearing takes up 100 ns, which equals to one tick of our 10 MHz clock
				
			// Other default values:
				event_good		= 1;
				
			// Initialize pointer:
				e0 = ScarletEvnt(h);
				e1 = e0[1];
				p = reinterpret_cast<int*>(e1.body());
			// ADC1 *******************************
				ReadADC1(&p,n_trig,&event_good,&n_bad_events);
				/*
				if (*p != 0xadc1adc1) {
					cout << "trig #" << n_trig << ", ADC1 marker not found where expected!" << endl;
					event_good = 0;
					n_bad_events++;
					break;
				}
				*p++; // move to ADC1 hit register
				x = int(*p & 0xffff ); // hit register, tells which channels were hit
				wordc = countbit(x);
				
				for (j=1; j<=wordc; j++) { // Loop over all ADC channels which have hits
					
					*p++;  // Increment pointer p to ADC channel with a hit
					x = int(*p & 0x0fff);
				
					// Get ADC channel of that hit then associate it with the data:
					adc_ch=int(*p & 0xf000);
					adc_ch=(adc_ch>>12)+1;
					//if (adc_ch == 1) {
					//	a_R_ge = x;
					//	ha_R_ge->Fill(x);
					//	na_R_ge++;
					//}
					if (adc_ch == 1) {
						a_T_mcpE = x;
						ha_T_mcpE->Fill(x);
						a_T_mcpE_corr = x - ped_T_mcpE + randgen->Rndm();
						ha_T_mcpE_corr->Fill(a_T_mcpE_corr);
						na_T_mcpE++;
					}
					if (adc_ch == 2) {
						a_B_dEa = x;
						ha_B_dEa->Fill(x);
						na_B_dEa++;
					}
					if (adc_ch == 3) {
						a_B_dEb = x;
						ha_B_dEb->Fill(x);
						na_B_dEb++;
					}
					if (adc_ch == 4) {
						a_B_E = x;
						ha_B_E->Fill(x);
						na_B_E++;
					}
					if (adc_ch == 5) {
						a_T_mcpA = x;
						ha_T_mcpA->Fill(x);
						a_T_mcpA_corr = x - ped_T_mcpA + randgen->Rndm();
						ha_T_mcpA_corr->Fill(a_T_mcpA_corr);
						na_T_mcpA++;
					}
					if (adc_ch == 6) {
						a_T_mcpB = x;
						ha_T_mcpB->Fill(x);
						a_T_mcpB_corr = x - ped_T_mcpB + randgen->Rndm();
						ha_T_mcpB_corr->Fill(a_T_mcpB_corr);
						na_T_mcpB++;
					}
					if (adc_ch == 7) {
						a_T_mcpC = x;
						ha_T_mcpC->Fill(x);
						a_T_mcpC_corr = x - ped_T_mcpC + randgen->Rndm();
						ha_T_mcpC_corr->Fill(a_T_mcpC_corr);
						na_T_mcpC++;
					}
					if (adc_ch == 8) {
						a_T_mcpD = x;
						ha_T_mcpD->Fill(x);
						a_T_mcpD_corr = x - ped_T_mcpD + randgen->Rndm();
						ha_T_mcpD_corr->Fill(a_T_mcpD_corr);
						na_T_mcpD++;
					}
					//if (adc_ch == 9) {
					//	a_T_ge = x;
					//	ha_T_ge->Fill(x);
					//	na_T_ge++;
					//}
					if (adc_ch == 9) {
						a_R_mcpE = x;
						ha_R_mcpE->Fill(x);
						a_R_mcpE_corr = x - ped_R_mcpE + randgen->Rndm();
						ha_R_mcpE_corr->Fill(a_R_mcpE_corr);
						na_R_mcpE++;
					}
					if (adc_ch == 10) {
						a_L_dEa = x;
						ha_L_dEa->Fill(x);
						na_L_dEa++;
					}
					if (adc_ch == 11) {
						a_L_dEb = x;
						ha_L_dEb->Fill(x);
						na_L_dEb++;
					}
					if (adc_ch == 12) {
						a_L_E = x;
						ha_L_E->Fill(x);
						na_L_E++;
					}
					if (adc_ch == 13) {
						a_R_mcpA = x;
						ha_R_mcpA->Fill(x);
						a_R_mcpA_corr = x - ped_R_mcpA + randgen->Rndm();
						ha_R_mcpA_corr->Fill(a_R_mcpA_corr);
						na_R_mcpA++;
					}
					if (adc_ch == 14) {
						a_R_mcpB = x;
						ha_R_mcpB->Fill(x);
						a_R_mcpB_corr = x - ped_R_mcpB + randgen->Rndm();
						ha_R_mcpB_corr->Fill(a_R_mcpB_corr);
						na_R_mcpB++;
					}
					if (adc_ch == 15) {
						a_R_mcpC = x;
						ha_R_mcpC->Fill(x);
						a_R_mcpC_corr = x - ped_R_mcpC + randgen->Rndm();
						ha_R_mcpC_corr->Fill(a_R_mcpC_corr);
						na_R_mcpC++;
					}
					if (adc_ch == 16) {
						a_R_mcpD = x;
						ha_R_mcpD->Fill(x);
						a_R_mcpD_corr = x - ped_R_mcpD + randgen->Rndm();
						ha_R_mcpD_corr->Fill(a_R_mcpD_corr);
						na_R_mcpD++;
					}
					
				} // for (wordc)
				*/
			// ADC2 *******************************
				ReadADC2(&p,n_trig,&event_good,&n_bad_events);
				/*
				*p++; // move pointer to ADC2 marker
				if (*p != 0xadc2adc2) {
					cout << "trig #" << n_trig << ", ADC2 marker not found where expected!" << endl;
					event_good = 0;
					n_bad_events++;
					break;
				}
				*p++; // move to ADC2 hit register
				
				x = int(*p & 0xffff ); // hit register, tells which channels were hit
				wordc = countbit(x);
				
				for (j=1; j<=wordc; j++){ // Loop over all ADC channels which have hits
					
					*p++;  // Increment pointer p to ADC channel with a hit
					x = int(*p & 0x0fff);
					
					// Get ADC channel of that hit then associate it with the data:
					adc_ch=int(*p & 0xf000);
					adc_ch=(adc_ch>>12)+1;
					
					if (adc_ch == 1) {
						y = x + randgen->Rndm();
						a_T_ge_highE = x;
						e_T_ge_highE = T_ge_highE_coeff[0] + y*T_ge_highE_coeff[1] + y*y*T_ge_highE_coeff[2];
						ha_T_ge_highE->Fill(a_T_ge_highE);
						//he_T_ge_highE->Fill(e_T_ge_highE);
						//he_ge_highE	 ->Fill(e_T_ge_highE);
						na_T_ge_highE++;
					}
					if (adc_ch == 2) {
						y = x + randgen->Rndm();
						a_R_ge_highE = x;
						e_R_ge_highE = R_ge_highE_coeff[0] + y*R_ge_highE_coeff[1] + y*y*R_ge_highE_coeff[2];
						ha_R_ge_highE->Fill(x);
						//he_R_ge_highE->Fill(e_R_ge_highE);
						//he_ge_highE	 ->Fill(e_R_ge_highE);
						na_R_ge_highE++;
					}
					//if (adc_ch == 5) {
					//	a_T_mcpE = x;
					//	ha_T_mcpE->Fill(x);
					//	na_T_mcpE++;
					//}
					//if (adc_ch == 6) {
					//	a_R_mcpE = x;
					//	ha_R_mcpE->Fill(x);
					//	na_R_mcpE++;
					//}
					if (n_run < 1682) {
						if (adc_ch == 9) {
							y = x + randgen->Rndm();
							a_T_ge = x;
							e_T_ge = T_ge_coeff[0] + y*T_ge_coeff[1] + y*y*T_ge_coeff[2];
							ha_T_ge	->Fill(x);
							he_T_ge	->Fill(e_T_ge);
							he_ge	->Fill(e_T_ge);
							na_T_ge++;
						}
					}
					else {
						if (adc_ch == 7) {
							y = x + randgen->Rndm();
							a_T_ge = x;
							e_T_ge = T_ge_coeff[0] + y*T_ge_coeff[1] + y*y*T_ge_coeff[2];
							ha_T_ge	->Fill(x);
							he_T_ge	->Fill(e_T_ge);
							he_ge	->Fill(e_T_ge);
							na_T_ge++;
						}
					}
					if (adc_ch == 8) {
						y = x + randgen->Rndm();
						a_R_ge = x;
						e_R_ge = R_ge_coeff[0] + y*R_ge_coeff[1] + y*y*R_ge_coeff[2];
						ha_R_ge	->Fill(x);
						he_R_ge	->Fill(e_R_ge);
						he_ge	->Fill(e_R_ge);
						na_R_ge++;
					}
				} // for (wordc)
				*/
				
			// TDC1 *******************************
				ReadTDC1(&p,n_trig,&event_good,&n_bad_events);
				/*
				*p++; // move pointer to TDC1 marker
				if (*p != 0x2dc12dc1) {
					cout << "trig #" << n_trig << ", TDC1 marker not found where expected!" << endl;
					event_good = 0;
					n_bad_events++;
					break;
				}
				*p++;
				while (*p != 0x2dc22dc2) {
					tdc_ch = *p;
					*p++;
					x = int(*p & 0x00ffffff); // take only the 24-bit data word
					if (x & 0x0080000) x -= 0x00ffffff; // test for neg value
						// if neg then you need to shift because the leading 1 in 24-bit
						// is not leading in 32-bit; the shift is by "-0x00ffffff"
					if (tdc_ch==1) {
						t_T_mcp = x;
						ht_T_mcp->Fill(x);
						nt_T_mcp++;
					}
					if (tdc_ch==2) {
						t_R_mcp = x;
						ht_R_mcp->Fill(x);
						nt_R_mcp++;
					}
					if (tdc_ch==3) {
						t_B_dEa = x;
						ht_B_dEa->Fill(x);
						nt_B_dEa++;
					}
					if (tdc_ch==4) {
						t_B_dEb = x;
						ht_B_dEb->Fill(x);
						nt_B_dEb++;
					}
					if (tdc_ch==5) {
						t_B_E = x;
						ht_B_E->Fill(x);
						nt_B_E++;
					}
					if (tdc_ch==6) {
						t_L_dEa = x;
						ht_L_dEa->Fill(x);
						nt_L_dEa++;
					}
					if (tdc_ch==7) {
						t_L_dEb = x;
						ht_L_dEb->Fill(x);
						nt_L_dEb++;
					}
					if (tdc_ch==8) {
						t_L_E = x;
						ht_L_E->Fill(x);
						nt_L_E++;
					}
				} // while
				*/
			// TDC2 *******************************				
				ReadTDC2(&p,n_trig,&event_good,&n_bad_events);
				/*
				if (*p != 0x2dc22dc2) {
					cout << "trig #" << n_trig << ", TDC2 marker not found where expected!" << endl;
					event_good = 0;
					n_bad_events++;
					break;
				}
				*p++;
				while (*p != 0x100cca1e) {
					tdc_ch = *p;
					*p++;
					x = int(*p & 0x00ffffff); // take only the 24-bit data word
					if (x & 0x0080000) x -= 0x00ffffff; // test for neg value
						// if neg then you need to shift because the leading 1 in 24-bit
						// is not leading in 32-bit; the shift is by "-0x00ffffff"
					if (tdc_ch==1) {
						t_rf = x;
						ht_rf->Fill(x);
					}
					if (tdc_ch==2) {
						t_T_ge = x;
						ht_T_ge->Fill(x);
						nt_T_ge++;
					}
					if (tdc_ch==3) {
						t_R_ge = x;
						ht_R_ge->Fill(x);
						nt_R_ge++;
					}
				} // end while
				*/
			// Scalers ****************************	
				ReadScalers(&p, n_trig, n_run, &all_trigs, &event_good, &n_bad_events);
				/*
				if (n_run < 1201) { // old scaler readout
				
				// Capt Scaler ************************
					if (*p != 0x100cca1e) {
						cout << "trig #" << n_trig << ", Capt Scaler marker not found where expected!" << endl;
						event_good = 0;
						n_bad_events++;
						break;
					}
					*p++; // p is at time since capture in ms
					s_ms_since_capt = int(*p & 0xffffff);
					*p++; // p is at trap state | 0 = trap full, 1 = trap empty
					s_capt_state = int(*p & 0xffffff);
					*p++; // move pointer to eject scaler
					
				// Eject Scaler ***********************
					if (*p != 0x100eca1e) {
						cout << "trig #" << n_trig << ", Eject Scaler marker not found where expected!" << endl;
						event_good = 0;
						n_bad_events++;
						break;
					}
					*p++; //p is at time since eject in ms
					s_ms_since_eject = int(*p & 0xffffff);
					*p++; //p is at # of capt since last eject
					s_capt = int(*p & 0xffffff);
					*p++; //p is at # of SiX4 hits since last eject
					s_SiX4 = int(*p & 0xffffff);
				
				} // end old scaler readout
				
				else {	// new scaler readout
				
				// Live Time Scaler ***********************
					if (*p != 0x100cca1e) {
						cout << "trig #" << n_trig << ", Livetime Scaler marker not found where expected!" << endl;
						event_good = 0;
						n_bad_events++;
						break;
					}
					
					*p++; // p is..  this is a bit iffy
					s_liveTime_us = int(*p & 0xffffff);
					*p++;
					all_trigs = int(*p & 0xffffff);
					*p++;
					s_runTime = int(*p & 0xffffff);
					//*p++; //p is at # of SiX4 hits since last eject
					//s_SiX4 = int(*p & 0xffffff);
					
					*p++; //this also iffy, check if it works?
				// Capt Scaler ************************
					if (*p != 0x100dca1e) {
						cout << "trig #" << n_trig << ", Capt Scaler marker not found where expected!" << endl;
						event_good = 0;
						n_bad_events++;
						break;
					}
					*p++; // p is at time since capture in ms
					s_ms_since_capt = int(*p & 0xffffff);
					*p++; // p is at trap state | 0 = trap full, 1 = trap empty
					s_capt_state = int(*p & 0xffffff);
					*p++; // move pointer to eject scaler
					
				// Eject Scaler ***********************
					if (*p != 0x100eca1e) {
						cout << "trig #" << n_trig << ", Eject Scaler marker not found where expected!" << endl;
						event_good = 0;
						n_bad_events++;
						break;
					}
					*p++; //p is at time since eject in ms
					s_ms_since_eject = int(*p & 0xffffff);
					*p++; //p is at # of capt since last eject
					s_capt = int(*p & 0xffffff);
					//*p++; //p is at # of SiX4 hits since last eject
					//s_SiX4 = int(*p & 0xffffff);
					
				} // end new scaler readout
				*/
				/*
			// MCP Pulse-heights correction for data above ADC range:
				//printf(argv[2]);
//				if (!strcmp(argv[2],"alpha")) {
//					if (a_R_mcpA < a_missing_mcp_post) { a_R_mcpA = 6251; na_R_mcpA++; na_R_mcpA_missing++; }
//					if (a_R_mcpB < a_missing_mcp_post) { a_R_mcpB = 6661; na_R_mcpB++; na_R_mcpB_missing++; }
//					if (a_R_mcpC < a_missing_mcp_post) { a_R_mcpC = 5280; na_R_mcpC++; na_R_mcpC_missing++; }
//					if (a_R_mcpD < a_missing_mcp_post) { a_R_mcpD = 6093; na_R_mcpD++; na_R_mcpD_missing++; }
//					if (a_T_mcpC < a_missing_mcp_post) { a_T_mcpC = 5518; na_T_mcpC++; na_T_mcpC_missing++; }
//				}
				*/
			// Reconstruct one missing post
				if (a_R_mcpA_corr < a_missing_mcp_post) {
					bdn.miss_R_mcpA = 1;
					na_R_mcpA_missing++;
					if (a_mcp_lo < a_R_mcpB_corr + a_R_mcpC_corr + a_R_mcpD_corr) a_R_mcpA_corr = a_R_mcpB_corr*a_R_mcpD_corr/a_R_mcpC_corr;
				}
				if (a_R_mcpB_corr < a_missing_mcp_post) {
					bdn.miss_R_mcpB = 1;
					na_R_mcpB_missing++;
					if (a_mcp_lo < a_R_mcpA_corr + a_R_mcpC_corr + a_R_mcpD_corr) a_R_mcpB_corr = a_R_mcpA_corr*a_R_mcpC_corr/a_R_mcpD_corr;
				}
				if (a_R_mcpC_corr < a_missing_mcp_post) {
					bdn.miss_R_mcpC = 1;
					na_R_mcpC_missing++;
					if (a_mcp_lo < a_R_mcpA_corr + a_R_mcpB_corr + a_R_mcpD_corr) a_R_mcpC_corr = a_R_mcpB_corr*a_R_mcpD_corr/a_R_mcpA_corr;
				}
				if (a_R_mcpD_corr < a_missing_mcp_post) {
					bdn.miss_R_mcpD = 1;
					na_R_mcpD_missing++;
					if (a_mcp_lo < a_R_mcpA_corr + a_R_mcpB_corr + a_R_mcpC_corr) a_R_mcpD_corr = a_R_mcpA_corr*a_R_mcpC_corr/a_R_mcpB_corr;
				}
				if (a_T_mcpA_corr < a_missing_mcp_post) {
					bdn.miss_T_mcpA = 1;
					na_T_mcpA_missing++;
					if (a_mcp_lo < a_T_mcpB_corr + a_T_mcpC_corr + a_T_mcpD_corr) a_T_mcpA_corr = a_T_mcpB_corr*a_T_mcpD_corr/a_T_mcpC_corr;
				}
				if (a_T_mcpB_corr < a_missing_mcp_post) {
					bdn.miss_T_mcpB = 1;
					na_T_mcpB_missing++;
					if (a_mcp_lo < a_T_mcpA_corr + a_T_mcpC_corr + a_T_mcpD_corr) a_T_mcpB_corr = a_T_mcpA_corr*a_T_mcpC_corr/a_T_mcpD_corr;
				}
				if (a_T_mcpC_corr < a_missing_mcp_post) {
					bdn.miss_T_mcpC = 1;
					na_T_mcpC_missing++;
					if (a_mcp_lo < a_T_mcpA_corr + a_T_mcpB_corr + a_T_mcpD_corr) a_T_mcpC_corr = a_T_mcpB_corr*a_T_mcpD_corr/a_T_mcpA_corr;
				}
				if (a_T_mcpD_corr < a_missing_mcp_post) {
					bdn.miss_T_mcpD = 1;
					na_T_mcpD_missing++;
					if (a_mcp_lo < a_T_mcpA_corr + a_T_mcpB_corr + a_T_mcpC_corr) a_T_mcpD_corr = a_T_mcpA_corr*a_T_mcpC_corr/a_T_mcpB_corr;
				}
				/*
			// Fill tree with ADC data
				bdn.a_R_ge			= a_R_ge;
				bdn.a_R_ge_highE	= a_R_ge_highE;
				bdn.a_B_dEa			= a_B_dEa;
				bdn.a_B_dEb			= a_B_dEb;
				bdn.a_B_E			= a_B_E;
				bdn.a_T_mcpA		= a_T_mcpA;
				bdn.a_T_mcpB		= a_T_mcpB;
				bdn.a_T_mcpC		= a_T_mcpC;
				bdn.a_T_mcpD		= a_T_mcpD;
				bdn.a_T_mcpE		= a_T_mcpE;
				bdn.a_T_ge			= a_T_ge;
				bdn.a_T_ge_highE	= a_T_ge_highE;
				bdn.a_L_dEa			= a_L_dEa;
				bdn.a_L_dEb			= a_L_dEb;
				bdn.a_L_E			= a_L_E;
				bdn.a_R_mcpA		= a_R_mcpA;
				bdn.a_R_mcpB		= a_R_mcpB;
				bdn.a_R_mcpC		= a_R_mcpC;
				bdn.a_R_mcpD		= a_R_mcpD;
				bdn.a_R_mcpE		= a_R_mcpE;
				*/
				/*
			// Fill tree with TDC data
				bdn.t_T_mcp			= t_T_mcp;
				bdn.t_R_mcp			= t_R_mcp;
				bdn.t_B_dEa			= t_B_dEa;
				bdn.t_B_dEb			= t_B_dEb;
				bdn.t_B_E			= t_B_E;
				bdn.t_L_dEa			= t_L_dEa;
				bdn.t_L_dEb			= t_L_dEb;
				bdn.t_L_E			= t_L_E;
				bdn.t_rf			= t_rf;
				bdn.t_T_ge			= t_T_ge;
				bdn.t_R_ge			= t_R_ge;
				*/
				
			// Fill tree with Scaler data
				bdn.s_ms_since_capt	= s_ms_since_capt;
				bdn.s_capt_state	= s_capt_state;
				bdn.s_ms_since_eject= s_ms_since_eject;
				bdn.s_capt			= s_capt;
				bdn.s_SiX4			= s_SiX4;
				
				tot_liveTime_us += s_liveTime_us;
				tot_runTime_us += s_runTime;
				tot_trigs += all_trigs;
				
				deadTime_us = (s_runTime-s_liveTime_us);
				bdn.deadTime_us	= deadTime_us;
				
				//cout << s_liveTime_us << "\t" << tot_liveTime_us << "\t" << "run time" << "\t";
				//cout << s_runTime << "\t" << tot_runTime_us << endl;
				//cout << "dead time" << "\t" << (s_runTime-s_liveTime_us) << endl;
				
			// Other data
				bdn.rf_phase	= (stBDNCase.dRFFrequencyHz/1000000000)*(1710.0-t_rf-randgen->Rndm());
				h_all_vs_rf_phase_observed->Fill(bdn.rf_phase);
//				ht_rf_phase_observed->Fill(bdn.rf_phase);
				bdn.a_B_dEsum	= a_B_dEa + a_B_dEb;
				bdn.a_L_dEsum	= a_L_dEa + a_L_dEb;
				bdn.t_B_dE		= 0.5*static_cast<double>(t_B_dEa + t_B_dEb);
				bdn.t_L_dE		= 0.5*static_cast<double>(t_L_dEa + t_L_dEb);
				bdn.a_T_mcpSum	= a_T_mcpA + a_T_mcpB + a_T_mcpC + a_T_mcpD;
				ha_T_mcpSum		->Fill(bdn.a_T_mcpSum);
				a_T_mcpSum_corr = a_T_mcpA_corr + a_T_mcpB_corr + a_T_mcpC_corr + a_T_mcpD_corr;
				ha_T_mcpSum_corr->Fill(a_T_mcpSum_corr);
				bdn.T_mcpX		= (a_T_mcpC_corr + a_T_mcpD_corr - a_T_mcpA_corr - a_T_mcpB_corr)/a_T_mcpSum_corr;
				bdn.T_mcpY		= (a_T_mcpA_corr + a_T_mcpD_corr - a_T_mcpB_corr - a_T_mcpC_corr)/a_T_mcpSum_corr;
				bdn.a_R_mcpSum	= a_R_mcpA + a_R_mcpB + a_R_mcpC + a_R_mcpD;
				ha_R_mcpSum		->Fill(bdn.a_R_mcpSum);
				a_R_mcpSum_corr = a_R_mcpA_corr + a_R_mcpB_corr + a_R_mcpC_corr + a_R_mcpD_corr;
				ha_R_mcpSum_corr->Fill(a_R_mcpSum_corr);
				bdn.R_mcpX		= (a_R_mcpC_corr + a_R_mcpD_corr - a_R_mcpA_corr - a_R_mcpB_corr)/a_R_mcpSum_corr;
				bdn.R_mcpY		= (a_R_mcpA_corr + a_R_mcpD_corr - a_R_mcpB_corr - a_R_mcpC_corr)/a_R_mcpSum_corr;
				bdn.event_good	=  event_good;
				bdn.event 		=  n_trig;
				bdn.run 		=  n_run;
				
		// MCP Maps -- Top
				if (a_mcp_lo < a_T_mcpSum_corr) {
					h_T_mcpX->Fill(bdn.T_mcpX);
					h_T_mcpY->Fill(bdn.T_mcpY);
				// 3-post events
					if (bdn.miss_T_mcpA == 1 || bdn.miss_T_mcpB == 1 || bdn.miss_T_mcpC == 1 || bdn.miss_T_mcpD == 1) {
						th = T_mcp_3post_theta;
						x0 = T_mcp_3post_x0;
						y0 = T_mcp_3post_y0;
						a0 = T_mcp_3post_a0;
						b0 = T_mcp_3post_b0;
						a2 = T_mcp_3post_a2;
						b2 = T_mcp_3post_b2;
						xx1 = bdn.T_mcpX;
						yy1 = bdn.T_mcpY;
						xx = (xx1-x0)*Cos(th+Pi()/4) - (yy1-y0)*Sin(th+Pi()/4);
						yy = (xx1-x0)*Sin(th+Pi()/4) + (yy1-y0)*Cos(th+Pi()/4);
						xx1 = xx;
						yy1 = yy;
						xx = (a0 + a2*xx1*xx1) * xx1;
						yy = (b0 + b2*yy1*yy1) * yy1;
						xx1 = xx;
						yy1 = yy;
						bdn.T_mcpPhysX = xx1*Cos(-Pi()/4) - yy1*Sin(-Pi()/4);
						bdn.T_mcpPhysY = xx1*Sin(-Pi()/4) + yy1*Cos(-Pi()/4);
					// Fiducial area cut for 3-post event, uses physical coords
						if ( fid_area_T_mcpPhysX_lo<bdn.T_mcpPhysX && bdn.T_mcpPhysX<fid_area_T_mcpPhysX_hi && fid_area_T_mcpPhysY_lo<bdn.T_mcpPhysY && bdn.T_mcpPhysY<fid_area_T_mcpPhysY_hi ) 
						{
							bdn.fid_area_hit_T_mcp = 1;
						}
						if (s_capt_state == 0) {
							h_T_mcpMapPhys_3post->Fill(bdn.T_mcpPhysX, bdn.T_mcpPhysY);
							if (!strcmp(argv[2],"posts"))		h_T_mcpMapPhys				->Fill(bdn.T_mcpPhysX, bdn.T_mcpPhysY);
							if (bdn.fid_area_hit_T_mcp == 1)	h_T_mcpMapPhysFidArea_3post	->Fill(bdn.T_mcpPhysX, bdn.T_mcpPhysY);
						}
						if (s_capt_state == 1) {
							h_bkgd_T_mcpMapPhys_3post							->Fill(bdn.T_mcpPhysX, bdn.T_mcpPhysY);
							if (!strcmp(argv[2],"posts")) h_bkgd_T_mcpMapPhys	->Fill(bdn.T_mcpPhysX, bdn.T_mcpPhysY);
						}
					}
				// 4-post events
					else {
						bdn.T_mcpPhysX	= T_mcp_a*((bdn.T_mcpX-T_mcp_x0)*Cos(T_mcp_theta)-(bdn.T_mcpY-T_mcp_y0)*Sin(T_mcp_theta));
						bdn.T_mcpPhysY	= T_mcp_a*((bdn.T_mcpX-T_mcp_x0)*Sin(T_mcp_theta)+(bdn.T_mcpY-T_mcp_y0)*Cos(T_mcp_theta));
					// Fiducial area cut for 4-post event, uses raw coords
						if ( fid_area_T_mcpX_lo<bdn.T_mcpX && bdn.T_mcpX<fid_area_T_mcpX_hi && fid_area_T_mcpY_lo<bdn.T_mcpY && bdn.T_mcpY<fid_area_T_mcpY_hi )
						{
							bdn.fid_area_hit_T_mcp = 1;
						}	
						if (s_capt_state == 0) {
							h_T_mcpMap		->Fill(bdn.T_mcpX,		bdn.T_mcpY);
							h_T_mcpMapPhys	->Fill(bdn.T_mcpPhysX,	bdn.T_mcpPhysY);
							if (bdn.fid_area_hit_T_mcp==1) h_T_mcpMapPhysFidArea->Fill(bdn.T_mcpPhysX, bdn.T_mcpPhysY);
						}
						if (s_capt_state == 1) {
							h_bkgd_T_mcpMap		->Fill(bdn.T_mcpX,		bdn.T_mcpY);
							h_bkgd_T_mcpMapPhys	->Fill(bdn.T_mcpPhysX,	bdn.T_mcpPhysY);
						}
					}
				}
		// MCP Maps -- Right
				if (a_mcp_lo < a_R_mcpSum_corr) {
					h_R_mcpX->Fill(bdn.R_mcpX);
					h_R_mcpY->Fill(bdn.R_mcpY);
				// 3-post events
					if (bdn.miss_R_mcpA == 1 || bdn.miss_R_mcpB == 1 || bdn.miss_R_mcpC == 1 || bdn.miss_R_mcpD == 1) {
						th = R_mcp_3post_theta;
						x0 = R_mcp_3post_x0;
						y0 = R_mcp_3post_y0;
						a0 = R_mcp_3post_a0;
						b0 = R_mcp_3post_b0;
						a2 = R_mcp_3post_a2;
						b2 = R_mcp_3post_b2;
						xx1 = bdn.R_mcpX;
						yy1 = bdn.R_mcpY;
						xx = (xx1-x0)*Cos(th+Pi()/4) - (yy1-y0)*Sin(th+Pi()/4);
						yy = (xx1-x0)*Sin(th+Pi()/4) + (yy1-y0)*Cos(th+Pi()/4);
						xx1 = xx;
						yy1 = yy;
						xx = (a0 + a2*xx1*xx1) * xx1;
						yy = (b0 + b2*yy1*yy1) * yy1;
						xx1 = xx;
						yy1 = yy;
						bdn.R_mcpPhysX = xx1*Cos(-Pi()/4) - yy1*Sin(-Pi()/4);
						bdn.R_mcpPhysY = xx1*Sin(-Pi()/4) + yy1*Cos(-Pi()/4);
					// Fiducial area cut for 3-post event, uses physical coords
						if ( fid_area_R_mcpPhysX_lo<bdn.R_mcpPhysX && bdn.R_mcpPhysX<fid_area_R_mcpPhysX_hi && fid_area_R_mcpPhysY_lo<bdn.R_mcpPhysY && bdn.R_mcpPhysY<fid_area_R_mcpPhysY_hi ) 
						{
							bdn.fid_area_hit_R_mcp = 1;
						}
						if (s_capt_state == 0) {
							h_R_mcpMapPhys_3post->Fill(bdn.R_mcpPhysX, bdn.R_mcpPhysY);
							if (!strcmp(argv[2],"posts"))		h_R_mcpMapPhys				->Fill(bdn.R_mcpPhysX, bdn.R_mcpPhysY);
							if (bdn.fid_area_hit_R_mcp == 1)	h_R_mcpMapPhysFidArea_3post	->Fill(bdn.R_mcpPhysX, bdn.R_mcpPhysY);
						}
						if (s_capt_state == 1) {
							h_bkgd_R_mcpMapPhys_3post							->Fill(bdn.R_mcpPhysX, bdn.R_mcpPhysY);
							if (!strcmp(argv[2],"posts")) h_bkgd_R_mcpMapPhys	->Fill(bdn.R_mcpPhysX, bdn.R_mcpPhysY);
						}
					}
				// 4-post events
					else {
						bdn.R_mcpPhysX	= R_mcp_a*((bdn.R_mcpX-R_mcp_x0)*Cos(R_mcp_theta)-(bdn.R_mcpY-R_mcp_y0)*Sin(R_mcp_theta));
						bdn.R_mcpPhysY	= R_mcp_a*((bdn.R_mcpX-R_mcp_x0)*Sin(R_mcp_theta)+(bdn.R_mcpY-R_mcp_y0)*Cos(R_mcp_theta));
						// Fiducial area cut for 4-post event, uses raw coords
						if ( fid_area_R_mcpX_lo<bdn.R_mcpX && bdn.R_mcpX<fid_area_R_mcpX_hi && fid_area_R_mcpY_lo<bdn.R_mcpY && bdn.R_mcpY<fid_area_R_mcpY_hi )
						{
							bdn.fid_area_hit_R_mcp = 1;
						}	
						if (s_capt_state == 0) {
							h_R_mcpMap		->Fill(bdn.R_mcpX,		bdn.R_mcpY);
							h_R_mcpMapPhys	->Fill(bdn.R_mcpPhysX,	bdn.R_mcpPhysY);
							if (bdn.fid_area_hit_R_mcp==1) h_R_mcpMapPhysFidArea->Fill(bdn.R_mcpPhysX, bdn.R_mcpPhysY);
						}
						if (s_capt_state == 1) {
							h_bkgd_R_mcpMap		->Fill(bdn.R_mcpX,		bdn.R_mcpY);
							h_bkgd_R_mcpMapPhys	->Fill(bdn.R_mcpPhysX,	bdn.R_mcpPhysY);
						}
					}
				}
				
				ha_B_dEsum->Fill(a_B_dEa + a_B_dEb);
				ha_L_dEsum->Fill(a_L_dEa + a_L_dEb);
				if (t_B_dEa > -100 && t_B_dEb > -100) {
					ht_B_dEdiff	->Fill(t_B_dEa-t_B_dEb);
					ht_B_dE		->Fill(bdn.t_B_dE);
					ht_B_dEmin	->Fill(TMath::Min(t_B_dEa,t_B_dEb));
					if (t_B_E > -100) ht_B_dE_E->Fill(t_B_E-bdn.t_B_dE);
				}
				if (t_L_dEa > -100 && t_L_dEb > -100) {
					ht_L_dEdiff	->Fill(t_L_dEa-t_L_dEb);
					ht_L_dE		->Fill(bdn.t_L_dE);
					ht_L_dEmin	->Fill(TMath::Min(t_L_dEa,t_L_dEb));
					if (t_L_E > -100) ht_L_dE_E->Fill(t_L_E-bdn.t_L_dE);
				}
				
				bdn_Tree->Fill();
				
//////////////////////////////////////////////////////////////////////////////////////////				
// Conditional filling:
	// EventLists, Trees, and Histograms
//////////////////////////////////////////////////////////////////////////////////////////		
				
// Beta-Recoil Events
				if ((event_good==1 && t_dE_lo<bdn.t_L_dE && t_mcp_lo<t_T_mcp) ||
					(event_good==1 && t_dE_lo<bdn.t_L_dE && t_mcp_lo<t_R_mcp) ||
					(event_good==1 && t_dE_lo<bdn.t_B_dE && t_mcp_lo<t_T_mcp) ||
					(event_good==1 && t_dE_lo<bdn.t_B_dE && t_mcp_lo<t_R_mcp))
				{
// LT
				  //if (event_good==1 &&      t_dE_lo<bdn.t_L_dE && bdn.t_L_dE<t_dE_hi      && a_dE_lo<bdn.a_L_dEsum &&     t_mcp_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum) {
					if (event_good==1 && t_trigger_lo<bdn.t_L_dE && bdn.t_L_dE<t_trigger_hi && a_dE_lo<bdn.a_L_dEsum && t_trigger_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum && bdn.fid_area_hit_T_mcp==1) {
						bdn.tof_LT	= t_T_mcp - bdn.t_L_dE - LT_zeroTime[0];
						t1			= 0.001 * tofToMCPGrid (stBDNCase, 'T', bdn.tof_LT); // need times in us
						z1			= stBDNCase.dTopGridDistance;
						t2			= 0.001 * (bdn.tof_LT - 0.5 + randgen->Rndm()); // need times in us
						x2			= bdn.T_mcpPhysX;
						y2			= bdn.T_mcpPhysY;
						s2			= Sqrt(x2*x2 + y2*y2);
						vs			= s2/t2;
						vz			= z1/t1;//stBDNCase.dTopMCPDistance/t2;//z1/t1;
						bdn.v_LT	= Sqrt(vs*vs + vz*vz);
						bdn.En_LT	= 0.5 * stBDNCase.dNeutronEnergyMassFactorKeV * bdn.v_LT * bdn.v_LT;
						if (s_capt_state == 0) {
							h_tof		->Fill(bdn.tof_LT);
							h_tof_LT	->Fill(bdn.tof_LT);
							h_v			->Fill(bdn.v_LT);
							h_v_LT		->Fill(bdn.v_LT);
							h_vInv		->Fill(1.0/(bdn.v_LT+0.000000001));
							h_vInv_LT	->Fill(1.0/(bdn.v_LT+0.000000001));
							h_En		->Fill(bdn.En_LT);
							h_En_LT		->Fill(bdn.En_LT);
						}
						if (s_capt_state == 1) {
							h_bkgd_tof		->Fill(bdn.tof_LT);
							h_bkgd_tof_LT	->Fill(bdn.tof_LT);
							h_bkgd_v		->Fill(bdn.v_LT);
							h_bkgd_v_LT		->Fill(bdn.v_LT);
							h_bkgd_vInv		->Fill(1.0/(bdn.v_LT+0.000000001));
							h_bkgd_vInv_LT	->Fill(1.0/(bdn.v_LT+0.000000001));
							h_bkgd_En		->Fill(bdn.En_LT);
							h_bkgd_En_LT	->Fill(bdn.En_LT);
						}
						if ( tof_zero_lo < bdn.tof_LT && bdn.tof_LT < tof_zero_hi) {
							h_LT_zero_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_zero_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) metadata.nZeroTOFCount[LT]		+= 1.0;
							if (s_capt_state == 1) metadata.nZeroTOFBkgdCount[LT]	+= 1.0;
						}
						if ( tof_lowTOF_lo < bdn.tof_LT && bdn.tof_LT < tof_lowTOF_hi) {
							h_LT_lowTOF_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_lowTOF_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) metadata.nLowTOFCount[LT]		+= 1.0;
							if (s_capt_state == 1) metadata.nLowTOFBkgdCount[LT]	+= 1.0;
						}
						if ( tof_T_fast_lo < bdn.tof_LT && bdn.tof_LT < tof_T_fast_hi) {
							h_LT_fast_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_fast_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) metadata.nFastCount[LT]		+= 1.0;
							if (s_capt_state == 1) metadata.nFastBkgdCount[LT]	+= 1.0;
						}
						if ( tof_T_slow_lo < bdn.tof_LT && bdn.tof_LT < tof_T_slow_hi) {
							h_LT_slow_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_slow_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) {
								   h_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								h_LT_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								metadata.nSlowCount[LT]		+= 1.0;
							}
							if (s_capt_state == 1) {
								   h_bkgd_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								h_bkgd_LT_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								metadata.nSlowBkgdCount[LT]	+= 1.0;
							}
						}
						if ( tof_oops_lo < bdn.tof_LT && bdn.tof_LT < tof_oops_hi) {
							h_LT_oops_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_oops_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
							 if (s_capt_state == 0) {
								   h_oops_vs_rf_phase_observed->Fill(bdn.rf_phase);
								h_LT_oops_vs_rf_phase_observed->Fill(bdn.rf_phase);
								metadata.nOopsCount[LT]		+= 1.0;
							}
							if (s_capt_state == 1)
								metadata.nOopsBkgdCount[LT]	+= 1.0;
						}
					}
// LR
				  //if (event_good==1 &&      t_dE_lo<bdn.t_L_dE && bdn.t_L_dE<t_dE_hi      && a_dE_lo<bdn.a_L_dEsum &&     t_mcp_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum) {
					if (event_good==1 && t_trigger_lo<bdn.t_L_dE && bdn.t_L_dE<t_trigger_hi && a_dE_lo<bdn.a_L_dEsum && t_trigger_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum && bdn.fid_area_hit_R_mcp==1) {
						bdn.tof_LR	= t_R_mcp - bdn.t_L_dE - LR_zeroTime[0];
						t1			= 0.001 * tofToMCPGrid (stBDNCase, 'R', bdn.tof_LR); // need times in us
						z1			= stBDNCase.dRightGridDistance;
						t2			= 0.001 * (bdn.tof_LR - 0.5 + randgen->Rndm()); // need times in us
						x2			= bdn.R_mcpPhysX;
						y2			= bdn.R_mcpPhysY;
						s2			= Sqrt(x2*x2 + y2*y2);
						vs			= s2/t2;
						vz			= z1/t1;//stBDNCase.dRightMCPDistance/t2;//z1/t1;
						bdn.v_LR	= Sqrt(vs*vs + vz*vz);
						bdn.En_LR	= 0.5 * stBDNCase.dNeutronEnergyMassFactorKeV * bdn.v_LR * bdn.v_LR;
						if (s_capt_state == 0) {
							h_tof		->Fill(bdn.tof_LR);
							h_tof_LR	->Fill(bdn.tof_LR);
							h_v			->Fill(bdn.v_LR);
							h_v_LR		->Fill(bdn.v_LR);
							h_vInv		->Fill(1.0/(bdn.v_LR+0.000000001));
							h_vInv_LR	->Fill(1.0/(bdn.v_LR+0.000000001));
							h_En		->Fill(bdn.En_LR);
							h_En_LR		->Fill(bdn.En_LR);
						}
						if (s_capt_state == 1) {
							h_bkgd_tof		->Fill(bdn.tof_LR);
							h_bkgd_tof_LR	->Fill(bdn.tof_LR);
							h_bkgd_v		->Fill(bdn.v_LR);
							h_bkgd_v_LR		->Fill(bdn.v_LR);
							h_bkgd_vInv		->Fill(1.0/(bdn.v_LR+0.000000001));
							h_bkgd_vInv_LR	->Fill(1.0/(bdn.v_LR+0.000000001));
							h_bkgd_En		->Fill(bdn.En_LR);
							h_bkgd_En_LR	->Fill(bdn.En_LR);
						}
						if ( tof_zero_lo < bdn.tof_LR && bdn.tof_LR < tof_zero_hi) {
							h_LR_zero_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_zero_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) metadata.nZeroTOFCount[LR]		+= 1.0;
							if (s_capt_state == 1) metadata.nZeroTOFBkgdCount[LR]	+= 1.0;
						}
						if ( tof_lowTOF_lo < bdn.tof_LR && bdn.tof_LR < tof_lowTOF_hi) {
							h_LR_lowTOF_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_lowTOF_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) metadata.nLowTOFCount[LR]		+= 1.0;
							if (s_capt_state == 1) metadata.nLowTOFBkgdCount[LR]	+= 1.0;
						}
						if ( tof_R_fast_lo < bdn.tof_LR && bdn.tof_LR < tof_R_fast_hi) {
							h_LR_fast_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_fast_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) metadata.nFastCount[LR]		+= 1.0;
							if (s_capt_state == 1) metadata.nFastBkgdCount[LR]	+= 1.0;
						}
						if ( tof_R_slow_lo < bdn.tof_LR && bdn.tof_LR < tof_R_slow_hi) {
							h_LR_slow_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_slow_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) {
								   h_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								h_LR_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								metadata.nSlowCount[LR]		+= 1.0;
							}
							if (s_capt_state == 1) {
								   h_bkgd_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								h_bkgd_LR_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								metadata.nSlowBkgdCount[LR]	+= 1.0;
							}
						}
						if ( tof_oops_lo < bdn.tof_LR && bdn.tof_LR < tof_oops_hi) {
							h_LR_oops_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_oops_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
							 if (s_capt_state == 0) {
								   h_oops_vs_rf_phase_observed->Fill(bdn.rf_phase);
								h_LR_oops_vs_rf_phase_observed->Fill(bdn.rf_phase);
								metadata.nOopsCount[LR]		+= 1.0;
							}
							if (s_capt_state == 1)
								metadata.nOopsBkgdCount[LR]	+= 1.0;
						}
					}
// BTbb
				  //if (event_good==1 &&      t_dE_lo<bdn.t_B_dE && bdn.t_B_dE<t_dE_hi      && a_dE_lo<bdn.a_B_dEsum &&     t_mcp_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum) {
					if (event_good==1 && t_trigger_lo<bdn.t_B_dE && bdn.t_B_dE<t_trigger_hi && a_dE_lo<bdn.a_B_dEsum && t_trigger_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum && bdn.fid_area_hit_T_mcp==1) {
						bdn.tof_BT	= t_T_mcp - bdn.t_B_dE - BT_zeroTime[0];
						t1			= 0.001 * tofToMCPGrid (stBDNCase, 'T', bdn.tof_BT); // need times in us
						z1			= stBDNCase.dTopGridDistance;
						t2			= 0.001 * (bdn.tof_BT - 0.5 + randgen->Rndm()); // need times in us
						x2			= bdn.T_mcpPhysX;
						y2			= bdn.T_mcpPhysY;
						s2			= Sqrt(x2*x2 + y2*y2);
						vs			= s2/t2;
						vz			= z1/t1;//stBDNCase.dTopMCPDistance/t2;//z1/t1;
						bdn.v_BT	= Sqrt(vs*vs + vz*vz);
						bdn.En_BT	= 0.5 * stBDNCase.dNeutronEnergyMassFactorKeV * bdn.v_BT * bdn.v_BT;
						if (s_capt_state == 0) {
							h_tof		->Fill(bdn.tof_BT);
							h_tof_BT	->Fill(bdn.tof_BT);
							h_v			->Fill(bdn.v_BT);
							h_v_BT		->Fill(bdn.v_BT);
							h_vInv		->Fill(1.0/(bdn.v_BT+0.000000001));
							h_vInv_BT	->Fill(1.0/(bdn.v_BT+0.000000001));
							h_En		->Fill(bdn.En_BT);
							h_En_BT		->Fill(bdn.En_BT);
						}
						if (s_capt_state == 1) {
							h_bkgd_tof		->Fill(bdn.tof_BT);
							h_bkgd_tof_BT	->Fill(bdn.tof_BT);
							h_bkgd_v		->Fill(bdn.v_BT);
							h_bkgd_v_BT		->Fill(bdn.v_BT);
							h_bkgd_vInv		->Fill(1.0/(bdn.v_BT+0.000000001));
							h_bkgd_vInv_BT	->Fill(1.0/(bdn.v_BT+0.000000001));
							h_bkgd_En		->Fill(bdn.En_BT);
							h_bkgd_En_BT	->Fill(bdn.En_BT);
						}
						if ( tof_zero_lo < bdn.tof_BT && bdn.tof_BT < tof_zero_hi) {
							h_BT_zero_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_zero_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) metadata.nZeroTOFCount[BT]		+= 1.0;
							if (s_capt_state == 1) metadata.nZeroTOFBkgdCount[BT]	+= 1.0;
						}
						if ( tof_lowTOF_lo < bdn.tof_BT && bdn.tof_BT < tof_lowTOF_hi) {
							h_BT_lowTOF_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_lowTOF_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) metadata.nLowTOFCount[BT]		+= 1.0;
							if (s_capt_state == 1) metadata.nLowTOFBkgdCount[BT]	+= 1.0;
						}
						if ( tof_T_fast_lo < bdn.tof_BT && bdn.tof_BT < tof_T_fast_hi) {
							h_BT_fast_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_fast_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) metadata.nFastCount[BT]		+= 1.0;
							if (s_capt_state == 1) metadata.nFastBkgdCount[BT]	+= 1.0;
						}
						if ( tof_T_slow_lo < bdn.tof_BT && bdn.tof_BT < tof_T_slow_hi) {
							h_BT_slow_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_slow_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) {
								   h_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								h_BT_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								metadata.nSlowCount[BT]		+= 1.0;
							}
							if (s_capt_state == 1) {
								   h_bkgd_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								h_bkgd_BT_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								metadata.nSlowBkgdCount[BT]	+= 1.0;
							}
						}
						if ( tof_oops_lo < bdn.tof_BT && bdn.tof_BT < tof_oops_hi) {
							h_BT_oops_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_oops_mcpMap->Fill(bdn.T_mcpX,bdn.T_mcpY);
							 h_T_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
							 if (s_capt_state == 0) {
								   h_oops_vs_rf_phase_observed->Fill(bdn.rf_phase);
								h_BT_oops_vs_rf_phase_observed->Fill(bdn.rf_phase);
								metadata.nOopsCount[BT]		+= 1.0;
							}
							if (s_capt_state == 1)
								metadata.nOopsBkgdCount[BT]	+= 1.0;
						}
					}
// BR
				  //if (event_good==1 &&      t_dE_lo<bdn.t_B_dE && bdn.t_B_dE<t_dE_hi      && a_dE_lo<bdn.a_B_dEsum &&     t_mcp_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum) {
					if (event_good==1 && t_trigger_lo<bdn.t_B_dE && bdn.t_B_dE<t_trigger_hi && a_dE_lo<bdn.a_B_dEsum && t_trigger_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum && bdn.fid_area_hit_R_mcp==1) {
//						if (s_capt_state == 0){      h_tof->Fill(bdn.tof_BR);      h_tof_BR->Fill(bdn.tof_BR); }
//						if (s_capt_state == 1){ h_bkgd_tof->Fill(bdn.tof_BR); h_bkgd_tof_BR->Fill(bdn.tof_BR); }
						bdn.tof_BR	= t_R_mcp - bdn.t_B_dE - BR_zeroTime[0];
						t1			= 0.001 * tofToMCPGrid (stBDNCase, 'R', bdn.tof_BR); // need times in us
						z1			= stBDNCase.dRightGridDistance;
						t2			= 0.001 * (bdn.tof_BR - 0.5 + randgen->Rndm()); // need times in us
						x2			= bdn.R_mcpPhysX;
						y2			= bdn.R_mcpPhysY;
						s2			= Sqrt(x2*x2 + y2*y2);
						vs			= s2/t2;
						vz			= z1/t1;//stBDNCase.dRightMCPDistance/t2;//z1/t1;
//						printf("miss(%d,%d,%d,%d); (x, y) = (%f, %f); s = %f; r = %f; t1 = %f; t2 = %f; (vs, vz) = (%f, %f); v = %f\n", bdn.miss_R_mcpA, bdn.miss_R_mcpB, bdn.miss_R_mcpC, bdn.miss_R_mcpD, x2, y2, s2, Sqrt(s2*s2+Power(stBDNCase.dRightMCPDistance,2)), t1, t2, vs, vz);
						bdn.v_BR	= Sqrt(vs*vs + vz*vz);
						bdn.En_BR	= 0.5 * stBDNCase.dNeutronEnergyMassFactorKeV * bdn.v_BR * bdn.v_BR;
						if (s_capt_state == 0) {
							h_tof		->Fill(bdn.tof_BR);
							h_tof_BR	->Fill(bdn.tof_BR);
							h_v			->Fill(bdn.v_BR);
							h_v_BR		->Fill(bdn.v_BR);
							h_vInv		->Fill(1.0/(bdn.v_BR+0.000000001));
							h_vInv_BR	->Fill(1.0/(bdn.v_BR+0.000000001));
							h_En		->Fill(bdn.En_BR);
							h_En_BR		->Fill(bdn.En_BR);
						}
						if (s_capt_state == 1) {
							h_bkgd_tof		->Fill(bdn.tof_BR);
							h_bkgd_tof_BR	->Fill(bdn.tof_BR);
							h_bkgd_v		->Fill(bdn.v_BR);
							h_bkgd_v_BR		->Fill(bdn.v_BR);
							h_bkgd_vInv		->Fill(1.0/(bdn.v_BR+0.000000001));
							h_bkgd_vInv_BR	->Fill(1.0/(bdn.v_BR+0.000000001));
							h_bkgd_En		->Fill(bdn.En_BR);
							h_bkgd_En_BR	->Fill(bdn.En_BR);
						}
						if ( tof_zero_lo < bdn.tof_BR && bdn.tof_BR < tof_zero_hi) {
							h_BR_zero_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_zero_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) metadata.nZeroTOFCount[BR]		+= 1.0;
							if (s_capt_state == 1) metadata.nZeroTOFBkgdCount[BR]	+= 1.0;
						}
						if ( tof_lowTOF_lo < bdn.tof_BR && bdn.tof_BR < tof_lowTOF_hi) {
							h_BR_lowTOF_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_lowTOF_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) metadata.nLowTOFCount[BR]		+= 1.0;
							if (s_capt_state == 1) metadata.nLowTOFBkgdCount[BR]	+= 1.0;
						}
						if ( tof_R_fast_lo < bdn.tof_BR && bdn.tof_BR < tof_R_fast_hi) {
							h_BR_fast_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_fast_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) metadata.nFastCount[BR]		+= 1.0;
							if (s_capt_state == 1) metadata.nFastBkgdCount[BR]	+= 1.0;
						}
						if ( tof_R_slow_lo < bdn.tof_BR && bdn.tof_BR < tof_R_slow_hi) {
							h_BR_slow_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_slow_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
							if (s_capt_state == 0) {
								   h_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								h_BR_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								metadata.nSlowCount[BR]		+= 1.0;
							}
							if (s_capt_state == 1) {
								   h_bkgd_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								h_bkgd_BR_slow_vs_rf_phase_observed->Fill(bdn.rf_phase);
								metadata.nSlowBkgdCount[BR]	+= 1.0;
							}
						}
						if ( tof_oops_lo < bdn.tof_BR && bdn.tof_BR < tof_oops_hi) {
							h_BR_oops_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_oops_mcpMap->Fill(bdn.R_mcpX,bdn.R_mcpY);
							 h_R_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
							   h_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
							 if (s_capt_state == 0) {
								   h_oops_vs_rf_phase_observed->Fill(bdn.rf_phase);
								h_BR_oops_vs_rf_phase_observed->Fill(bdn.rf_phase);
								metadata.nOopsCount[BR]		+= 1.0;
							}
							if (s_capt_state == 1)
								metadata.nOopsBkgdCount[BR]	+= 1.0;
						}
					}
					
					beta_recoil_tree->Fill();
					
				} // end Beta-Recoil events
				
// vs Cycle Time (others that require tof are below in the TOF section)
			h_state_vs_cycle_time->Fill(s_ms_since_eject,s_capt_state);
			h_all_vs_cycle_time_observed->Fill(s_ms_since_eject);
			if (t_dE_lo<bdn.t_L_dE && a_dE_lo<bdn.a_L_dEsum) {
				h_L_betas_vs_cycle_time_observed->Fill(s_ms_since_eject);
				h_betas_vs_cycle_time_observed->Fill(s_ms_since_eject);
			}
			if (t_dE_lo<bdn.t_B_dE && a_dE_lo<bdn.a_B_dEsum) {
				h_B_betas_vs_cycle_time_observed->Fill(s_ms_since_eject);
				h_betas_vs_cycle_time_observed->Fill(s_ms_since_eject);
			}

// Ge singles: with ions
			if (event_good==1 && s_capt_state == 0) {
				ha_sgnl_R_ge->Fill(a_R_ge);
				ha_sgnl_T_ge->Fill(a_T_ge);
				
				he_sgnl_R_ge->Fill(e_R_ge);
				he_sgnl_T_ge->Fill(e_T_ge);
				he_sgnl_ge	->Fill(e_R_ge);
				he_sgnl_ge	->Fill(e_T_ge);
				//he_sgnl_R_ge_highE->Fill(e_R_ge_highE);
				//he_sgnl_T_ge_highE->Fill(e_T_ge_highE);
				//he_sgnl_ge_highE	->Fill(e_R_ge_highE);
				//he_sgnl_ge_highE	->Fill(e_T_ge_highE);
			}
// Ge singles: no ions
			if (event_good==1 && s_capt_state == 1) {
				ha_bkgd_R_ge->Fill(a_R_ge);
				ha_bkgd_T_ge->Fill(a_T_ge);
				
				he_bkgd_R_ge->Fill(e_R_ge);
				he_bkgd_T_ge->Fill(e_T_ge);
				he_bkgd_ge	->Fill(e_R_ge);
				he_bkgd_ge	->Fill(e_T_ge);
				
				//he_bkgd_R_ge_highE->Fill(e_R_ge_highE);
				//he_bkgd_T_ge_highE->Fill(e_T_ge_highE);
				//he_bkgd_ge_highE	->Fill(e_R_ge_highE);
				//he_bkgd_ge_highE	->Fill(e_T_ge_highE);
			}
			
// Beta-Gamma
			if ((event_good==1 && t_trigger_lo<bdn.t_L_dE && t_trigger_lo<t_T_ge) ||
				(event_good==1 && t_trigger_lo<bdn.t_L_dE && t_trigger_lo<t_R_ge) ||
				(event_good==1 && t_trigger_lo<bdn.t_B_dE && t_trigger_lo<t_T_ge) ||
				(event_good==1 && t_trigger_lo<bdn.t_B_dE && t_trigger_lo<t_R_ge))
			{
				beta_gamma_tree->Fill();
			}
			// old cut 2013-12-02:
			//if (event_good==1 && s_capt_state==0 && t_E_lo<t_L_E && t_dE_lo<bdn.t_L_dE && a_dE_lo<bdn.a_L_dEsum && a_E_lo<a_L_E && 0 < (t_T_ge-bdn.t_L_dE) && (t_T_ge-bdn.t_L_dE) < 1000) {
			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_L_dE && a_dE_lo<bdn.a_L_dEsum && 0 < (t_T_ge-bdn.t_L_dE) && (t_T_ge-bdn.t_L_dE) < 1000) {
				ha_bg_LT		->Fill(a_T_ge); // ADC
				he_bg_LT		->Fill(e_T_ge); // keV
				he_bg			->Fill(e_T_ge); // keV
				ha_sgnl_bg_LT	->Fill(a_T_ge); // ADC
				he_sgnl_bg_LT	->Fill(e_T_ge); // keV
				he_sgnl_bg		->Fill(e_T_ge); // keV
			}
			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_L_dE && a_dE_lo<bdn.a_L_dEsum && 0 < (t_R_ge-bdn.t_L_dE) && (t_R_ge-bdn.t_L_dE) < 1000) {
				ha_bg_LR		->Fill(a_R_ge); // ADC
				he_bg_LR		->Fill(e_R_ge); // keV
				he_bg			->Fill(e_R_ge); // keV
				ha_sgnl_bg_LR	->Fill(a_R_ge); // ADC
				he_sgnl_bg_LR	->Fill(e_R_ge); // keV
				he_sgnl_bg		->Fill(e_R_ge); // keV
			}
			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_B_dE && a_dE_lo<bdn.a_B_dEsum && 0 < (t_T_ge-bdn.t_B_dE) && (t_T_ge-bdn.t_B_dE) < 1000) {
				ha_bg_BT		->Fill(a_T_ge); // ADC
				he_bg_BT		->Fill(e_T_ge); // keV
				he_bg			->Fill(e_T_ge); // keV
				ha_sgnl_bg_BT	->Fill(a_T_ge); // ADC
				he_sgnl_bg_BT	->Fill(e_T_ge); // keV
				he_sgnl_bg		->Fill(e_T_ge); // keV
			}
			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_B_dE && a_dE_lo<bdn.a_B_dEsum && 0 < (t_R_ge-bdn.t_B_dE) && (t_R_ge-bdn.t_B_dE) < 1000) {
				ha_bg_BR		->Fill(a_R_ge); // ADC
				he_bg_BR		->Fill(e_R_ge); // keV
				he_bg			->Fill(e_R_ge); // keV
				ha_sgnl_bg_BR	->Fill(a_R_ge); // ADC
				he_sgnl_bg_BR	->Fill(e_R_ge); // keV
				he_sgnl_bg		->Fill(e_R_ge); // keV
			}
			if (event_good==1 && s_capt_state==1 && t_dE_lo<bdn.t_L_dE && a_dE_lo<bdn.a_L_dEsum && 0 < (t_T_ge-bdn.t_L_dE) && (t_T_ge-bdn.t_L_dE) < 1000) {
				ha_bg_LT		->Fill(a_T_ge); // ADC
				he_bg_LT		->Fill(e_T_ge); // keV
				he_bg			->Fill(e_T_ge); // keV
				ha_bkgd_bg_LT	->Fill(a_T_ge); // ADC
				he_bkgd_bg_LT	->Fill(e_T_ge); // keV
				he_bkgd_bg		->Fill(e_T_ge); // keV
			}
			if (event_good==1 && s_capt_state==1 && t_dE_lo<bdn.t_L_dE && a_dE_lo<bdn.a_L_dEsum && 0 < (t_R_ge-bdn.t_L_dE) && (t_R_ge-bdn.t_L_dE) < 1000) {
				ha_bg_LR		->Fill(a_R_ge); // ADC
				he_bg_LR		->Fill(e_R_ge); // keV
				he_bg			->Fill(e_R_ge); // keV
				ha_bkgd_bg_LR	->Fill(a_R_ge); // ADC
				he_bkgd_bg_LR	->Fill(e_R_ge); // keV
				he_bkgd_bg		->Fill(e_R_ge); // keV
			}
			if (event_good==1 && s_capt_state==1 && t_dE_lo<bdn.t_B_dE && a_dE_lo<bdn.a_B_dEsum && 0 < (t_T_ge-bdn.t_B_dE) && (t_T_ge-bdn.t_B_dE) < 1000) {
				ha_bg_BT		->Fill(a_T_ge); // ADC
				he_bg_BT		->Fill(e_T_ge); // keV
				he_bg			->Fill(e_T_ge); // keV
				ha_bkgd_bg_BT	->Fill(a_T_ge); // ADC
				he_bkgd_bg_BT	->Fill(e_T_ge); // keV
				he_bkgd_bg		->Fill(e_T_ge); // keV
			}
			if (event_good==1 && s_capt_state==1 && t_dE_lo<bdn.t_B_dE && a_dE_lo<bdn.a_B_dEsum && 0 < (t_R_ge-bdn.t_B_dE) && (t_R_ge-bdn.t_B_dE) < 1000) {
				ha_bg_BR		->Fill(a_R_ge); // ADC
				he_bg_BR		->Fill(e_R_ge); // keV
				he_bg			->Fill(e_R_ge); // keV
				ha_bkgd_bg_BR	->Fill(a_R_ge); // ADC
				he_bkgd_bg_BR	->Fill(e_R_ge); // keV
				he_bkgd_bg		->Fill(e_R_ge); // keV
			}
/*			
//// LT-TOF
//			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_L_dE && bdn.t_L_dE<t_dE_hi && a_dE_lo<bdn.a_L_dEsum && t_mcp_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum) {
//				tof = t_T_mcp - bdn.t_L_dE - LT_zeroTime[0];
//				h_tof_LT	-> Fill(tof);
//				h_tof		-> Fill(tof);
//				if ( tof_zero_lo < tof && tof < tof_zero_hi) {
//					h_T_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_lowTOF_lo < tof && tof < tof_lowTOF_hi) {
//					h_T_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_fast_lo < tof && tof < tof_fast_hi) {
//					h_T_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_slow_lo < tof && tof < tof_slow_hi) {
//					h_T_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_oops_lo < tof && tof < tof_oops_hi) {
//					h_T_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//			}
//// LR-TOF
//			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_L_dE && bdn.t_L_dE<t_dE_hi && a_dE_lo<bdn.a_L_dEsum && t_mcp_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum) {
//				tof = t_R_mcp - bdn.t_L_dE - LR_zeroTime[0];
//				h_tof_LR	-> Fill(tof);
//				h_tof		-> Fill(tof);
//				if ( tof_zero_lo < tof && tof < tof_zero_hi) {
//					h_R_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_lowTOF_lo < tof && tof < tof_lowTOF_hi) {
//					h_R_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_fast_lo < tof && tof < tof_fast_hi) {
//					h_R_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_slow_lo < tof && tof < tof_slow_hi) {
//					h_R_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_oops_lo < tof && tof < tof_oops_hi) {
//					h_R_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//			}
//// BT-TOF
//			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_B_dE && bdn.t_B_dE<t_dE_hi && a_dE_lo<bdn.a_B_dEsum && t_mcp_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum) {
//				tof = t_T_mcp - bdn.t_B_dE - BT_zeroTime[0];
//				h_tof_BT	-> Fill(tof);
//				h_tof		-> Fill(tof);
//				if ( tof_zero_lo < tof && tof < tof_zero_hi) {
//					h_T_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_lowTOF_lo < tof && tof < tof_lowTOF_hi) {
//					h_T_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_fast_lo < tof && tof < tof_fast_hi) {
//					h_T_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_slow_lo < tof && tof < tof_slow_hi) {
//					h_T_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_oops_lo < tof && tof < tof_oops_hi) {
//					h_T_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//			}
//// BR-TOF
//			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_B_dE && bdn.t_B_dE<t_dE_hi && a_dE_lo<bdn.a_B_dEsum && t_mcp_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum) {
//				tof = t_R_mcp - bdn.t_B_dE - BR_zeroTime[0];
//				h_tof_BR	-> Fill(tof);
//				h_tof		-> Fill(tof);
//				if ( tof_zero_lo < tof && tof < tof_zero_hi) {
//					h_R_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_zero_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_lowTOF_lo < tof && tof < tof_lowTOF_hi) {
//					h_R_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_lowTOF_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_fast_lo < tof && tof < tof_fast_hi) {
//					h_R_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_fast_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_slow_lo < tof && tof < tof_slow_hi) {
//					h_R_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_slow_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//				if ( tof_oops_lo < tof && tof < tof_oops_hi) {
//					h_R_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
//					h_oops_vs_cycle_time_observed->Fill(s_ms_since_eject);
//				}
//			}
//// Background LT-TOF
//			if (event_good==1 && s_capt_state==1 && t_dE_lo<bdn.t_L_dE && bdn.t_L_dE<t_dE_hi && a_dE_lo<bdn.a_L_dEsum && t_mcp_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum) {
//				tof = t_T_mcp - bdn.t_L_dE - LT_zeroTime[0];
//				h_bkgd_tof_LT	-> Fill(tof);
//				h_bkgd_tof		-> Fill(tof);
//			}
//// Background LR-TOF
//			if (event_good==1 && s_capt_state==1 && t_dE_lo<bdn.t_L_dE && bdn.t_L_dE<t_dE_hi && a_dE_lo<bdn.a_L_dEsum && t_mcp_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum) {
//				tof = t_R_mcp - bdn.t_L_dE - LR_zeroTime[0];
//				h_bkgd_tof_LR	-> Fill(tof);
//				h_bkgd_tof		-> Fill(tof);
//			}
//// Background BT-TOF
//			if (event_good==1 && s_capt_state==1 && t_dE_lo<bdn.t_B_dE && bdn.t_B_dE<t_dE_hi && a_dE_lo<bdn.a_B_dEsum && t_mcp_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum) {
//				tof = t_T_mcp - bdn.t_B_dE - BT_zeroTime[0];
//				h_bkgd_tof_BT	-> Fill(tof);
//				h_bkgd_tof		-> Fill(tof);
//			// Was used for looking at coinc. gammas during 2013 setup:
//				h_bkgd_tof_dEmin-> Fill(t_T_mcp-TMath::Min(t_B_dEa,t_B_dEb));
//				//h_bkgd_tof_dEmin-> Fill(t_T_mcp-t_B_dEa);
//				if ((30 < t_T_mcp - bdn.t_B_dE) && (t_T_mcp - bdn.t_B_dE < 50)) {
//					ht_B_dE_zero_time_singles	->Fill(bdn.t_B_dE);
//					ht_B_dEa_zero_time_singles	->Fill(t_B_dEa);
//					ht_B_dEb_zero_time_singles	->Fill(t_B_dEb);
//					ht_T_mcp_zero_time_singles	->Fill(t_T_mcp);
//					zero_time_Tree->Fill();
//				}
//			}
//// Background BR-TOF
//			if (event_good==1 && s_capt_state==1 && t_dE_lo<bdn.t_B_dE && bdn.t_B_dE<t_dE_hi && a_dE_lo<bdn.a_B_dEsum && t_mcp_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum) {
//				tof = t_R_mcp - bdn.t_B_dE - BR_zeroTime[0];
//				h_bkgd_tof_BR	-> Fill(tof);
//				h_bkgd_tof		-> Fill(tof);
//			}
*/			
// LT-E_TOF
			if (event_good==1 && s_capt_state==0 && t_E_lo<t_L_E && t_L_E<t_E_hi && a_E_lo<a_L_E && t_mcp_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum) {
				h_E_tof_LT	-> Fill(t_T_mcp - t_L_E - LT_zeroTime_E[0]);
				h_E_tof		-> Fill(t_T_mcp - t_L_E - LT_zeroTime_E[0]);
			}
// LR-E_TOF
			if (event_good==1 && s_capt_state==0 && t_E_lo<t_L_E && t_L_E<t_E_hi && a_E_lo<a_L_E && t_mcp_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum) {
				h_E_tof_LR	-> Fill(t_R_mcp - t_L_E - LR_zeroTime_E[0]);
				h_E_tof		-> Fill(t_R_mcp - t_L_E - LR_zeroTime_E[0]);
			}
// BT-E_TOF
			if (event_good==1 && s_capt_state==0 && t_E_lo<t_B_E && t_B_E<t_E_hi && a_E_lo<a_B_E && t_mcp_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum) {
				h_E_tof_BT	-> Fill(t_T_mcp - t_B_E - BT_zeroTime_E[0]);
				h_E_tof		-> Fill(t_T_mcp - t_B_E - BT_zeroTime_E[0]);
			}
// BR-E_TOF
			if (event_good==1 && s_capt_state==0 && t_E_lo<t_B_E && t_B_E<t_E_hi && a_E_lo<a_B_E && t_mcp_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum) {
				h_E_tof_BR	-> Fill(t_R_mcp - t_B_E - BR_zeroTime_E[0]);
				h_E_tof		-> Fill(t_R_mcp - t_B_E - BR_zeroTime_E[0]);
			}
// Background LT-E_TOF
			if (event_good==1 && s_capt_state==1 && t_E_lo<t_L_E && t_L_E<t_E_hi && a_E_lo<a_L_E && t_mcp_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum) {
				h_bkgd_E_tof_LT	-> Fill(t_T_mcp - t_L_E - LT_zeroTime_E[0]);
				h_bkgd_E_tof	-> Fill(t_T_mcp - t_L_E - LT_zeroTime_E[0]);
			}
// Background LR-E_TOF
			if (event_good==1 && s_capt_state==1 && t_E_lo<t_L_E && t_L_E<t_E_hi && a_E_lo<a_L_E && t_mcp_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum) {
				h_bkgd_E_tof_LR	-> Fill(t_R_mcp - t_L_E - LR_zeroTime_E[0]);
				h_bkgd_E_tof	-> Fill(t_R_mcp - t_L_E - LR_zeroTime_E[0]);
			}
// Background BT-E_TOF
			if (event_good==1 && s_capt_state==1 && t_E_lo<t_B_E && t_B_E<t_E_hi && a_E_lo<a_B_E && t_mcp_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum) {
				h_bkgd_E_tof_BT	-> Fill(t_T_mcp - t_B_E - BT_zeroTime_E[0]);
				h_bkgd_E_tof	-> Fill(t_T_mcp - t_B_E - BT_zeroTime_E[0]);
			}
// Background BR-E_TOF
			if (event_good==1 && s_capt_state==1 && t_E_lo<t_B_E && t_B_E<t_E_hi && a_E_lo<a_B_E && t_mcp_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum) {
				h_bkgd_E_tof_BR	-> Fill(t_R_mcp - t_B_E - BR_zeroTime_E[0]);
				h_bkgd_E_tof	-> Fill(t_R_mcp - t_B_E - BR_zeroTime_E[0]);
			}

// RT-ge_TOF
			//if (t_ge_lo<t_R_ge && t_R_ge<t_ge_hi && a_ge_lo<a_R_ge && t_mcp_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum) {
			if (event_good==1 && s_capt_state==0 && t_ge_lo<t_R_ge && t_R_ge<t_ge_hi && a_ge_lo<a_R_ge && t_mcp_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum) {
				h_ge_tof_RT	-> Fill((t_T_mcp - t_R_ge));
				h_ge_tof	-> Fill((t_T_mcp - t_R_ge));
			}
// RR-ge_TOF
			if (event_good==1 && s_capt_state==0 && t_ge_lo<t_R_ge && t_R_ge<t_ge_hi && a_ge_lo<a_R_ge && t_mcp_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum) {
				h_ge_tof_RR	-> Fill((t_R_mcp - t_R_ge));
				h_ge_tof	-> Fill((t_R_mcp - t_R_ge));
			}
// TT-ge_TOF
			if (event_good==1 && s_capt_state==0 && t_ge_lo<t_T_ge && t_T_ge<t_ge_hi && a_ge_lo<a_T_ge && t_mcp_lo<t_T_mcp && a_mcp_lo<bdn.a_T_mcpSum) {
				h_ge_tof_TT	-> Fill((t_T_mcp - t_T_ge));
				h_ge_tof	-> Fill((t_T_mcp - t_T_ge));
			}
// TR-ge_TOF
			if (event_good==1 && s_capt_state==0 && t_ge_lo<t_T_ge && t_T_ge<t_ge_hi && a_ge_lo<a_T_ge && t_mcp_lo<t_R_mcp && a_mcp_lo<bdn.a_R_mcpSum) {
				h_ge_tof_TR	-> Fill((t_R_mcp - t_T_ge));
				h_ge_tof	-> Fill((t_R_mcp - t_T_ge));
			}
			
			//////////////////////////////////////////////////////////////////////////////////////////			
			
// LT-dE-ge_TOF
			//if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_L_dE && bdn.t_L_dE<t_dE_hi && a_dE_lo<bdn.a_L_dEsum && t_ge_lo<t_T_ge && a_ge_lo<a_T_ge) {
			if (t_dE_lo<bdn.t_L_dE && bdn.t_L_dE<t_dE_hi && a_dE_lo<bdn.a_L_dEsum && t_ge_lo<t_T_ge) { // && a_ge_lo<a_T_ge) {
				h_dE_ge_tof_LT	-> Fill((t_T_ge - bdn.t_L_dE));
				h_dE_ge_tof		-> Fill((t_T_ge - bdn.t_L_dE));
			}
// LR-dE-ge_TOF
			//if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_L_dE && bdn.t_L_dE<t_dE_hi && a_dE_lo<bdn.a_L_dEsum && t_ge_lo<t_R_ge && a_ge_lo<a_R_ge) {
			if (event_good==1 && t_dE_lo<bdn.t_L_dE && bdn.t_L_dE<t_dE_hi && a_dE_lo<bdn.a_L_dEsum && t_ge_lo<t_R_ge) { // && a_ge_lo<a_R_ge) {
				h_dE_ge_tof_LR	-> Fill((t_R_ge - bdn.t_L_dE));
				h_dE_ge_tof		-> Fill((t_R_ge - bdn.t_L_dE));
			}
// BT-dE-ge_TOF
			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_B_dE && bdn.t_B_dE<t_dE_hi && a_dE_lo<bdn.a_B_dEsum && t_ge_lo<t_T_ge) { // && a_ge_lo<a_T_ge) {
				h_dE_ge_tof_BT	-> Fill((t_T_ge - bdn.t_B_dE));
				h_dE_ge_tof		-> Fill((t_T_ge - bdn.t_B_dE));
			}
// BR-dE-ge_TOF
			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_B_dE && bdn.t_B_dE<t_dE_hi && a_dE_lo<bdn.a_B_dEsum && t_ge_lo<t_R_ge) { // && a_ge_lo<a_R_ge) {
				h_dE_ge_tof_BR	-> Fill((t_R_ge - bdn.t_B_dE));
				h_dE_ge_tof		-> Fill((t_R_ge - bdn.t_B_dE));
			}
			
			//////////////////////////////////////////////////////////////////////////////////////////			
			
// LL-dE-E_TOF
			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_L_dE && bdn.t_L_dE<t_dE_hi && a_dE_lo<bdn.a_L_dEsum && t_E_lo<t_L_E && a_E_lo<a_L_E) {
				h_dE_E_tof_LL	-> Fill((t_L_E - bdn.t_L_dE));
				h_dE_E_tof		-> Fill((t_L_E - bdn.t_L_dE));
			}
// LB-dE-E_TOF
			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_L_dE && bdn.t_L_dE<t_dE_hi && a_dE_lo<bdn.a_L_dEsum && t_E_lo<t_B_E && a_E_lo<a_B_E) {
				h_dE_E_tof_LB	-> Fill((t_B_E - bdn.t_L_dE));
				h_dE_E_tof		-> Fill((t_B_E - bdn.t_L_dE));
			}
// BL-dE-E_TOF
			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_B_dE && bdn.t_B_dE<t_dE_hi && a_dE_lo<bdn.a_B_dEsum && t_E_lo<t_L_E && a_E_lo<a_L_E) {
				h_dE_E_tof_BL	-> Fill((t_L_E - bdn.t_B_dE));
				h_dE_E_tof		-> Fill((t_L_E - bdn.t_B_dE));
			}
// BB-dE-E_TOF
			if (event_good==1 && s_capt_state==0 && t_dE_lo<bdn.t_B_dE && bdn.t_B_dE<t_dE_hi && a_dE_lo<bdn.a_B_dEsum && t_E_lo<t_B_E && a_E_lo<a_B_E) {
				h_dE_E_tof_BB	-> Fill((t_B_E - bdn.t_B_dE));
				h_dE_E_tof		-> Fill((t_B_E - bdn.t_B_dE));
			}
			
// TRAP FULL
			if (event_good==1 && s_capt_state==0) {
/*// Left-Top
//				if (-100<t_L_dEa && -100<t_L_dEb && -100<t_T_mcp) {
//					 tree_LT->Fill();
//					 list_LT->Enter(n_trig-1);
//					 if (t_L_dEa<0 && t_L_dEb<0 && bdn.a_T_mcpSum>200) {
//						 tof 	= 0.001*((double)t_T_mcp - (t_L_dEa + t_L_dEb)/2.0); // convert ns to us
//						 if (fast_tof1 <= tof && tof < fast_tof2) n_fast_LT++;
//						 if (slow_tof1 <= tof && tof < slow_tof2) n_slow_LT++;
//						 mcp_x		= 25.0 * bdn.T_mcpX; // mm
//						 mcp_y		= 25.0 * bdn.T_mcpY; // mm
//						 mcp_r		= sqrt(pow(mcp_dist,2) + pow(mcp_x,2) + pow(mcp_y,2));	// for KE=(1/2)*m*d^2/t^2
//						 vInv	= tof/mcp_r;
//						 En		= mc2_factor/vInv/vInv/c/c;
//						 //h_tof_LT	->Fill(tof);
//						 //h_tof		->Fill(tof);
//						 h_vInv_LT	->Fill(vInv);
//						 h_vInv		->Fill(vInv);
//						 h_En_LT	->Fill(En);
//						 h_En		->Fill(En);
//					}
//				}
// Left-Right
//				if (-100<t_L_dEa && -100<t_L_dEb && -100<t_R_mcp) {
//					 tree_LR->Fill();
//					 list_LR->Enter(n_trig-1);
//					 if (t_L_dEa<0 && t_L_dEb<0 && a_R_mcpA>100 && a_R_mcpB>100 && a_R_mcpC>100 && a_R_mcpD>100) {
//						 tof 	= 0.001*((double)t_R_mcp - (t_L_dEa + t_L_dEb)/2.0); // convert ns to us
//						 if (fast_tof1 <= tof && tof < fast_tof2) n_fast_LR++;
//						 if (slow_tof1 <= tof && tof < slow_tof2) n_slow_LR++;
//						 mcp_x		= 25.0 * bdn.R_mcpX;
//						 mcp_y		= 25.0 * bdn.R_mcpY;
//						 mcp_r		= sqrt(pow(mcp_dist,2) + pow(mcp_x,2) + pow(mcp_y,2));	// for KE=(1/2)*m*d^2/t^2
//						 vInv	= tof/mcp_r;
//						 En		= mc2_factor/vInv/vInv/c/c;
//						 //h_tof_LR	->Fill(tof);
//						 //h_tof		->Fill(tof);
//						 h_vInv_LR	->Fill(vInv);
//						 h_vInv		->Fill(vInv);
//						 h_En_LR	->Fill(En);
//						 h_En		->Fill(En);
//					}
//				}
// Bottom-Top
//				if (-100<t_B_dEa && -100<t_B_dEb && -100<t_T_mcp) {
//					 tree_BT->Fill();
//					 list_BT->Enter(n_trig-1);
//					 if (t_B_dEa<0 && t_B_dEb<0 && a_T_mcpA>100 && a_T_mcpB>100 && a_T_mcpC>100 && a_T_mcpD>100) {
//						 tof 	= 0.001*((double)t_T_mcp - (t_B_dEa + t_B_dEb)/2.0); // convert ns to us
//						 if (fast_tof1 <= tof && tof < fast_tof2) n_fast_BT++;
//						 if (slow_tof1 <= tof && tof < slow_tof2) n_slow_BT++;
//						 mcp_x		= 25.0 * bdn.T_mcpX;
//						 mcp_y		= 25.0 * bdn.T_mcpY;
//						 mcp_r		= sqrt(pow(mcp_dist,2) + pow(mcp_x,2) + pow(mcp_y,2));	// for KE=(1/2)*m*d^2/t^2
//						 vInv	= tof/mcp_r;
//						 En		= mc2_factor/vInv/vInv/c/c;
//						 //h_tof_BT	->Fill(tof);
//						 //h_tof		->Fill(tof);
//						 h_vInv_BT	->Fill(vInv);
//						 h_vInv		->Fill(vInv);
//						 h_En_BT	->Fill(En);
//						 h_En		->Fill(En);
//					}
//				}
// Bottom-Right
//				if (-100<t_B_dEa && -100<t_B_dEb && -100<t_R_mcp) {
//					 tree_BR->Fill();
//					 list_BR->Enter(n_trig-1);
//					 if (t_B_dEa<0 && t_B_dEb<0 && a_R_mcpA>100 && a_R_mcpB>100 && a_R_mcpC>100 && a_R_mcpD>100) {
//						 tof 	= 0.001*((double)t_R_mcp - (t_B_dEa + t_B_dEb)/2.0); // convert ns to us
//						 if (fast_tof1 <= tof && tof < fast_tof2) n_fast_BR++;
//						 if (slow_tof1 <= tof && tof < slow_tof2) n_slow_BR++;
//						 mcp_x		= 25.0 * bdn.R_mcpX;
//						 mcp_y		= 25.0 * bdn.R_mcpY;
//						 mcp_r		= sqrt(pow(mcp_dist,2) + pow(mcp_x,2) + pow(mcp_y,2));	// for KE=(1/2)*m*d^2/t^2
//						 vInv	= tof/mcp_r;
//						 En		= mc2_factor/vInv/vInv/c/c;
//						 //h_tof_BR	->Fill(tof);
//						 //h_tof		->Fill(tof);
//						 h_vInv_BR	->Fill(vInv);
//						 h_vInv		->Fill(vInv);
//						 h_En_BR	->Fill(En);
//						 h_En		->Fill(En);
//					}
//				}
// Top MCP maps */
//				if (a_T_mcpA>0 && a_T_mcpB>0 && a_T_mcpC>0 && a_T_mcpD>0) {
				if (bdn.a_T_mcpSum > a_missing_mcp_post) {
					h_T_mcpMap_all->Fill(bdn.T_mcpX,bdn.T_mcpY);
					if (bdn.a_T_mcpSum > 50) { //if (a_T_mcpA>50 && a_T_mcpB>50 && a_T_mcpC>50 && a_T_mcpD>50) {
						h_T_mcpMap_post50->Fill(bdn.T_mcpX,bdn.T_mcpY);
						if (bdn.a_T_mcpSum > 100) { //if (a_T_mcpA>100 && a_T_mcpB>100 && a_T_mcpC>100 && a_T_mcpD>100) {
							h_T_mcpMap_post100->Fill(bdn.T_mcpX,bdn.T_mcpY);
							if (bdn.a_T_mcpSum > 200) {
								h_T_mcpMap_post200->Fill(bdn.T_mcpX,bdn.T_mcpY);
								if (bdn.a_T_mcpSum > 250) { //if (a_T_mcpA>250 && a_T_mcpB>250 && a_T_mcpC>250 && a_T_mcpD>250) {
									h_T_mcpMap_post250->Fill(bdn.T_mcpX,bdn.T_mcpY);
								}
							}
						}
					}
				}
// Right MCP maps
//				if (a_R_mcpA>0 && a_R_mcpB>0 && a_R_mcpC>0 && a_R_mcpD>0) {
				if (bdn.a_R_mcpSum > a_missing_mcp_post) {
					h_R_mcpMap_all->Fill(bdn.R_mcpX,bdn.R_mcpY);
					if (bdn.a_R_mcpSum > 50) { //if (a_R_mcpA>50 && a_R_mcpB>50 && a_R_mcpC>50 && a_R_mcpD>50) {
						h_R_mcpMap_post50->Fill(bdn.R_mcpX,bdn.R_mcpY);
						if (bdn.a_R_mcpSum > 100) { //if (a_R_mcpA>100 && a_R_mcpB>100 && a_R_mcpC>100 && a_R_mcpD>100) {
							h_R_mcpMap_post100->Fill(bdn.R_mcpX,bdn.R_mcpY);
							if (bdn.a_R_mcpSum > 200) {
								h_R_mcpMap_post200->Fill(bdn.R_mcpX,bdn.R_mcpY);
								if (bdn.a_R_mcpSum > 250) { //if (a_R_mcpA>250 && a_R_mcpB>250 && a_R_mcpC>250 && a_R_mcpD>250) {
									h_R_mcpMap_post250->Fill(bdn.R_mcpX,bdn.R_mcpY);
								}
							}
						}
					}
				}
			} // end TRAP FULL
			
// TRAP EMPTY
			if (event_good==1 && s_capt_state==1) { // conditions for all histograms
/*				
// Left-Top
//				if (-100<t_L_dEa && -100<t_L_dEb && -100<t_T_mcp) {
//					 tree_bkgd_LT->Fill();
//					 list_bkgd_LT->Enter(n_trig-1);
//					 if (t_L_dEa<0 && t_L_dEb<0 && a_T_mcpA>100 && a_T_mcpB>100 && a_T_mcpC>100 && a_T_mcpD>100) {
//						 tof 	= 0.001*((double)t_T_mcp - (t_L_dEa + t_L_dEb)/2.0); // convert ns to us
//						 mcp_x		= 25.0 * bdn.T_mcpX;
//						 mcp_y		= 25.0 * bdn.T_mcpY;
//						 mcp_r		= sqrt(pow(mcp_dist,2) + pow(mcp_x,2) + pow(mcp_y,2));	// for KE=(1/2)*m*d^2/t^2
//						 vInv	= tof/mcp_r;
//						 En		= mc2_factor/vInv/vInv/c/c;
//						 //h_bkgd_tof_LT	->Fill(tof);
//						 //h_bkgd_tof		->Fill(tof);
//						 h_bkgd_vInv_LT	->Fill(vInv);
//						 h_bkgd_vInv	->Fill(vInv);
//						 h_bkgd_En_LT	->Fill(En);
//						 h_bkgd_En		->Fill(En);
//					}
//				}
// Left-Right
//				if (-100<t_L_dEa && -100<t_L_dEb && -100<t_R_mcp) {
//					 tree_bkgd_LR->Fill();
//					 list_bkgd_LR->Enter(n_trig-1);
//					 if (t_L_dEa<0 && t_L_dEb<0 && a_R_mcpA>100 && a_R_mcpB>100 && a_R_mcpC>100 && a_R_mcpD>100) {
//						 tof 	= 0.001*((double)t_R_mcp - (t_L_dEa + t_L_dEb)/2.0); // convert ns to us
//						 mcp_x		= 25.0 * bdn.R_mcpX;
//						 mcp_y		= 25.0 * bdn.R_mcpY;
//						 mcp_r		= sqrt(pow(mcp_dist,2) + pow(mcp_x,2) + pow(mcp_y,2));	// for KE=(1/2)*m*d^2/t^2
//						 vInv	= tof/mcp_r;
//						 En		= mc2_factor/vInv/vInv/c/c;
//						 //h_bkgd_tof_LR	->Fill(tof);
//						 //h_bkgd_tof		->Fill(tof);
//						 h_bkgd_vInv_LR	->Fill(vInv);
//						 h_bkgd_vInv	->Fill(vInv);
//						 h_bkgd_En_LR	->Fill(En);
//						 h_bkgd_En		->Fill(En);
//					}
//				}
// Bottom-Top
//				if (-100<t_B_dEa && -100<t_B_dEb && -100<t_T_mcp) {
//					 tree_bkgd_BT->Fill();
//					 list_bkgd_BT->Enter(n_trig-1);
//					 if (t_B_dEa<0 && t_B_dEb<0 && a_T_mcpA>100 && a_T_mcpB>100 && a_T_mcpC>100 && a_T_mcpD>100) {
//						 tof 	= 0.001*((double)t_T_mcp - (t_B_dEa + t_B_dEb)/2.0); // convert ns to us
//						 mcp_x		= 25.0 * bdn.T_mcpX;
//						 mcp_y		= 25.0 * bdn.T_mcpY;
//						 mcp_r		= sqrt(pow(mcp_dist,2) + pow(mcp_x,2) + pow(mcp_y,2));	// for KE=(1/2)*m*d^2/t^2
//						 vInv	= tof/mcp_r;
//						 En		= mc2_factor/vInv/vInv/c/c;
//						 //h_bkgd_tof_BT	->Fill(tof);
//						 //h_bkgd_tof		->Fill(tof);
//						 h_bkgd_vInv_BT	->Fill(vInv);
//						 h_bkgd_vInv	->Fill(vInv);
//						 h_bkgd_En_BT	->Fill(En);
//						 h_bkgd_En		->Fill(En);
//					}
//				}
// Bottom-Right
//				if (-100<t_B_dEa && -100<t_B_dEb && -100<t_R_mcp) {
//					 tree_bkgd_BR->Fill();
//					 list_bkgd_BR->Enter(n_trig-1);
//					 if (t_B_dEa<0 && t_B_dEb<0 && a_R_mcpA>100 && a_R_mcpB>100 && a_R_mcpC>100 && a_R_mcpD>100) {
//						 tof 	= 0.001*((double)t_R_mcp - (t_B_dEa + t_B_dEb)/2.0); // convert ns to us
//						 mcp_x		= 25.0 * bdn.R_mcpX;
//						 mcp_y		= 25.0 * bdn.R_mcpY;
//						 mcp_r		= sqrt(pow(mcp_dist,2) + pow(mcp_x,2) + pow(mcp_y,2));	// for KE=(1/2)*m*d^2/t^2
//						 vInv	= tof/mcp_r;
//						 En		= mc2_factor/vInv/vInv/c/c;
//						 //h_bkgd_tof_BR	->Fill(tof);
//						 //h_bkgd_tof		->Fill(tof);
//						 h_bkgd_vInv_BR	->Fill(vInv);
//						 h_bkgd_vInv	->Fill(vInv);
//						 h_bkgd_En_BR	->Fill(En);
//						 h_bkgd_En		->Fill(En);
//					}
//				}
*/
// Top MCP maps
//				if (a_T_mcpA>0 && a_T_mcpB>0 && a_T_mcpC>0 && a_T_mcpD>0) {
				if (bdn.a_T_mcpSum > a_missing_mcp_post) {
					h_bkgd_T_mcpMap_all->Fill(bdn.T_mcpX,bdn.T_mcpY);
					if (bdn.a_T_mcpSum > 50) { //if (a_T_mcpA>50 && a_T_mcpB>50 && a_T_mcpC>50 && a_T_mcpD>50) {
						h_bkgd_T_mcpMap_post50->Fill(bdn.T_mcpX,bdn.T_mcpY);
						if (bdn.a_T_mcpSum > 100) { //if (a_T_mcpA>100 && a_T_mcpB>100 && a_T_mcpC>100 && a_T_mcpD>100) {
							h_bkgd_T_mcpMap_post100->Fill(bdn.T_mcpX,bdn.T_mcpY);
							if (bdn.a_T_mcpSum > 200) {
								h_bkgd_T_mcpMap_post200->Fill(bdn.T_mcpX,bdn.T_mcpY);
								if (bdn.a_T_mcpSum > 250) { //if (a_T_mcpA>250 && a_T_mcpB>250 && a_T_mcpC>250 && a_T_mcpD>250) {
									h_bkgd_T_mcpMap_post250->Fill(bdn.T_mcpX,bdn.T_mcpY);
								}
							}
						}
					}
				}
// Right MCP maps
//				if (a_R_mcpA>0 && a_R_mcpB>0 && a_R_mcpC>0 && a_R_mcpD>0) {
				if (bdn.a_R_mcpSum > a_missing_mcp_post) {
					h_bkgd_R_mcpMap_all->Fill(bdn.R_mcpX,bdn.R_mcpY);
					if (bdn.a_R_mcpSum > 50) { //if (a_R_mcpA>50 && a_R_mcpB>50 && a_R_mcpC>50 && a_R_mcpD>50) {
						h_bkgd_R_mcpMap_post50->Fill(bdn.R_mcpX,bdn.R_mcpY);
						if (bdn.a_R_mcpSum > 100) { //if (a_R_mcpA>100 && a_R_mcpB>100 && a_R_mcpC>100 && a_R_mcpD>100) {
							h_bkgd_R_mcpMap_post100->Fill(bdn.R_mcpX,bdn.R_mcpY);
							if (bdn.a_R_mcpSum > 200) { 
								h_bkgd_R_mcpMap_post200->Fill(bdn.R_mcpX,bdn.R_mcpY);
								if (bdn.a_R_mcpSum > 250) { //if (a_R_mcpA>250 && a_R_mcpB>250 && a_R_mcpC>250 && a_R_mcpD>250) {
									h_bkgd_R_mcpMap_post250->Fill(bdn.R_mcpX,bdn.R_mcpY);
								}
							}
						}
					}
				}
			} // end TRAP EMPTY
				/*
				//if (bdn.T_mcpX != 0) cout<<"x= "<<bdn.T_mcpX<<endl;
				
//				// Move pointer to EjectScaler data:
//				while (*p++ != 0x100eca1e) {
//					// skip a place
//				};
//				
//				// Read EjectScaler data:
//				s_ms_since_eject	= (*p++ & 0xffffff);
//				s_capt		= (*p++ & 0xffffff);
//				s_SiX4		= (*p++	& 0xffffff);
//				
//				// Every t_cycle, update histos and reset cycle counts:
//				//if (s_ms_since_eject < s_ms_since_eject_last) { // new trap cycle
//				if (s_ms_since_eject/1000 < t_cycle/2 && t_cycle/2 < s_ms_since_eject_last/1000) { // new trap cycle
//					
//					trig_bin++;
//					if (trig_bin == 1) first_cycle_flag = 1;
//					
//					// Write last cycle's values to histos:
//					h_SiX4	->	SetBinContent(trig_bin, s_SiX4_last/t_cycle);
//					
//					// Add last cycle's values to totals:
//					sTot_SiX4	+= s_SiX4_last;
//					
//					// Reset nCyc values for next cycle:
//					
//				} // if new trap cycle
				*/
				break;
				
			case SE_TYPE_SYNC:
				
				n_sync++;
				sync_day_last	= sync_day;
				
				// Initialize pointer:
				e0 = ScarletEvnt(h);
				e1 = e0[1];
				p = reinterpret_cast<int*>(e1.body());
				
				// Move pointer to TrigSyncScaler data:
				while (*p++ != 0x1002ca1e) {
					// skip a place
				};
				
				// Read scaler data; these are counts/sync:
				s_T_mcp		= (*p++ & 0xffffff);
				s_R_mcp		= (*p++ & 0xffffff);
				s_B_dEa		= (*p++ & 0xffffff);
				s_B_E		= (*p++ & 0xffffff);
				s_L_dEa		= (*p++ & 0xffffff);
				s_L_dEb		= (*p++ & 0xffffff);
				s_L_E		= (*p++ & 0xffffff);
				s_T_ge		= (*p++	& 0xffffff);
				s_R_ge		= (*p++ & 0xffffff);
				s_B_dEb		= (*p++ & 0xffffff);
				s_SiX4_ts	= (*p++ & 0xffffff);
				
				// Move pointer to timestamp:
				while (*p++ != 0x0000abcd) {
					// skip a place
				};
				
				sync_day	= *p++;
				sync_hour	= *p++;
				sync_min	= *p++;
				sync_sec	= *p++;
				if (sync_day_last > sync_day)	fake_day = sync_day_last + 1; // this happens once when day turns over...
				if (sync_day == 0)				sync_day = fake_day; // ... after which we always replace it with the fake day (eg. Oct 32)
				sync_time_sec = time_in_seconds(sync_day, sync_hour, sync_min, sync_sec) - start_time_sec;
				//printf("\n%d",sync_time_sec);
				
				// Add latest counts to totals for current cycle:
				sCyc_B_dEa	+= s_B_dEa;
				sCyc_B_dEb	+= s_B_dEb;
				sCyc_B_E	+= s_B_E;
				sCyc_L_dEa	+= s_L_dEa;
				sCyc_L_dEb	+= s_L_dEb;
				sCyc_L_E	+= s_L_E;
				sCyc_R_mcp	+= s_R_mcp;
				sCyc_R_ge	+= s_R_ge;
				sCyc_T_mcp	+= s_T_mcp;
				sCyc_T_ge	+= s_T_ge;
				sCyc_SiX4_ts+= s_SiX4_ts;
				
				// Every n_sync_update syncs, update histos and reset cycle counts:
				if (n_sync % n_sync_update == 0) {
					
					sync_bin++;
					//printf("Last sync:      d:%02d, h:%02d, m:%02d, s:%02d\n", sync_day, sync_hour, sync_min, sync_sec);
					//bin = t_per_sync*n_sync/t_cycle; // up to experimenter to make this an int
					/*
//					// Write last cycle's values to histos:
//					hs_B_dEa-> SetBinContent(sync_bin, sCyc_B_dEa	/t_cycle);
//					hs_B_dEb-> SetBinContent(sync_bin, sCyc_B_dEb	/t_cycle);
//					hs_B_E	-> SetBinContent(sync_bin, sCyc_B_E		/t_cycle);
//					hs_L_dEa-> SetBinContent(sync_bin, sCyc_L_dEa	/t_cycle);
//					hs_L_dEb-> SetBinContent(sync_bin, sCyc_L_dEb	/t_cycle);
//					hs_L_E	-> SetBinContent(sync_bin, sCyc_L_E		/t_cycle);
//					hs_R_mcp-> SetBinContent(sync_bin, sCyc_R_mcp	/t_cycle);
//					hs_R_ge	-> SetBinContent(sync_bin, sCyc_R_ge	/t_cycle);
//					hs_T_mcp-> SetBinContent(sync_bin, sCyc_T_mcp	/t_cycle);
//					hs_T_ge	-> SetBinContent(sync_bin, sCyc_T_ge	/t_cycle);
//					hs_SiX4	-> SetBinContent(sync_bin, sCyc_SiX4_ts	/t_cycle);
					*/
					// Add last cycle's values to totals:
					sTot_B_dEa	+= sCyc_B_dEa;
					sTot_B_dEb	+= sCyc_B_dEb;
					sTot_B_E	+= sCyc_B_E;
					sTot_L_dEa	+= sCyc_L_dEa;
					sTot_L_dEb	+= sCyc_L_dEb;
					sTot_L_E	+= sCyc_L_E;
					sTot_R_mcp	+= sCyc_R_mcp;
					sTot_R_ge	+= sCyc_R_ge;
					sTot_T_mcp	+= sCyc_T_mcp;
					sTot_T_ge	+= sCyc_T_ge;
					sTot_SiX4_ts+= sCyc_SiX4_ts;
					/*
//					printf("\n sCyc_B_dEa %d", sCyc_B_dEa);
//					printf("\n sCyc_B_dEb %d", sCyc_B_dEb);
//					printf("\n sCyc_B_E %d"  , sCyc_B_E);
//					printf("\n sCyc_L_dEa %d", sCyc_L_dEa);
//					printf("\n sCyc_L_dEb %d", sCyc_L_dEb);
//					printf("\n sCyc_L_E %d"  , sCyc_L_E);
//					printf("\n sCyc_R_ge %d" , sCyc_R_ge);
//					printf("\n sCyc_T_ge %d" , sCyc_T_ge);
//					printf("\n sCyc_R_mcp %d", sCyc_R_mcp);
//					printf("\n sCyc_T_mcp %d", sCyc_T_mcp);
//					printf("\n sCyc_SiX4 %d" , sCyc_SiX4_ts);
					*/
					// Reset nCyc values for next cycle:
					sCyc_B_dEa	= 0;
					sCyc_B_dEb	= 0;
					sCyc_B_E	= 0;
					sCyc_L_dEa	= 0;
					sCyc_L_dEb	= 0;
					sCyc_L_E	= 0;
					sCyc_R_mcp	= 0;
					sCyc_R_ge	= 0;
					sCyc_T_mcp	= 0;
					sCyc_T_ge	= 0;
					sCyc_SiX4_ts= 0;
										
				} // if new trap cycle
				
				break;
				
			case SE_TYPE_ACQUIRE:
				
				// Initialize pointer:
				e0 = ScarletEvnt(h);
				e1 = e0[1];
				p = reinterpret_cast<int*>(e1.body());
				
				// Move pointer to timestamp:
				while (*p++ != 0x0000abcd) {
					// skip a place
				};
				
				start_day	= *p++;
				start_hour	= *p++;
				start_min	= *p++;
				start_sec	= *p++;
				start_time_sec	= time_in_seconds(start_day, start_hour, start_min, start_sec);
				
				break;
				
			case SE_TYPE_STOP:
				
				stop_flag = 1;
				
				// Initialize pointer:
				e0 = ScarletEvnt(h);
				e1 = e0[1];
				p = reinterpret_cast<int*>(e1.body());
				
				// Move pointer to timestamp:
				while (*p++ != 0x0000abcd) {
					// skip a place
				};
				
				stop_day	= *p++;
				stop_hour	= *p++;
				stop_min	= *p++;
				stop_sec	= *p++;
				stop_time_sec	= time_in_seconds(stop_day, stop_hour, stop_min, stop_sec);
				
				break;
				
			default:
				// Ignore other types of events
				// See ScarletEvntHdr.h for other event types
				break;
				
		} // switch
		
	} //while (getevent()!=0)
	
	//**** Now have data for all events in data tree ****
	//**** Now fill metadata line in metadata tree ****
	
	if (stop_flag == 0)	\
		run_time_sec = sync_time_sec;
	else \
		run_time_sec = stop_time_sec - start_time_sec;
	sTot_all	= sTot_B_dEa + sTot_B_dEb + sTot_B_E + sTot_L_dEa + sTot_L_dEb + sTot_L_E + sTot_R_mcp + sTot_R_ge + sTot_T_mcp + sTot_T_ge;
	nt_all		= nt_B_dEa + nt_B_dEb + nt_B_E + nt_L_dEa + nt_L_dEb + nt_L_E + nt_R_mcp + nt_R_ge + nt_T_mcp + nt_T_ge;
	
//~~~~~~~~ Fill metadata (tree) ~~~~~~~~//
	
	metadata.n_run			= n_run;
	metadata.n_trigs		= n_trig;
	metadata.tot_trigs		= tot_trigs;
	metadata.n_syncs		= n_sync;
	metadata.n_treeEntries	= n_trig;
	metadata.n_bad_events	= n_bad_events;
	metadata.bkgd_good		= bkgd_good;
	metadata.start_month	= 0;//start_month;
	metadata.start_day		= start_day;
	metadata.start_hour		= start_hour;
	metadata.start_min		= start_min;
	metadata.start_sec		= start_sec;
	metadata.start_time_sec	= start_time_sec;
	metadata.stop_month		= 0;//stop_month;
	metadata.stop_day		= stop_day;
	metadata.stop_hour		= stop_hour;
	metadata.stop_min		= stop_min;
	metadata.stop_sec		= stop_sec;
	metadata.stop_time_sec	= stop_time_sec;
	metadata.run_time_sec	= run_time_sec;
	
	metadata.n_scaler_hits_B_dEa	= sTot_B_dEa;
	metadata.n_scaler_hits_B_dEb	= sTot_B_dEb;
	metadata.n_scaler_hits_B_E		= sTot_B_E;
	metadata.n_scaler_hits_L_dEa	= sTot_L_dEa;
	metadata.n_scaler_hits_L_dEb	= sTot_L_dEb;
	metadata.n_scaler_hits_L_E		= sTot_L_E;
	metadata.n_scaler_hits_R_mcp	= sTot_R_mcp;
	metadata.n_scaler_hits_R_ge		= sTot_R_ge;
	metadata.n_scaler_hits_T_mcp	= sTot_T_mcp;
	metadata.n_scaler_hits_T_ge		= sTot_T_ge;
	
	metadata.n_tdc_hits_B_dEa		= nt_B_dEa;
	metadata.n_tdc_hits_B_dEb		= nt_B_dEb;
	metadata.n_tdc_hits_B_E			= nt_B_E;
	metadata.n_tdc_hits_L_dEa		= nt_L_dEa;
	metadata.n_tdc_hits_L_dEb		= nt_L_dEb;
	metadata.n_tdc_hits_L_E			= nt_L_E;
	metadata.n_tdc_hits_R_mcp		= nt_R_mcp;
	metadata.n_tdc_hits_R_ge		= nt_R_ge;
	metadata.n_tdc_hits_T_mcp		= nt_T_mcp;
	metadata.n_tdc_hits_T_ge		= nt_T_ge;
	
	metadata.n_adc_hits_B_dEa		= na_B_dEa;
	metadata.n_adc_hits_B_dEb		= na_B_dEb;
	metadata.n_adc_hits_B_E			= na_B_E;
	metadata.n_adc_hits_L_dEa		= na_L_dEa;
	metadata.n_adc_hits_L_dEb		= na_L_dEb;
	metadata.n_adc_hits_L_E			= na_L_E;
	metadata.n_adc_hits_R_mcpA		= na_R_mcpA;
	metadata.n_adc_hits_R_mcpB		= na_R_mcpB;
	metadata.n_adc_hits_R_mcpC		= na_R_mcpC;
	metadata.n_adc_hits_R_mcpD		= na_R_mcpD;
	metadata.n_adc_hits_R_mcpE		= na_R_mcpE;
	metadata.n_adc_hits_R_ge		= na_R_ge;
	metadata.n_adc_hits_R_ge_highE	= na_R_ge_highE;
	metadata.n_adc_hits_T_mcpA		= na_T_mcpA;
	metadata.n_adc_hits_T_mcpB		= na_T_mcpB;
	metadata.n_adc_hits_T_mcpC		= na_T_mcpC;
	metadata.n_adc_hits_T_mcpD		= na_T_mcpD;
	metadata.n_adc_hits_T_mcpE		= na_T_mcpE;
	metadata.n_adc_hits_T_ge		= na_T_ge;
	metadata.n_adc_hits_T_ge_highE	= na_T_ge_highE;
	metadata.n_missing_adc_hits_R_mcpA	= na_R_mcpA_missing;
	metadata.n_missing_adc_hits_R_mcpB	= na_R_mcpB_missing;
	metadata.n_missing_adc_hits_R_mcpC	= na_R_mcpC_missing;
	metadata.n_missing_adc_hits_R_mcpD	= na_R_mcpD_missing;
	metadata.n_missing_adc_hits_R_mcpE	= na_R_mcpE_missing;
	metadata.n_missing_adc_hits_T_mcpA	= na_T_mcpA_missing;
	metadata.n_missing_adc_hits_T_mcpB	= na_T_mcpB_missing;
	metadata.n_missing_adc_hits_T_mcpC	= na_T_mcpC_missing;
	metadata.n_missing_adc_hits_T_mcpD	= na_T_mcpD_missing;
	metadata.n_missing_adc_hits_T_mcpE	= na_T_mcpE_missing;
	/*
//	metadata.n_fast_LT	= n_fast_LT;
//	metadata.n_fast_LR	= n_fast_LR;
//	metadata.n_fast_BT	= n_fast_BT;
//	metadata.n_fast_BR	= n_fast_BR;
//	metadata.n_slow_LT	= n_slow_LT;
//	metadata.n_slow_LR	= n_slow_LR;
//	metadata.n_slow_BT	= n_slow_BT;
//	metadata.n_slow_BR	= n_slow_BR;
	*/
	metadata.tot_liveTime_us = tot_liveTime_us;
	metadata.tot_runTime_us = tot_runTime_us;
	
// Count Fast and Slow recoils, and other TOF regions
	Float_t oopsPerNs, oopsPerNsBkgd;
	for (Int_t iCombo; iCombo < (nDetectors_dE * nDetectors_MCP); iCombo++) {
		
	// By counting
		// Trap full
		metadata.nNetZeroTOFCount[iCombo]	= metadata.nZeroTOFCount[iCombo]	- metadata.nOopsCount[iCombo] * (tof_zero_hi - tof_zero_lo)     / (tof_oops_hi - tof_oops_lo);
		metadata.nNetLowTOFCount[iCombo]	= metadata.nLowTOFCount[iCombo]		- metadata.nOopsCount[iCombo] * (tof_lowTOF_hi - tof_lowTOF_lo) / (tof_oops_hi - tof_oops_lo);
//		metadata.nNetFastCount[iCombo]		= metadata.nFastCount[iCombo]		- metadata.nOopsCount[iCombo] * (tof_fast_hi - tof_fast_lo)     / (tof_oops_hi - tof_oops_lo); // these have MCP-specific ion ranges and have to be treated in the switch
//		metadata.nNetSlowCount[iCombo]		= metadata.nSlowCount[iCombo]		- metadata.nOopsCount[iCombo] * (tof_slow_hi - tof_slow_lo)     / (tof_oops_hi - tof_oops_lo);
		// Trap empty
		metadata.nNetZeroTOFBkgdCount[iCombo]	= metadata.nZeroTOFBkgdCount[iCombo]	- metadata.nOopsBkgdCount[iCombo] * (tof_zero_hi - tof_zero_lo)     / (tof_oops_hi - tof_oops_lo);
		metadata.nNetLowTOFBkgdCount[iCombo]	= metadata.nLowTOFBkgdCount[iCombo]		- metadata.nOopsBkgdCount[iCombo] * (tof_lowTOF_hi - tof_lowTOF_lo) / (tof_oops_hi - tof_oops_lo);
//		metadata.nNetFastBkgdCount[iCombo]		= metadata.nFastBkgdCount[iCombo]		- metadata.nOopsBkgdCount[iCombo] * (tof_fast_hi - tof_fast_lo)     / (tof_oops_hi - tof_oops_lo);
//		metadata.nNetSlowBkgdCount[iCombo]		= metadata.nSlowBkgdCount[iCombo]		- metadata.nOopsBkgdCount[iCombo] * (tof_slow_hi - tof_slow_lo)     / (tof_oops_hi - tof_oops_lo);
		
	// By integrals
		switch (iCombo) { // switch is needed because the h_tof histos are not in an array... this could be amended but would affect many other codes
		// In these calculations, the Integral method takes integer bin #'s as arguments (as of now one bin is 0.5ns)
		// while the background subtraction is done in counts per ns
		// To avoid this ugly switching, you would need to put everything that has a B/L/R/T into arrays that are enumerated according to the combos LT, LR, etc.
		// This would require updates to pretty much every piece of code I have, so I won't be doing it. -SC
			case LT:
				// Trap Full
				metadata.nZeroTOFIntegral[iCombo]	= h_tof_LT->Integral(tofBin_zero_lo,   tofBin_zero_hi);
				metadata.nLowTOFIntegral[iCombo]	= h_tof_LT->Integral(tofBin_lowTOF_lo, tofBin_lowTOF_hi);
				metadata.nFastIntegral[iCombo]		= h_tof_LT->Integral(tofBin_T_fast_lo, tofBin_T_fast_hi);
				metadata.nSlowIntegral[iCombo]		= h_tof_LT->Integral(tofBin_T_slow_lo, tofBin_T_slow_hi);
				metadata.nOopsIntegral[iCombo]		= h_tof_LT->Integral(tofBin_oops_lo,   tofBin_oops_hi);
				oopsPerNs							= metadata.nOopsIntegral[iCombo] / (tof_oops_hi - tof_oops_lo);
				metadata.nNetFastIntegral[iCombo]	= metadata.nFastIntegral[iCombo]	- oopsPerNs * (tof_T_fast_hi - tof_T_fast_lo);
				metadata.nNetSlowIntegral[iCombo]	= metadata.nSlowIntegral[iCombo]	- oopsPerNs * (tof_T_slow_hi - tof_T_slow_lo);
				metadata.nNetFastCount[iCombo]		= metadata.nFastCount[iCombo]		- metadata.nOopsCount[iCombo] * (tof_T_fast_hi - tof_T_fast_lo)     / (tof_oops_hi - tof_oops_lo);
				metadata.nNetSlowCount[iCombo]		= metadata.nSlowCount[iCombo]		- metadata.nOopsCount[iCombo] * (tof_T_slow_hi - tof_T_slow_lo)     / (tof_oops_hi - tof_oops_lo);
				// Trap Empty
				metadata.nZeroTOFBkgdIntegral[iCombo]	= h_bkgd_tof_LT->Integral(tofBin_zero_lo,   tofBin_zero_hi);
				metadata.nLowTOFBkgdIntegral[iCombo]	= h_bkgd_tof_LT->Integral(tofBin_lowTOF_lo, tofBin_lowTOF_hi);
				metadata.nFastBkgdIntegral[iCombo]		= h_bkgd_tof_LT->Integral(tofBin_T_fast_lo, tofBin_T_fast_hi);
				metadata.nSlowBkgdIntegral[iCombo]		= h_bkgd_tof_LT->Integral(tofBin_T_slow_lo, tofBin_T_slow_hi);
				metadata.nOopsBkgdIntegral[iCombo]		= h_bkgd_tof_LT->Integral(tofBin_oops_lo,   tofBin_oops_hi);
				oopsPerNsBkgd							= metadata.nOopsBkgdIntegral[iCombo] / (tof_oops_hi - tof_oops_lo);
				metadata.nNetFastBkgdIntegral[iCombo]	= metadata.nFastBkgdIntegral[iCombo]	- oopsPerNsBkgd * (tof_T_fast_hi - tof_T_fast_lo);
				metadata.nNetSlowBkgdIntegral[iCombo]	= metadata.nSlowBkgdIntegral[iCombo]	- oopsPerNsBkgd * (tof_T_slow_hi - tof_T_slow_lo);
				metadata.nNetFastBkgdCount[iCombo]		= metadata.nFastBkgdCount[iCombo]		- metadata.nOopsBkgdCount[iCombo] * (tof_T_fast_hi - tof_T_fast_lo)     / (tof_oops_hi - tof_oops_lo);
				metadata.nNetSlowBkgdCount[iCombo]		= metadata.nSlowBkgdCount[iCombo]		- metadata.nOopsBkgdCount[iCombo] * (tof_T_slow_hi - tof_T_slow_lo)     / (tof_oops_hi - tof_oops_lo);
				break;
			case LR:
				// Trap Full
				metadata.nZeroTOFIntegral[iCombo]	= h_tof_LR->Integral(tofBin_zero_lo,   tofBin_zero_hi);
				metadata.nLowTOFIntegral[iCombo]	= h_tof_LR->Integral(tofBin_lowTOF_lo, tofBin_lowTOF_hi);
				metadata.nFastIntegral[iCombo]		= h_tof_LR->Integral(tofBin_R_fast_lo, tofBin_R_fast_hi);
				metadata.nSlowIntegral[iCombo]		= h_tof_LR->Integral(tofBin_R_slow_lo, tofBin_R_slow_hi);
				metadata.nOopsIntegral[iCombo]		= h_tof_LR->Integral(tofBin_oops_lo,   tofBin_oops_hi);
				oopsPerNs							= metadata.nOopsIntegral[iCombo] / (tof_oops_hi - tof_oops_lo);
				metadata.nNetFastIntegral[iCombo]	= metadata.nFastIntegral[iCombo]	- oopsPerNs * (tof_R_fast_hi - tof_R_fast_lo);
				metadata.nNetSlowIntegral[iCombo]	= metadata.nSlowIntegral[iCombo]	- oopsPerNs * (tof_R_slow_hi - tof_R_slow_lo);
				metadata.nNetFastCount[iCombo]		= metadata.nFastCount[iCombo]		- metadata.nOopsCount[iCombo] * (tof_R_fast_hi - tof_R_fast_lo)     / (tof_oops_hi - tof_oops_lo);
				metadata.nNetSlowCount[iCombo]		= metadata.nSlowCount[iCombo]		- metadata.nOopsCount[iCombo] * (tof_R_slow_hi - tof_R_slow_lo)     / (tof_oops_hi - tof_oops_lo);
				// Trap Empty
				metadata.nZeroTOFBkgdIntegral[iCombo]	= h_bkgd_tof_LR->Integral(tofBin_zero_lo,   tofBin_zero_hi);
				metadata.nLowTOFBkgdIntegral[iCombo]	= h_bkgd_tof_LR->Integral(tofBin_lowTOF_lo, tofBin_lowTOF_hi);
				metadata.nFastBkgdIntegral[iCombo]		= h_bkgd_tof_LR->Integral(tofBin_R_fast_lo, tofBin_R_fast_hi);
				metadata.nSlowBkgdIntegral[iCombo]		= h_bkgd_tof_LR->Integral(tofBin_R_slow_lo, tofBin_R_slow_hi);
				metadata.nOopsBkgdIntegral[iCombo]		= h_bkgd_tof_LR->Integral(tofBin_oops_lo,   tofBin_oops_hi);
				oopsPerNsBkgd							= metadata.nOopsBkgdIntegral[iCombo] / (tof_oops_hi - tof_oops_lo);
				metadata.nNetFastBkgdIntegral[iCombo]	= metadata.nFastBkgdIntegral[iCombo]	- oopsPerNsBkgd * (tof_R_fast_hi - tof_R_fast_lo);
				metadata.nNetSlowBkgdIntegral[iCombo]	= metadata.nSlowBkgdIntegral[iCombo]	- oopsPerNsBkgd * (tof_R_slow_hi - tof_R_slow_lo);
				metadata.nNetFastBkgdCount[iCombo]		= metadata.nFastBkgdCount[iCombo]		- metadata.nOopsBkgdCount[iCombo] * (tof_R_fast_hi - tof_R_fast_lo)     / (tof_oops_hi - tof_oops_lo);
				metadata.nNetSlowBkgdCount[iCombo]		= metadata.nSlowBkgdCount[iCombo]		- metadata.nOopsBkgdCount[iCombo] * (tof_R_slow_hi - tof_R_slow_lo)     / (tof_oops_hi - tof_oops_lo);
				break;
			case BT:
				// Trap Full
				metadata.nZeroTOFIntegral[iCombo]	= h_tof_BT->Integral(tofBin_zero_lo,   tofBin_zero_hi);
				metadata.nLowTOFIntegral[iCombo]	= h_tof_BT->Integral(tofBin_lowTOF_lo, tofBin_lowTOF_hi);
				metadata.nFastIntegral[iCombo]		= h_tof_BT->Integral(tofBin_T_fast_lo, tofBin_T_fast_hi);
				metadata.nSlowIntegral[iCombo]		= h_tof_BT->Integral(tofBin_T_slow_lo, tofBin_T_slow_hi);
				metadata.nOopsIntegral[iCombo]		= h_tof_BT->Integral(tofBin_oops_lo,   tofBin_oops_hi);
				oopsPerNs							= metadata.nOopsIntegral[iCombo] / (tof_oops_hi - tof_oops_lo);
				metadata.nNetFastIntegral[iCombo]	= metadata.nFastIntegral[iCombo]	- oopsPerNs * (tof_T_fast_hi - tof_T_fast_lo);
				metadata.nNetSlowIntegral[iCombo]	= metadata.nSlowIntegral[iCombo]	- oopsPerNs * (tof_T_slow_hi - tof_T_slow_lo);
				metadata.nNetFastCount[iCombo]		= metadata.nFastCount[iCombo]		- metadata.nOopsCount[iCombo] * (tof_T_fast_hi - tof_T_fast_lo)     / (tof_oops_hi - tof_oops_lo);
				metadata.nNetSlowCount[iCombo]		= metadata.nSlowCount[iCombo]		- metadata.nOopsCount[iCombo] * (tof_T_slow_hi - tof_T_slow_lo)     / (tof_oops_hi - tof_oops_lo);
				// Trap Empty
				metadata.nZeroTOFBkgdIntegral[iCombo]	= h_bkgd_tof_BT->Integral(tofBin_zero_lo,   tofBin_zero_hi);
				metadata.nLowTOFBkgdIntegral[iCombo]	= h_bkgd_tof_BT->Integral(tofBin_lowTOF_lo, tofBin_lowTOF_hi);
				metadata.nFastBkgdIntegral[iCombo]		= h_bkgd_tof_BT->Integral(tofBin_T_fast_lo, tofBin_T_fast_hi);
				metadata.nSlowBkgdIntegral[iCombo]		= h_bkgd_tof_BT->Integral(tofBin_T_slow_lo, tofBin_T_slow_hi);
				metadata.nOopsBkgdIntegral[iCombo]		= h_bkgd_tof_BT->Integral(tofBin_oops_lo,   tofBin_oops_hi);
				oopsPerNsBkgd							= metadata.nOopsBkgdIntegral[iCombo] / (tof_oops_hi - tof_oops_lo);
				metadata.nNetFastBkgdIntegral[iCombo]	= metadata.nFastBkgdIntegral[iCombo]	- oopsPerNsBkgd * (tof_T_fast_hi - tof_T_fast_lo);
				metadata.nNetSlowBkgdIntegral[iCombo]	= metadata.nSlowBkgdIntegral[iCombo]	- oopsPerNsBkgd * (tof_T_slow_hi - tof_T_slow_lo);
				metadata.nNetFastBkgdCount[iCombo]		= metadata.nFastBkgdCount[iCombo]		- metadata.nOopsBkgdCount[iCombo] * (tof_T_fast_hi - tof_T_fast_lo)     / (tof_oops_hi - tof_oops_lo);
				metadata.nNetSlowBkgdCount[iCombo]		= metadata.nSlowBkgdCount[iCombo]		- metadata.nOopsBkgdCount[iCombo] * (tof_T_slow_hi - tof_T_slow_lo)     / (tof_oops_hi - tof_oops_lo);
				break;
			case BR:
				// Trap Full
				metadata.nZeroTOFIntegral[iCombo]	= h_tof_BR->Integral(tofBin_zero_lo,   tofBin_zero_hi);
				metadata.nLowTOFIntegral[iCombo]	= h_tof_BR->Integral(tofBin_lowTOF_lo, tofBin_lowTOF_hi);
				metadata.nFastIntegral[iCombo]		= h_tof_BR->Integral(tofBin_R_fast_lo, tofBin_R_fast_hi);
				metadata.nSlowIntegral[iCombo]		= h_tof_BR->Integral(tofBin_R_slow_lo, tofBin_R_slow_hi);
				metadata.nOopsIntegral[iCombo]		= h_tof_BR->Integral(tofBin_oops_lo,   tofBin_oops_hi);
				oopsPerNs							= metadata.nOopsIntegral[iCombo] / (tof_oops_hi - tof_oops_lo);
				metadata.nNetFastIntegral[iCombo]	= metadata.nFastIntegral[iCombo]	- oopsPerNs * (tof_R_fast_hi - tof_R_fast_lo);
				metadata.nNetSlowIntegral[iCombo]	= metadata.nSlowIntegral[iCombo]	- oopsPerNs * (tof_R_slow_hi - tof_R_slow_lo);
				metadata.nNetFastCount[iCombo]		= metadata.nFastCount[iCombo]		- metadata.nOopsCount[iCombo] * (tof_R_fast_hi - tof_R_fast_lo)     / (tof_oops_hi - tof_oops_lo);
				metadata.nNetSlowCount[iCombo]		= metadata.nSlowCount[iCombo]		- metadata.nOopsCount[iCombo] * (tof_R_slow_hi - tof_R_slow_lo)     / (tof_oops_hi - tof_oops_lo);
				// Trap Empty
				metadata.nZeroTOFBkgdIntegral[iCombo]	= h_bkgd_tof_BR->Integral(tofBin_zero_lo,   tofBin_zero_hi);
				metadata.nLowTOFBkgdIntegral[iCombo]	= h_bkgd_tof_BR->Integral(tofBin_lowTOF_lo, tofBin_lowTOF_hi);
				metadata.nFastBkgdIntegral[iCombo]		= h_bkgd_tof_BR->Integral(tofBin_R_fast_lo, tofBin_R_fast_hi);
				metadata.nSlowBkgdIntegral[iCombo]		= h_bkgd_tof_BR->Integral(tofBin_R_slow_lo, tofBin_R_slow_hi);
				metadata.nOopsBkgdIntegral[iCombo]		= h_bkgd_tof_BR->Integral(tofBin_oops_lo,   tofBin_oops_hi);
				oopsPerNsBkgd							= metadata.nOopsBkgdIntegral[iCombo] / (tof_oops_hi - tof_oops_lo);
				metadata.nNetFastBkgdIntegral[iCombo]	= metadata.nFastBkgdIntegral[iCombo]	- oopsPerNsBkgd * (tof_R_fast_hi - tof_R_fast_lo);
				metadata.nNetSlowBkgdIntegral[iCombo]	= metadata.nSlowBkgdIntegral[iCombo]	- oopsPerNsBkgd * (tof_R_slow_hi - tof_R_slow_lo);
				metadata.nNetFastBkgdCount[iCombo]		= metadata.nFastBkgdCount[iCombo]		- metadata.nOopsBkgdCount[iCombo] * (tof_R_fast_hi - tof_R_fast_lo)     / (tof_oops_hi - tof_oops_lo);
				metadata.nNetSlowBkgdCount[iCombo]		= metadata.nSlowBkgdCount[iCombo]		- metadata.nOopsBkgdCount[iCombo] * (tof_R_slow_hi - tof_R_slow_lo)     / (tof_oops_hi - tof_oops_lo);
				break;
		} // end switch over detector combos
		
		// Trap full
//		oopsPerNs								= metadata.nOopsIntegral[iCombo] / (tof_oops_hi - tof_oops_lo);
		metadata.nNetZeroTOFIntegral[iCombo]	= metadata.nZeroTOFIntegral[iCombo] - oopsPerNs * (tof_zero_hi   - tof_zero_lo);
		metadata.nNetLowTOFIntegral[iCombo]		= metadata.nLowTOFIntegral[iCombo]	- oopsPerNs * (tof_lowTOF_hi - tof_lowTOF_lo);
//		metadata.nNetFastIntegral[iCombo]		= metadata.nFastIntegral[iCombo]	- oopsPerNs * (tof_T_fast_hi - tof_T_fast_lo);
//		metadata.nNetSlowIntegral[iCombo]		= metadata.nSlowIntegral[iCombo]	- oopsPerNs * (tof_T_slow_hi - tof_T_slow_lo);
		// Trap empty
//		oopsPerNsBkgd								= metadata.nOopsBkgdIntegral[iCombo] / (tof_oops_hi - tof_oops_lo);
		metadata.nNetZeroTOFBkgdIntegral[iCombo]	= metadata.nZeroTOFBkgdIntegral[iCombo] - oopsPerNsBkgd * (tof_zero_hi   - tof_zero_lo);
		metadata.nNetLowTOFBkgdIntegral[iCombo]		= metadata.nLowTOFBkgdIntegral[iCombo]	- oopsPerNsBkgd * (tof_lowTOF_hi - tof_lowTOF_lo);
//		metadata.nNetFastBkgdIntegral[iCombo]		= metadata.nFastBkgdIntegral[iCombo]	- oopsPerNsBkgd * (tof_T_fast_hi - tof_T_fast_lo);
//		metadata.nNetSlowBkgdIntegral[iCombo]		= metadata.nSlowBkgdIntegral[iCombo]	- oopsPerNsBkgd * (tof_T_slow_hi - tof_T_slow_lo);
	} // for loop over detector combos
/*	
//	For some reason this doesn't work.  The same code in a .c macro gives sensible answers, but this gives 0 or ~1 billion
// Fast counts
//	int tof1_bin = floor(TOFBins*(fast_tof1 - TOFMin)/(TOFMax - TOFMin));
//	int tof2_bin = floor(TOFBins*(fast_tof2 - TOFMin)/(TOFMax - TOFMin));
//	metadata.n_fast_LT	= h_tof_LT->Integral(tof1_bin,tof2_bin) - h_bkgd_tof_LT->Integral(tof1_bin,tof2_bin);
//	metadata.n_fast_LR	= h_tof_LR->Integral(tof1_bin,tof2_bin) - h_bkgd_tof_LR->Integral(tof1_bin,tof2_bin);
//	metadata.n_fast_BT	= h_tof_BT->Integral(tof1_bin,tof2_bin) - h_bkgd_tof_BT->Integral(tof1_bin,tof2_bin);
//	metadata.n_fast_BR	= h_tof_BR->Integral(tof1_bin,tof2_bin) - h_bkgd_tof_BR->Integral(tof1_bin,tof2_bin);
//	
//	// Slow counts
//	tof1_bin = ceil(TOFBins*(slow_tof1 - TOFMin)/(TOFMax - TOFMin));
//	tof2_bin = ceil(TOFBins*(slow_tof2 - TOFMin)/(TOFMax - TOFMin));
//	metadata.n_slow_LT	= h_tof_LT->Integral(tof1_bin,tof2_bin) - h_bkgd_tof_LT->Integral(tof1_bin,tof2_bin);
//	metadata.n_slow_LR	= h_tof_LR->Integral(tof1_bin,tof2_bin) - h_bkgd_tof_LR->Integral(tof1_bin,tof2_bin);
//	metadata.n_slow_BT	= h_tof_BT->Integral(tof1_bin,tof2_bin) - h_bkgd_tof_BT->Integral(tof1_bin,tof2_bin);
//	metadata.n_slow_BR	= h_tof_BR->Integral(tof1_bin,tof2_bin) - h_bkgd_tof_BR->Integral(tof1_bin,tof2_bin);
*/	
	metadata_Tree->Fill();
	
//~~~~~~~~ PRINT-OUT ~~~~~~~//
	
	printf("\nRun %d", n_run);
	cout << "\nSaved ROOT file: " << rootFileName;
	printf("\nRun start time: d:%02d, h:%02d, m:%02d, s:%02d, t = %d",start_day,start_hour,start_min,start_sec,start_time_sec);
	if (stop_flag == 0)	\
		printf("\nRun is ongoing");
	else \
		printf("\nRun stop time:  d:%02d, h:%02d, m:%02d, s:%02d, t = %d", stop_day, stop_hour, stop_min, stop_sec, stop_time_sec);
	printf("\nRun time:       %d seconds",run_time_sec);
	cout<<endl<< "# trigs:"<< setw(12) << right << n_trig;
	cout<<endl<< "# syncs:"<< setw(12) << right << n_sync;
	cout<<endl<< "# bad events:"<< setw(12) << right << n_bad_events << endl;
//	cout<<endl<< "# non-vetoed trigs:"<< setw(12) << tot_trigs;
	printf("\n~~~~~~~~~~~~~ RATES (Hz) ~~~~~~~~~~~~~~");
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	printf("\n                 Scaler     TDC     ADC");
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	printf("\nBottom dEa   : %8.2f%8.2f%8.2f", sTot_B_dEa	/ (float)sync_time_sec, nt_B_dEa	/ (float)run_time_sec, na_B_dEa	/ (float)run_time_sec);
	printf("\nBottom dEb   : %8.2f%8.2f%8.2f", sTot_B_dEb	/ (float)sync_time_sec, nt_B_dEb	/ (float)run_time_sec, na_B_dEb	/ (float)run_time_sec);
	printf("\nBottom E     : %8.2f%8.2f%8.2f", sTot_B_E		/ (float)sync_time_sec, nt_B_E		/ (float)run_time_sec, na_B_E	/ (float)run_time_sec);
	printf("\nLeft   dEa   : %8.2f%8.2f%8.2f", sTot_L_dEa	/ (float)sync_time_sec, nt_L_dEa	/ (float)run_time_sec, na_L_dEa	/ (float)run_time_sec);
	printf("\nLeft   dEb   : %8.2f%8.2f%8.2f", sTot_L_dEb	/ (float)sync_time_sec, nt_L_dEb	/ (float)run_time_sec, na_L_dEb	/ (float)run_time_sec);
	printf("\nLeft   E     : %8.2f%8.2f%8.2f", sTot_L_E		/ (float)sync_time_sec, nt_L_E		/ (float)run_time_sec, na_L_E	/ (float)run_time_sec);
	printf("\nRight  Ge    : %8.2f%8.2f%8.2f", sTot_R_ge	/ (float)sync_time_sec, nt_R_ge		/ (float)run_time_sec, na_R_ge	/ (float)run_time_sec);
	printf("\nTop    Ge    : %8.2f%8.2f%8.2f", sTot_T_ge	/ (float)sync_time_sec, nt_T_ge		/ (float)run_time_sec, na_T_ge	/ (float)run_time_sec);
	printf("\nRight  MCP   : %8.2f%8.2f"	 , sTot_R_mcp	/ (float)sync_time_sec, nt_R_mcp	/ (float)run_time_sec);
	printf("\nRight  MCP A :                 %8.2f", na_R_mcpA	/ (float)run_time_sec);
	printf("\nRight  MCP B :                 %8.2f", na_R_mcpB	/ (float)run_time_sec);
	printf("\nRight  MCP C :                 %8.2f", na_R_mcpC	/ (float)run_time_sec);
	printf("\nRight  MCP D :                 %8.2f", na_R_mcpD	/ (float)run_time_sec);
	printf("\nTop    MCP   : %8.2f%8.2f"	 , sTot_T_mcp	/ (float)sync_time_sec, nt_T_mcp	/ (float)run_time_sec);
	printf("\nTop    MCP A :                 %8.2f", na_T_mcpA	/ (float)run_time_sec);
	printf("\nTop    MCP B :                 %8.2f", na_T_mcpB	/ (float)run_time_sec);
	printf("\nTop    MCP C :                 %8.2f", na_T_mcpC	/ (float)run_time_sec);
	printf("\nTop    MCP D :                 %8.2f", na_T_mcpD	/ (float)run_time_sec);
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	printf("\nSums         : %8.2f%8.2f"	 , sTot_all		/ (float)sync_time_sec, nt_all		/ (float)run_time_sec);
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	printf("\nScarlet      :         %8.2f", n_trig	/ (float)run_time_sec);
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	printf("\nSiX4         : %8.2f"			 , sTot_SiX4_ts	/ (float)run_time_sec);
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	printf("\nAverage Dead Time per Event (us)  : %8.2f", ((float)(tot_runTime_us - tot_liveTime_us))/(1.0*n_trig));
	printf("\nLive Time (s)                     : %8.2f", ((float)tot_liveTime_us)/1000000.0);
	printf("\nTotal run time (s)                : %8.2f", ((float)tot_runTime_us)/1000000.0);
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	printf("\ndE-MCP Coincidences by counting events:");
	printf("\n  LT-- Zero-TOF: %6d; \"Low-TOF\": %6d; Fast ions: %6d; Slow ions: %6d; Accidentals: %6d", Nint(metadata.nZeroTOFCount[LT]), Nint(metadata.nLowTOFCount[LT]), Nint(metadata.nFastCount[LT]), Nint(metadata.nSlowCount[LT]), Nint(metadata.nOopsCount[LT]));
	printf("\n  LR-- Zero-TOF: %6d; \"Low-TOF\": %6d; Fast ions: %6d; Slow ions: %6d; Accidentals: %6d", Nint(metadata.nZeroTOFCount[LR]), Nint(metadata.nLowTOFCount[LR]), Nint(metadata.nFastCount[LR]), Nint(metadata.nSlowCount[LR]), Nint(metadata.nOopsCount[LR]));
	printf("\n  BT-- Zero-TOF: %6d; \"Low-TOF\": %6d; Fast ions: %6d; Slow ions: %6d; Accidentals: %6d", Nint(metadata.nZeroTOFCount[BT]), Nint(metadata.nLowTOFCount[BT]), Nint(metadata.nFastCount[BT]), Nint(metadata.nSlowCount[BT]), Nint(metadata.nOopsCount[BT]));
	printf("\n  BR-- Zero-TOF: %6d; \"Low-TOF\": %6d; Fast ions: %6d; Slow ions: %6d; Accidentals: %6d", Nint(metadata.nZeroTOFCount[BR]), Nint(metadata.nLowTOFCount[BR]), Nint(metadata.nFastCount[BR]), Nint(metadata.nSlowCount[BR]), Nint(metadata.nOopsCount[BR]));
	printf("\ndE-MCP Coincidences by h_tof integrals:");
	printf("\n  LT-- Zero-TOF: %6d; \"Low-TOF\": %6d; Fast ions: %6d; Slow ions: %6d; Accidentals: %6d", Nint(metadata.nZeroTOFIntegral[LT]), Nint(metadata.nLowTOFIntegral[LT]), Nint(metadata.nFastIntegral[LT]), Nint(metadata.nSlowIntegral[LT]), Nint(metadata.nOopsIntegral[LT]));
	printf("\n  LR-- Zero-TOF: %6d; \"Low-TOF\": %6d; Fast ions: %6d; Slow ions: %6d; Accidentals: %6d", Nint(metadata.nZeroTOFIntegral[LR]), Nint(metadata.nLowTOFIntegral[LR]), Nint(metadata.nFastIntegral[LR]), Nint(metadata.nSlowIntegral[LR]), Nint(metadata.nOopsIntegral[LR]));
	printf("\n  BT-- Zero-TOF: %6d; \"Low-TOF\": %6d; Fast ions: %6d; Slow ions: %6d; Accidentals: %6d", Nint(metadata.nZeroTOFIntegral[BT]), Nint(metadata.nLowTOFIntegral[BT]), Nint(metadata.nFastIntegral[BT]), Nint(metadata.nSlowIntegral[BT]), Nint(metadata.nOopsIntegral[BT]));
	printf("\n  BR-- Zero-TOF: %6d; \"Low-TOF\": %6d; Fast ions: %6d; Slow ions: %6d; Accidentals: %6d", Nint(metadata.nZeroTOFIntegral[BR]), Nint(metadata.nLowTOFIntegral[BR]), Nint(metadata.nFastIntegral[BR]), Nint(metadata.nSlowIntegral[BR]), Nint(metadata.nOopsIntegral[BR]));
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	printf("\nNet fast and slow ions:");
	printf("\n  LT-- Net fast by counting: %8.2f; Net fast by integrals: %8.2f; Net slow by counting: %8.2f; Net slow by integrals: %8.2f", metadata.nNetFastCount[LT],  metadata.nNetFastIntegral[LT], metadata.nNetSlowCount[LT], metadata.nNetSlowIntegral[LT]);
	printf("\n  LR-- Net fast by counting: %8.2f; Net fast by integrals: %8.2f; Net slow by counting: %8.2f; Net slow by integrals: %8.2f", metadata.nNetFastCount[LR],  metadata.nNetFastIntegral[LR], metadata.nNetSlowCount[LR], metadata.nNetSlowIntegral[LR]);
	printf("\n  BT-- Net fast by counting: %8.2f; Net fast by integrals: %8.2f; Net slow by counting: %8.2f; Net slow by integrals: %8.2f", metadata.nNetFastCount[BT],  metadata.nNetFastIntegral[BT], metadata.nNetSlowCount[BT], metadata.nNetSlowIntegral[BT]);
	printf("\n  BR-- Net fast by counting: %8.2f; Net fast by integrals: %8.2f; Net slow by counting: %8.2f; Net slow by integrals: %8.2f", metadata.nNetFastCount[BR],  metadata.nNetFastIntegral[BR], metadata.nNetSlowCount[BR], metadata.nNetSlowIntegral[BR]);
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	printf("\nDone sorting run %d", n_run);
	
	//cout<<endl<< "# cycles by syncs:     "<< sync_bin;
	//cout<<endl<< "# cycles by triggers:  "<< trig_bin;
	//cout<<	"# cycles by triggers:  "<< trig_bin	<< endl;
	
	cout<<endl<<endl;
	
	f->Write();
	f->Close();
	delete esrc;
	return 0;
	
} // main

// Function to convert dhms time to seconds:
// Could be improved by allowing the smallest unit (sec) to be a float
int time_in_seconds(int day, int hour, int min, int sec) {
	int t = 0;
	t += 60*60*24*day;
	t += 60*60*hour;
	t += 60*min;
	t += sec;
	return t;
}

int countbit(int x) {
	int n=0;
	if (x==0) return 0;
	else {
		if (x & 1) n++;
		if (x &	2) n++;
		if (x & 4) n++;
		if (x & 8) n++;
		if (x & 16) n++;
		if (x & 32) n++;
		if (x & 64) n++;
		if (x & 128) n++;
		if (x & 256) n++;
		if (x & 512) n++;
		if (x & 1024) n++;
		if (x & 2048) n++;
		if (x & 0x1000) n++;
		if (x & 0x2000) n++;
		if (x & 0x4000) n++;
		if (x & 0x8000) n++;
	}
	return n;
}

int ReadADC1(int **p,int n_trig, int *event_good,int *n_bad_events) {
// ADC1 *******************************
	int x, adc_ch, wordc;
	if (int(*(*p)) != 0xadc1adc1) {
		cout << "trig #" << n_trig << ", ADC1 marker not found where expected!" << endl;
		(*event_good) = 0;
		(*n_bad_events)++;
		return -1;
	}
	(*p)++; // move to ADC1 hit register
	x = int(int(*(*p)) & 0xffff ); // hit register, tells which channels were hit
	wordc = countbit(x);
	
	for (int j=1; j<=wordc; j++) { // Loop over all ADC channels which have hits
		(*p)++;  // Increment pointer p to ADC channel with a hit
		x = int((*(*p)) & 0x0fff);
	
		// Get ADC channel of that hit then associate it with the data:
		adc_ch=int((*(*p)) & 0xf000);
		adc_ch=(adc_ch>>12)+1;
		//if (adc_ch == 1) {
		//	a_R_ge = x;
		//	ha_R_ge->Fill(x);
		//	na_R_ge++;
		//}
		switch(adc_ch) {
		case 1:
			bdn.a_T_mcpE = x;
			ha_T_mcpE->Fill(x);
			ha_T_mcpE_corr->Fill(x - ped_T_mcpE + randgen->Rndm());
			metadata.n_adc_hits_T_mcpE++;
		break;
		case 2: 
			bdn.a_B_dEa = x;
			ha_B_dEa->Fill(x);
			metadata.n_adc_hits_B_dEa++;
		break;
		case 3:
			bdn.a_B_dEb = x;
			ha_B_dEb->Fill(x);
			metadata.n_adc_hits_B_dEb++;
		break;
		case 4: 
			bdn.a_B_E = x;
			ha_B_E->Fill(x);
			metadata.n_adc_hits_B_E++;
		break;
		case 5:
			bdn.a_T_mcpA = x;
			ha_T_mcpA->Fill(x);
			ha_T_mcpA_corr->Fill(x - ped_T_mcpA + randgen->Rndm());
			metadata.n_adc_hits_T_mcpA++;
		break;
		case 6:
			bdn.a_T_mcpB = x;
			ha_T_mcpB->Fill(x);
			ha_T_mcpB_corr->Fill(x - ped_T_mcpB + randgen->Rndm());
			metadata.n_adc_hits_T_mcpB++;
		break;
		case 7:
			bdn.a_T_mcpC = x;
			ha_T_mcpC->Fill(x);
			ha_T_mcpC_corr->Fill(x - ped_T_mcpC + randgen->Rndm());
			metadata.n_adc_hits_T_mcpC++;
		break;
		case 8:
			bdn.a_T_mcpD = x;
			ha_T_mcpD->Fill(x);
			ha_T_mcpD_corr->Fill(x - ped_T_mcpD + randgen->Rndm());
			metadata.n_adc_hits_T_mcpD++;
		break;
		//if (adc_ch == 9) {
		//	a_T_ge = x;
		//	ha_T_ge->Fill(x);
		//	na_T_ge++;
		//}
		case 9:
			bdn.a_R_mcpE = x;
			ha_R_mcpE->Fill(x);
			ha_R_mcpE_corr->Fill(x - ped_R_mcpE + randgen->Rndm());
			metadata.n_adc_hits_R_mcpE++;
		break;
		case 10:
			bdn.a_L_dEa = x;
			ha_L_dEa->Fill(x);
			metadata.n_adc_hits_L_dEa++;
		break;
		case 11:
			bdn.a_L_dEb = x;
			ha_L_dEb->Fill(x);
			metadata.n_adc_hits_L_dEb++;
		break;
		case 12:
			bdn.a_L_E = x;
			ha_L_E->Fill(x);
			metadata.n_adc_hits_L_E++;
		break;
		case 13:
			bdn.a_R_mcpA = x;
			ha_R_mcpA->Fill(x);
			ha_R_mcpA_corr->Fill(x - ped_R_mcpA + randgen->Rndm());
			metadata.n_adc_hits_R_mcpA++;
		break;
		case 14:
			bdn.a_R_mcpB = x;
			ha_R_mcpB->Fill(x);
			ha_R_mcpB_corr->Fill(x - ped_R_mcpB + randgen->Rndm());
			metadata.n_adc_hits_R_mcpB++;
		break;
		case 15:
			bdn.a_R_mcpC = x;
			ha_R_mcpC->Fill(x);
			ha_R_mcpC_corr->Fill(x - ped_R_mcpC + randgen->Rndm());
			metadata.n_adc_hits_R_mcpC++;
		break;
		case 16:
			bdn.a_R_mcpD = x;
			ha_R_mcpD->Fill(x);
			ha_R_mcpD_corr->Fill( x - ped_R_mcpD + randgen->Rndm());
			metadata.n_adc_hits_R_mcpD++;
		break;
		default:
		break;
		}
	} // for (wordc)
	return 0;
}


int ReadADC2(int **p,int n_trig, int *event_good,int *n_bad_events) {
// ADC1 *******************************
	int x,y, adc_ch, wordc;
	if (int(*(*p)) != 0xadc2adc2) {
		cout << "trig #" << n_trig << ", ADC2 marker not found where expected!" << endl;
		(*event_good) = 0;
		(*n_bad_events)++;
		return -1;
	}
	(*p)++; // move to ADC1 hit register
	x = int(int(*(*p)) & 0xffff ); // hit register, tells which channels were hit
	wordc = countbit(x);
	
	for (int j=1; j<=wordc; j++) { // Loop over all ADC channels which have hits
		(*p)++;  // Increment pointer p to ADC channel with a hit
		x = int((*(*p)) & 0x0fff);
	
		// Get ADC channel of that hit then associate it with the data:
		adc_ch=int((*(*p)) & 0xf000);
		adc_ch=(adc_ch>>12)+1;
		//if (adc_ch == 1) {
		//	a_R_ge = x;
		//	ha_R_ge->Fill(x);
		//	na_R_ge++;
		//}
		switch(adc_ch) {
		case 1:
			y = x + randgen->Rndm();
			bdn.a_T_ge_highE = x;
			ha_T_ge_highE->Fill(T_ge_highE_coeff[0] + y*T_ge_highE_coeff[1] + y*y*T_ge_highE_coeff[2]);
			//he_T_ge_highE->Fill(e_T_ge_highE);
			//he_ge_highE	 ->Fill(e_T_ge_highE);
			metadata.n_adc_hits_T_ge_highE++;
		break;
		case 2: 
			y = x + randgen->Rndm();
			bdn.a_R_ge_highE = x;
			ha_R_ge_highE->Fill(R_ge_highE_coeff[0] + y*R_ge_highE_coeff[1] + y*y*R_ge_highE_coeff[2]);
			//he_T_ge_highE->Fill(e_T_ge_highE);
			//he_ge_highE	 ->Fill(e_T_ge_highE);
			metadata.n_adc_hits_R_ge_highE++;
		break;
		case 7:
			y = x + randgen->Rndm();
			bdn.a_T_ge = x;
			ha_T_ge->Fill(T_ge_coeff[0] + y*T_ge_coeff[1] + y*y*T_ge_coeff[2]);
			//he_T_ge_highE->Fill(e_T_ge_highE);
			//he_ge_highE	 ->Fill(e_T_ge_highE);
			metadata.n_adc_hits_T_ge++;
		break;
		case 8:
			y = x + randgen->Rndm();
			bdn.a_R_ge = x;
			ha_R_ge->Fill(R_ge_coeff[0] + y*R_ge_coeff[1] + y*y*R_ge_coeff[2]);
			//he_T_ge_highE->Fill(e_T_ge_highE);
			//he_ge_highE	 ->Fill(e_T_ge_highE);
			metadata.n_adc_hits_R_ge++;
		break;
		//if (adc_ch == 9) {
		//	a_T_ge = x;
		//	ha_T_ge->Fill(x);
		//	na_T_ge++;
		//}
		case 9:
			y = x + randgen->Rndm();
			bdn.a_T_ge = x;
			ha_T_ge->Fill(T_ge_coeff[0] + y*T_ge_coeff[1] + y*y*T_ge_coeff[2]);
			//he_T_ge_highE->Fill(e_T_ge_highE);
			//he_ge_highE	 ->Fill(e_T_ge_highE);
			metadata.n_adc_hits_T_ge++;
		break;
		default:
		break;
		}
	} // for (wordc)
	return 0;
}

int ReadTDC1(int **p, int n_trig, int *event_good, int *n_bad_events){
	int x, tdc_ch;
	(*p)++;
	if(**p != 0x2dc12dc1) {
		cout << "trig #" << n_trig << ", TDC1 marker not found where expected!" << endl;
		(*event_good) = 0;
		(*n_bad_events)++;
		return -1;
	}
	(*p)++;
	while(**p != 0x2dc22dc2) {
		tdc_ch=**p;
		(*p)++;
		x = int(**p & 0x00ffffff); // take only the 24-bit data word
		if (x & 0x0080000) x -= 0x00ffffff; // test for neg value
			// if neg then you need to shift because the leading 1 in 24-bit
			// is not leading in 32-bit; the shift is by "-0x00ffffff"
		switch(tdc_ch) {
		case 1:
			bdn.t_T_mcp = x;
			ht_T_mcp->Fill(x);
			metadata.n_tdc_hits_T_mcp++;
		break;
		case 2:
			bdn.t_R_mcp = x;
			ht_R_mcp->Fill(x);
			metadata.n_tdc_hits_R_mcp++;
		break;
		case 3:
			bdn.t_B_dEa = x;
			ht_B_dEa->Fill(x);
			metadata.n_tdc_hits_B_dEa++;
		break;
		case 4:
			bdn.t_B_dEb = x;
			ht_B_dEb->Fill(x);
			metadata.n_tdc_hits_B_dEb++;
		break;
		case 5:
			bdn.t_B_E = x;
			ht_B_E->Fill(x);
			metadata.n_tdc_hits_B_E++;
		break;
		case 6:
			bdn.t_L_dEa = x;
			ht_L_dEa->Fill(x);
			metadata.n_tdc_hits_L_dEa++;
		break;
		case 7:
			bdn.t_L_dEb = x;
			ht_L_dEb->Fill(x);
			metadata.n_tdc_hits_L_dEb++;
		break;
		case 8:
			bdn.t_L_E = x;
			ht_L_E->Fill(x);
			metadata.n_tdc_hits_L_E++;
		break;
		default:
		break;
		}
	}
}

int ReadTDC2(int **p, int n_trig, int *event_good, int *n_bad_events){
	int x, tdc_ch;
	if(**p != 0x2dc22dc2) {
		cout << "trig #" << n_trig << ", TDC2 marker not found where expected!" << endl;
		(*event_good) = 0;
		(*n_bad_events)++;
		return -1;
	}
	(*p)++;
	while(**p != 0x100cca1e) {
		tdc_ch=**p;
		(*p)++;
		x = int(**p & 0x00ffffff); // take only the 24-bit data word
		if (x & 0x0080000) x -= 0x00ffffff; // test for neg value
			// if neg then you need to shift because the leading 1 in 24-bit
			// is not leading in 32-bit; the shift is by "-0x00ffffff"
		switch(tdc_ch) {
		case 1:
			bdn.t_rf = x;
			ht_rf->Fill(x);
		break;
		case 2:
			bdn.t_T_ge = x;
			ht_T_ge->Fill(x);
			metadata.n_tdc_hits_T_ge++;
		break;
		case 3:
			bdn.t_R_ge = x;
			ht_R_ge->Fill(x);
			metadata.n_tdc_hits_R_ge++;
		break;
		default:
		break;
		}
	}
}

int ReadScalers(int **p, int n_trig, int n_run,int *all_trigs, int *event_good, int *n_bad_events){
	if (n_run < 1201) { // old scaler readout
	// Capt Scaler ************************
		if (**p != 0x100cca1e) {
			cout << "trig #" << n_trig << ", Capt Scaler marker not found where expected!" << endl;
			event_good = 0;
			n_bad_events++;
			return -1;
		}
		(*p)++; // p is at time since capture in ms
		bdn.s_ms_since_capt = int(**p & 0xffffff);
		(*p)++; // p is at trap state | 0 = trap full, 1 = trap empty
		bdn.s_capt_state = int(**p & 0xffffff);
		(*p)++; // move pointer to eject scaler
		
	// Eject Scaler ***********************
		if (**p != 0x100eca1e) {
			cout << "trig #" << n_trig << ", Eject Scaler marker not found where expected!" << endl;
			event_good = 0;
			n_bad_events++;
			return -1;
		}
		(*p)++; //p is at time since eject in ms
		bdn.s_ms_since_eject = int(**p & 0xffffff);
		(*p)++; //p is at # of capt since last eject
		bdn.s_capt = int(**p & 0xffffff);
		(*p)++; //p is at # of SiX4 hits since last eject
		bdn.s_SiX4 = int(**p & 0xffffff);

	} // end old scaler readout

	else {	// new scaler readout

	// Live Time Scaler ***********************
		if (**p != 0x100cca1e) {
			cout << "trig #" << n_trig << ", Livetime Scaler marker not found where expected!" << endl;
			event_good = 0;
			n_bad_events++;
			return -1;
		}
		
		(*p)++; // p is..  this is a bit iffy
		bdn.s_liveTime_us = int(**p & 0xffffff);
		(*p)++;
		(*all_trigs) = int(**p & 0xffffff);
		(*p)++;
		bdn.s_runTime = int(**p & 0xffffff);
		//*p++; //p is at # of SiX4 hits since last eject
		//s_SiX4 = int(*p & 0xffffff);
		
		(*p)++; //this also iffy, check if it works?
	// Capt Scaler ************************
		if (**p != 0x100dca1e) {
			cout << "trig #" << n_trig << ", Capt Scaler marker not found where expected!" << endl;
			event_good = 0;
			n_bad_events++;
			return -1;
		}
		(*p)++; // p is at time since capture in ms
		bdn.s_ms_since_capt = int(**p & 0xffffff);
		(*p)++; // p is at trap state | 0 = trap full, 1 = trap empty
		bdn.s_capt_state = int(**p & 0xffffff);
		(*p)++; // move pointer to eject scaler
		
	// Eject Scaler ***********************
		if (**p != 0x100eca1e) {
			cout << "trig #" << n_trig << ", Eject Scaler marker not found where expected!" << endl;
			event_good = 0;
			n_bad_events++;
			return -1;
		}
		(*p)++; //p is at time since eject in ms
		bdn.s_ms_since_eject = int(**p & 0xffffff);
		(*p)++; //p is at # of capt since last eject
		bdn.s_capt = int(**p & 0xffffff);
		//*p++; //p is at # of SiX4 hits since last eject
		//s_SiX4 = int(*p & 0xffffff);
		
	} // end new scaler readout
	return 0;
}
/*
// void sync_sort(const struct ScarletEvntHdr *e) {
	
	// e0 = ScarletEvnt(e);
	// e1 = e0[1];
	// p = reinterpret_cast<int*>(e1.body());
	// //cout << "first word: " << p << endl;
	// *p++;
	
	// s_T_mcp = (*p++ & 0xffffff);
	// s_R_mcp = (*p++ & 0xffffff);
	// s_T_ge = (*p++ & 0xffffff);
	// s_R_ge = (*p++ & 0xffffff);
	// s_L_dEa = (*p++ & 0xffffff);
	// s_L_dEb = (*p++ & 0xffffff);
	// s_L_E = (*p++ & 0xffffff);
	// s_B_dEa = (*p++ & 0xffffff);
	// s_B_dEb = (*p++ & 0xffffff);
	// s_B_E = (*p++ & 0xffffff);
	
	// sCyc_T_mcp += s_T_mcp;
	// sCyc_R_mcp += s_R_mcp;
	// sCyc_T_ge += s_T_ge;
	// sCyc_R_ge += s_R_ge;
	// sCyc_L_dEa += s_L_dEa;
	// sCyc_L_dEb += s_L_dEb;
	// sCyc_L_E += s_L_E;
	// sCyc_B_dEa += s_B_dEa;
	// sCyc_B_dEb += s_B_dEb;
	// sCyc_B_E += s_B_E;
	
// } // sync_sort

// void event_sort(const struct ScarletEvtHdr *e) {
	
	// e0 = ScarletEvnt(e);
	// e1 = e0[1];
	// p = reinterpret_cast<int*>(e1.body());

	// n_trig++; //# of events
	
// } // event_sort
*/