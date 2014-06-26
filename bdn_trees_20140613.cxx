#define _bdn_trees_20140613_cxx "bdn_trees_20140613.cxx"
#include "bdn_trees_20140613.h"

void book_trees()
{

list_LT = new TEventList("list_LT");
list_LR = new TEventList("list_LR");
list_BT = new TEventList("list_BT");
list_BR = new TEventList("list_BR");
list_bkgd_LT = new TEventList("list_bkgd_LT");
list_bkgd_LR = new TEventList("list_bkgd_LR");
list_bkgd_BT = new TEventList("list_bkgd_BT");
list_bkgd_BR = new TEventList("list_bkgd_BR");

tree_LT = new TTree("tree_LT","all data that looks like a coincidence");
tree_LR = new TTree("tree_LR","all data that looks like a coincidence");
tree_BT = new TTree("tree_BT","all data that looks like a coincidence");
tree_BR = new TTree("tree_BR","all data that looks like a coincidence");
tree_bkgd_LT = new TTree("tree_bkgd_LT","all data that looks like a coincidence (BACKGROUND)");
tree_bkgd_LR = new TTree("tree_bkgd_LR","all data that looks like a coincidence (BACKGROUND)");
tree_bkgd_BT = new TTree("tree_bkgd_BT","all data that looks like a coincidence (BACKGROUND)");
tree_bkgd_BR = new TTree("tree_bkgd_BR","all data that looks like a coincidence (BACKGROUND)");

bdn_Tree 	= new TTree("bdn_Tree", "beta delayed neutron data");
metadata_Tree = new TTree("metadata_Tree","data about each file");
zero_time_Tree = new TTree("zero_time_Tree", "Events in zero-time TOF peaks");
beta_recoil_tree = new TTree("beta_recoil_tree","doubles with plastic and mcp data");
beta_gamma_tree = new TTree("beta_gamma_tree","doubles with plastic and hpge data");

metadata_Tree->Branch("metadata", &metadata, \
	"n_run/I:n_trigs:tot_trigs:n_syncs:n_treeEntries:n_bad_events:bkgd_good:"\
	"start_month:start_day:start_hour:start_min:start_sec:start_time_sec:"\
	"stop_month:stop_day:stop_hour:stop_min:stop_sec:stop_time_sec:run_time_sec:tot_liveTime_us:tot_runTime_us:"\
	"n_scaler_hits_B_dEa:n_scaler_hits_B_dEb:n_scaler_hits_B_E:"\
	"n_scaler_hits_L_dEa:n_scaler_hits_L_dEb:n_scaler_hits_L_E:"\
	"n_scaler_hits_R_mcp:n_scaler_hits_R_ge:"\
	"n_scaler_hits_T_mcp:n_scaler_hits_T_ge:"\
	"n_tdc_hits_B_dEa:n_tdc_hits_B_dEb:n_tdc_hits_B_E:"\
	"n_tdc_hits_L_dEa:n_tdc_hits_L_dEb:n_tdc_hits_L_E:"\
	"n_tdc_hits_R_mcp:n_tdc_hits_R_ge:"\
	"n_tdc_hits_T_mcp:n_tdc_hits_T_ge:"\
	"n_adc_hits_B_dEa:n_adc_hits_B_dEb:n_adc_hits_B_E:"\
	"n_adc_hits_L_dEa:n_adc_hits_L_dEb:n_adc_hits_L_E:"\
	"n_adc_hits_R_mcpA:n_adc_hits_R_mcpB:n_adc_hits_R_mcpC:n_adc_hits_R_mcpD:n_adc_hits_R_mcpE:n_adc_hits_R_ge:n_adc_hits_R_ge_highE:"\
	"n_adc_hits_T_mcpA:n_adc_hits_T_mcpB:n_adc_hits_T_mcpC:n_adc_hits_T_mcpD:n_adc_hits_T_mcpE:n_adc_hits_T_ge:n_adc_hits_T_ge_highE:"\
	"n_missing_adc_hits_R_mcpA:n_missing_adc_hits_R_mcpB:n_missing_adc_hits_R_mcpC:n_missing_adc_hits_R_mcpD:n_missing_adc_hits_R_mcpE:"\
	"n_missing_adc_hits_T_mcpA:n_missing_adc_hits_T_mcpB:n_missing_adc_hits_T_mcpC:n_missing_adc_hits_T_mcpD:n_missing_adc_hits_T_mcpE:"\
	"nZeroTOFCount[4]/F:nLowTOFCount[4]:nFastCount[4]:nSlowCount[4]:nOopsCount[4]:"\
	"nNetZeroTOFCount[4]:nNetLowTOFCount[4]:nNetFastCount[4]:nNetSlowCount[4]:"\
	"nZeroTOFBkgdCount[4]:nLowTOFBkgdCount[4]:nFastBkgdCount[4]:nSlowBkgdCount[4]:nOopsBkgdCount[4]:"\
	"nNetZeroTOFBkgdCount[4]:nNetLowTOFBkgdCount[4]:nNetFastBkgdCount[4]:nNetSlowBkgdCount[4]:"\
	"nZeroTOFIntegral[4]:nLowTOFIntegral[4]:nFastIntegral[4]:nSlowIntegral[4]:nOopsIntegral[4]:"\
	"nNetZeroTOFIntegral[4]:nNetLowTOFIntegral[4]:nNetFastIntegral[4]:nNetSlowIntegral[4]:"\
	"nZeroTOFBkgdIntegral[4]:nLowTOFBkgdIntegral[4]:nFastBkgdIntegral[4]:nSlowBkgdIntegral[4]:nOopsBkgdIntegral[4]:"\
	"nNetZeroTOFBkgdIntegral[4]:nNetLowTOFBkgdIntegral[4]:nNetFastBkgdIntegral[4]:nNetSlowBkgdIntegral[4]");
	
bdn_Tree->Branch("bdn", &bdn, "miss_R_mcpA/O:miss_R_mcpB:miss_R_mcpC:miss_R_mcpD:miss_T_mcpA:miss_T_mcpB:miss_T_mcpC:miss_T_mcpD:fid_area_hit_R_mcp:fid_area_hit_T_mcp:"\
    "a_R_ge/I:a_T_ge:a_R_ge_highE:a_T_ge_highE:"\
    "a_L_dEa:a_L_dEb:a_L_dEsum:a_L_E:"\
    "a_B_dEa:a_B_dEb:a_B_dEsum:a_B_E:"\
    "a_T_mcpA:a_T_mcpB:a_T_mcpC:a_T_mcpD:a_T_mcpE:a_T_mcpSum:"\
    "a_R_mcpA:a_R_mcpB:a_R_mcpC:a_R_mcpD:a_R_mcpE:a_R_mcpSum:"\
    "t_T_mcp:t_R_mcp:"\
    "t_B_dEa:t_B_dEb:t_B_E:"\
    "t_L_dEa:t_L_dEb:t_L_E:"\
    "t_rf:t_T_ge:t_R_ge:"\
    "s_ms_since_capt:s_capt_state:s_ms_since_eject:s_capt:deadTime_us:s_SiX4:event_good:event:run:"\
    "T_mcpX/D:T_mcpY:R_mcpX:R_mcpY:T_mcpPhysX:T_mcpPhysY:R_mcpPhysX:R_mcpPhysY:t_B_dE:t_L_dE:tof_LT:tof_LR:tof_BT:tof_BR:"\
    "v_LT:v_LR:v_BT:v_BR:En_LT:En_LR:En_BT:En_BR:rf_phase");
	
tree_LT->Branch("bdn", &bdn, "miss_R_mcpA/O:miss_R_mcpB:miss_R_mcpC:miss_R_mcpD:miss_T_mcpA:miss_T_mcpB:miss_T_mcpC:miss_T_mcpD:fid_area_hit_R_mcp:fid_area_hit_T_mcp:"\
    "a_R_ge/I:a_T_ge:a_R_ge_highE:a_T_ge_highE:"\
    "a_L_dEa:a_L_dEb:a_L_dEsum:a_L_E:"\
    "a_B_dEa:a_B_dEb:a_B_dEsum:a_B_E:"\
    "a_T_mcpA:a_T_mcpB:a_T_mcpC:a_T_mcpD:a_T_mcpE:a_T_mcpSum:"\
    "a_R_mcpA:a_R_mcpB:a_R_mcpC:a_R_mcpD:a_R_mcpE:a_R_mcpSum:"\
    "t_T_mcp:t_R_mcp:"\
    "t_B_dEa:t_B_dEb:t_B_E:"\
    "t_L_dEa:t_L_dEb:t_L_E:"\
    "t_rf:t_T_ge:t_R_ge:"\
    "s_ms_since_capt:s_capt_state:s_ms_since_eject:s_capt:deadTime_us:s_SiX4:event_good:event:run:"\
    "T_mcpX/D:T_mcpY:R_mcpX:R_mcpY:T_mcpPhysX:T_mcpPhysY:R_mcpPhysX:R_mcpPhysY:t_B_dE:t_L_dE:tof_LT:tof_LR:tof_BT:tof_BR:"\
    "v_LT:v_LR:v_BT:v_BR:En_LT:En_LR:En_BT:En_BR:rf_phase");
	
tree_LR->Branch("bdn", &bdn, "miss_R_mcpA/O:miss_R_mcpB:miss_R_mcpC:miss_R_mcpD:miss_T_mcpA:miss_T_mcpB:miss_T_mcpC:miss_T_mcpD:fid_area_hit_R_mcp:fid_area_hit_T_mcp:"\
    "a_R_ge/I:a_T_ge:a_R_ge_highE:a_T_ge_highE:"\
    "a_L_dEa:a_L_dEb:a_L_dEsum:a_L_E:"\
    "a_B_dEa:a_B_dEb:a_B_dEsum:a_B_E:"\
    "a_T_mcpA:a_T_mcpB:a_T_mcpC:a_T_mcpD:a_T_mcpE:a_T_mcpSum:"\
    "a_R_mcpA:a_R_mcpB:a_R_mcpC:a_R_mcpD:a_R_mcpE:a_R_mcpSum:"\
    "t_T_mcp:t_R_mcp:"\
    "t_B_dEa:t_B_dEb:t_B_E:"\
    "t_L_dEa:t_L_dEb:t_L_E:"\
    "t_rf:t_T_ge:t_R_ge:"\
    "s_ms_since_capt:s_capt_state:s_ms_since_eject:s_capt:deadTime_us:s_SiX4:event_good:event:run:"\
    "T_mcpX/D:T_mcpY:R_mcpX:R_mcpY:T_mcpPhysX:T_mcpPhysY:R_mcpPhysX:R_mcpPhysY:t_B_dE:t_L_dE:tof_LT:tof_LR:tof_BT:tof_BR:"\
    "v_LT:v_LR:v_BT:v_BR:En_LT:En_LR:En_BT:En_BR:rf_phase");
	
tree_BT->Branch("bdn", &bdn, "miss_R_mcpA/O:miss_R_mcpB:miss_R_mcpC:miss_R_mcpD:miss_T_mcpA:miss_T_mcpB:miss_T_mcpC:miss_T_mcpD:fid_area_hit_R_mcp:fid_area_hit_T_mcp:"\
    "a_R_ge/I:a_T_ge:a_R_ge_highE:a_T_ge_highE:"\
    "a_L_dEa:a_L_dEb:a_L_dEsum:a_L_E:"\
    "a_B_dEa:a_B_dEb:a_B_dEsum:a_B_E:"\
    "a_T_mcpA:a_T_mcpB:a_T_mcpC:a_T_mcpD:a_T_mcpE:a_T_mcpSum:"\
    "a_R_mcpA:a_R_mcpB:a_R_mcpC:a_R_mcpD:a_R_mcpE:a_R_mcpSum:"\
    "t_T_mcp:t_R_mcp:"\
    "t_B_dEa:t_B_dEb:t_B_E:"\
    "t_L_dEa:t_L_dEb:t_L_E:"\
    "t_rf:t_T_ge:t_R_ge:"\
    "s_ms_since_capt:s_capt_state:s_ms_since_eject:s_capt:deadTime_us:s_SiX4:event_good:event:run:"\
    "T_mcpX/D:T_mcpY:R_mcpX:R_mcpY:T_mcpPhysX:T_mcpPhysY:R_mcpPhysX:R_mcpPhysY:t_B_dE:t_L_dE:tof_LT:tof_LR:tof_BT:tof_BR:"\
    "v_LT:v_LR:v_BT:v_BR:En_LT:En_LR:En_BT:En_BR:rf_phase");
	
tree_BR->Branch("bdn", &bdn, "miss_R_mcpA/O:miss_R_mcpB:miss_R_mcpC:miss_R_mcpD:miss_T_mcpA:miss_T_mcpB:miss_T_mcpC:miss_T_mcpD:fid_area_hit_R_mcp:fid_area_hit_T_mcp:"\
    "a_R_ge/I:a_T_ge:a_R_ge_highE:a_T_ge_highE:"\
    "a_L_dEa:a_L_dEb:a_L_dEsum:a_L_E:"\
    "a_B_dEa:a_B_dEb:a_B_dEsum:a_B_E:"\
    "a_T_mcpA:a_T_mcpB:a_T_mcpC:a_T_mcpD:a_T_mcpE:a_T_mcpSum:"\
    "a_R_mcpA:a_R_mcpB:a_R_mcpC:a_R_mcpD:a_R_mcpE:a_R_mcpSum:"\
    "t_T_mcp:t_R_mcp:"\
    "t_B_dEa:t_B_dEb:t_B_E:"\
    "t_L_dEa:t_L_dEb:t_L_E:"\
    "t_rf:t_T_ge:t_R_ge:"\
    "s_ms_since_capt:s_capt_state:s_ms_since_eject:s_capt:deadTime_us:s_SiX4:event_good:event:run:"\
    "T_mcpX/D:T_mcpY:R_mcpX:R_mcpY:T_mcpPhysX:T_mcpPhysY:R_mcpPhysX:R_mcpPhysY:t_B_dE:t_L_dE:tof_LT:tof_LR:tof_BT:tof_BR:"\
    "v_LT:v_LR:v_BT:v_BR:En_LT:En_LR:En_BT:En_BR:rf_phase");
	
tree_bkgd_LT->Branch("bdn", &bdn, "miss_R_mcpA/O:miss_R_mcpB:miss_R_mcpC:miss_R_mcpD:miss_T_mcpA:miss_T_mcpB:miss_T_mcpC:miss_T_mcpD:fid_area_hit_R_mcp:fid_area_hit_T_mcp:"\
    "a_R_ge/I:a_T_ge:a_R_ge_highE:a_T_ge_highE:"\
    "a_L_dEa:a_L_dEb:a_L_dEsum:a_L_E:"\
    "a_B_dEa:a_B_dEb:a_B_dEsum:a_B_E:"\
    "a_T_mcpA:a_T_mcpB:a_T_mcpC:a_T_mcpD:a_T_mcpE:a_T_mcpSum:"\
    "a_R_mcpA:a_R_mcpB:a_R_mcpC:a_R_mcpD:a_R_mcpE:a_R_mcpSum:"\
    "t_T_mcp:t_R_mcp:"\
    "t_B_dEa:t_B_dEb:t_B_E:"\
    "t_L_dEa:t_L_dEb:t_L_E:"\
    "t_rf:t_T_ge:t_R_ge:"\
    "s_ms_since_capt:s_capt_state:s_ms_since_eject:s_capt:deadTime_us:s_SiX4:event_good:event:run:"\
    "T_mcpX/D:T_mcpY:R_mcpX:R_mcpY:T_mcpPhysX:T_mcpPhysY:R_mcpPhysX:R_mcpPhysY:t_B_dE:t_L_dE:tof_LT:tof_LR:tof_BT:tof_BR:"\
    "v_LT:v_LR:v_BT:v_BR:En_LT:En_LR:En_BT:En_BR:rf_phase");

tree_bkgd_LR->Branch("bdn", &bdn, "miss_R_mcpA/O:miss_R_mcpB:miss_R_mcpC:miss_R_mcpD:miss_T_mcpA:miss_T_mcpB:miss_T_mcpC:miss_T_mcpD:fid_area_hit_R_mcp:fid_area_hit_T_mcp:"\
    "a_R_ge/I:a_T_ge:a_R_ge_highE:a_T_ge_highE:"\
    "a_L_dEa:a_L_dEb:a_L_dEsum:a_L_E:"\
    "a_B_dEa:a_B_dEb:a_B_dEsum:a_B_E:"\
    "a_T_mcpA:a_T_mcpB:a_T_mcpC:a_T_mcpD:a_T_mcpE:a_T_mcpSum:"\
    "a_R_mcpA:a_R_mcpB:a_R_mcpC:a_R_mcpD:a_R_mcpE:a_R_mcpSum:"\
    "t_T_mcp:t_R_mcp:"\
    "t_B_dEa:t_B_dEb:t_B_E:"\
    "t_L_dEa:t_L_dEb:t_L_E:"\
    "t_rf:t_T_ge:t_R_ge:"\
    "s_ms_since_capt:s_capt_state:s_ms_since_eject:s_capt:deadTime_us:s_SiX4:event_good:event:run:"\
    "T_mcpX/D:T_mcpY:R_mcpX:R_mcpY:T_mcpPhysX:T_mcpPhysY:R_mcpPhysX:R_mcpPhysY:t_B_dE:t_L_dE:tof_LT:tof_LR:tof_BT:tof_BR:"\
    "v_LT:v_LR:v_BT:v_BR:En_LT:En_LR:En_BT:En_BR:rf_phase");

tree_bkgd_BT->Branch("bdn", &bdn, "miss_R_mcpA/O:miss_R_mcpB:miss_R_mcpC:miss_R_mcpD:miss_T_mcpA:miss_T_mcpB:miss_T_mcpC:miss_T_mcpD:fid_area_hit_R_mcp:fid_area_hit_T_mcp:"\
    "a_R_ge/I:a_T_ge:a_R_ge_highE:a_T_ge_highE:"\
    "a_L_dEa:a_L_dEb:a_L_dEsum:a_L_E:"\
    "a_B_dEa:a_B_dEb:a_B_dEsum:a_B_E:"\
    "a_T_mcpA:a_T_mcpB:a_T_mcpC:a_T_mcpD:a_T_mcpE:a_T_mcpSum:"\
    "a_R_mcpA:a_R_mcpB:a_R_mcpC:a_R_mcpD:a_R_mcpE:a_R_mcpSum:"\
    "t_T_mcp:t_R_mcp:"\
    "t_B_dEa:t_B_dEb:t_B_E:"\
    "t_L_dEa:t_L_dEb:t_L_E:"\
    "t_rf:t_T_ge:t_R_ge:"\
    "s_ms_since_capt:s_capt_state:s_ms_since_eject:s_capt:deadTime_us:s_SiX4:event_good:event:run:"\
    "T_mcpX/D:T_mcpY:R_mcpX:R_mcpY:T_mcpPhysX:T_mcpPhysY:R_mcpPhysX:R_mcpPhysY:t_B_dE:t_L_dE:tof_LT:tof_LR:tof_BT:tof_BR:"\
    "v_LT:v_LR:v_BT:v_BR:En_LT:En_LR:En_BT:En_BR:rf_phase");

tree_BR->Branch("bdn", &bdn, "miss_R_mcpA/O:miss_R_mcpB:miss_R_mcpC:miss_R_mcpD:miss_T_mcpA:miss_T_mcpB:miss_T_mcpC:miss_T_mcpD:fid_area_hit_R_mcp:fid_area_hit_T_mcp:"\
    "a_R_ge/I:a_T_ge:a_R_ge_highE:a_T_ge_highE:"\
    "a_L_dEa:a_L_dEb:a_L_dEsum:a_L_E:"\
    "a_B_dEa:a_B_dEb:a_B_dEsum:a_B_E:"\
    "a_T_mcpA:a_T_mcpB:a_T_mcpC:a_T_mcpD:a_T_mcpE:a_T_mcpSum:"\
    "a_R_mcpA:a_R_mcpB:a_R_mcpC:a_R_mcpD:a_R_mcpE:a_R_mcpSum:"\
    "t_T_mcp:t_R_mcp:"\
    "t_B_dEa:t_B_dEb:t_B_E:"\
    "t_L_dEa:t_L_dEb:t_L_E:"\
    "t_rf:t_T_ge:t_R_ge:"\
    "s_ms_since_capt:s_capt_state:s_ms_since_eject:s_capt:deadTime_us:s_SiX4:event_good:event:run:"\
    "T_mcpX/D:T_mcpY:R_mcpX:R_mcpY:T_mcpPhysX:T_mcpPhysY:R_mcpPhysX:R_mcpPhysY:t_B_dE:t_L_dE:tof_LT:tof_LR:tof_BT:tof_BR:"\
    "v_LT:v_LR:v_BT:v_BR:En_LT:En_LR:En_BT:En_BR:rf_phase");

zero_time_Tree->Branch("bdn", &bdn, "miss_R_mcpA/O:miss_R_mcpB:miss_R_mcpC:miss_R_mcpD:miss_T_mcpA:miss_T_mcpB:miss_T_mcpC:miss_T_mcpD:fid_area_hit_R_mcp:fid_area_hit_T_mcp:"\
    "a_R_ge/I:a_T_ge:a_R_ge_highE:a_T_ge_highE:"\
    "a_L_dEa:a_L_dEb:a_L_dEsum:a_L_E:"\
    "a_B_dEa:a_B_dEb:a_B_dEsum:a_B_E:"\
    "a_T_mcpA:a_T_mcpB:a_T_mcpC:a_T_mcpD:a_T_mcpE:a_T_mcpSum:"\
    "a_R_mcpA:a_R_mcpB:a_R_mcpC:a_R_mcpD:a_R_mcpE:a_R_mcpSum:"\
    "t_T_mcp:t_R_mcp:"\
    "t_B_dEa:t_B_dEb:t_B_E:"\
    "t_L_dEa:t_L_dEb:t_L_E:"\
    "t_rf:t_T_ge:t_R_ge:"\
    "s_ms_since_capt:s_capt_state:s_ms_since_eject:s_capt:deadTime_us:s_SiX4:event_good:event:run:"\
    "T_mcpX/D:T_mcpY:R_mcpX:R_mcpY:T_mcpPhysX:T_mcpPhysY:R_mcpPhysX:R_mcpPhysY:t_B_dE:t_L_dE:tof_LT:tof_LR:tof_BT:tof_BR:"\
    "v_LT:v_LR:v_BT:v_BR:En_LT:En_LR:En_BT:En_BR:rf_phase");

beta_recoil_tree->Branch("bdn", &bdn, "miss_R_mcpA/O:miss_R_mcpB:miss_R_mcpC:miss_R_mcpD:miss_T_mcpA:miss_T_mcpB:miss_T_mcpC:miss_T_mcpD:fid_area_hit_R_mcp:fid_area_hit_T_mcp:"\
    "a_R_ge/I:a_T_ge:a_R_ge_highE:a_T_ge_highE:"\
    "a_L_dEa:a_L_dEb:a_L_dEsum:a_L_E:"\
    "a_B_dEa:a_B_dEb:a_B_dEsum:a_B_E:"\
    "a_T_mcpA:a_T_mcpB:a_T_mcpC:a_T_mcpD:a_T_mcpE:a_T_mcpSum:"\
    "a_R_mcpA:a_R_mcpB:a_R_mcpC:a_R_mcpD:a_R_mcpE:a_R_mcpSum:"\
    "t_T_mcp:t_R_mcp:"\
    "t_B_dEa:t_B_dEb:t_B_E:"\
    "t_L_dEa:t_L_dEb:t_L_E:"\
    "t_rf:t_T_ge:t_R_ge:"\
    "s_ms_since_capt:s_capt_state:s_ms_since_eject:s_capt:deadTime_us:s_SiX4:event_good:event:run:"\
    "T_mcpX/D:T_mcpY:R_mcpX:R_mcpY:T_mcpPhysX:T_mcpPhysY:R_mcpPhysX:R_mcpPhysY:t_B_dE:t_L_dE:tof_LT:tof_LR:tof_BT:tof_BR:"\
    "v_LT:v_LR:v_BT:v_BR:En_LT:En_LR:En_BT:En_BR:rf_phase");

beta_gamma_tree->Branch("bdn", &bdn, "miss_R_mcpA/O:miss_R_mcpB:miss_R_mcpC:miss_R_mcpD:miss_T_mcpA:miss_T_mcpB:miss_T_mcpC:miss_T_mcpD:fid_area_hit_R_mcp:fid_area_hit_T_mcp:"\
    "a_R_ge/I:a_T_ge:a_R_ge_highE:a_T_ge_highE:"\
    "a_L_dEa:a_L_dEb:a_L_dEsum:a_L_E:"\
    "a_B_dEa:a_B_dEb:a_B_dEsum:a_B_E:"\
    "a_T_mcpA:a_T_mcpB:a_T_mcpC:a_T_mcpD:a_T_mcpE:a_T_mcpSum:"\
    "a_R_mcpA:a_R_mcpB:a_R_mcpC:a_R_mcpD:a_R_mcpE:a_R_mcpSum:"\
    "t_T_mcp:t_R_mcp:"\
    "t_B_dEa:t_B_dEb:t_B_E:"\
    "t_L_dEa:t_L_dEb:t_L_E:"\
    "t_rf:t_T_ge:t_R_ge:"\
    "s_ms_since_capt:s_capt_state:s_ms_since_eject:s_capt:deadTime_us:s_SiX4:event_good:event:run:"\
    "T_mcpX/D:T_mcpY:R_mcpX:R_mcpY:T_mcpPhysX:T_mcpPhysY:R_mcpPhysX:R_mcpPhysY:t_B_dE:t_L_dE:tof_LT:tof_LR:tof_BT:tof_BR:"\
    "v_LT:v_LR:v_BT:v_BR:En_LT:En_LR:En_BT:En_BR:rf_phase");

}
