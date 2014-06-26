// File to include experiment-specific data
// Ideally this data will be the same for all BDN runs in a campaign

#ifndef _bdn_h
#define _bdn_h

enum BRCombos_t {LT, LR, BT, BR};
static const Int_t nDetectors_dE  = 2;
static const Int_t nDetectors_MCP = 2;
static const Int_t nDetectors_Ge  = 2;

//static const Double_t rf_freq	= 310000.0;// Hz

//static const Double_t R_mcpX_width		= 46.0;// dist from trap center to MCP grid in mm
//static const Double_t R_mcpY_width		= 46.0;// dist from trap center to MCP grid in mm
//static const Double_t T_mcpX_width		= 46.0;// dist from trap center to MCP grid in mm
//static const Double_t T_mcpY_width		= 46.0;// dist from trap center to MCP grid in mm
//static const Double_t mcp_dist			= 52.5;	// dist from trap center to MCP face in mm
//static const Double_t mcp_grid_dist		= 48.0;// dist from trap center to MCP grid in mm
//static const Double_t mcp_r_avg			= 53.298471; // geometric average of straight path lengths (r) to MCP
//static const Double_t bias_T_mcp		= 2.5;	// kV
//static const Double_t bias_R_mcp		= 2.5;	// kV

// selfTime is a window: {low,high}
//static const Double_t B_dEa_selfTime[2] = {-24.0,-21.0};
//static const Double_t B_dEb_selfTime[2] = {-25.0,-21.0};
//static const Double_t B_E_selfTime[2]	 = {-22.0,-21.0};
//static const Double_t L_dEa_selfTime[2] = {-24.0,-21.0};
//static const Double_t L_dEb_selfTime[2] = {-25.0,-21.0};
//static const Double_t L_E_selfTime[2] 	 = {-22.0,-21.0};
//static const Double_t T_mcp_selfTime[2] = {-24.0,-22.0};
//static const Double_t R_mcp_selfTime[2] = {-24.0,-22.0};
//static const Double_t T_ge_selfTime[2];
//static const Double_t R_ge_selfTime[2];

// zeroTime is a value and sigma: {t0,sigma}
static const Double_t LT_zeroTime[2] = {-72.244,	1.2}; // left dE, top mcp
static const Double_t LR_zeroTime[2] = {-7.255,	1.0}; // left dE, right mcp
static const Double_t BT_zeroTime[2] = {-63.855,	1.0}; // bottom dE, top mcp
static const Double_t BR_zeroTime[2] = {1.930,	1.1}; // bottom dE, right mcp

static const Double_t LT_zeroTime_E[2] = {-81.99,	0.91}; // left E, top mcp
static const Double_t LR_zeroTime_E[2] = {-19.777,	1.433}; // left E, right mcp
static const Double_t BT_zeroTime_E[2] = {-81.018,	1.210}; // bottom E, top mcp
static const Double_t BR_zeroTime_E[2] = {-15.808,	1.363}; // bottom E, right mcp

// HPGe calibration
static const Double_t R_ge_coeff[3]		= {-10.64667, 0.8468166, 0.0000002009597};//0.8471;
static const Double_t T_ge_coeff[3]		= {-7.943132, 0.8370388, 0.0000010671100};//-10.674;
static const Double_t R_ge_highE_coeff[3]= {-10.64667, 2.5000000, 0.0000000000000};//0.8471;
static const Double_t T_ge_highE_coeff[3]= {-7.943132, 2.5000000, 0.0000000000000};//-10.674;

// MCP Map parameters
// 4-post map
static const Double_t R_mcp_theta	=  0.007302;
static const Double_t R_mcp_x0		= -0.004353;
static const Double_t R_mcp_y0		= -0.010621;
static const Double_t R_mcp_a		= 25.728395;
static const Double_t T_mcp_theta	=  0.002577;
static const Double_t T_mcp_x0		=  0.006452;
static const Double_t T_mcp_y0		=  0.024603;
static const Double_t T_mcp_a		= 25.645184;
// 3-post map
static const Double_t R_mcp_3post_theta	= -0.012313;
static const Double_t R_mcp_3post_x0	=  0.026255;
static const Double_t R_mcp_3post_y0	= -0.025641;
static const Double_t R_mcp_3post_a0	= 24.232701;
static const Double_t R_mcp_3post_b0	= 25.897749;
static const Double_t R_mcp_3post_a2	=  4.1616254;
static const Double_t R_mcp_3post_b2	=  3.9205337;
static const Double_t T_mcp_3post_theta	=  0.017888;
static const Double_t T_mcp_3post_x0	=  0.062248;
static const Double_t T_mcp_3post_y0	= -0.022670;
static const Double_t T_mcp_3post_a0	= 29.076854;
static const Double_t T_mcp_3post_b0	= 24.637936;
static const Double_t T_mcp_3post_a2	=  2.3815485;
static const Double_t T_mcp_3post_b2	=  6.8938076;

// ADC pedestals
static const Double_t ped_R_mcpA = 19.17;//18.6412;
static const Double_t ped_R_mcpB = 16.01;//15.8493;
static const Double_t ped_R_mcpC = 18.23;//17.4128;
static const Double_t ped_R_mcpD = 19.56;//19.1173;
static const Double_t ped_R_mcpE = 22.93;
static const Double_t ped_T_mcpA = 19.89;//20.0549;//15.30;
static const Double_t ped_T_mcpB = 15.18;//14.8766;//15.70;
static const Double_t ped_T_mcpC = 20.15;//20.2492;//16.66;
static const Double_t ped_T_mcpD = 20.11;//19.8245;//16.47;
static const Double_t ped_T_mcpE = 19.97;

// Thresholds for cuts
static const Double_t t_trigger_lo	= -100.0;
static const Double_t t_trigger_hi	=    0.0;
static const Double_t t_E_lo		= -100.0;
static const Double_t t_E_hi		=    0.0;
static const Double_t t_dE_lo		= -100.0;
static const Double_t t_dE_hi		=    0.0;
static const Double_t t_ge_lo		= -100.0;
static const Double_t t_ge_hi		=    0.0;
static const Double_t t_mcp_lo		= -100.0;
static const Double_t a_E_lo		=   50.0;
static const Double_t a_dE_lo		=  100.0;
static const Double_t a_ge_lo		=  400.0;
static const Double_t a_mcp_lo		=  200.0;// sum of posts cut

static const Double_t a_R_mcpA_lo	=  5.0;
static const Double_t a_R_mcpB_lo	=  5.0;
static const Double_t a_R_mcpC_lo	= 10.0;
static const Double_t a_R_mcpD_lo	=  5.0;
static const Double_t a_T_mcpA_lo	= 34.0;
static const Double_t a_T_mcpB_lo	= 30.0;
static const Double_t a_T_mcpC_lo	= 34.0;
static const Double_t a_T_mcpD_lo	= 36.0;
static const Double_t a_missing_mcp_post = -1000.0;

// dE-MCP TOF ranges for cuts
static const Double_t tof_zero_lo	=-5.0;
static const Double_t tof_zero_hi	= 5.0;
static const Double_t tof_lowTOF_lo	= 10.0;
static const Double_t tof_lowTOF_hi	= 230.0;
static const Double_t tof_fast_lo	= 230.0;
static const Double_t tof_fast_hi	= 1550.0;
static const Double_t tof_slow_lo	= 1550.0;
static const Double_t tof_slow_hi	= 10000.0;
static const Double_t tof_oops_lo	= 15000.0;
static const Double_t tof_oops_hi	= 20000.0;

// MCP Fiducial Areas -- physical values in mm
static const Double_t fid_area_R_mcpPhysX_lo	= -23.00;
static const Double_t fid_area_R_mcpPhysX_hi	=  23.00;
static const Double_t fid_area_R_mcpPhysY_lo	= -23.00;
static const Double_t fid_area_R_mcpPhysY_hi	=  23.00;
static const Double_t fid_area_T_mcpPhysX_lo	= -23.00;
static const Double_t fid_area_T_mcpPhysX_hi	=  23.00;
static const Double_t fid_area_T_mcpPhysY_lo	= -23.00;
static const Double_t fid_area_T_mcpPhysY_hi	=  23.00;
// MCP Fiducial Areas -- raw coord values
static const Double_t fid_area_R_mcpX_lo		= -0.862;
static const Double_t fid_area_R_mcpX_hi		=  0.860;
static const Double_t fid_area_R_mcpY_lo		= -0.866;
static const Double_t fid_area_R_mcpY_hi		=  0.835;
static const Double_t fid_area_T_mcpX_lo		= -0.876;
static const Double_t fid_area_T_mcpX_hi		=  0.861;
static const Double_t fid_area_T_mcpY_lo		= -0.845;
static const Double_t fid_area_T_mcpY_hi		=  0.865;

#endif
