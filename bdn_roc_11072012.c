// 2012-02-20 Adapted from MGS's Li_ROC1.c by SC
// 2012-05-07 Adding TDC1
// 2012-05-16 Adding Delay Generator to test accuracy of TDC
// 2012-05-18 Making the DG cycle through active tables to send each of 
// ... 8 delays in turn; for testing TDC linearity with only one 'acquire'
// 2012-09-07 Adding lower threshold to ADC
// 2012-09-18 Remove clear command from TrigSyncScaler in surf_user_acquire
// 2012-09-19 Changed ADC1 from 7 to 9
// 2012-09-19 Uncommented the NAF(ADC1,0,9) and NAF(ADC1,0,23) commands in surf_user_acquire
//            ... why were these commented out??
// 2012-09-20 Added ADC2 to hold low-gain/high-full-scale Ge channels and MCP_E channels
//            ... requires new sort code
// 2012-09-27 Commented out all ADC2 code.  ADC1 lowthr back to 30 from 120
// ... (had been using 120 when diagnosing ADC1 fishiness last weekend)
// 2012-09-27 Rearranged scalers as described on Sept 20 in online log.
// ... In surf_user_acquire: clear TrigSyncScaler; do not clear EjectScaler 
// ... and CaptScaler
// ... Read TrigSyncScaler only in surf_user_trigsync
// 2012-10-22 Added day to trigsync timestamp
// 2012-10-23 Changed ADC1 and ADC channels from 7 and 9 to 6 and 8
// 2012-11-07 Scaler channels have been messed up since c. 2012-09-20 (see 2012-11 log book)
//	Adding several lines of readout to Eject scaler and TrigSync scaler and changing comments in scaler readouts to reflect reality
//	Reducing ADC lower threshold from 30 to 0 to allow us to see low level noise

#include <stdio.h>
#include <unistd.h>
#include <time.h>   //to get the system time
#include "surf_cc32.h"

// Names for CAMAC channels
#define ADC1		6		// CAMAC channel of ADC 1 (Phillips 7164)
//#define ADC2		8		// CAMAC channel of ADC 2 (Phillips 7164)
#define TDC1		12	// CAMAC channel of TDC 1 (LeCroy 4208)
#define TDC2    14  // CAMAC channel of TDC 2 (LeCroy 4208)
#define CaptScaler	17	// Scaler that stops at capt pulse (LeCroy 2551)
#define EjectScaler	19	// Scaler that stops at eject pulse (LeCroy 2551)
#define TrigSyncScaler	21	// Scaler that stops at trigsync function (CAEN C257)
//#define DG			15	// CAMAC channel of delay generator (LeCroy DG11A)

time_t now;
struct tm *timenow;
int i,j,n,x;

	// Set ADC pedestals & thresholds
	int ped1=0;
	int lowthr1=0;
	int hithr1=4095;
//	int ped2=0;
//	int lowthr2=120;
//	int hithr2=4095;

void *surf_user_acquire(void* start) // Called when Scarlet is told to 'acquire'
{
	// Clear and reset all modules
	NAF(ADC1,0,9);		// clear and reset the module
	NAF(ADC1,0,23);		// usleep(10); // wait for 10us// reset module control register
//	NAF(ADC2,0,9);		// clear and reset the module
//	NAF(ADC2,0,23);		// usleep(10); // wait for 10us// reset module control register

/* HERE IS WHAT I'VE BEEN USING BEFORE DOING MANUAL INIT 
*/
	//No! NAF(CaptScaler,0,9);		// clear module
	//No! NAF(EjectScaler,0,9);		// clear module
	// Want the above to keep their states when we start a file
	NAF(TrigSyncScaler,0,9);		// clear module
	// This scaler now counts hits on all detectors

	NAF(TDC1,0,9);		// clear and reset the module
	NAF(TDC2,0,9);		// clear and reset the module

	NAF(ADC1,3,11);    //Reset hit regester,LAM and data memory.
	NAF(ADC1,0,26);    //Enable LAM
	NAF(ADC1,0,19)=7;  // 7 = 00...00111  Enable UT,LT,ped

//	NAF(ADC2,3,11);    //Reset hit regester,LAM and data memory.
//	NAF(ADC2,0,26);    //Enable LAM
//	NAF(ADC2,0,19)=7;  // 7 = 00...00111  Enable UT,LT,ped
	
	for (i=0; i<16; i++) {
		NAF(ADC1,0,17) = 0;
		NAF(ADC1,i,20) = ped1;
		NAF(ADC1,1,17) = 0;
		NAF(ADC1,i,20) = lowthr1;
		NAF(ADC1,2,17) = 0;
		NAF(ADC1,i,20) = hithr1;

//		NAF(ADC2,0,17) = 0;
//		NAF(ADC2,i,20) = ped2;
//		NAF(ADC2,1,17) = 0;
//		NAF(ADC2,i,20) = lowthr2;
//		NAF(ADC2,2,17) = 0;
//		NAF(ADC2,i,20) = hithr2;
	}
	
    int *p = (int*)start;
	
    time(&now);
    timenow = localtime(&now);
    
    *p++ = 0xABCD;
    *p++ = timenow->tm_mday;
    *p++ = timenow->tm_hour;
    *p++ = timenow->tm_min;
    *p++ = timenow->tm_sec;
    *p++ = 0xABCD;
	
    return p;
}

void *surf_user_event(void *start)
{
  int *p = (int*)start;
	
/***** ADC READ-OUT *****/

	int j;
  int maxtries = 20;
	for (j=1; (x=NAF(ADC1,0,8)) > 0 && j < maxtries; ++j){
		// Loops until ADC1 is ready
		// Each F(8)A(0) call takes 2 mic-sec
	}
	
	x=NAF(ADC1,0,8); // Extra call to make sure all ADC's ready
	
	*p++ = 0xadc1adc1; // Marker for ADC1 data
	*p++ = NAF(ADC1,1,6);	//Read hit register
	while( (x=NAF(ADC1,0,4)) < 0 ) *p++ = x;
	// Read sparse data: Valid data < 0
	NAF(ADC1,3,11);	// Reset Hit Registr,LAM and data memory
	
//	*p++ = 0xadc2adc2; // Marker for ADC2 data
//	*p++ = NAF(ADC2,1,6);	//Read hit register
//	while( (x=NAF(ADC2,0,4)) < 0 ) *p++ = x;
//	// Read sparse data: Valid data < 0
//	NAF(ADC2,3,11);	// Reset Hit Registr,LAM and data memory

// End ADC Read-out

/***** TDC READ-OUT *****/

	//if (NAF(TDC1,0,8) < 0){
	*p++ = 0x2dc12dc1; // Marker for TDC1 data
	for(j=0; j<8; j++){
		x = NAF(TDC1,j,0); // Read all TDC channels
		if (x & 0x80000000){
			*p++ = j+1;
			*p++ = x;
		}
	}
	//usleep(100); 	// wait for this number of microseconds
	NAF(TDC1,0,9);	// Clear TDC and reset LAM
	//}
	//else *p++ = 0xdead2dc1;

	*p++ = 0x2dc22dc2; // Marker for TDC2 data
	for(j=0; j<8; j++){
		x = NAF(TDC2,j,0);
		if (x & 0x80000000){
			*p++ = j+1;
			*p++ = x;
		}
	}
	NAF(TDC2,0,9);
// End TDC Read-out

// Read out scalers
	*p++ = 0x100cca1e;  // Capture scaler
	*p++ = NAF(CaptScaler,0,0); // Read Channel 0 (time since capture 1 kHz clock)
	*p++ = NAF(CaptScaler,1,0); // Read Channel 1 (trap state: 0 = trap mode, 1 = waiting)
	// scaler is hardware-cleared	by BPT capture pulse
	
	*p++ = 0x100eca1e;  // Eject scaler
	*p++ = NAF(EjectScaler,0,0); // Read Channel 0 (time since ejection 1 kHz clock)
	*p++ = NAF(EjectScaler,1,0); // Read Channel 1 (# of captures since last eject)
	*p++ = NAF(EjectScaler,2,0); // Read Channel 2 (time since ejection 2 Hz clock) <-- This is what we've really been getting
	*p++ = NAF(EjectScaler,3,0); // Read Channel 3 (# of SiX4 counts) <-- Added 2012-11-07
	// scaler is hardware-cleared	by BPT ejection pulse

	return p;
}

void *surf_user_trigsync(void *start)
{
    int *p = (int*)start;

	*p++ = 0x1002ca1e;

// 2012-11-07 Comments are rearranged to reflect true arrangement of these signals
// See Nov log book for discussion
	*p++ = NAF(TrigSyncScaler,0,0); // Top MCP hits
	*p++ = NAF(TrigSyncScaler,1,0); // Right MCP hits
	*p++ = NAF(TrigSyncScaler,2,0); // Bottom dEa hits
	*p++ = NAF(TrigSyncScaler,3,0); // Bottom E hits
	*p++ = NAF(TrigSyncScaler,4,0); // Left dEa hits
	*p++ = NAF(TrigSyncScaler,5,0); // Left dEb hits
	*p++ = NAF(TrigSyncScaler,6,0); // Left E hits
	*p++ = NAF(TrigSyncScaler,7,0); // Top Ge hits
	*p++ = NAF(TrigSyncScaler,8,0); // Right Ge hits
	*p++ = NAF(TrigSyncScaler,9,0); // Bottom dEb hits
	*p++ = NAF(TrigSyncScaler,10,0); // SiX4 counts
	NAF(TrigSyncScaler,0,9);	//Clear & Reset Scalar

	time(&now);
	timenow = localtime(&now);
    
	*p++ = 0xABCD;
	*p++ = timenow->tm_mday;
	*p++ = timenow->tm_hour;
	*p++ = timenow->tm_min;
	*p++ = timenow->tm_sec;
	*p++ = 0xABCD;
    return p;
}

void *surf_user_stop(void *start)
{

	int *p = (int*)start;
    
	time(&now);
	timenow = localtime(&now);
    
    *p++ = 0xABCD;
    *p++ = timenow->tm_mday;
    *p++ = timenow->tm_hour;
    *p++ = timenow->tm_min;
    *p++ = timenow->tm_sec;
    *p++ = 0xABCD;
    return p; 
}
