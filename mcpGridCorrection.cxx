#include "Rtypes.h"
#include "CSVtoStruct.h"
#include "TMath.h"

Double_t tofToMCPGrid (BDNCase_t stBDNCase, char whichMCP, Double_t t2) {
	
	using namespace TMath;
	Double_t a, z1, z2;
	if (whichMCP == 'R') {
		a		= stBDNCase.dRightGridAcceleration;
		z1		= stBDNCase.dRightGridDistance;
		z2		= stBDNCase.dRightMCPDistance;
	}
	if (whichMCP == 'T') {
		a		= stBDNCase.dTopGridAcceleration;
		z1		= stBDNCase.dTopGridDistance;
		z2		= stBDNCase.dTopMCPDistance;
	}
	
	t2 = t2/1000.0; // convert ns to us
	
	Double_t A		= - 9.0*Power(a,2.0)*(3.0*z1-2.0*z2)*t2 - Power(a,3.0)*Power(t2,3.0);
	Double_t B		= 27.0*Power(a,3.0)*(8.0*Power(z2,3.0) - a*(27.0*Power(z1,2.0)-36.0*z1*z2+8.0*Power(z2,2.0))*Power(t2,2.0) + 2.0*Power(a,2.0)*(z2-z1)*Power(t2,4.0));
	Double_t theta	= ATan(Sqrt(B)/A)+Pi();
	Double_t t1		= (1.0/6.0)*(4.0*t2-((6*z2+a*t2*t2)/Power(A*A+B,1.0/6.0)+Power(A*A+B,1.0/6.0)/a)*(Cos(theta/3)-Sqrt(3)*Sin(theta/3)));
	
// Diagnostics
//	printf("%s MCP\n", &whichMCP);
//	printf("a = %f\n", a);
//	printf("z1 = %f\n", z1);
//	printf("z2 = %f\n", z2);
//	printf("A = %f\n", A);
//	printf("B = %f\n", B);
//	printf("theta = %f\n", theta);
//	printf("t1 = %f\n", t1);
	
	return 1000.0 * t1; // return in ns
}
