#include "calibLineCam.h"


int main()
{
	CcalibrateLineCamera calib;
	calib.Calibrate("C:\\Users\\YYKJ05\\Desktop\\LittleCamera", 21, "D:\\hg-450.cpd", "D:\\hg-450.ps", "C:\\Users\\YYKJ05\\Desktop\\CalibResult.csv");

	return 0;
}

