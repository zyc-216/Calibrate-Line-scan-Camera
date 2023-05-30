#include "calibLineCam.h"
#include "HalconCpp.h"
#include "HDevThread.h"
//#include <ctime>
#include <random>
#include <cmath>
//#include <pcl/sample_consensus/ransac.h>
//#include <pcl/sample_consensus/sac_model_plane.h>
//#include <pcl/visualization/pcl_visualizer.h>
//typedef pcl::PointXYZ PointT;
//typedef pcl::PointCloud<PointT> PointCloudT;

CcalibrateLineCamera::CcalibrateLineCamera() {
	m_vecLaserPlaneCoeffs.clear();
	m_vecLaserQuadricSurfaceCoeffs.clear();
}


bool CcalibrateLineCamera::fnMatchCornersAndWorldPoints()
{
	bool bRet = false;

	HalconCpp::HObject Img;
	HalconCpp::HTuple RowPoint, ColPoint, Index, Pose;
	//创建标定板文件，设置相机参数
	HalconCpp::HTuple markRow, markCol;
	markRow.Clear();
	markCol.Clear();
	int temp1 = 0;
	int temp2 = 0;
	for (auto i = begin(m_calibObject.m_iMarkRow); i != end(m_calibObject.m_iMarkRow); ++i) {
		markRow[temp1] = *i;
		temp1++;
	}
	for (auto j = begin(m_calibObject.m_iMarkCol); j != end(m_calibObject.m_iMarkCol); ++j) {
		markCol[temp2] = *j;
		temp2++;
	}

	HalconCpp::CreateCaltab(m_calibObject.m_iCalTabRow, m_calibObject.m_iCalTabCol, m_calibObject.m_cornerDiameter, markRow, markCol
		, m_calibObject.m_strPolarity.c_str(),m_strCalTabCpdPath.c_str(), m_strCalTabPsPath.c_str());
	HalconCpp::HTuple CalibDataID;
	HalconCpp::CreateCalibData("calibration_object", 1, 1, &CalibDataID);
	HalconCpp::HTuple startCamPar;
	startCamPar.Clear();
	int temp3 = 0;
	for (auto k = begin(m_calibObject.m_startCamPar); k != end(m_calibObject.m_startCamPar); ++k) {
		startCamPar[temp3] = *k;
		temp3++;
	}
	HalconCpp::SetCalibDataCamParam(CalibDataID, 0, HalconCpp::HTuple(), startCamPar);
	HalconCpp::SetCalibDataCalibObject(CalibDataID, 0, m_strCalTabCpdPath.c_str());

	vector<string> imgNames;
	cv::glob(m_strImgDir, imgNames);
	int imageNums = imgNames.size();

	for (int i = 0; i < m_iImgNum; ++i) {
		string imgName = m_strImgDir + "\\" + to_string(i+1) + ".bmp";
		//HalconCpp::ReadImage(&Img, imgNames[i].c_str());
		HalconCpp::ReadImage(&Img, imgName.c_str());
		HalconCpp::HTuple imgWidth = m_calibObject.m_iImgWidth;
		HalconCpp::HTuple imgHeight = m_calibObject.m_iImgHeight;
		HalconCpp::GetImageSize(Img, &imgWidth, &imgHeight);
		m_calibObject.m_iImgWidth = imgWidth.I();
		m_calibObject.m_iImgHeight = imgHeight.I();
		HalconCpp::FindCalibObject(Img, CalibDataID, 0, 0, i + 1, HalconCpp::HTuple(), HalconCpp::HTuple());
		HalconCpp::GetCalibDataObservPoints(CalibDataID, 0, 0, i + 1, &ColPoint, &RowPoint, &Index, &Pose);

		for (int j = 0; j < RowPoint.Length(); ++j) {
			m_cornersPerImg.push_back(cv::Point2f(RowPoint[j], ColPoint[j]));
			m_worldPointsPerImg.push_back(cv::Point3f(m_calibObject.m_XofWorldPt[Index[j]], m_calibObject.m_YofWorldPt[Index[j]], 0.0));
		}
		m_cornersAllImg.push_back(m_cornersPerImg);
		m_cornersPerImg.clear();
		m_worldPointsAllImg.push_back(m_worldPointsPerImg);
		m_worldPointsPerImg.clear();
	}
	bRet = true;
	return bRet;
}


bool CcalibrateLineCamera::CalibrateUsingCV(const char* outputCamParamCsv) 
{
	bool bRet = false;

	this->CcalibrateLineCamera::fnMatchCornersAndWorldPoints();
	cv::Size m_imageSize = cv::Size((int)m_calibObject.m_iImgWidth, (int)m_calibObject.m_iImgHeight);

	cv::calibrateCamera(m_worldPointsAllImg, m_cornersAllImg, m_imageSize, m_cameraInstricMatrix, m_distCoeffs, m_rVecs, m_tVecs,
		cv::CALIB_TILTED_MODEL | cv::CALIB_RATIONAL_MODEL| cv::CALIB_THIN_PRISM_MODEL);
	/*cout << m_cameraInstricMatrix << cameraInstricMatrix << endl;
	cout << m_distCoeffs << distCoeffs << endl;*/
	FILE* fp = nullptr;
	fopen_s(&fp, outputCamParamCsv, "w+");
	if (nullptr != fp) {
		int cnt = 0;
		fprintf_s(fp, "%s\n", "InstricMatrix:");
		for (int i = 0; i < m_cameraInstricMatrix.cols; ++i) {
			for (int j = 0; j < m_cameraInstricMatrix.rows; ++j) {
				fprintf_s(fp,"%.10f,", m_cameraInstricMatrix.at<double>(i, j));
				cnt++;
				if (cnt % 3 == 0) {
					fprintf_s(fp,"\n");
					cnt = 0;
				}
			}
		}
		fprintf_s(fp, "%s\n", "distCoeffs:");
		for (int i = 0; i < m_distCoeffs.cols; ++i) {
			fprintf_s(fp, "%.18f,", m_distCoeffs.at<double>(0,i));
		}
		fclose(fp);
	}
	PrintCamParam();

	bRet = true;
	return bRet;
}

bool CcalibrateLineCamera::planeFitting(const vector<cv::Point3d>& input)
{
	bool bRet = false;

	cv::Mat dst = cv::Mat(3, 3, CV_64F, cv::Scalar(0));//初始化系数矩阵A
	cv::Mat out = cv::Mat(3, 1, CV_64F, cv::Scalar(0));//初始化矩阵b
	for (int i = 0; i < input.size(); ++i)
	{
		//计算3*3的系数矩阵
		dst.at<double>(0, 0) = dst.at<double>(0, 0) + pow((double)input[i].x, 2);
		dst.at<double>(0, 1) = dst.at<double>(0, 1) + (double)input[i].x * (double)input[i].y;
		dst.at<double>(0, 2) = dst.at<double>(0, 2) + (double)input[i].x;
		dst.at<double>(1, 0) = dst.at<double>(1, 0) + (double)input[i].x * (double)input[i].y;
		dst.at<double>(1, 1) = dst.at<double>(1, 1) + pow((double)input[i].y, 2);
		dst.at<double>(1, 2) = dst.at<double>(1, 2) + (double)input[i].y;
		dst.at<double>(2, 0) = dst.at<double>(2, 0) + (double)input[i].x;
		dst.at<double>(2, 1) = dst.at<double>(2, 1) + (double)input[i].y;
		dst.at<double>(2, 2) = (double)input.size();
		//计算3*1的结果矩阵
		out.at<double>(0, 0) = out.at<double>(0, 0) + (double)input[i].x * (double)input[i].z;
		out.at<double>(1, 0) = out.at<double>(1, 0) + (double)input[i].y * (double)input[i].z;
		out.at<double>(2, 0) = out.at<double>(2, 0) + (double)input[i].z;
	}
	//判断矩阵是否奇异
	//double determ = determinant(dst);
	//if (abs(determ) < 0.001) {
	//	cout << "矩阵奇异" << endl;
	//	return false;
	//}
	//Mat inv;
	//invert(dst, inv);//求矩阵的逆
	//Mat output = inv*out;//计算输出
	//coeffient.clear();//把结果输出
	//coeffient.push_back(output.at<float>(0, 0));
	//coeffient.push_back(output.at<float>(1, 0));
	//coeffient.push_back(output.at<float>(2, 0));

	// 修改此处代码为SVD分解
	cv::Mat result = cv::Mat(3, 1, CV_64F, cv::Scalar(0));
	cout << dst << out << endl;
	solve(dst, out, result, cv::DECOMP_SVD);
	cout << "系数为：" << result << endl;
	m_vecLaserPlaneCoeffs.push_back(result.at<double>(0, 0));
	m_vecLaserPlaneCoeffs.push_back(result.at<double>(1, 0));
	m_vecLaserPlaneCoeffs.push_back(result.at<double>(2, 0));

	//TestAcurracyUsingTrueData("C:\\Users\\YYKJ05\\Desktop\\trueData.csv", coeffient);
	bRet = true;
	return bRet;
}

bool CcalibrateLineCamera::GetProjectionError(const char* fileName)
{
	bool bRet = false;

	vector<cv::Point2f> ProjectCorners;
	FILE* fp = nullptr;
	fopen_s(&fp, fileName, "w+");
	if (nullptr != fp) {
		for (int i = 0; i < m_worldPointsAllImg.size(); i++){
			projectPoints(m_worldPointsAllImg[i], m_rVecs[i], m_tVecs[i], m_cameraInstricMatrix, m_distCoeffs, ProjectCorners, cv::noArray());
			cv::Mat UsubRealU = cv::Mat(1, ProjectCorners.size(), CV_64FC1);  //改动
			cv::Mat VsubRealV = cv::Mat(1, ProjectCorners.size(), CV_64FC1);	//改动
			for (int j = 0; j < ProjectCorners.size(); ++j) {
				UsubRealU.at<double>(0, j) = m_cornersAllImg[i][j].x - ProjectCorners[j].x;
				VsubRealV.at<double>(0, j) = m_cornersAllImg[i][j].y - ProjectCorners[j].y;
			}
			cv::Mat Final; //每个点的误差组成的mat
			cv::sqrt(UsubRealU.mul(UsubRealU) + VsubRealV.mul(VsubRealV), Final);
			cv::Scalar mean;  //均值
			cv::Scalar stdDev;  //标准差
			cv::meanStdDev(Final, mean, stdDev);
			fprintf_s(fp, "%s%d\n%s,%.5f,%s,%.5f\n", "img", i + 1,"mean:",mean[0],"stdDev:",stdDev[0]);
			for (int j = 0; j < ProjectCorners.size(); ++j) {
				fprintf_s(fp, "%.5f,%.5f,%.5f,%.5f,%.5f\n", Final.at<double>(0, j), m_cornersAllImg[i][j].x, m_cornersAllImg[i][j].y   //改动
							, ProjectCorners[j].x, ProjectCorners[j].y);
			}
			ProjectCorners.clear();
		}
		fclose(fp);
	}

	bRet = true;
	return bRet;
}

bool CcalibrateLineCamera::GetUndistortedCorners(const char* outputfileName)
{
	bool bRet = false;
	//CTimeDiffer	time;
	//time.StartTiming();
	//for (int i = 0; i < 1000; i++);
	FILE* fp = nullptr;
	fopen_s(&fp, outputfileName, "w+");
	if (nullptr != fp) {
		fprintf_s(fp, "%s,,%s\n", "畸变矫正前","畸变校正后");
		for (int i = 0; i < m_cornersAllImg.size(); ++i) {
			fprintf_s(fp, "%s%d%s\n", "imp", i + 1,":");
			cv::undistortPoints(m_cornersAllImg[i], m_undistortCornersPerImg, m_cameraInstricMatrix, m_distCoeffs, cv::noArray(), m_cameraInstricMatrix);
			for (int j = 0; j < m_undistortCornersPerImg.size(); ++j) {
				fprintf_s(fp, "%.5f,%.5f,%.5f,%.5f\n", m_cornersAllImg[i][j].x, m_cornersAllImg[i][j].y,
					m_undistortCornersPerImg[j].x, m_undistortCornersPerImg[j].y);
			}
		}
		fclose(fp);
	}

	bRet = true;
	return bRet;
}

bool CcalibrateLineCamera::CalUndistortedCorners(vector<cv::Point2f>& srcPointArray, vector<cv::Point2f>& dstPointArray) {
	bool bRet = false;

	cv::undistortPoints(srcPointArray, dstPointArray, m_cameraInstricMatrix, m_distCoeffs, cv::noArray(), m_cameraInstricMatrix);

	bRet = true;
	return bRet;
}

bool CcalibrateLineCamera::GetUndistortedCornersFromCsv(const char* outputFileName) {
	bool bRet = false;

	vector<cv::Point2f> distortCorners, undistortCorners;
	double x, y;  //改动
	FILE* fp = nullptr;
	fopen_s(&fp, outputFileName, "r+");
	char strTemp[1024] = { 0 };
	if (nullptr != fp) {
		fgets(strTemp, 1024, fp);
		while(!feof(fp))
		{
			fscanf_s(fp, "%f,%f\n", &y, &x);
			distortCorners.push_back(cv::Point2f(x,y));
		}
		cout << distortCorners.size() << endl;
		cv::undistortPoints(distortCorners, undistortCorners, m_cameraInstricMatrix, m_distCoeffs, cv::noArray(), m_cameraInstricMatrix);
		rewind(fp);
		fprintf_s(fp, "%s,,,,%s\n","畸变矫正前", "畸变矫正后：");
		for (int i = 0; i < undistortCorners.size(); ++i) {
			fprintf_s(fp, "%.5f,%.5f,,,%.5f,%.5f\n", distortCorners[i].y, distortCorners[i].x,undistortCorners[i].y, undistortCorners[i].x);
		}
		fclose(fp);
	}

	bRet = true;
	return bRet;
}

bool CcalibrateLineCamera::GetUndistortedCornersFromBmp(const char* bmpPath, const char* outputFileName) {
	bool bRet = false;

	HalconCpp::HObject Img;
	HalconCpp::HTuple RowPoint, ColPoint, Index, Pose, imgWidth, imgHeight;
	vector<cv::Point2f> ptCorners,undistortCorners; 
	vector<cv::Point3f> worldPoints;
	HalconCpp::HTuple CalibDataID;
	HalconCpp::CreateCalibData("calibration_object", 1, 1, &CalibDataID);
	HalconCpp::HTuple startCamPar;
	startCamPar.Clear();
	int temp3 = 0;
	for (auto k = begin(m_calibObject.m_startCamPar); k != end(m_calibObject.m_startCamPar); ++k) {
		startCamPar[temp3] = *k;
		temp3++;
	}
	HalconCpp::SetCalibDataCamParam(CalibDataID, 0, HalconCpp::HTuple(), startCamPar);
	HalconCpp::SetCalibDataCalibObject(CalibDataID, 0,m_strCalTabCpdPath.c_str());

	//读取图片，找点
	HalconCpp::ReadImage(&Img, bmpPath);
	HalconCpp::GetImageSize(Img, &imgWidth, &imgHeight);
	HalconCpp::FindCalibObject(Img, CalibDataID, 0, 0,  1, HalconCpp::HTuple(), HalconCpp::HTuple());
	HalconCpp::GetCalibDataObservPoints(CalibDataID, 0, 0,  1, &ColPoint, &RowPoint, &Index, &Pose);

	for (int j = 0; j < RowPoint.Length(); ++j) {
		ptCorners.push_back(cv::Point2f(RowPoint[j], ColPoint[j]));
		worldPoints.push_back(cv::Point3f(m_calibObject.m_XofWorldPt[Index[j]], m_calibObject.m_YofWorldPt[Index[j]], 0.0));
	}
	cv::undistortPoints(ptCorners, undistortCorners, m_cameraInstricMatrix, m_distCoeffs, cv::noArray(), m_cameraInstricMatrix);

	FILE* fp = nullptr;
	fopen_s(&fp, outputFileName, "w+");
	if (nullptr != fp) {
		fprintf_s(fp, "%s,,%s\n", "畸变矫正前：","畸变矫正后");
		for (int i = 0; i < undistortCorners.size(); ++i) {
			fprintf_s(fp, "%.5f,%.5f,%.5f,%.5f\n", ptCorners[i].x, ptCorners[i].y,undistortCorners[i].x, undistortCorners[i].y);
		}
		fclose(fp);
	}
	ptCorners.clear();
	undistortCorners.clear();
	worldPoints.clear();

	bRet = true;
	return bRet;
}

//bool CcalibrateLineCamera::Clone(CParammeter* pParam) {
//	bool bRet = false;
//
//	m_strImgDir = pParam->m_strImgDir;
//	m_imgNum = pParam->m_imgNum;
//	m_strCalTabCpdPath = pParam->m_strCalTabCpdPath;
//	m_strCalTabPsPath = pParam->m_strCalTabPsPath;
//
//	bRet = true;
//	return bRet;
//}

bool CcalibrateLineCamera::SetParam(CPathParam* pParam) {
	bool bRet = false;

	m_strImgDir = pParam->m_strImgDir;
	m_iImgNum = pParam->m_iImgNum;
	m_strCalTabCpdPath = pParam->m_strCalTabCpdPath;
	m_strCalTabPsPath = pParam->m_strCalTabPsPath;
	bRet = true;
	return	bRet;
}



bool CPathParam::SetParamValue(string strImgDir, int imgNum, string strCalTabCpdPath, string strCalTabPsPath) {
	bool bRet = false;
	//参与标定的照片文件夹路径和照片数量
	m_strImgDir = strImgDir;
	m_iImgNum = imgNum;
	//标定板文件生成位置，重复运行程序只会覆盖已生成文件，不会产生新文件
	m_strCalTabCpdPath = strCalTabCpdPath;
	m_strCalTabPsPath = strCalTabPsPath;

	bRet = true;
	return bRet;
}




CCalibObject::CCalibObject() {
	//kXw,kYw是HG-450标定板Corners的世界坐标，Zw=0，使用以下默认值即可
	const double kXw[195] = { -0.21, -0.18, -0.15, -0.12, -0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, -0.225, -0.195, -0.165, -0.135, -0.105, -0.075, -0.045, -0.015, 0.015, 0.045, 0.075, 0.105, 0.135, 0.165, 0.195, -0.21, -0.18, -0.15, -0.12, -0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, -0.225, -0.195, -0.165, -0.135, -0.105, -0.075, -0.045, -0.015, 0.015, 0.045, 0.075, 0.105, 0.135, 0.165, 0.195, -0.21, -0.18, -0.15, -0.12, -0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, -0.225, -0.195, -0.165, -0.135, -0.105, -0.075, -0.045, -0.015, 0.015, 0.045, 0.075, 0.105, 0.135, 0.165, 0.195, -0.21, -0.18, -0.15, -0.12, -0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, -0.225, -0.195, -0.165, -0.135, -0.105, -0.075, -0.045, -0.015, 0.015, 0.045, 0.075, 0.105, 0.135, 0.165, 0.195, -0.21, -0.18, -0.15, -0.12, -0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, -0.225, -0.195, -0.165, -0.135, -0.105, -0.075, -0.045, -0.015, 0.015, 0.045, 0.075, 0.105, 0.135, 0.165, 0.195, -0.21, -0.18, -0.15, -0.12, -0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, -0.225, -0.195, -0.165, -0.135, -0.105, -0.075, -0.045, -0.015, 0.015, 0.045, 0.075, 0.105, 0.135, 0.165, 0.195, -0.21, -0.18, -0.15, -0.12, -0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21 };
	const double kYw[195] = { -0.155885, -0.155885, -0.155885, -0.155885, -0.155885, -0.155885, -0.155885, -0.155885, -0.155885, -0.155885, -0.155885, -0.155885, -0.155885, -0.155885, -0.155885, -0.129904, -0.129904, -0.129904, -0.129904, -0.129904, -0.129904, -0.129904, -0.129904, -0.129904, -0.129904, -0.129904, -0.129904, -0.129904, -0.129904, -0.129904, -0.103923, -0.103923, -0.103923, -0.103923, -0.103923, -0.103923, -0.103923, -0.103923, -0.103923, -0.103923, -0.103923, -0.103923, -0.103923, -0.103923, -0.103923, -0.0779423, -0.0779423, -0.0779423, -0.0779423, -0.0779423, -0.0779423, -0.0779423, -0.0779423, -0.0779423, -0.0779423, -0.0779423, -0.0779423, -0.0779423, -0.0779423, -0.0779423, -0.0519615, -0.0519615, -0.0519615, -0.0519615, -0.0519615, -0.0519615, -0.0519615, -0.0519615, -0.0519615, -0.0519615, -0.0519615, -0.0519615, -0.0519615, -0.0519615, -0.0519615, -0.0259808, -0.0259808, -0.0259808, -0.0259808, -0.0259808, -0.0259808, -0.0259808, -0.0259808, -0.0259808, -0.0259808, -0.0259808, -0.0259808, -0.0259808, -0.0259808, -0.0259808, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0259808, 0.0259808, 0.0259808, 0.0259808, 0.0259808, 0.0259808, 0.0259808, 0.0259808, 0.0259808, 0.0259808, 0.0259808, 0.0259808, 0.0259808, 0.0259808, 0.0259808, 0.0519615, 0.0519615, 0.0519615, 0.0519615, 0.0519615, 0.0519615, 0.0519615, 0.0519615, 0.0519615, 0.0519615, 0.0519615, 0.0519615, 0.0519615, 0.0519615, 0.0519615, 0.0779423, 0.0779423, 0.0779423, 0.0779423, 0.0779423, 0.0779423, 0.0779423, 0.0779423, 0.0779423, 0.0779423, 0.0779423, 0.0779423, 0.0779423, 0.0779423, 0.0779423, 0.103923, 0.103923, 0.103923, 0.103923, 0.103923, 0.103923, 0.103923, 0.103923, 0.103923, 0.103923, 0.103923, 0.103923, 0.103923, 0.103923, 0.103923, 0.129904, 0.129904, 0.129904, 0.129904, 0.129904, 0.129904, 0.129904, 0.129904, 0.129904, 0.129904, 0.129904, 0.129904, 0.129904, 0.129904, 0.129904, 0.155885, 0.155885, 0.155885, 0.155885, 0.155885, 0.155885, 0.155885, 0.155885, 0.155885, 0.155885, 0.155885, 0.155885, 0.155885, 0.155885, 0.155885 };
	for (int i = 0; i < 195; ++i) {
		m_XofWorldPt[i] = kXw[i];
		m_YofWorldPt[i] = kYw[i];
	}
	//使用以下默认值即可，halcon需要设置相机参数才能提取点，但标定不需要halcon
	//m_StartCamPar.Clear();
	m_startCamPar[0] = 0.55;
	m_startCamPar[1] = 0;
	m_startCamPar[2] = 0.45;
	m_startCamPar[3] = 0.45;
	m_startCamPar[4] = 968;
	m_startCamPar[5] = 732;
	m_startCamPar[6] = 1936;
	m_startCamPar[7] = 1464;
	//HG-450标定板行列数、点直径、mark区域，使用默认值即可
	m_cornerDiameter = 15 * 0.001;  //mm
	m_iCalTabRow = 13;
	m_iCalTabCol = 15;
	//m_markRow.Clear();
	m_iMarkRow[0] = 6;
	m_iMarkRow[1] = 2;
	m_iMarkRow[2] = 2;
	m_iMarkRow[3] = 10;
	m_iMarkRow[4] = 10;
	//m_markCol.Clear();
	m_iMarkCol[0] = 7;
	m_iMarkCol[1] = 2;
	m_iMarkCol[2] = 12;
	m_iMarkCol[3] = 2;
	m_iMarkCol[4] = 12;
	//图片长宽使用默认值0即可，提取过程中会由提取图像的函数改为正确值
	m_iImgWidth = 0;
	m_iImgHeight = 0;
	//标定板黑底白点
	m_strPolarity = "light_on_dark";

}

bool CcalibrateLineCamera::PrintCamParam() {
	bool bRet = false;
	int cnt = 0;
	cout << " Camera intrinsic: " /*<< cameraMatrix.rows << "x" << cameraMatrix.cols*/ << endl;
	cout << m_cameraInstricMatrix << endl;
	/*for (int i = 0; i < m_cameraInstricMatrix.cols; ++i) {
		for (int j = 0; j < m_cameraInstricMatrix.rows; ++j) {
			printf("%f,", m_cameraInstricMatrix.at<double>(i, j));
			cnt++;
			if (cnt % 3 == 0) {
				printf("\n");
				cnt = 0;
			}
		}
	}*/
	cout << "distCoeffs:" << endl;
	cout << m_distCoeffs << endl;
	//cout << m_distCoeffs.rows << endl;
	//printf("%.10f\t", m_distCoeffs.at<double>(1, 0));
	bRet = true;
	return bRet;
}
CcalibrateLineCamera::CcalibrateLineCamera(const char* inputCamParamCsv) {
	m_cameraInstricMatrix = cv::Mat::zeros(3, 3, CV_64F);
	m_distCoeffs = cv::Mat::zeros(14, 1, CV_64F);
	FILE* fp = nullptr;
	fopen_s(&fp, inputCamParamCsv, "r+");
	if (nullptr != fp) {
		char strTemp[1024] = { 0 };
		char c[1024];
		char strTemp2[1024] = { 0 };
		fgets(strTemp, 1024, fp);
		//cout << strTemp << endl;
		for (int i = 0; i < 3; ++i) {
			fscanf_s(fp, "%s\n", c, 1024);
			//cout << c << endl;
			double fTemp1, fTemp2, fTemp3;  //改动
			sscanf_s(c, "%f,%f,%f,", &fTemp1, &fTemp2, &fTemp3);
			//fscanf_s(fp, "%f,%f,%f\n", &fTemp1, &fTemp2, &fTemp3);
			//cout << fTemp1 << " " << fTemp2 << " " << fTemp3 << endl;
			m_cameraInstricMatrix.at<double>(i, 0) = (double)fTemp1;
			m_cameraInstricMatrix.at<double>(i, 1) = (double)fTemp2;
			m_cameraInstricMatrix.at<double>(i, 2) = (double)fTemp3;
			////fscanf_s(fp, "%f,%f,%f\n", &m_cameraInstricMatrix.at<float>(i,0), &m_cameraInstricMatrix.at<float>(i, 1), &m_cameraInstricMatrix.at<float>(i, 2));
			//cout << m_cameraInstricMatrix.at<double>(i, 0) << " " << m_cameraInstricMatrix.at<double>(i, 1) << " " << m_cameraInstricMatrix.at<double>(i, 2) << endl;

		}
		fgets(strTemp2, 1024, fp);
		//fscanf_s(fp, "%s\n",sTemp);
		for (int i = 0; i < 14; ++i) {
			double fTemp;     //改动
			fscanf_s(fp, "%f,", &fTemp);
			//printf("%.18f\n", fTemp);
			m_distCoeffs.at<double>(i, 0) = fTemp;
			//printf("%.18f\n", m_distCoeffs.at<double>(i, 0));
			//cout << m_distCoeffs.at<double>(i, 0) << " ";
		}
		fclose(fp);
	}
}

bool CcalibrateLineCamera::Calibrate(string strImgDir, int imgNum, string strCalTabCpdPath, string strCalTabPsPath, string strOutputCsvOfCamParam) {
	bool bRet = false;

	CPathParam param;
	param.SetParamValue(strImgDir, imgNum, strCalTabCpdPath,strCalTabPsPath);
	SetParam(&param);
	CalibrateUsingCV(strOutputCsvOfCamParam.c_str());

	bRet = true;
	return bRet;
}

void CcalibrateLineCamera::GetLaserLinePlaneFunc(const char* sArr[5]) {
	double u, v;
	/*NOTE!!!!!! 照片数量改动，这里也需要改动！！*/
	vector<int> numOfLaserImg = { 16,17,18,19,20 };
	int numOfLaserImg_ = 5;
	vector<vector<cv::Point2d>> laserLinePointsWithDistortAll, laserLinePointsWithoutDistortAll;
	vector<cv::Point2d> laserLinePointsWithDistortPer, laserLinePoints2WithDistort, laserLinePointsWithoutDistortPer, laserLinePoints2WithoutDistort;

	for (int i = 0; i < numOfLaserImg_; ++i) {
		FILE* fp = nullptr;
		fopen_s(&fp, sArr[i], "r+");
		char strTemp[1024] = { 0 };
		if (nullptr != fp) {
			fgets(strTemp, 1024, fp);
			while (!feof(fp))
			{
				fscanf_s(fp, "%lf,%lf\n", &u, &v);
				laserLinePointsWithDistortPer.push_back(cv::Point2d(u, v));
			}
			cv::undistortPoints(laserLinePointsWithDistortPer, laserLinePointsWithoutDistortPer, m_cameraInstricMatrix, m_distCoeffs, cv::noArray(), m_cameraInstricMatrix);
			laserLinePointsWithDistortAll.push_back(laserLinePointsWithDistortPer);
			laserLinePointsWithoutDistortAll.push_back(laserLinePointsWithoutDistortPer);
			laserLinePointsWithDistortPer.clear();
			laserLinePointsWithoutDistortPer.clear();
			fclose(fp);
		}
	}

	vector<cv::Point3d> xcyczc;

	for (int i=0;i<laserLinePointsWithoutDistortAll.size();++i)
	{	
		cv::Mat rVec_64F,r, r_64F,rt,rt_inv,rt_inv_64F;
		cv::Mat mat0001 = (cv::Mat_<double>(1, 4) << 0,0,0,1);//直接赋初始值的方法
		m_rVecs[numOfLaserImg[i]].convertTo(rVec_64F, CV_64FC1);
		cv::Rodrigues(rVec_64F, r);
		r.convertTo(r_64F, CV_64FC1);
		cv::hconcat(r_64F, m_tVecs[numOfLaserImg[i]],rt);
		cv::vconcat(rt, mat0001, rt);
		cout << rt << endl;
		cv::invert(rt, rt_inv);
		rt_inv.convertTo(rt_inv_64F, CV_64FC1);
		for (int j = 0; j < laserLinePointsWithoutDistortAll[i].size(); ++j) {
			cv::Mat A(3, 3, CV_64FC1, cv::Scalar(0));
			cv::Mat x(3, 1, CV_64FC1, cv::Scalar(0));
			cv::Mat b(3, 1, CV_64FC1, cv::Scalar(0));
			A.at<double>(0, 0) = m_cameraInstricMatrix.at<double>(0, 0);
			A.at<double>(1, 1) = m_cameraInstricMatrix.at<double>(1, 1);
			A.at<double>(0, 2) = m_cameraInstricMatrix.at<double>(0, 2) - laserLinePointsWithoutDistortAll[i][j].x;
			A.at<double>(1, 2) = m_cameraInstricMatrix.at<double>(1, 2) - laserLinePointsWithoutDistortAll[i][j].y;
			A.at<double>(2, 0) = rt_inv_64F.at<double>(2, 0);
			A.at<double>(2, 1) = rt_inv_64F.at<double>(2, 1);
			A.at<double>(2, 2) = rt_inv_64F.at<double>(2, 2);
			b.at<double>(2, 0) = -1 * rt_inv_64F.at<double>(2, 3);; //double float?????????????????????????????
			cv::solve(A, b, x);
			xcyczc.push_back(cv::Point3d(x.at<double>(0, 0), x.at<double>(1, 0), x.at<double>(2, 0)));
		}
	}
	vector<double> PlaneCoffs;
	//planeFittingUsingPclRanSac(m_vecLaserPlaneCoeffs, xcyczc);
	planeFitting(m_vecLaserPlaneCoeffs,xcyczc);
	//cout << PlaneCoffs[0] << PlaneCoffs[1] << PlaneCoffs[02] << endl;


}
//bool CcalibrateLineCamera::GetLaserLinePlaneFunc(const char* sArr[5])
//{
//	bool bRet = false;
//
//	
//	//从csv文件中读取五张光条中心像素坐标并去畸变//
//	double u, v;
//	vector<int> numOfLaserImg = { 16,17,18,19,20 };
//	int numOfLaserImg_ = 5;
//	vector<vector<cv::Point2d>> laserLinePointsWithDistortAll, laserLinePointsWithoutDistortAll;
//	vector<cv::Point2d> laserLinePointsWithDistortPer, laserLinePoints2WithDistort, laserLinePointsWithoutDistortPer, laserLinePoints2WithoutDistort;
//
//	for (int i = 0; i < numOfLaserImg_; ++i) {
//		FILE* fp = nullptr;
//		fopen_s(&fp, sArr[i], "r+");
//		char strTemp[1024] = { 0 };
//		if (nullptr != fp) {
//			fgets(strTemp, 1024, fp);
//			while (!feof(fp))
//			{
//				fscanf_s(fp, "%lf,%lf\n", &u, &v);
//				laserLinePointsWithDistortPer.push_back(cv::Point2d(u, v));
//			}
//			cv::undistortPoints(laserLinePointsWithDistortPer, laserLinePointsWithoutDistortPer, m_cameraInstricMatrix, m_distCoeffs, cv::noArray(), m_cameraInstricMatrix);
//			laserLinePointsWithDistortAll.push_back(laserLinePointsWithDistortPer);
//			laserLinePointsWithoutDistortAll.push_back(laserLinePointsWithoutDistortPer);
//			laserLinePointsWithDistortPer.clear();
//			laserLinePointsWithoutDistortPer.clear();
//			fclose(fp);
//		}
//	}
//
//	//计算五张图片中光条中心点的xcyczc坐标//
//
//	vector<cv::Mat> rotationMatrix;
//	vector<cv::Mat> RT, Homo34, Homo;
//	for (int i = 0; i <numOfLaserImg_; ++i) {
//		cv::Mat tempRotationMatrix, tempRT, tempH_12, tempH_4, tempHomo;
//		cv::Rodrigues(m_rVecs[numOfLaserImg[i]], tempRotationMatrix);
//		rotationMatrix.push_back(tempRotationMatrix);
//		cv::hconcat(rotationMatrix[i], m_tVecs[numOfLaserImg[i]], tempRT);
//		RT.push_back(tempRT);
//		cv::Mat cameraInstricMatrix_64F(3, 3, CV_64FC1), RT_64F(3, 4, CV_64FC1);
//		m_cameraInstricMatrix.convertTo(cameraInstricMatrix_64F, CV_64FC1);
//		RT[i].convertTo(RT_64F, CV_64FC1);
//		Homo34.push_back(cameraInstricMatrix_64F * RT_64F); //!!!NOTE: CV_64F !!!!
//		tempH_12 = Homo34[i](cv::Rect(0, 0, 2, 3));
//		tempH_4 = Homo34[i](cv::Rect(3, 0, 1, 3));
//		cv::hconcat(tempH_12, tempH_4, tempHomo);
//		Homo.push_back(tempHomo);
//		}
//
//	try {
//		vector<cv::Mat> uv1_T, XY, XcYcZc;
//		for (int i = 0; i < numOfLaserImg_; i++) {
//			cv::Mat temp_uv1_T, temp_H_inverse(3, 3, CV_64FC1), temp_H_inverse_64F(3, 3, CV_64FC1);
//			cv::Mat uv(laserLinePointsWithoutDistortAll[i]);
//			cout << laserLinePointsWithoutDistortAll[i].size() << endl;
//			cv::Mat dest = uv.reshape(1, laserLinePointsWithoutDistortAll[i].size()).clone();
//			cv::Mat ones = cv::Mat::ones(cv::Size(1, laserLinePointsWithoutDistortAll[i].size()), CV_64FC1);
//			cv::Mat dest_64F(dest.size(), CV_64FC1);
//			dest.convertTo(dest_64F, CV_64FC1);
//			cv::hconcat(dest_64F, ones, temp_uv1_T);
//			uv1_T.push_back(temp_uv1_T);
//			cv::invert(Homo[i], temp_H_inverse);
//			temp_H_inverse.convertTo(temp_H_inverse_64F, CV_64FC1);
//			XY.push_back(temp_H_inverse_64F * (uv1_T[i].t()));
//			XY[i].row(0) /= XY[i].row(2);
//			XY[i].row(1) /= XY[i].row(2);
//			XY[i].row(2) /= XY[i].row(2);
//			XY[i].row(2) = 0;
//			XY[i].push_back(cv::Mat::ones(cv::Size(laserLinePointsWithoutDistortAll[i].size(), 1), CV_64FC1));
//
//			cv::Mat temp_RT_64F;
//			RT[i].convertTo(temp_RT_64F, CV_64FC1);
//			XcYcZc.push_back(temp_RT_64F * XY[i]);
//			}
//		cv::Mat XcYcZcAll,XcYcZcAll_64F;
//		XcYcZcAll = XcYcZc[0];
//		for (int i = 0; i < numOfLaserImg_-1; ++i) {
//				
//				cv::hconcat(XcYcZcAll, XcYcZc[i+1], XcYcZcAll);
//			}
//		cout << XcYcZcAll.size() << endl;
//		XcYcZcAll.convertTo(XcYcZcAll_64F, CV_64FC1);
//		vector<cv::Point3d> xcyczc;
//		for (int i = 0; i < XcYcZcAll_64F.cols; ++i) {
//				xcyczc.push_back(cv::Point3d(XcYcZcAll_64F.at<double>(0, i), XcYcZcAll_64F.at<double>(1, i), XcYcZcAll_64F.at<double>(2, i)));
//			}
//		cout << xcyczc.size() << endl;
//		vector<double> coeffient,coeffient2;
//	//	planeFittingUsingPclRanSac(m_vecLaserPlaneCoeffs, xcyczc);
//	//	planeFitting(m_vecLaserPlaneCoeffs,xcyczc);
//		//cout << coeffient2[0] << " " << coeffient2[1] << " " << coeffient2[2] << " " << coeffient2[3] << endl;
//		FILE* fp1 = nullptr;
//		fopen_s(&fp1, "C:\\Users\\YYKJ05\\Desktop\\planeBias.csv", "w+");
//		char strTemp[1024] = { 0 };
//		if (nullptr != fp1) {
//			fgets(strTemp, 1024, fp1);
//			for (auto &i : xcyczc)
//			{
//				double bias1;
//				bias1 = abs(m_vecLaserPlaneCoeffs[0] * i.x + m_vecLaserPlaneCoeffs[1] * i.y - i.z + m_vecLaserPlaneCoeffs[2])
//					/ sqrt(pow(m_vecLaserPlaneCoeffs[0], 2) + pow(m_vecLaserPlaneCoeffs[1], 2) + 1);
//// 				bias2 = abs(m_vecLaserPlaneCoeffs[0] * i.x + m_vecLaserPlaneCoeffs[1] * i.y + m_vecLaserPlaneCoeffs[2] * i.z + m_vecLaserPlaneCoeffs[3])
//// 					/ sqrt(pow(m_vecLaserPlaneCoeffs[0], 2) + pow(m_vecLaserPlaneCoeffs[1], 2) + pow(m_vecLaserPlaneCoeffs[2], 2));
//				fprintf_s(fp1, "%f\n", bias1);
//			}
//		}
//		//二次曲面拟合激光平面//
//	/*	quadricSurfaceFitting(m_vecLaserQuadricSurfaceCoeffs, xcyczc);
//		cout << m_vecLaserQuadricSurfaceCoeffs[0] << m_vecLaserQuadricSurfaceCoeffs[1] << m_vecLaserQuadricSurfaceCoeffs[2] << m_vecLaserQuadricSurfaceCoeffs[3]
//			<< m_vecLaserQuadricSurfaceCoeffs[4] << m_vecLaserQuadricSurfaceCoeffs[5] << endl;
//
//		cv::Mat x2 = XcYcZcAll_64F.row(0).mul(XcYcZcAll_64F.row(0));
//		cv::Mat xy = XcYcZcAll_64F.row(0).mul(XcYcZcAll_64F.row(1));
//		cv::Mat y2 = XcYcZcAll_64F.row(1).mul(XcYcZcAll_64F.row(1));
//		cv::Mat x = XcYcZcAll_64F.row(0);
//		cv::Mat y = XcYcZcAll_64F.row(1);
//		cv::Mat shouldBeZero = m_vecLaserQuadricSurfaceCoeffs[0] * x2 + m_vecLaserQuadricSurfaceCoeffs[1] * xy + m_vecLaserQuadricSurfaceCoeffs[2] * y2 + m_vecLaserQuadricSurfaceCoeffs[3] * x
//			+ m_vecLaserQuadricSurfaceCoeffs[4] * y + m_vecLaserQuadricSurfaceCoeffs[5];*/
//	//	cv::Mat shouldBeZero = XcYcZcAll_64F.row(0) * coeffient[0] + XcYcZcAll_64F.row(1) * coeffient[1] + coeffient[2] - XcYcZcAll_64F.row(2);
//	//	cv::Mat shouldBeZero_ = XcYcZcAll_64F.row(0) * coeffient2[0] + XcYcZcAll_64F.row(1) * coeffient2[1] + XcYcZcAll_64F.row(2)*coeffient2[2]+coeffient2[3];
//
//		}
//		catch (cv::Exception& e) { cout << e.what() << endl; }
//		bRet = true;
//		return bRet;
//	
//}




bool CcalibrateLineCamera::planeFitting(vector<double>& coeffient,const vector<cv::Point3d>& input)
{	
	bool bRet = false;
	
	cv::Mat dst = cv::Mat(3, 3, CV_64F, cv::Scalar(0));//初始化系数矩阵A
	cv::Mat out = cv::Mat(3, 1, CV_64F, cv::Scalar(0));//初始化矩阵b
	for (int i = 0; i < input.size(); ++i)
	{
		//计算3*3的系数矩阵
		dst.at<double>(0, 0) = dst.at<double>(0, 0) + pow((double)input[i].x, 2);
		dst.at<double>(0, 1) = dst.at<double>(0, 1) + (double)input[i].x * (double)input[i].y;
		dst.at<double>(0, 2) = dst.at<double>(0, 2) + (double)input[i].x;
		dst.at<double>(1, 0) = dst.at<double>(1, 0) + (double)input[i].x * (double)input[i].y;
		dst.at<double>(1, 1) = dst.at<double>(1, 1) + pow((double)input[i].y, 2);
		dst.at<double>(1, 2) = dst.at<double>(1, 2) + (double)input[i].y;
		dst.at<double>(2, 0) = dst.at<double>(2, 0) + (double)input[i].x;
		dst.at<double>(2, 1) = dst.at<double>(2, 1) + (double)input[i].y;
		dst.at<double>(2, 2) = (double)input.size();
		//计算3*1的结果矩阵
		out.at<double>(0, 0) = out.at<double>(0, 0) + (double)input[i].x * (double)input[i].z;
		out.at<double>(1, 0) = out.at<double>(1, 0) + (double)input[i].y * (double)input[i].z;
		out.at<double>(2, 0) = out.at<double>(2, 0) + (double)input[i].z;
	}
	//判断矩阵是否奇异
	//double determ = determinant(dst);
	//if (abs(determ) < 0.001) {
	//	cout << "矩阵奇异" << endl;
	//	return false;
	//}
	//Mat inv;
	//invert(dst, inv);//求矩阵的逆
	//Mat output = inv*out;//计算输出
	//coeffient.clear();//把结果输出
	//coeffient.push_back(output.at<float>(0, 0));
	//coeffient.push_back(output.at<float>(1, 0));
	//coeffient.push_back(output.at<float>(2, 0));

	// 修改此处代码为SVD分解
	cv::Mat result = cv::Mat(3, 1, CV_64F,cv::Scalar(0));
	//cout << dst << out << endl;
	solve(dst, out, result, cv::DECOMP_SVD);
	cout <<"系数为：" << result << endl;
	coeffient.clear();//把结果输出
	coeffient.push_back(result.at<double>(0, 0));
	coeffient.push_back(result.at<double>(1, 0));
	coeffient.push_back(result.at<double>(2, 0));
	//TestAcurracyUsingTrueData("C:\\Users\\YYKJ05\\Desktop\\trueData.csv", coeffient);
	bRet = true;
	return bRet;
}

bool CcalibrateLineCamera::quadricSurfaceFitting(vector<double>& coeffs, const vector<cv::Point3d>& xyz)
{	
	bool bRet = false;
	//z=ax*2+bxy+cy*2+dx+ey+f,最小二乘法求偏导为0得超定方程组Ax=b
	cv::Mat A = cv::Mat(6, 6, CV_64FC1, cv::Scalar(0));
	cv::Mat b = cv::Mat(6, 1, CV_64FC1, cv::Scalar(0));
	cv::Mat x = cv::Mat(6, 1, CV_64FC1, cv::Scalar(0));

	for (auto point : xyz) {
		A.at<double>(0, 0) += pow(point.x, 4);
		A.at<double>(1, 1) += pow(point.x * point.y, 2);
		A.at<double>(2, 2)+= pow(point.y, 4);
		A.at<double>(3, 3) += pow(point.x, 2);
		A.at<double>(4, 4) += pow(point.y, 2);
		A.at<double>(5, 5) += 1;

		A.at<double>(0, 1) += (pow(point.x, 3) * point.y);
		A.at<double>(0, 2) += (pow(point.x, 2) * pow(point.y, 2));
		A.at<double>(0, 3) += pow(point.x, 3);
		A.at<double>(0, 1) += (pow(point.x, 2) * point.y);
		A.at<double>(0, 3) += pow(point.x, 2);

		A.at<double>(1, 2) += (point.x * pow(point.y, 3));
		A.at<double>(1, 3) += (pow(point.x, 2) * pow(point.y, 1));
		A.at<double>(1, 4) += (pow(point.x, 1) * pow(point.y, 2));
		A.at<double>(1, 5) += (pow(point.x, 1) * pow(point.y, 1));

		A.at<double>(2, 3) += (pow(point.x, 1) * pow(point.y, 2));
		A.at<double>(2, 4) += (1 * pow(point.y, 3));
		A.at<double>(2, 5) += (1 * pow(point.y, 2));

		A.at<double>(3, 4) += (pow(point.x, 1) * pow(point.y, 1));
		A.at<double>(3, 5) += (pow(point.x, 1) * 1);

		A.at<double>(4, 5) += (1 * pow(point.y, 1));
	
		b.at<double>(0, 0) += (pow(point.z, 1) * pow(point.x, 2));
		b.at<double>(0, 1) += (pow(point.z, 1) * pow(point.x, 1)* pow(point.y, 1));
		b.at<double>(0, 2) += (pow(point.z, 1) * pow(point.y, 2));
		b.at<double>(0, 3) += (pow(point.z, 1) * pow(point.x, 1));
		b.at<double>(0, 4) += (pow(point.z, 1) * pow(point.y, 1));
		b.at<double>(0, 5) += (pow(point.z, 1) * 1);

	}
	A.at<double>(1, 0) = A.at<double>(0, 1);
	A.at<double>(2, 0) = A.at<double>(0, 2);
	A.at<double>(2, 1) = A.at<double>(1, 2);
	A.at<double>(3, 0) = A.at<double>(0, 3);
	A.at<double>(3, 1) = A.at<double>(1, 3);
	A.at<double>(3, 2) = A.at<double>(2, 3);
	A.at<double>(4, 0) = A.at<double>(0, 4);
	A.at<double>(4, 1) = A.at<double>(1, 4);
	A.at<double>(4, 2) = A.at<double>(2, 4);
	A.at<double>(4, 3) = A.at<double>(3, 4);
	A.at<double>(5, 0) = A.at<double>(0, 5);
	A.at<double>(5, 1) = A.at<double>(1, 5);
	A.at<double>(5, 2) = A.at<double>(2, 5);
	A.at<double>(5, 3) = A.at<double>(3, 5);
	A.at<double>(5, 4) = A.at<double>(4, 5);

	//求解Ax=b,用SVD

	cv::solve(A, b, x, cv::DECOMP_SVD);

	for (int i = 0; i < x.rows; ++i) {
		coeffs.push_back(x.at<double>(i, 0));
	}

	bRet = true;
	return bRet;
}

//
//bool CcalibrateLineCamera::planeFittingUsingPclRanSac(vector<double>& coefficient, const vector<cv::Point3d>& xcyczc) {
//	try
//	{
//		bool bRet = false;
//
//		PointCloudT::Ptr ptrXcyczcPointCloud(new PointCloudT);
//
//		for (int i = 0; i < xcyczc.size(); ++i) {
//			ptrXcyczcPointCloud->push_back(PointT(xcyczc[i].x, xcyczc[i].y, xcyczc[i].z));
//		}
//		pcl::SampleConsensusModelPlane<PointT>::Ptr ptrXcyczcModelPlane(new pcl::SampleConsensusModelPlane<PointT>(ptrXcyczcPointCloud));
//		pcl::RandomSampleConsensus<PointT> ransac(ptrXcyczcModelPlane);
//
//		ransac.setDistanceThreshold(0.000035);	//设置距离阈值，与平面距离小于0.01的点作为内点
//		ransac.computeModel();				//执行模型估计
//
//		PointCloudT::Ptr cloud_plane(new PointCloudT);
//
//		//---------- 根据索引提取内点 ----------
//		///方法1
//		vector<int> inliers;				//存储内点索引的向量
//		ransac.getInliers(inliers);			//提取内点对应的索引
//		pcl::copyPointCloud<PointT>(*ptrXcyczcPointCloud, inliers, *cloud_plane);
//
//		/// 输出模型参数Ax+By+Cz+D=0
//		Eigen::VectorXf coefficient_;
//		ransac.getModelCoefficients(coefficient_);
//		cout << coefficient_ << endl;
//		coefficient.push_back(-1*coefficient_[0]/ coefficient_[2]);
//		coefficient.push_back(-1*coefficient_[1]/ coefficient_[2]);
//		coefficient.push_back(-1*coefficient_[3]/ coefficient_[2]);
//		//cout << "平面方程为：\n"
//		//	<< coefficient_[0] << "x + "
//		//	<< coefficient_[1] << "y + "
//		//	<< coefficient_[2] << "z + "
//		//	<< coefficient_[3] << " = 0"
//		//	<< endl;
//		//cout << "系数为\n" 
//		//	<< -coefficient_[0] / coefficient_[2] 
//		//	<< " " << -coefficient_[1] / coefficient_[2] 
//		//	<< " " << -coefficient_[3] / coefficient_[2] 
//		//	<< endl;
//		//========================== 模型估计 ==========================
//
//		//----------------------- 可视化拟合结果 ----------------------
//		//new pcl::visualization::PCLVisualizer("result");
//		//pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("result"));
//
//		//viewer->addPointCloud<pcl::PointXYZ>(ptrXcyczcPointCloud, "cloud");													//添加原始点云
//		//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 1, 1, "cloud");	//颜色
//		//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud");	//点的大小
//
//		//viewer->addPointCloud<pcl::PointXYZ>(cloud_plane, "plane");												//添加平面点云
//		//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1, 0, "plane");	//颜色
//		//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "plane");	//点的大小
//
//		//while (!viewer->wasStopped())
//		//{
//		//	viewer->spinOnce(100);
//		//	/*	system("pause");
//		//		boost::this_thread::sleep(boost::posix_time::microseconds(100000));*/
//		//}
//
//		bRet = true;
//		return bRet;
//	}
//	catch (pcl::PCLException(&e))
//	{
//		cout << e.what() << endl;
//	}
//}


//bool CcalibrateLineCamera::TestAcurracyUsingTrueData(const char* trueDataFromCsv,vector<double>& coefficient) {
//	bool bRet = false;
//
//
//	vector<cv::Point2f> uv, xy,uv_withoutDistortion;
//
//	FILE* fp = nullptr;
//	fopen_s(&fp, trueDataFromCsv, "r+");
//	if (nullptr != fp) {
//		char temp[1024] = { 0 };
//		fgets(temp, 1024, fp);
//		double u, v, x, y;
//		while (!feof(fp)) {
//
//			fscanf_s(fp, "%lf,%lf,,,%lf,%lf\n", &u, &v, &x, &y);
//			uv.push_back(cv::Point2f(u, v));
//			xy.push_back(cv::Point2f(0.001*x, 0.001*y));
//			fgets(temp, 1024, fp);	
//		}
//		//cout << xy << endl;
//		fclose(fp);
//	}
//
//	cv::undistortPoints(uv, uv_withoutDistortion, m_cameraInstricMatrix, m_distCoeffs, cv::noArray(), m_cameraInstricMatrix);
//	 
//	vector<cv::Point3d> XY0_real;
//	vector<cv::Point2d> uv_real_without,uv_real_with;
//	for (int i = 0; i < uv_withoutDistortion.size(); ++i) {
//		uv_real_without.push_back(cv::Point2d((double)uv_withoutDistortion[i].x, (double)uv_withoutDistortion[i].y));
//		uv_real_with.push_back(cv::Point2d((double)uv[i].x, (double)uv[i].y));
//	}
//	for (int i = 0; i <xy.size(); ++i) {
//		XY0_real.push_back(cv::Point3d((double)xy[i].x, (double)xy[i].y, 0));
//	}
//
//
//
//	//cv::Mat	rVecs_real, tVecs_real;
//	//cv::solvePnP(XY0_real, uv_real_with, m_cameraInstricMatrix, m_distCoeffs, rVecs_real, tVecs_real);
//	//cv::Mat R_real;
//	//cv::Rodrigues(rVecs_real, R_real);
//	//cv::Mat H1;
//	//cv::hconcat(R_real, tVecs_real,H1);
//	//cv::Mat H1_12, H1_4;
//	//H1_12 = H1(cv::Rect(0, 0, 2, 3));
//	//H1_4 = H1(cv::Rect(3, 0, 1, 3));
//	//cv::Mat RT_final,RT_final_64F;
//	//cv::hconcat(H1_12, H1_4, RT_final);
//	//RT_final.convertTo(RT_final_64F, CV_64F);
//	//
//	//cv::Mat temp(xy);
//	//cv::Mat temp1=temp.reshape(1, xy.size()).clone().t();
//	//cv::Mat xy1;
//	//temp1.convertTo(xy1, CV_64F);
//	//cout << xy1.type() << endl;
//	////cout << cv::Mat::ones(cv::Size(xy.size(), 1), CV_64F) << endl;
//	//xy1.push_back(cv::Mat::ones(cv::Size(xy.size(), 1), CV_64F));
//	//cv::Mat xy1_64F,RT_64F;
//	//xy1.convertTo(xy1_64F, CV_64F);
//	//cv::Mat xcyczc=RT_final_64F* xy1_64F;
//	//cv::Mat xcyczc_cal(xcyczc.rows, xcyczc.cols, CV_64F);
//	//cv::Mat coeffi_mat=cv::Mat(3, 3, CV_64F,cv::Scalar(0));
//	//coeffi_mat.at<double>(0, 0) = -coefficient[0];
//	//coeffi_mat.at<double>(0, 1) = -coefficient[1];
//	//coeffi_mat.at<double>(0, 2) = 1;
//	//coeffi_mat.at<double>(1, 0) = m_cameraInstricMatrix.at<double>(0,0);
//	//coeffi_mat.at<double>(2, 1) = m_cameraInstricMatrix.at<double>(1, 1);
//	//cv::Mat B = cv::Mat(3, 1,CV_64F, cv::Scalar(0));
//	//B.at<double>(0, 0) = coefficient[2];
//
//	//for (int i = 0; i < xcyczc.cols; ++i) {
//	//	coeffi_mat.at<double>(1, 2) = m_cameraInstricMatrix.at<double>(0, 2)-uv_withoutDistortion[i].x;
//	//	coeffi_mat.at<double>(2, 2) = m_cameraInstricMatrix.at<double>(1, 2) - uv_withoutDistortion[i].y;
//	////	cout << coeffi_mat << B << endl;
//	//	double temp=1 + coefficient[0] * coeffi_mat.at<double>(1, 2) / coeffi_mat.at<double>(1, 0) + coefficient[1] * coeffi_mat.at<double>(2, 2) / coeffi_mat.at<double>(2, 1);
//	//	double zc = coefficient[2] / temp;
//	//	double xc = -coeffi_mat.at<double>(1, 2) / coeffi_mat.at<double>(1, 0) * zc;
//	//	double yc = -coeffi_mat.at<double>(2, 2) / coeffi_mat.at<double>(2,1) * zc;
//	//	xcyczc_cal.at<double>(2, i) = zc;
//	//	xcyczc_cal.at<double>(0, i) = xc;
//	//	xcyczc_cal.at<double>(1, i) = yc;
//
//	//	//cout << zc << endl;
//	////	cout << coeffi_mat * xcyczc_cal.col(i) << B << endl;
//	//	cv::solve(coeffi_mat, B, xcyczc_cal.col(i),cv::DECOMP_SVD);
//	//	//cout << coeffi_mat * xcyczc << endl;                                                                                                                                                                        _cal.col(i) << B << endl;
//	//	//cout << xcyczc_cal.col(i) << endl;
//	//	//cout << xcyczc.col(i) << endl;
//	//}
//	////cout << xcyczc.type() << xcyczc_cal.type() << endl;
//
//
//	vector<cv::Point3d> XY0_real_train, XY0_real_test;
//	vector<cv::Point2d> uv_real_train, uv_real_test;
//	cv::Mat xcyczc_train,xcyczc_test;
//	separateDataSet(XY0_real, uv_real_without, XY0_real_train, uv_real_train, XY0_real_test, uv_real_test);
//	xcyczc_train=getXcyczcFromUv(uv_real_train, m_cameraInstricMatrix);
//	xcyczc_test = getXcyczcFromUv(uv_real_test, m_cameraInstricMatrix);
//
//	cv::Mat rt= CorrectLaserPlane(xcyczc_train, XY0_real_train, uv_real_train);
//	//vector<cv::Mat> rt_test= CorrectLaserPlane(xcyczc_test, XY0_real_test, uv_real_test);
//	//cout << rt[0].size() << xcyczc_test.t().size() << rt[1].size() << endl;
//	//cout << rt[0].type() << xcyczc_test.type() << rt[1].type() << endl;
//	//cout << xcyczc_test.t().col(0).size() << endl;
//	//cout << xcyczc_test.size().height << endl;
//	for (int i = 0; i < xcyczc_test.size().width; ++i) {
//		xcyczc_test.col(i)-=rt.col(3);
//	}
//	cv::Mat r_inv, r_inv_64F;
//	cv::invert(rt(cv::Rect(0,0,3,3)), r_inv);
//	r_inv.convertTo(r_inv_64F, CV_64FC1);
//	cv::Mat xwywzw_cal_test = r_inv_64F * xcyczc_test;
//	/*cout << xwywzw_cal_test.size() << xwywzw_cal_test.row(0) << endl;
//	for (int i = 0; i < xcyczc_test.size().width; ++i) {
//		xwywzw_cal_test.col(i) += rt[1];
//	}*/
//	cv::Mat XY0_real_test_mat(3,XY0_real_test.size(),CV_64FC1,cv::Scalar(0));
//	for (int i = 0; i < XY0_real_test.size(); ++i) {
//		XY0_real_test_mat.at<double>(0, i) = XY0_real_test[i].x;
//		XY0_real_test_mat.at<double>(1, i) = XY0_real_test[i].y;
//		XY0_real_test_mat.at<double>(2, i) = XY0_real_test[i].z;
//	}
//
//	cout << "计算值与真值之差：" << endl;
//	cout << xwywzw_cal_test << endl;
//	cout << XY0_real_test_mat << endl;
//	cv::Mat sub_final = xwywzw_cal_test - XY0_real_test_mat;
//	FILE* fp4 = nullptr;
//	fopen_s(&fp4, "C:\\Users\\YYKJ05\\Desktop\\testFinalmethod2.csv", "w+");
//	if (nullptr != fp4) {
//		char temp[1024] = { 0 };
//		fgets(temp, 1024, fp4);
//		double delta;
//		for (int i = 0; i < sub_final.cols; ++i) {
//			delta = sub_final.at<double>(0, i) * sub_final.at<double>(0, i) + sub_final.at<double>(1, i) * sub_final.at<double>(1, i);
//			delta = sqrt(delta);
//			fprintf_s(fp4, "%f\n", delta);
//		}
//		fclose(fp4);
//	}
//
//	/*for (int i = 0; i < XY0_real_test.size(); ++i) {
//		double pow_ = pow((xwywzw_cal_test.at<double>(i,0) - XY0_real_test_mat.at<double>(0, i)), 2) 
//					+ pow((xwywzw_cal_test.at<double>(i, 1) - XY0_real_test_mat.at<double>(1, i)), 2);
//		double sqrt_ = sqrt(pow_);
//		cout << sqrt_ << endl;
//	}*/
//
//	//cv::Mat sub = xcyczc - xcyczc_cal;
//	//cv::Mat RT_final_inverse(3, 3, CV_64F), RT_final_inverse_64F;;
//	//cv::invert(RT_final_64F, RT_final_inverse,cv::DECOMP_SVD);
//	//RT_final_inverse.convertTo(RT_final_inverse_64F, CV_64F);
//	//cv::Mat xy_real = RT_final_inverse_64F * xcyczc_cal;
//	////vector<double> coefficient;
//	//vector<cv::Point3d> xwywzw;
//	//for (int i = 0; i < xwywzw.size(); ++i) {
//	//	xwywzw.push_back(cv::Point3d(xy_real.at<double>(0, i), xy_real.at<double>(1, i), xy_real.at<double>(2, i)));
//	//}
//	////planeFitting(coefficient, xwywzw);
//
//	//
//
//
//	//cout << xy_real.size() << endl;
//	//cout << xy_real - xy1_64F << endl;
//	//cout << "计算与真之差" << endl;
//	//cout<<sub << endl;
//	//FILE* fp3 = nullptr;
//	//fopen_s(&fp3, "C:\\Users\\YYKJ05\\Desktop\\testWithout.csv", "w+");
//	//if (nullptr != fp3) {
//	//	char temp[1024] = { 0 };
//	//	fgets(temp, 1024, fp3);
//	//	double delta;
//	//	for (int i = 0; i < sub.cols; ++i) {
//	//		delta = sub.at<double>(0, i) * sub.at<double>(0, i) + sub.at<double>(1, i) * sub.at<double>(1, i) + sub.at<double>(2, i) * sub.at<double>(2, i);
//	//		delta = sqrt(delta);
//	//		fprintf_s(fp3, "%f\n", delta);
//	//	}
//	//	fclose(fp3);
//	//}
//
//
//
//	bRet = true;
//	return bRet;
//}


cv::Mat CcalibrateLineCamera::CorrectLaserPlane(const cv::Mat& xcyczc_mat, vector<cv::Point3d> xwywzw_vecPoint3d, const vector<cv::Point2d>& uv)  //xwywzw_vecPoint3d uv xcyczc_mat是部分数据
{
	
	vector<cv::Point3d> xwyw0_vecPoint3d(xwywzw_vecPoint3d);
	cv::Mat rt;
	vector<double> coef_final;




	//迭代
	int iter = 200000;
	for (int i = 0; i < iter; ++i) {
		vector<double> coef;
		cv::Mat	rVec, tVec,rVec_64F,tVec_64F;
		cv::solvePnP(xwywzw_vecPoint3d, uv, m_cameraInstricMatrix, m_distCoeffs, rVec, tVec);
		rVec.convertTo(rVec_64F, CV_64FC1);
		tVec.convertTo(tVec_64F, CV_64FC1);
		cv::Mat R, R_inv, R_inv_64F,R_64F;
		cv::Rodrigues(rVec_64F, R);
		R.convertTo(R_64F, CV_64FC1);
		//cv::Mat H, H_64F;
		//cv::hconcat(R, tVec, H);
		//H.convertTo(H_64F, CV_64FC1);

		//R_invert*(xcyczc_mat-tvec)=>xwywzw_mat
		cv::Mat xwywzw_mat(xcyczc_mat.rows, xcyczc_mat.cols, CV_64FC1);
		//cout << tVec.size() << endl;
		for (int i = 0; i < xcyczc_mat.cols; ++i) {
			xwywzw_mat.col(i) = xcyczc_mat.col(i) - tVec_64F;
		}
		cv::invert(R_64F, R_inv);
		R_inv.convertTo(R_inv_64F, CV_64FC1);
		xwywzw_mat = R_inv_64F * xwywzw_mat;

		xwywzw_vecPoint3d.clear();
		for (int i = 0; i < xwywzw_mat.cols; ++i) {
			xwywzw_vecPoint3d.push_back(cv::Point3d(xwywzw_mat.at<double>(0, i), xwywzw_mat.at<double>(1, i), xwywzw_mat.at<double>(2, i)));
		}
	//	cout << xwywzw_vecPoint3d.size() << endl;
		planeFitting(coef, xwywzw_vecPoint3d);
		if (iter == (i + 1)) {
			for (int i = 0; i < xwywzw_vecPoint3d.size(); ++i) {
				cout << xwywzw_vecPoint3d[i].x - xwyw0_vecPoint3d[i].x
					<< xwywzw_vecPoint3d[i].y - xwyw0_vecPoint3d[i].y << endl;
			}
		}

		for (int i = 0; i < xwywzw_vecPoint3d.size(); ++i) {
			xwywzw_vecPoint3d[i].x = xwyw0_vecPoint3d[i].x;
			xwywzw_vecPoint3d[i].y = xwyw0_vecPoint3d[i].y;
			xwywzw_vecPoint3d[i].z = coef[0] * xwywzw_vecPoint3d[i].x + coef[1] * xwywzw_vecPoint3d[i].y + coef[2];
			//cout << xwywzw_vecPoint3d[i].x << " " << xwywzw_vecPoint3d[i].y << " " << xwywzw_vecPoint3d[i].z << endl;
		}

		if (iter == (i+1)) {
			cv::hconcat(R_64F, tVec_64F,rt);
			coef_final = coef;
		}
		
	}


	return rt;
}



//input:uv像素坐标值和xwywzw真实世界坐标值//
//output：uv和xwywzw的train、test//
bool CcalibrateLineCamera::separateDataSet( vector<cv::Point3d>& xwywzw_vecPoint3d_all, vector<cv::Point2d>& uv_all,  vector<cv::Point3d>& xwywzw_vecPoint3d_train, vector<cv::Point2d>& uv_train, vector<cv::Point3d>& xwywzw_vecPoint3d_test, vector<cv::Point2d>& uv_test)
{
	bool bRet = false;

	std::default_random_engine e;
	std::bernoulli_distribution u(0.8); // 生成1的概率为0.9
	e.seed(time(0));
	int flag;
	for (int i = 0; i < xwywzw_vecPoint3d_all.size(); i++) {
		flag = u(e);
		if (1 == flag) {
			xwywzw_vecPoint3d_train.push_back(xwywzw_vecPoint3d_all[i]);
			uv_train.push_back(uv_all[i]);
		}
		else {
			xwywzw_vecPoint3d_test.push_back(xwywzw_vecPoint3d_all[i]);
			uv_test.push_back(uv_all[i]);
		}
	}

	bRet = true;
	return bRet;
}

cv::Mat CcalibrateLineCamera::getXcyczcFromUv(vector<cv::Point2d> uv, cv::Mat cameraInstricMattrix)
{	
	vector<cv::Point2d> uv_withoutDistortion(uv);
	
	//cv::undistortPoints(uv, uv_withoutDistortion, m_cameraInstricMatrix, m_distCoeffs, cv::noArray(), m_cameraInstricMatrix);

	cv::Mat xcyczc_cal(3, uv.size(), CV_64F);
	cv::Mat coeffi_mat = cv::Mat(3, 3, CV_64F, cv::Scalar(0));
	coeffi_mat.at<double>(0, 0) = -m_vecLaserPlaneCoeffs[0];
	coeffi_mat.at<double>(0, 1) = -m_vecLaserPlaneCoeffs[1];
	coeffi_mat.at<double>(0, 2) = 1;
	coeffi_mat.at<double>(1, 0) = m_cameraInstricMatrix.at<double>(0, 0);
	coeffi_mat.at<double>(2, 1) = m_cameraInstricMatrix.at<double>(1, 1);
	cv::Mat B = cv::Mat(3, 1, CV_64F, cv::Scalar(0));
	B.at<double>(0, 0) = m_vecLaserPlaneCoeffs[2];

	for (int i = 0; i < uv.size(); ++i) {
		coeffi_mat.at<double>(1, 2) = m_cameraInstricMatrix.at<double>(0, 2) - uv_withoutDistortion[i].x;
		coeffi_mat.at<double>(2, 2) = m_cameraInstricMatrix.at<double>(1, 2) - uv_withoutDistortion[i].y;
		double temp = 1 + m_vecLaserPlaneCoeffs[0] * coeffi_mat.at<double>(1, 2) / coeffi_mat.at<double>(1, 0) + m_vecLaserPlaneCoeffs[1] * coeffi_mat.at<double>(2, 2) / coeffi_mat.at<double>(2, 1);
		double zc = m_vecLaserPlaneCoeffs[2] / temp;
		double xc = -coeffi_mat.at<double>(1, 2) / coeffi_mat.at<double>(1, 0) * zc;
		double yc = -coeffi_mat.at<double>(2, 2) / coeffi_mat.at<double>(2, 1) * zc;
		xcyczc_cal.at<double>(2, i) = zc;
		xcyczc_cal.at<double>(0, i) = xc;
		xcyczc_cal.at<double>(1, i) = yc;

		cv::solve(coeffi_mat, B, xcyczc_cal.col(i), cv::DECOMP_SVD);
	//	cout << xcyczc_cal.col(i) << endl;
	//	cout << xcyczc_cal.size() << endl;
	}

	return xcyczc_cal;
}






bool CcalibrateLineCamera::planeFittingUsingVecLinearLS(vector<double>&coefficient, const vector<cv::Point3d>& xcyczc) {
	bool bRet = false;
	cv::Mat A=cv::Mat::ones(xcyczc.size(), 3, CV_64FC1);
	cv::Mat b(xcyczc.size(), 1, CV_64FC1);


	for (int i = 0; i < xcyczc.size(); ++i) {
		A.at<double>(i, 0) = xcyczc[i].x;
		A.at<double>(i, 1) = xcyczc[i].y;
		b.at<double>(i, 0) = xcyczc[i].z;
	}
	/*以下在求解超定线性方程组时，使用最小二乘法+直接求导后取逆，后续可改成SVD或QR，避免求逆矩阵*/
	cv::Mat coef(3, 1, CV_64FC1);
	cv::Mat temp1 = A.t() * A;
	cv::Mat temp2;
	cv::invert(A.t() * A, temp1);
	temp1.convertTo(temp2, CV_64FC1);
	coef = temp2 * A.t() * b;
	for (int i = 0; i < 3; ++i) {
		coefficient.push_back(coef.at<double>(i, 0));
	//	cout << coef.at<double>(i, 0) << endl;
	}

	bRet = true;
	return bRet;
}