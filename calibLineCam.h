#pragma once
#ifndef _CALIB_CAM_H_
#define _CALIB_CAM_H_

#include <vector>
#include <opencv2/opencv.hpp>
using namespace std;

//标定图片和标定文件路径类
class CPathParam
{
public:
	CPathParam() = default;
	bool SetParamValue(string imgDir,int imgNum, string calTabCpdPath,string calTabPsPath);

public:
	string				m_strImgDir;
	int					m_iImgNum;
	string				m_strCalTabCpdPath;
	string				m_strCalTabPsPath;
};


//HG-450标定板类
class CCalibObject 
{
public:
	CCalibObject();
public:
	int				    m_iCalTabRow;
	int					m_iCalTabCol;
	double			    m_cornerDiameter;
	double			    m_XofWorldPt[195];
	double				m_YofWorldPt[195];		
	int					m_iMarkRow[5];
	int					m_iMarkCol[5];
	double				m_startCamPar[8];
	int					m_iImgWidth;
	int					m_iImgHeight;
	string			    m_strPolarity;
};



class CcalibrateLineCamera
{
public:
	CcalibrateLineCamera();
	CcalibrateLineCamera(const char* inputCamParamCsv);

	/*Calibrate函数第一个参数为参与标定图片的路径，参数二为参与标定图片数量，后两个参数为halcon生成的标定文件。
	  注意！！参与标定的图片必须能够被halcon检测出标定点，有的图片因为拍摄角度以及光线等无法检测成功，后续可添加一个筛选可参与标定图片的函数*/
	bool Calibrate(string strImgDir, int imgNum, string strCalTabCpdPath, string strCalTabPsPath, string strOutputCsvOfCamParam);

	bool CalUndistortedCorners(vector<cv::Point2f>& srcPointArray, vector<cv::Point2f>& dstPointArray);  
	bool GetUndistortedCorners(const char* outputFileName) ;							//参与标定的点进行畸变矫正					
	bool GetUndistortedCornersFromCsv(const char* csvName);							//csv文件中的点进行畸变矫正
	bool GetUndistortedCornersFromBmp(const char* bmpPath, const char* outputFileName);	//未参与标定的单张图片上的点进行畸变矫正，参数为图片的路径
	bool GetProjectionError(const char* outputFilename);	//参与标定的点的反投影误差
	bool PrintCamParam();          //打印输出相机参数
	void GetLaserLinePlaneFunc(const char* sArr[5]);//激光线上的点
	//bool TestAcurracyUsingTrueData(const char*, vector<double>&);
	//bool SetCalObjParam(CCalibObject* pCalObj);
	//bool Clone(CParammeter* pParam);
	//bool GetUndistortedCornersToCsv(vector<cv::Point2f>& srcPointArray,const char* csvName);	
	//bool GetLaserLinePlaneFunc(const char *sArr[5]); 


protected:
	bool SetParam(CPathParam* pParam);
	bool fnMatchCornersAndWorldPoints();
	bool CalibrateUsingCV(const char* outputCamParamCsv);
	bool planeFitting(const vector<cv::Point3d>&);
	bool planeFitting(vector<double>& , const vector<cv::Point3d>& );  //z=ax+by+c, coeffient为a,b,c
	bool quadricSurfaceFitting(vector<double>&, const vector<cv::Point3d>&); //z=ax*2+bxy+cy*2+dx+ey+f, coeffient为a b c d e f
	bool planeFittingUsingVecLinearLS(vector<double>&,const vector<cv::Point3d>&);
//	bool planeFittingUsingPclRanSac(vector<double>&, const vector<cv::Point3d>&);
	cv::Mat CorrectLaserPlane(const cv::Mat&,  vector<cv::Point3d>,const vector<cv::Point2d>&);  //返回一个真正的相机坐标系和光学平台坐标系之间的R矩阵(矩阵形式）和T矩阵（向量形式）
	bool separateDataSet( vector<cv::Point3d> &xwywzw_vecPoint3d_all, vector<cv::Point2d> &uv_all, vector<cv::Point3d> &xwywzw_vecPoint3d_train, vector<cv::Point2d> &uv_train,
							 vector<cv::Point3d> &xwywzw_vecPoint3d_test, vector<cv::Point2d>& uv_test);
	//需要用在 GetLaserLinePlaneFunc之后
	cv::Mat getXcyczcFromUv(vector<cv::Point2d> uv, cv::Mat cameraInstricMattrix);

private:
	cv::Mat								m_cameraInstricMatrix;      // = cv::Mat::zeros(3, 3, CV_64F);
	cv::Mat								m_distCoeffs;               // = cv::Mat::zeros(14, 1, CV_64F);
	vector<cv::Mat>						m_rVecs;
	vector<cv::Mat>						m_tVecs;
	vector<cv::Point2f>					m_cornersPerImg;
	vector<vector<cv::Point2f>>			m_cornersAllImg;
	vector<cv::Point3f>					m_worldPointsPerImg;
	vector<vector<cv::Point3f>>			m_worldPointsAllImg;
	vector<cv::Point2f>					m_undistortCornersPerImg;
	CCalibObject						m_calibObject;
	string								m_strImgDir;
	int									m_iImgNum;
	string								m_strCalTabCpdPath, m_strCalTabPsPath;
	vector<double>						m_vecLaserPlaneCoeffs;         // z=ax+by+c, coeffient为a,b,c
	vector<double>						m_vecLaserQuadricSurfaceCoeffs; //z=ax*2+bxy+cy*2+dx+ey+f, coeffient为a b c d e f
};



#endif // !_CALIB_CAM_H_
