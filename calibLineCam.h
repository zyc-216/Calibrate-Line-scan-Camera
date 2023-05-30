#pragma once
#ifndef _CALIB_CAM_H_
#define _CALIB_CAM_H_

#include <vector>
#include <opencv2/opencv.hpp>
using namespace std;

//�궨ͼƬ�ͱ궨�ļ�·����
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


//HG-450�궨����
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

	/*Calibrate������һ������Ϊ����궨ͼƬ��·����������Ϊ����궨ͼƬ����������������Ϊhalcon���ɵı궨�ļ���
	  ע�⣡������궨��ͼƬ�����ܹ���halcon�����궨�㣬�е�ͼƬ��Ϊ����Ƕ��Լ����ߵ��޷����ɹ������������һ��ɸѡ�ɲ���궨ͼƬ�ĺ���*/
	bool Calibrate(string strImgDir, int imgNum, string strCalTabCpdPath, string strCalTabPsPath, string strOutputCsvOfCamParam);

	bool CalUndistortedCorners(vector<cv::Point2f>& srcPointArray, vector<cv::Point2f>& dstPointArray);  
	bool GetUndistortedCorners(const char* outputFileName) ;							//����궨�ĵ���л������					
	bool GetUndistortedCornersFromCsv(const char* csvName);							//csv�ļ��еĵ���л������
	bool GetUndistortedCornersFromBmp(const char* bmpPath, const char* outputFileName);	//δ����궨�ĵ���ͼƬ�ϵĵ���л������������ΪͼƬ��·��
	bool GetProjectionError(const char* outputFilename);	//����궨�ĵ�ķ�ͶӰ���
	bool PrintCamParam();          //��ӡ����������
	void GetLaserLinePlaneFunc(const char* sArr[5]);//�������ϵĵ�
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
	bool planeFitting(vector<double>& , const vector<cv::Point3d>& );  //z=ax+by+c, coeffientΪa,b,c
	bool quadricSurfaceFitting(vector<double>&, const vector<cv::Point3d>&); //z=ax*2+bxy+cy*2+dx+ey+f, coeffientΪa b c d e f
	bool planeFittingUsingVecLinearLS(vector<double>&,const vector<cv::Point3d>&);
//	bool planeFittingUsingPclRanSac(vector<double>&, const vector<cv::Point3d>&);
	cv::Mat CorrectLaserPlane(const cv::Mat&,  vector<cv::Point3d>,const vector<cv::Point2d>&);  //����һ���������������ϵ�͹�ѧƽ̨����ϵ֮���R����(������ʽ����T����������ʽ��
	bool separateDataSet( vector<cv::Point3d> &xwywzw_vecPoint3d_all, vector<cv::Point2d> &uv_all, vector<cv::Point3d> &xwywzw_vecPoint3d_train, vector<cv::Point2d> &uv_train,
							 vector<cv::Point3d> &xwywzw_vecPoint3d_test, vector<cv::Point2d>& uv_test);
	//��Ҫ���� GetLaserLinePlaneFunc֮��
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
	vector<double>						m_vecLaserPlaneCoeffs;         // z=ax+by+c, coeffientΪa,b,c
	vector<double>						m_vecLaserQuadricSurfaceCoeffs; //z=ax*2+bxy+cy*2+dx+ey+f, coeffientΪa b c d e f
};



#endif // !_CALIB_CAM_H_
