#pragma once

#include <opencv2/opencv.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Target.h"
#include "STCResultType.h"
#include "STCResult.h"

using namespace cv;
using namespace std;

class STCTracker
{
//��Ա����
private:
	
	int mFrameIndex;				//��ǰ����֡������
	int mIntervalNumber;	//ÿ������֡ͼƬ����һ�γ߶Ȳ���
	double sigma;			//������������ʹ�ʽ�еĳ߶Ȳ���������ʽ5
	double alpha;			//������ʣ�������ͼ����ʽ�еĳ߶����Ų���������ʽ6
	double beta;			//������ʣ�������ͼ����ʽ�еĿ���ͼ�μ���ȵĲ���������ʽ6
	double mLearnRate;		//��ʽ12�е�ѧϰ����

	int mWidth;			//��Ƶ֡����
	int mHeight;		//��Ƶ֡�߶�

	
	Size2d mRectSize;		//���پ��ο��һ��ߴ�
	Size mCxtRegionSize;	//����������ȫ�ߴ�
	double mScale;			//�߶����ű���
	double mScaleLearnRate;	//�߶�����ѧϰ����

	int mFrameCount;				//��Ƶ����֡��
	double *mEachFrameConfd;	//ÿһ֡�����Ŷȣ����ڸ��³߶Ȳ���

	//�ֱ�洢���Ƕ�Ӧ�ľֲ��ռ�������������ͼ�����������
	int *mCxtX;
	int *mCxtY;

	Mat mCxtMat;			//�����������Ӧ�����ع淶����ľ��󣬼�����ǿ�ȣ���ʽ4�е�I(z)

	Mat mCxtPriProMat;			//������������ʾ���

	Mat mCxtPostProMat;		//�����ĺ�����ʾ��󣬼�����ͼ
	Mat mCxtPostProFFTMat;	//��������Ҷת����ĺ�����ʾ���

	Mat mFilterWindowMat;	//��ǰ֡���˲����ڣ�����ʽ4�л�δ��I(z)ʱ�ľ������ݣ���Ϊһ���˲�����������

	Mat mSCModelFFT;			//�õ���ÿһ֡�Ŀռ�������ģ�;���FFTת�����
	Mat mSTCModelFFT;			//�ۼ�ȫ��֡��ʱ�ռ�������ģ�;���FFTת�����

	Mat mHammingWinMat;		//����������������ͼ���Ե������Ƶ��Ӱ��
	Mat mEucliDistMat;		//�洢������������ÿһ�����ص㵽���ٶ������ĵ��ŷ����þ���ƽ����������������ʾ���Ҳ���ں�����ʾ���

	Mat mNextCxtPostProMat;		//�洢����һ֡Ԥ��ĺ�����ʣ�������ͼ��

	double mMinConfiPosi;

public:
	STCTracker();
	virtual ~STCTracker();

public:
	//��ʼ������
	bool init(const Mat& mFrame, const Target& mTarget, int mFrameCount, ofstream& out);

	double track(const Mat& mFrame, Target& mTarget, ofstream& out);

	void updateWithSTC(const Mat& mFrame, Target& mTarget, ofstream& out);

	void updateWithParti(const Mat& mFrame, const Target& mTarget, ofstream& out);

	void updateWithPartiByDeny(const Mat& mFrame, const Target& mTarget, ofstream& out);

	bool checkConfidenceLegal(const double confid, ofstream& out, const int mFrameNum, const double mThreshold);

	int getFrameIndex();

private:
	//���캺����
	void createHammingWindow();
	//��ʼ��ŷ����þ������ͺ�����ʾ���
	void createEucliDistMatAndCxtPostPro();
	//�����˷�
	void complexMultiply(const Mat& src1, const Mat& src2, Mat& dst);
	//��������
	void complexDivide(const Mat& src1, const Mat& src2, Mat& dst);

	//�õ������Ķ�Ӧ������ǿ�ȹ淶������
	void getContext(const Mat& mFrame, const Target& mTarget);

	void getFilterWindow();

	//ѧϰ��ǰ֡�Ŀռ�������ģ��
	void learnSCModol();

	//����ʱ�ռ�ģ��
	void updateSTCModol();

};
