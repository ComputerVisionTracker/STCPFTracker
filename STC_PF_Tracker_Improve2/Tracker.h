#include <opencv2/opencv.hpp>
#include "Target.h"
#include "STCTracker.h"
#include "ParticleFilter.h"
#include "STCResultType.h"
#include "STCResult.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace cv;

class Tracker
{
private:
	double mConfiThreshold;;

private:
	ofstream out;
	STCResult mSTCResult;
	STCTracker mSTCTracker;
	ParticleFilter mParticleFilter;

public:
	Tracker();
	virtual ~Tracker();


public:
	/*
	 * ��ʼ������
	 * @param mFrame ��ɫͼ��֡
	 * @param mGrayFrame �Ҷ�ͼͼ��֡
	 * @param mTarget ����Ŀ��
	 * @param mFrameCount ͼ������֡��
	 * @return true��ʾ��ʼ���ɹ���false��ʾ��ʼ��ʧ��
	 */
	bool init(const Mat& mFrame, const Mat& mGrayFrame, const Target& mTarget, int mFrameCount);

	/*
	 * @param mFrame ��ɫͼ��֡
	 * @param mGrayFrame �Ҷ�ͼͼ��֡
	 * @param mTarget ����Ŀ��
	 */
	void track(const Mat& mFrame, const Mat& mGrayFrame, Target& mTarget);

};