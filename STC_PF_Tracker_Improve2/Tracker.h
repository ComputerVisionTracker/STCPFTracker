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
	 * 初始化方法
	 * @param mFrame 彩色图像帧
	 * @param mGrayFrame 灰度图图像帧
	 * @param mTarget 跟踪目标
	 * @param mFrameCount 图像序列帧数
	 * @return true表示初始化成功，false表示初始化失败
	 */
	bool init(const Mat& mFrame, const Mat& mGrayFrame, const Target& mTarget, int mFrameCount);

	/*
	 * @param mFrame 彩色图像帧
	 * @param mGrayFrame 灰度图图像帧
	 * @param mTarget 跟踪目标
	 */
	void track(const Mat& mFrame, const Mat& mGrayFrame, Target& mTarget);

};