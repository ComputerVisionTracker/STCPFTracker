#pragma once

#include <opencv2/opencv.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include "Target.h"

using namespace std;
using namespace cv;

/************************************************************************/
/* 粒子                                                                 */
/************************************************************************/
class Particle{

private:
	static const float POSITION_DISTURB_RATE;
	static const float VELOCITY_DISTURB;
	static const float WINDOW_DISTURB;
	static const float SCALE_DISTURB;

public:
	int x;					/* x坐标位置 */
	int y;					/* x坐标位置 */
	float mXVelocity;		/* x方向运动速度 */
	float mYVelocity;		/* y方向运动速度 */
	int mHalfWidth;			/* x方向半窗宽 */
	int mHalfHeight;		/* y方向半窗宽 */
	float mScaleVelocity;	/* 尺度变换速度 */

	//单独更新
	float weight;			//此粒子权重

	//构造函数和析构函数
public:
	Particle();
	Particle(const Particle& r);
	virtual ~Particle();

	//重载运算符
public:
	Particle& operator=(const Particle &r);

	//公有接口
public:
	void setWeight(float weight);
	//根据指定的高斯随机数生成器更新当前粒子状态
	void update(const gsl_rng* rng, const double sigma);
	//给定图像帧，绘制粒子
	void drawParticle(Mat& mMat) const;

	/*
	* 初始化成员变量，两种初始化方法
	*/
	void init(const Target& mTarget, const float weight);
	void init(const Target& mTarget, const gsl_rng* rng, const double sigma, const float weight);

};