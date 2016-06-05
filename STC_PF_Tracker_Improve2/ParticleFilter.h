#pragma once

#include <opencv2/opencv.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Particle.h"
#include "Target.h"
#include<time.h>

using namespace std;
using namespace cv;

# define R_BIN      8  /* 红色分量的直方图条数 */
# define G_BIN      8  /* 绿色分量的直方图条数 */
# define B_BIN      8  /* 兰色分量的直方图条数 */ 

# define R_SHIFT    5  /* 与上述直方图条数对应 */
# define G_SHIFT    5  /* 的R、G、B分量左移位数 */
# define B_SHIFT    5  /* log2( 256/8 )为移动位数 */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
//#define MASK 123459876

class ParticleFilter
{

private:
	const int mParticleNum = 100;			//粒子总数
	const float mUpdateWeightThres = (float)0.9;	//更新模型时的权重阈值
	const float mLegalWeightThresD = (float)0.0001;	//合法权重阈值，用于确定这个最大权重粒子是否可以算作跟踪结果
	const float mGaussianSigma = float(0.6);		//随机的高斯分布的方差
	const float mLearnRate = 0.3f;			//目标模型更新权重取0.1-0.3
	const float mBha2SigmaSquare = 0.02f;	//根据巴氏距离计算权重所用的Sigma的平方的2倍，一般sigma=0.1，所以sigma平方的2倍=0.02;

private:
	long seed;			//生成随机数的种子

	int mHistogramNum;		//直方图条数
	float* mHistogram;		//彩色直方图
	float* mTmpHistogram;	//彩色临时直方图
	Particle* mParticles;	//粒子
	Particle* mTmpParticles;

	int mWidth;			//视频帧宽度
	int mHeight;		//视频帧高度

	gsl_rng* rng;		//高斯分布随机数生成器

	float* mAddWeiArr;	//累计权重数组
	int* mResamIndArr;	//重采样索引数组

public:
	ParticleFilter();
	virtual ~ParticleFilter();

//公有接口
public:
	/*
	* 初始化方法
	* @return  true表示初始化成功；false表示初始化失败
	*/
	bool init(const Mat& mFrame, const Target& mTarget, int mFrameCount);
	//进行每一帧的跟踪
	void track(const Mat& mFrame, Target& mTarget);

	/*
	* 更新模型
	* 即更新彩色直方图
	*/
	void updateModel(const Mat& mMat, const Target& mTarget);

	void updateModelWithTarget(const Mat& mMat, const Target& mTarget);

//私有接口
private:
	/*
	* 获得一个[0,1]之间均匀分布的随机数
	*/
	float random0_1();

	/*
	* 计算一个图像帧中指定区域的彩色直方图分布
	* 返回结果true表示结果合法，false表示结果非法
	*/
	void calcuColorHistogram(const Mat& mMat, const Target& mTarget, float* histogram, int histogramNum);

	/*
	* 样本选择，从N个输入样本中根据权重重新挑选出N个
	*/
	void reselectSample();

	/*
	* 重要性采样
	*/
	void importanceSampling();

	/*
	* 规范化累计权重数组
	*/
	void normalizeAddWeightArr();

	/*
	* 二叉查找，从累计权重数组中找出接近指定权重的最邻近索引
	*/
	int binearySearchAddWeight(float v);

	/*
	* 传播：根据系统状态方程求取状态预测量
	* 状态方程为： S(t) = A S(t-1) + W(t-1)
	* W(t-1)为高斯噪声
	* 根据成员变量中的粒子集作为输入并完成更新
	*/
	void propagate();

	/*
	* 观测：对状态量进行更新
	*/
	void observe(const Mat& mMat);

	/*
	* 计算彩色直方图和临时彩色直方图的巴氏系数
	*/
	float calcuBhattacharyya();

	/*
	* 根据巴氏系数计算此粒子的权重值
	*/
	float calcuWeightedWithBha(float mBha);

	/*
	* 估计：对状态量进行估计，提取位置量
	*/
	void estimation(Particle& mPar);

};