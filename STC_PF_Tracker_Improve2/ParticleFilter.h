#pragma once

#include <opencv2/opencv.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Particle.h"
#include "Target.h"
#include<time.h>

using namespace std;
using namespace cv;

# define R_BIN      8  /* ��ɫ������ֱ��ͼ���� */
# define G_BIN      8  /* ��ɫ������ֱ��ͼ���� */
# define B_BIN      8  /* ��ɫ������ֱ��ͼ���� */ 

# define R_SHIFT    5  /* ������ֱ��ͼ������Ӧ */
# define G_SHIFT    5  /* ��R��G��B��������λ�� */
# define B_SHIFT    5  /* log2( 256/8 )Ϊ�ƶ�λ�� */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
//#define MASK 123459876

class ParticleFilter
{

private:
	const int mParticleNum = 100;			//��������
	const float mUpdateWeightThres = (float)0.9;	//����ģ��ʱ��Ȩ����ֵ
	const float mLegalWeightThresD = (float)0.0001;	//�Ϸ�Ȩ����ֵ������ȷ��������Ȩ�������Ƿ�����������ٽ��
	const float mGaussianSigma = float(0.6);		//����ĸ�˹�ֲ��ķ���
	const float mLearnRate = 0.3f;			//Ŀ��ģ�͸���Ȩ��ȡ0.1-0.3
	const float mBha2SigmaSquare = 0.02f;	//���ݰ��Ͼ������Ȩ�����õ�Sigma��ƽ����2����һ��sigma=0.1������sigmaƽ����2��=0.02;

private:
	long seed;			//���������������

	int mHistogramNum;		//ֱ��ͼ����
	float* mHistogram;		//��ɫֱ��ͼ
	float* mTmpHistogram;	//��ɫ��ʱֱ��ͼ
	Particle* mParticles;	//����
	Particle* mTmpParticles;

	int mWidth;			//��Ƶ֡���
	int mHeight;		//��Ƶ֡�߶�

	gsl_rng* rng;		//��˹�ֲ������������

	float* mAddWeiArr;	//�ۼ�Ȩ������
	int* mResamIndArr;	//�ز�����������

public:
	ParticleFilter();
	virtual ~ParticleFilter();

//���нӿ�
public:
	/*
	* ��ʼ������
	* @return  true��ʾ��ʼ���ɹ���false��ʾ��ʼ��ʧ��
	*/
	bool init(const Mat& mFrame, const Target& mTarget, int mFrameCount);
	//����ÿһ֡�ĸ���
	void track(const Mat& mFrame, Target& mTarget);

	/*
	* ����ģ��
	* �����²�ɫֱ��ͼ
	*/
	void updateModel(const Mat& mMat, const Target& mTarget);

	void updateModelWithTarget(const Mat& mMat, const Target& mTarget);

//˽�нӿ�
private:
	/*
	* ���һ��[0,1]֮����ȷֲ��������
	*/
	float random0_1();

	/*
	* ����һ��ͼ��֡��ָ������Ĳ�ɫֱ��ͼ�ֲ�
	* ���ؽ��true��ʾ����Ϸ���false��ʾ����Ƿ�
	*/
	void calcuColorHistogram(const Mat& mMat, const Target& mTarget, float* histogram, int histogramNum);

	/*
	* ����ѡ�񣬴�N�����������и���Ȩ��������ѡ��N��
	*/
	void reselectSample();

	/*
	* ��Ҫ�Բ���
	*/
	void importanceSampling();

	/*
	* �淶���ۼ�Ȩ������
	*/
	void normalizeAddWeightArr();

	/*
	* ������ң����ۼ�Ȩ���������ҳ��ӽ�ָ��Ȩ�ص����ڽ�����
	*/
	int binearySearchAddWeight(float v);

	/*
	* ����������ϵͳ״̬������ȡ״̬Ԥ����
	* ״̬����Ϊ�� S(t) = A S(t-1) + W(t-1)
	* W(t-1)Ϊ��˹����
	* ���ݳ�Ա�����е����Ӽ���Ϊ���벢��ɸ���
	*/
	void propagate();

	/*
	* �۲⣺��״̬�����и���
	*/
	void observe(const Mat& mMat);

	/*
	* �����ɫֱ��ͼ����ʱ��ɫֱ��ͼ�İ���ϵ��
	*/
	float calcuBhattacharyya();

	/*
	* ���ݰ���ϵ����������ӵ�Ȩ��ֵ
	*/
	float calcuWeightedWithBha(float mBha);

	/*
	* ���ƣ���״̬�����й��ƣ���ȡλ����
	*/
	void estimation(Particle& mPar);

};