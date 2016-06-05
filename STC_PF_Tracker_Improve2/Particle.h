#pragma once

#include <opencv2/opencv.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include "Target.h"

using namespace std;
using namespace cv;

/************************************************************************/
/* ����                                                                 */
/************************************************************************/
class Particle{

private:
	static const float POSITION_DISTURB_RATE;
	static const float VELOCITY_DISTURB;
	static const float WINDOW_DISTURB;
	static const float SCALE_DISTURB;

public:
	int x;					/* x����λ�� */
	int y;					/* x����λ�� */
	float mXVelocity;		/* x�����˶��ٶ� */
	float mYVelocity;		/* y�����˶��ٶ� */
	int mHalfWidth;			/* x����봰�� */
	int mHalfHeight;		/* y����봰�� */
	float mScaleVelocity;	/* �߶ȱ任�ٶ� */

	//��������
	float weight;			//������Ȩ��

	//���캯������������
public:
	Particle();
	Particle(const Particle& r);
	virtual ~Particle();

	//���������
public:
	Particle& operator=(const Particle &r);

	//���нӿ�
public:
	void setWeight(float weight);
	//����ָ���ĸ�˹��������������µ�ǰ����״̬
	void update(const gsl_rng* rng, const double sigma);
	//����ͼ��֡����������
	void drawParticle(Mat& mMat) const;

	/*
	* ��ʼ����Ա���������ֳ�ʼ������
	*/
	void init(const Target& mTarget, const float weight);
	void init(const Target& mTarget, const gsl_rng* rng, const double sigma, const float weight);

};