#include "ParticleFilter.h"

ParticleFilter::ParticleFilter()
{

	seed = (long)time(NULL);	//ʹ��ϵͳʱ����Ϊ����

	mParticles = new Particle[mParticleNum];
	mTmpParticles = new Particle[mParticleNum];

	mHistogramNum = R_BIN * G_BIN * B_BIN;
	mHistogram = new float[mHistogramNum];
	mTmpHistogram = new float[mHistogramNum];

	//��ʼ��GSL��ѧ����⻷�������ڻ�ø�˹�ֲ������
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, seed);

	mAddWeiArr = new float[mParticleNum + 1];
	mResamIndArr = new int[mParticleNum];
}


ParticleFilter::~ParticleFilter()
{
	delete[] mParticles;
	delete[] mTmpParticles;
	delete[] mHistogram;
	delete[] mTmpHistogram;

	gsl_rng_free(rng);

	delete[] mAddWeiArr;
	delete[] mResamIndArr;
}

bool ParticleFilter::init(const Mat& mFrame, const Target& mTarget, int mFrameCount)
{
	mWidth = mFrame.cols;
	mHeight = mFrame.rows;

	if (mWidth == 0 || mHeight == 0)
	{
		return false;
	}

	calcuColorHistogram(mFrame, mTarget, mHistogram, mHistogramNum);
	float weight = (float)(1.0 / mParticleNum);
	mParticles[0].init(mTarget, weight);
	for (int i = 1; i < mParticleNum; i++)
	{
		mParticles[i].init(mTarget, rng, mGaussianSigma, weight);
	}
	return true;
}

void ParticleFilter::track(const Mat& mFrame, Target& mTarget)
{
	Particle mParticle;
	reselectSample();
	//����
	propagate();
	//�۲����
	observe(mFrame);
	//����
	estimation(mParticle);
	mTarget.set(mParticle);
	updateModel(mFrame, mTarget);
}

void ParticleFilter::updateModel(const Mat& mMat, const Target& mTarget)
{
	calcuColorHistogram(mMat, mTarget, mTmpHistogram, mHistogramNum);
	float mBha = calcuBhattacharyya();
	float mWeight = calcuWeightedWithBha(mBha);
	if (mWeight > mUpdateWeightThres)
	{
		for (int i = 0; i < mHistogramNum; i++)
		{
			mHistogram[i] = (float)((1 - mLearnRate)*mHistogram[i] + mLearnRate*mTmpHistogram[i]);
		}
	}

}

void ParticleFilter::updateModelWithTarget(const Mat& mMat, const Target& mTarget)
{
	calcuColorHistogram(mMat, mTarget, mTmpHistogram, mHistogramNum);
	float mBha = calcuBhattacharyya();
	float mWeight = calcuWeightedWithBha(mBha);
	if (mWeight > mUpdateWeightThres)
	{
		for (int i = 0; i < mHistogramNum; i++)
		{
			mHistogram[i] = (float)((1 - mLearnRate)*mHistogram[i] + mLearnRate*mTmpHistogram[i]);
		}
	}

	reselectSample();
	//����
	propagate();
	//�۲����
	observe(mMat);
}

/*
����Park and Miller��������[0,1]֮����ȷֲ���α�����
�㷨��ϸ��������
[1] NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING.Cambridge University Press. 1992. pp.278-279.
[2] Park, S.K., and Miller, K.W. 1988, Communications of the ACM,vol. 31, pp. 1192�C1201.
*/
float ParticleFilter::random0_1()
{
	long k;
	float ans;

	/* *idum ^= MASK;*/      /* XORing with MASK allows use of zero and other */
	k = seed / IQ;            /* simple bit patterns for idum.                 */
	seed = IA*(seed - k*IQ) - IR*k;  /* Compute idum=(IA*idum) % IM without over- */
	if (seed < 0) seed += IM;  /* flows by Schrage��s method.               */
	ans = (float)(AM*seed);          /* Convert idum to a floating result.            */
	/* *idum ^= MASK;*/      /* Unmask before return.                         */
	return ans;
}

void ParticleFilter::calcuColorHistogram(const Mat& mMat, const Target& mTarget, float* histogram, int histogramNum)
{
	for (int i = 0; i < histogramNum; i++)
	{
		histogram[i] = 0.0;
	}
	//���Ŀ�����򲻺Ϸ������������أ�ʹ��ɫֱ��ͼΪ0
	if (!mTarget.isValid(mWidth, mHeight))
	{
		return;
	}

	int x_min = mTarget.x - mTarget.wid, x_max = mTarget.x + mTarget.wid;
	int y_min = mTarget.y - mTarget.hei, y_max = mTarget.y + mTarget.hei;
	if (x_min < 0)
	{
		x_min = 0;
	}
	if (x_max >= mWidth)
	{
		x_max = mWidth - 1;
	}
	if (y_min < 0)
	{
		y_min = 0;
	}
	if (y_max >= mHeight)
	{
		y_max = mHeight - 1;
	}
	int index;		//����ֵ����ǵ�ǰ������ص�������ɫֱ��ͼ������
	int r, g, b;	//��ǰ���ص��RGBֵ
	float mKernel, mTotal = 0;	//�ֱ�Ϊ��ǰ���ص�˺���ֵ���ܹ��ĺ˺���ֵ(���ھ�һ��)�����ƽ��ֵ
	const int mSquare = mTarget.wid*mTarget.wid + mTarget.hei*mTarget.hei;

	uchar* temp;
	for (int y = y_min; y < y_max; y++)
	{
		temp = mMat.data + y * mMat.step + x_min * 3;
		for (int x = x_min; x < x_max; x++, temp++)
		{
			r = (*temp) >> R_SHIFT;
			temp++;
			g = (*temp) >> G_SHIFT;
			temp++;
			b = (*temp) >> B_SHIFT;
			index = r * G_BIN * B_BIN + g * B_BIN + b;
			mKernel = 1 - (float)(((y - mTarget.y)*(y - mTarget.y) + (x - mTarget.x)*(x - mTarget.x))*1.0 / mSquare); /* ����뾶ƽ��r^2  �˺���k(r) = 1-r^2, |r| < 1; ����ֵ k(r) = 0 */
			mTotal += mKernel;
			histogram[index] += mKernel;
		}
	}

	for (int i = 0; i < histogramNum; i++)
	{
		histogram[i] /= mTotal;
	}

}

void ParticleFilter::reselectSample()
{
	importanceSampling();
	for (int i = 0; i < mParticleNum; i++)
	{
		mTmpParticles[i] = mParticles[mResamIndArr[i]];
	}
	memcpy(mParticles, mTmpParticles, sizeof(Particle)*mParticleNum);
}

void ParticleFilter::importanceSampling()
{
	normalizeAddWeightArr();
	int x;
	for (int i = 0; i < mParticleNum; i++)
	{
		x = binearySearchAddWeight(random0_1());
		if (x == mParticleNum)
		{
			x--;
		}
		mResamIndArr[i] = x;
	}
}

void ParticleFilter::normalizeAddWeightArr()
{
	mAddWeiArr[0] = 0;
	for (int i = 0; i < mParticleNum; i++)
	{
		mAddWeiArr[i + 1] = mParticles[i].weight + mAddWeiArr[i];
	}
	for (int i = 0; i <= mParticleNum; i++)
	{
		mAddWeiArr[i] /= mAddWeiArr[mParticleNum];
	}
}

int ParticleFilter::binearySearchAddWeight(float v)
{
	int l, r, m;
	l = 0; 	r = mParticleNum;
	while (r >= l)
	{
		m = (l + r) / 2;
		if (v >= mAddWeiArr[m] && v < mAddWeiArr[m + 1])
			return(m);
		if (v < mAddWeiArr[m])
			r = m - 1;
		else
			l = m + 1;
	}
	return 0;
}


void ParticleFilter::propagate()
{
	for (int i = 0; i < mParticleNum; i++)
	{
		mParticles[i].update(rng, mGaussianSigma);
	}
}

void ParticleFilter::observe(const Mat& mMat)
{
	Target* mTarget = new Target();
	float mBha = 0;
	for (int i = 0; i < mParticleNum; i++)
	{
		(*mTarget).set(mParticles[i]);
		calcuColorHistogram(mMat, *mTarget, mTmpHistogram, mHistogramNum);
		mBha = calcuBhattacharyya();
		mParticles[i].setWeight(calcuWeightedWithBha(mBha));

	}
	delete mTarget;
}

float ParticleFilter::calcuBhattacharyya()
{
	float rho = 0;
	for (int i = 0; i < mHistogramNum; i++)
	{
		rho = rho + sqrt(mHistogram[i] * mTmpHistogram[i]);
	}
	return rho;
}


float ParticleFilter::calcuWeightedWithBha(float mBha)
{
	return (float)(exp((mBha - 1) / mBha2SigmaSquare));
}

void ParticleFilter::estimation(Particle& mPar)
{
	float x = 0, y = 0, xv = 0, yv = 0, wid = 0, hei = 0, sca = 0, weightSum = 0;
	for (int i = 0; i < mParticleNum; i++)
	{
		x += mParticles[i].x*mParticles[i].weight;
		y += mParticles[i].y*mParticles[i].weight;
		xv += mParticles[i].mXVelocity*mParticles[i].weight;
		yv += mParticles[i].mYVelocity*mParticles[i].weight;
		wid += mParticles[i].mHalfWidth*mParticles[i].weight;
		hei += mParticles[i].mHalfHeight*mParticles[i].weight;
		sca += mParticles[i].mScaleVelocity*mParticles[i].weight;
		weightSum += mParticles[i].weight;
	}
	if (weightSum <= 0)
	{
		weightSum = 1;
	}
	mPar.x = (int)(x / weightSum + 0.5);
	mPar.y = (int)(y / weightSum + 0.5);
	mPar.mXVelocity = xv / weightSum;
	mPar.mYVelocity = yv / weightSum;
	mPar.mHalfWidth = (int)(wid / weightSum + 0.5);
	mPar.mHalfHeight = (int)(hei / weightSum + 0.5);
	mPar.mScaleVelocity = sca / weightSum;
}