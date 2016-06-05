#include "STCTracker.h"


STCTracker::STCTracker()
{
	
	mFrameIndex = 0;
	mIntervalNumber = 5;
	alpha = 2.25;
	beta = 1;
	mLearnRate = 0.075;
	mScale = 1.0;
	mScaleLearnRate = 0.25;

	mEachFrameConfd = NULL;
	mCxtX = NULL;
	mCxtY = NULL;


	mMinConfiPosi = 1;
}


STCTracker::~STCTracker()
{
	if (mEachFrameConfd != NULL)
	{
		delete[] mEachFrameConfd;
	}
	if (mCxtX != NULL)
	{
		delete[] mCxtX;
	}
	if (mCxtY != NULL)
	{
		delete[] mCxtY;
	}
}


void STCTracker::createHammingWindow()
{
	int mCols = mHammingWinMat.cols;
	int mRows = mHammingWinMat.rows;

	double *mHammingValue = new double[mRows];
	double *mHanningValue = new double[mCols];

	for (int i = 0; i < mRows; i++)
	{
		mHammingValue[i] = 0.54 - 0.46*cos(2 * CV_PI*i / (mRows - 1));
	}

	for (int i = 0; i < mCols; i++)
	{
		mHanningValue[i] = 0.5 - 0.5*cos(2 * CV_PI*i / (mCols - 1));
	}

	double *mData = NULL;

	for (int i = 0; i < mRows; i++)
	{
		mData = (double *)(mHammingWinMat.data + i * mHammingWinMat.step);
		for (int j = 0; j < mCols; j++)
		{
			mData[j] = mHanningValue[j] * mHammingValue[i];
		}

	}

	delete[] mHammingValue;
	delete[] mHanningValue;
}

void STCTracker::createEucliDistMatAndCxtPostPro()
{
	double *mXArray, *mYArray, *mXTempPtr, *mYTempPtr;

	mXArray = new double[mCxtRegionSize.width];
	mYArray = new double[mCxtRegionSize.height];

	double mXMiddle = floor((mCxtRegionSize.width - 1 + 0.000001) / 2);
	double mYMiddle = floor((mCxtRegionSize.height - 1 + 0.000001) / 2);

	mXTempPtr = mXArray;
	for (int i = 0; i < mCxtRegionSize.width; i++)
	{
		*mXTempPtr = (i - mXMiddle);
		mXTempPtr++;
	}

	mYTempPtr = mYArray;
	for (int i = 0; i < mCxtRegionSize.height; i++)
	{
		*mYTempPtr = (i - mYMiddle);
		mYTempPtr++;
	}
	double *mEucData = NULL;
	double *mCxtPostProData = NULL;
	double mPostSum = 0;
	for (int i = 0; i < mCxtRegionSize.height; i++)
	{
		mEucData = (double *)(mEucliDistMat.data + mEucliDistMat.step * i);
		mCxtPostProData = (double *)(mCxtPostProMat.data + mCxtPostProMat.step * i);

		mXTempPtr = mXArray;
		mYTempPtr = mYArray + i;
		for (int j = 0; j < mCxtRegionSize.width; j++)
		{
			mEucData[j] = *mYTempPtr * *mYTempPtr + *mXTempPtr * *mXTempPtr;
			mXTempPtr++;

			mCxtPostProData[j] = exp(-pow(0.5 * sqrt(mEucData[j]) / alpha, beta));
			mPostSum += mCxtPostProData[j];
		}
	}
	mCxtPostProMat.convertTo(mCxtPostProMat, -1, 1.0 / mPostSum);

	//Ԥ�ȼ����������ʣ������Ŷ�
	Mat planes1[] = { mCxtPostProMat, Mat::zeros(mCxtPostProMat.size(), CV_64FC1) };
	merge(planes1, 2, mCxtPostProFFTMat);
	dft(mCxtPostProFFTMat, mCxtPostProFFTMat);

	delete[] mXArray;
	delete[] mYArray;
}

void STCTracker::complexMultiply(const Mat& src1, const Mat& src2, Mat& dst)
{
	Mat *mAReal, *mAImag, *mBReal, *mBImag, *mRReal, *mRImag;

	vector<Mat> mAPlanes;
	split(src1, mAPlanes);
	mAReal = &mAPlanes[0];
	mAImag = &mAPlanes[1];

	vector<Mat> mBPlanes;
	split(src2, mBPlanes);
	mBReal = &mBPlanes[0];
	mBImag = &mBPlanes[1];

	dst.create(src1.rows, src1.cols, CV_64FC2);
	vector<Mat> mRPlanes;
	split(dst, mRPlanes);
	mRReal = &mRPlanes[0];
	mRImag = &mRPlanes[1];

	double *mARealData, *mAImagData, *mBRealData, *mBImagData, *mRRealData, *mRImagData;

	for (int i = 0; i < src1.rows; i++)
	{
		mARealData = (double *)(mAReal->data + mAReal->step * i);
		mAImagData = (double *)(mAImag->data + mAImag->step * i);
		mBRealData = (double *)(mBReal->data + mBReal->step * i);
		mBImagData = (double *)(mBImag->data + mBImag->step * i);
		mRRealData = (double *)(mRReal->data + mRReal->step * i);
		mRImagData = (double *)(mRImag->data + mRImag->step * i);

		for (int j = 0; j < src1.cols; j++)
		{
			mRRealData[j] = mARealData[j] * mBRealData[j] - mAImagData[j] * mBImagData[j];
			mRImagData[j] = mARealData[j] * mBImagData[j] + mAImagData[j] * mBRealData[j];
		}

	}
	merge(mRPlanes, dst);
}


void STCTracker::complexDivide(const Mat& src1, const Mat& src2, Mat& dst)
{
	Mat *mAReal, *mAImag, *mBReal, *mBImag, *mRReal, *mRImag;

	vector<Mat> mAPlanes;
	split(src1, mAPlanes);
	mAReal = &mAPlanes[0];
	mAImag = &mAPlanes[1];

	vector<Mat> mBPlanes;
	split(src2, mBPlanes);
	mBReal = &mBPlanes[0];
	mBImag = &mBPlanes[1];

	dst.create(src1.rows, src1.cols, CV_64FC2);
	vector<Mat> mRPlanes;
	split(dst, mRPlanes);
	mRReal = &mRPlanes[0];
	mRImag = &mRPlanes[1];

	double *mARealData, *mAImagData, *mBRealData, *mBImagData, *mRRealData, *mRImagData;

	for (int i = 0; i < src1.rows; i++)
	{
		mARealData = (double *)(mAReal->data + mAReal->step * i);
		mAImagData = (double *)(mAImag->data + mAImag->step * i);
		mBRealData = (double *)(mBReal->data + mBReal->step * i);
		mBImagData = (double *)(mBImag->data + mBImag->step * i);
		mRRealData = (double *)(mRReal->data + mRReal->step * i);
		mRImagData = (double *)(mRImag->data + mRImag->step * i);

		for (int j = 0; j < src1.cols; j++)
		{
			mRRealData[j] = (mARealData[j] * mBRealData[j] + mAImagData[j] * mBImagData[j]) / (mBRealData[j] * mBRealData[j] + mBImagData[j] * mBImagData[j]);
			mRImagData[j] = (mAImagData[j] * mBRealData[j] - mARealData[j] * mBImagData[j]) / (mBRealData[j] * mBRealData[j] + mBImagData[j] * mBImagData[j]);
		}

	}
	merge(mRPlanes, dst);
}

/*
* ��Ϊsigmaÿ֡���ڱ仯�������˲���ÿһ֡�������
*/
void STCTracker::getFilterWindow()
{
	double *mFilterData = NULL, *mEucliData = NULL, *mHammData = NULL, mSum = 0;

	for (int i = 0; i < mCxtRegionSize.height; i++)
	{
		mFilterData = (double *)(mFilterWindowMat.data + mFilterWindowMat.step * i);
		mEucliData = (double *)(mEucliDistMat.data + mEucliDistMat.step * i);
		mHammData = (double *)(mHammingWinMat.data + mHammingWinMat.step * i);
		for (int j = 0; j < mCxtRegionSize.width; j++)
		{
			mFilterData[j] = mHammData[j] * exp(-mEucliData[j] / (2 * sigma * sigma));
			mSum += mFilterData[j];
		}

	}
	mFilterWindowMat.convertTo(mFilterWindowMat, -1, 1 / mSum);

}

/************************************************************************/
/*	������Ҫ���Ǳ߽����⣬�������Ҫ�ü��ľֲ��ռ������Ķ�Ӧ��������ֵ��Ȼ��ü����淶��
/************************************************************************/
void STCTracker::getContext(const Mat& mFrame, const Target& mTarget)
{
	int *temp = NULL;
	temp = mCxtX;
	for (int i = 0; i < mCxtRegionSize.width; i++)
	{
		*temp = (int)floor(mTarget.x + (i - mCxtRegionSize.width / 2.0));
		if (*temp < 0)
		{
			*temp = 0;
		}
		if (*temp >= mFrame.cols)
		{
			*temp = mFrame.cols - 1;
		}
		temp++;
	}

	temp = mCxtY;
	for (int i = 0; i < mCxtRegionSize.height; i++)
	{
		*temp = (int)floor(mTarget.y + (i - mCxtRegionSize.height / 2.0));
		if (*temp < 0)
		{
			*temp = 0;
		}
		if (*temp >= mFrame.rows)
		{
			*temp = mFrame.rows - 1;
		}
		temp++;
	}

	double *mData = NULL;
	uchar *mFrameData = NULL;
	double mPxSum = 0, mPxAverage = 0;

	for (int i = 0; i < mCxtRegionSize.height; i++)
	{

		mData = (double *)(mCxtMat.data + mCxtMat.step * i);
		mFrameData = (uchar *)(mFrame.data + mFrame.step * mCxtY[i]);
		for (int j = 0; j < mCxtRegionSize.width; j++)
		{
			mData[j] = mFrameData[mCxtX[j]];
			mPxSum += mData[j];
		}

	}
	mPxAverage = mPxSum / (mCxtRegionSize.width * mCxtRegionSize.height);

	mCxtMat.convertTo(mCxtMat, -1, 1, -mPxAverage);

	mCxtPriProMat = mFilterWindowMat.mul(mCxtMat);

}


void STCTracker::learnSCModol()
{
	//��������ʽ���FFTת��
	Mat priorFourier;
	Mat planes1[] = { mCxtPriProMat, Mat::zeros(mCxtPriProMat.size(), CV_64FC1) };
	merge(planes1, 2, priorFourier);
	dft(priorFourier, priorFourier);

	complexDivide(mCxtPostProFFTMat, priorFourier, mSCModelFFT);

	//��FFTת���õ���ǰ֡�Ŀռ�������ģ��
	//dft(conditionalFourier, mSCModelFFT, DFT_INVERSE | DFT_REAL_OUTPUT | DFT_SCALE);

}

void STCTracker::updateSTCModol()
{
	addWeighted(mSTCModelFFT, 1.0 - mLearnRate, mSCModelFFT, mLearnRate, 0.0, mSTCModelFFT);
}



/************************************************************************/
/* ���ȳ�ʼ��һ���������������ڴ洢ÿһ֡��Ӧ�����Ŷȣ����ڴ���߶�����
/************************************************************************/
bool STCTracker::init(const Mat& mFrame, const Target& mTarget, int mFrameCount, ofstream& out)
{
	out << "��ʼ���ٿ�";
	out << mTarget << endl;
	mWidth = mFrame.cols;
	mHeight = mFrame.rows;
	if (mWidth == 0 || mHeight == 0)
	{
		return false;
	}

	mFrameIndex = 0;

	if (mFrameCount <= 0)
	{
		this->mFrameCount = 1000;
	}
	else
	{
		this->mFrameCount = mFrameCount + 5;
	}

	mEachFrameConfd = new double[this->mFrameCount];
	mEachFrameConfd[0] = 0.000001;
	for (int i = 1; i < this->mFrameCount; i++)
	{
		mEachFrameConfd[i] = 0;
	}
	//memset(mEachFrameConfd, 0, this->mFrameCount * sizeof(double));

	sigma = 0.5 * (mTarget.wid + mTarget.hei) * 2;

	mRectSize.width = mTarget.wid;
	mRectSize.height = mTarget.hei;

	mCxtRegionSize.width = mTarget.wid * 4;
	mCxtRegionSize.height = mTarget.hei * 4;

	if (mCxtRegionSize.width > mWidth)
	{
		mCxtRegionSize.width = mWidth;
	}

	if (mCxtRegionSize.height > mHeight)
	{
		mCxtRegionSize.height = mHeight;
	}


	mCxtX = new int[mCxtRegionSize.width];
	mCxtY = new int[mCxtRegionSize.height];

	mCxtPriProMat = Mat::zeros(mCxtRegionSize.height, mCxtRegionSize.width, CV_64FC1);
	mCxtPostProMat = Mat::zeros(mCxtRegionSize.height, mCxtRegionSize.width, CV_64FC1);
	mSCModelFFT = Mat::zeros(mCxtRegionSize.height, mCxtRegionSize.width, CV_64FC2);
	mSTCModelFFT = Mat::zeros(mCxtRegionSize.height, mCxtRegionSize.width, CV_64FC2);
	mFilterWindowMat = Mat::zeros(mCxtRegionSize.height, mCxtRegionSize.width, CV_64FC1);
	mEucliDistMat = Mat::zeros(mCxtRegionSize.height, mCxtRegionSize.width, CV_64FC1);
	createEucliDistMatAndCxtPostPro();

	mHammingWinMat = Mat::zeros(mCxtRegionSize.height, mCxtRegionSize.width, CV_64FC1);
	createHammingWindow();

	mCxtMat = Mat::zeros(mCxtRegionSize.height, mCxtRegionSize.width, CV_64FC1);

	getFilterWindow();
	getContext(mFrame, mTarget);
	learnSCModol();
	mSCModelFFT.copyTo(mSTCModelFFT);

	out << endl;
	return true;
}


double STCTracker::track(const Mat& mFrame, Target& mTarget, ofstream& out)
{
	mFrameIndex++;
	out << "��" << mFrameIndex << "֡" << endl;
	out << "��ʼ���ٿ�:";
	out << mTarget << endl;

	sigma = sigma * mScale;

	getFilterWindow();
	getContext(mFrame, mTarget);

	Mat priorFourier;
	Mat planes1[] = { mCxtPriProMat, Mat::zeros(mCxtPriProMat.size(), CV_64FC1) };
	merge(planes1, 2, priorFourier);
	dft(priorFourier, priorFourier);

	Mat conditionalFourier;
	complexMultiply(mSTCModelFFT, priorFourier, conditionalFourier);
	dft(conditionalFourier, mNextCxtPostProMat, DFT_INVERSE | DFT_REAL_OUTPUT | DFT_SCALE);

	Point point;
	minMaxLoc(mNextCxtPostProMat, 0, 0, 0, &point);
	double mConfidence = *(((double *)(mNextCxtPostProMat.data + point.y * mNextCxtPostProMat.step)) + point.x);
	out << "STC����������Ŷ�����:" << point << ";�������Ŷ�:" << mConfidence << endl;;

	//����Ҫ��һ�����ص�ƫ����������Ϊpoint�õ����Ǿ�����Ԫ�ص����꣬�߶ȷ�ΧΪ0��(size-1)����mCxtRegionSize�Ŀ�߳����õ�����������ĳ߶ȷ�ΧΪ1��size
	//���������һ�����ص�ƫ������ʹ���������䵽һ���߶ȷ�Χ�ڣ���ɵõ����������ƶ���ƫ����
	mTarget.x = mTarget.x + point.x + 1 - mCxtRegionSize.width / 2;
	mTarget.y = mTarget.y + point.y + 1 - mCxtRegionSize.height / 2;
	out << "STCԤ�����ٿ�:";
	out << mTarget << endl;

	return mConfidence;
}


void STCTracker::updateWithSTC(const Mat& mFrame, Target& mTarget, ofstream& out)
{
	getFilterWindow();
	getContext(mFrame, mTarget);

	Mat priorFourier;
	Mat planes2[] = { mCxtPriProMat, Mat::zeros(mCxtPriProMat.size(), CV_64FC1) };
	merge(planes2, 2, priorFourier);
	dft(priorFourier, priorFourier);

	Mat conditionalFourier2;
	complexMultiply(mSTCModelFFT, priorFourier, conditionalFourier2);
	dft(conditionalFourier2, mNextCxtPostProMat, DFT_INVERSE | DFT_REAL_OUTPUT | DFT_SCALE);

	Point point;
	minMaxLoc(mNextCxtPostProMat, 0, 0, 0, &point);
	double mConfidence = mEachFrameConfd[mFrameIndex] = *(((double *)(mNextCxtPostProMat.data + point.y * mNextCxtPostProMat.step)) + point.x);
	if (mConfidence > 0.001)
	{
		if (mConfidence < mMinConfiPosi)
		{
			mMinConfiPosi = mConfidence;
		}
	}
	else
	{
		mEachFrameConfd[mFrameIndex] = mMinConfiPosi;
	}
	out << "����STC������º������Ŷ�:" << mEachFrameConfd[mFrameIndex] << ";��ʵ���Ŷ�:" << mConfidence << endl;

	if (mFrameIndex % (mIntervalNumber + 1) == 0)
	{
		double mSum = 0;
		for (int i = 0; i < mIntervalNumber; i++)
		{
			mSum = mSum + sqrt(mEachFrameConfd[mFrameIndex - i] / mEachFrameConfd[mFrameIndex - i - 1]);
		}
		mScale = (1 - mScaleLearnRate)*mScale + mScaleLearnRate*(mSum / mIntervalNumber);
	}

	mRectSize.width = mRectSize.width * mScale;
	mRectSize.height = mRectSize.height * mScale;
	mTarget.wid = (int)mRectSize.width;
	mTarget.hei = (int)mRectSize.height;

	out << "Scale-->" << mScale << endl;
	out << "ʹ��STC���ٽ��:";
	out << mTarget << endl << endl;

	learnSCModol();
	updateSTCModol();

}

void STCTracker::updateWithParti(const Mat& mFrame, const Target& mTarget, ofstream& out)
{
	getFilterWindow();
	getContext(mFrame, mTarget);

	Mat priorFourier;
	Mat planes2[] = { mCxtPriProMat, Mat::zeros(mCxtPriProMat.size(), CV_64FC1) };
	merge(planes2, 2, priorFourier);
	dft(priorFourier, priorFourier);

	Mat conditionalFourier2;
	complexMultiply(mSTCModelFFT, priorFourier, conditionalFourier2);
	dft(conditionalFourier2, mNextCxtPostProMat, DFT_INVERSE | DFT_REAL_OUTPUT | DFT_SCALE);

	Point point;
	minMaxLoc(mNextCxtPostProMat, 0, 0, 0, &point);
	double mConfidence = mEachFrameConfd[mFrameIndex] = *(((double *)(mNextCxtPostProMat.data + point.y * mNextCxtPostProMat.step)) + point.x);
	
	if (mConfidence > 0.001)
	{
		if (mConfidence < mMinConfiPosi)
		{
			mMinConfiPosi = mConfidence;
		}
	}
	else
	{
		mEachFrameConfd[mFrameIndex] = mMinConfiPosi;
	}
	out << "����PF������º������Ŷ�:" << mEachFrameConfd[mFrameIndex] << ";��ʵ���Ŷ�:" << mConfidence << endl;

	if (mFrameIndex % (mIntervalNumber + 1) == 0)
	{
		double mSum = 0;
		for (int i = 0; i < mIntervalNumber; i++)
		{
			mSum = mSum + sqrt(mEachFrameConfd[mFrameIndex - i] / mEachFrameConfd[mFrameIndex - i - 1]);
		}
		mScale = (1 - mScaleLearnRate)*mScale + mScaleLearnRate*(mSum / mIntervalNumber);
	}

	mRectSize.width = mTarget.wid;
	mRectSize.height = mTarget.hei;

	out << "Scale-->" << mScale << endl;
	out << "ʹ��PF���ٽ��:";
	out << mTarget << endl << endl;

	learnSCModol();
	updateSTCModol();

}


void STCTracker::updateWithPartiByDeny(const Mat& mFrame, const Target& mTarget, ofstream& out)
{
	mEachFrameConfd[mFrameIndex] = mEachFrameConfd[mFrameIndex - 1];

	if (mFrameIndex % (mIntervalNumber + 1) == 0)
	{
		double mSum = 0;
		for (int i = 0; i < mIntervalNumber; i++)
		{
			mSum = mSum + sqrt(mEachFrameConfd[mFrameIndex - i] / mEachFrameConfd[mFrameIndex - i - 1]);
		}
		mScale = (1 - mScaleLearnRate)*mScale + mScaleLearnRate*(mSum / mIntervalNumber);
	}

	mRectSize.width = mTarget.wid;
	mRectSize.height = mTarget.hei;

	out << "ʹ��PF���ٽ��:";
	out << mTarget << endl << endl;
}


bool STCTracker::checkConfidenceLegal(const double confid, ofstream& out, const int mFrameNum, const double mThreshold)
{
	if (mFrameIndex <= mFrameNum)
	{
		return true;
	}

	double sum = 0;
	for (int i = 1; i <= mFrameNum; i++)
	{
		sum += mEachFrameConfd[mFrameIndex - i];
	}
	double mean = sum / mFrameNum;

	if ((mean - confid) / mean > mThreshold)
	{
		out << "�������ŶȽ�������" << mThreshold << "�����������˲����" << endl;
		return false;
	}

	return true;
}

int STCTracker::getFrameIndex()
{
	return mFrameIndex;
}
