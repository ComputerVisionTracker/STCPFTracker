#include "Tracker.h"


Tracker::Tracker()
{
	mConfiThreshold = 0.25;

	out.open("STC_PF_Improve2_Result.txt",ios::out);
}

Tracker::~Tracker()
{
	out.flush();
	out.close();
}


bool Tracker::init(const Mat& mFrame, const Mat& mGrayFrame, const Target& mTarget, int mFrameCount)
{
	if (out.fail())
	{
		return false;
	}
	bool result;
	result = mSTCTracker.init(mGrayFrame, mTarget, mFrameCount, out);
	if (!result)
	{
		return result;
	}
	result = mParticleFilter.init(mFrame, mTarget, mFrameCount);
	return result;
}

void Tracker::track(const Mat& mFrame, const Mat& mGrayFrame, Target& mTarget)
{

	double result = mSTCTracker.track(mGrayFrame, mTarget, out);

	if (mSTCTracker.checkConfidenceLegal(result, out, 3, 0.30))
	{
		mSTCTracker.updateWithSTC(mGrayFrame, mTarget, out);
		mParticleFilter.updateModelWithTarget(mFrame, mTarget);
	}
	else
	{
		mParticleFilter.track(mFrame, mTarget);
		mSTCTracker.updateWithParti(mGrayFrame, mTarget, out);
		//mSTCTracker.updateWithPartiByDeny(mGrayFrame, mTarget, out);
	}
	

}

