#include "STCResult.h"

STCResult::STCResult(){}

STCResult::STCResult(const STCResult& r)
{
	mType = r.mType;
	mResult = r.mResult;
}

STCResult::~STCResult(){}

STCResult& STCResult::operator=(const STCResult& r)
{
	if (this == &r)
	{
		return *this;
	}
	this->mType = r.mType;
	this->mResult = r.mResult;
	return *this;
}

void STCResult::setResultType(STCResultType mType)
{
	this->mType = mType;
}

void STCResult::setResult(double mResult)
{
	this->mResult = mResult;
}

void STCResult::setPreResult(double mPreResult)
{
	this->mPreResult = mPreResult;
}

STCResultType STCResult::getResultType() const
{
	return mType;
}

double STCResult::getResult() const
{
	return mResult;
}

double STCResult::getPreResult() const
{
	return mPreResult;
}

bool STCResult::isLegal(double mThreshold) const
{
	return (mPreResult - mResult) / mPreResult < mThreshold;
}

ofstream& operator<<(ofstream&out, STCResult& mSTCResult)
{
	if (mSTCResult.mType == PRIOR)
	{
		out << "先验置信度：" << mSTCResult.mResult;
	}
	else
	{
		out << "后验置信度：" << mSTCResult.mResult;
	}
	return out;
}

