#include "Target.h"
#include "Particle.h"

Target::Target(){}

Target::Target(const Target& mTarget)
{
	x = mTarget.x;
	y = mTarget.y;
	wid = mTarget.wid;
	hei = mTarget.hei;
}

Target::~Target(){}

Target& Target::operator=(const Target& r)
{
	if (&r == this)
	{
		return *this;
	}
	x = r.x;
	y = r.y;
	wid = r.wid;
	hei = r.hei;
	return *this;
}

void Target::set(const Rect& mRect)
{
	wid = mRect.width / 2;
	hei = mRect.height / 2;
	x = mRect.x + wid;
	y = mRect.y + hei;
}

void Target::set(const Particle& mPar)
{
	x = mPar.x;
	y = mPar.y;
	wid = mPar.mHalfWidth;
	hei = mPar.mHalfHeight;
}

void Target::set(int x, int y, int wid, int hei)
{
	this->x = x;
	this->y = y;
	this->wid = wid;
	this->hei = hei;
}

bool Target::isValid(int mWidth, int mHeight) const
{
	if (x < 0 || x >= mWidth || y < 0 || y >= mHeight || wid <= 0 || hei <= 0)
	{
		return false;
	}
	return true;
}

ofstream& operator<<(ofstream& out, const Target& mTarget)
{
	out << "X-" << mTarget.x << ",Y-" << mTarget.y << ",Wid-" << mTarget.wid << ",Hei-" << mTarget.hei;
	return out;
}