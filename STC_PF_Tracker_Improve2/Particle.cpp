#include "Particle.h"

const float Particle::POSITION_DISTURB_RATE = (float)0.05;
const float Particle::VELOCITY_DISTURB = (float)40;
const float Particle::WINDOW_DISTURB = (float)0;
const float Particle::SCALE_DISTURB = (float)0.001;

Particle::Particle()
{

}

Particle::Particle(const Particle& r)
{
	x = r.x;
	y = r.y;
	mXVelocity = r.mXVelocity;
	mYVelocity = r.mYVelocity;
	mHalfWidth = r.mHalfWidth;
	mHalfHeight = r.mHalfHeight;
	mScaleVelocity = r.mScaleVelocity;
	weight = r.weight;
}


Particle::~Particle(){}


Particle& Particle::operator=(const Particle& r){
	if (this == &r)
	{
		return *this;
	}
	x = r.x;
	y = r.y;
	mXVelocity = r.mXVelocity;
	mYVelocity = r.mYVelocity;
	mHalfWidth = r.mHalfWidth;
	mHalfHeight = r.mHalfHeight;
	mScaleVelocity = r.mScaleVelocity;
	weight = r.weight;
	return *this;
}

void Particle::setWeight(float weight)
{
	this->weight = weight;
}

void Particle::update(const gsl_rng* rng, const double sigma)
{
	x = (int)(x + mXVelocity*POSITION_DISTURB_RATE + gsl_ran_gaussian(rng, sigma)*mHalfWidth + 0.5);
	y = (int)(y + mYVelocity*POSITION_DISTURB_RATE + gsl_ran_gaussian(rng, sigma)*mHalfHeight + 0.5);
	mXVelocity = (float)(mXVelocity + gsl_ran_gaussian(rng, sigma)*VELOCITY_DISTURB);
	mYVelocity = (float)(mYVelocity + gsl_ran_gaussian(rng, sigma)*VELOCITY_DISTURB);
	double t = gsl_ran_gaussian(rng, sigma);
	mHalfWidth = (int)(mHalfWidth + mHalfWidth*mScaleVelocity + t*WINDOW_DISTURB + 0.5);
	mHalfHeight = (int)(mHalfHeight + mHalfHeight*mScaleVelocity + t*WINDOW_DISTURB + 0.5);
	/*mHalfWidth = (int)(mHalfWidth + mHalfWidth*mScaleVelocity + gsl_ran_gaussian(rng, sigma)*WINDOW_DISTURB + 0.5);
	mHalfHeight = (int)(mHalfHeight + mHalfHeight*mScaleVelocity + gsl_ran_gaussian(rng, sigma)*WINDOW_DISTURB + 0.5);*/
	mScaleVelocity = (float)(mScaleVelocity + gsl_ran_gaussian(rng, sigma)*SCALE_DISTURB);
}

void Particle::drawParticle(Mat& mMat) const
{
	circle(mMat, Point(x, y), 3, CV_RGB(0, 255, 0), 1, 8, 3);
}


void Particle::init(const Target& mTarget, const float weight)
{
	x = mTarget.x;
	y = mTarget.y;
	mXVelocity = 0;
	mYVelocity = 0;
	mHalfWidth = mTarget.wid;
	mHalfHeight = mTarget.hei;
	mScaleVelocity = 0;
	this->weight = weight;
}

void Particle::init(const Target& mTarget, const gsl_rng* rng, const double sigma, const float weight)
{
	x = int(mTarget.x + gsl_ran_gaussian(rng, sigma)*mTarget.wid);
	y = int(mTarget.y + gsl_ran_gaussian(rng, sigma)*mTarget.hei);
	mXVelocity = float(0 + gsl_ran_gaussian(rng, sigma)*VELOCITY_DISTURB);
	mYVelocity = float(0 + gsl_ran_gaussian(rng, sigma)*VELOCITY_DISTURB);
	mHalfWidth = int(mTarget.wid + gsl_ran_gaussian(rng, sigma)*WINDOW_DISTURB);
	mHalfHeight = int(mTarget.hei + gsl_ran_gaussian(rng, sigma)*WINDOW_DISTURB);
	mScaleVelocity = float(0 + gsl_ran_gaussian(rng, sigma)*SCALE_DISTURB);
	this->weight = weight;
}