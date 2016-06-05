#pragma once

#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

class Particle;

class Target
{

public:
	int x;		//中心点X坐标
	int y;		//中心点Y坐标
	int wid;	//半长
	int hei;	//半宽

public:
	Target();
	Target(const Target& mTarget);
	virtual ~Target();

public:
	Target& operator=(const Target& r);

public:
	void set(const Rect& r);
	void set(const Particle& t);
	void set(int x, int y, int wid, int hei);

	/*
	* 根据指定的图像大小，判断当前这个目标区域的有效性
	* true为合法有效，false为无效
	*/
	bool isValid(int mWidth, int mHeight) const;

public:
	friend ofstream& operator<<(ofstream& out, const Target& mTarget);

};