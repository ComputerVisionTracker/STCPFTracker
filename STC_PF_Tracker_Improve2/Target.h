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
	int x;		//���ĵ�X����
	int y;		//���ĵ�Y����
	int wid;	//�볤
	int hei;	//���

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
	* ����ָ����ͼ���С���жϵ�ǰ���Ŀ���������Ч��
	* trueΪ�Ϸ���Ч��falseΪ��Ч
	*/
	bool isValid(int mWidth, int mHeight) const;

public:
	friend ofstream& operator<<(ofstream& out, const Target& mTarget);

};