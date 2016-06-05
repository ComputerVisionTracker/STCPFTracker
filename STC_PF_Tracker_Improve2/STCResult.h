#pragma once

#include <iostream>
#include <fstream>
#include "STCResultType.h"

using namespace std;

class STCResult
{
private:
	STCResultType mType;
	double mResult;
	double mPreResult;

public:
	STCResult();
	virtual ~STCResult();
	STCResult(const STCResult& r);

public:
	STCResult& operator=(const STCResult& r);
	

public:
	void setResultType(STCResultType mType);
	void setResult(double mResult);
	void setPreResult(double mPreResult);
	STCResultType getResultType() const;
	double getResult() const;
	double getPreResult() const;
	bool isLegal(double mThreshold) const;

public:
	friend ofstream& operator<<(ofstream& out, STCResult& mSTCResult);

};