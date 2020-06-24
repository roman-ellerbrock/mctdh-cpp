#pragma once
#include "DVR/Potential.h"

class CDVRModelV :
	public Potential
{
public:
	CDVRModelV(size_t f_ = 4, bool coupling = true);
	~CDVRModelV() = default;

	double Evaluate(const Vectord& Xv, size_t part)const;

protected:
	size_t f;
	bool coupling_;
};

