#pragma once
#include "TreeOperators/Potential.h"

class NOClPotential
	:public Potential
{
public:
	explicit NOClPotential(bool dissociation = true);
	~NOClPotential()override = default;
	double Evaluate(const Vectord& Xv, size_t part) const override ;

private:
    double EvaluateGS(const Vectord & Xv, size_t part);
    double EvaluateS1(const Vectord & Xv, size_t part);
    bool dissociation;
};

