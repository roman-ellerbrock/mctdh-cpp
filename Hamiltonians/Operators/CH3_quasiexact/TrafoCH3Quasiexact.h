#pragma once
#include "CoordinateTransformation.h"


class TrafoCH3Quasiexact :
	public CoordinateTransformation
{
  public:
	TrafoCH3Quasiexact(Vectord mass_);

	Vectord Transform(const Vectord& q)const override;

private:
	Vectord mass;
};


