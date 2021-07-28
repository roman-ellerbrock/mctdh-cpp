//
// Created by Roman Ellerbrock on 7/27/21.
//

#include "JordanWigner.h"
#include "TreeOperators/LeafMatrix.h"
#include "TreeOperators/SumOfProductsOperator.h"
#include "Util/QMConstants.h"

namespace JordanWigner {
	Matrixcd sigmaX() {
		Matrixcd x(2,2);
		x(1,0)=1.;
		x(0,1)=1.;
		return x;
	}

	Matrixcd sigmaY() {
		Matrixcd y(2,2);
		y(1,0)=QM::im;
		y(0,1)=-QM::im;
		return y;
	}

	Matrixcd sigmaZ() {
		Matrixcd z(2,2);
		z(0,0)=1.;
		z(1,1)=1.;
		return z;
	}

	Matrixcd sigmaPlus() {
		/// Ref. [1] Eq. (19)
		return 0.5 * (sigmaX() - QM::im * sigmaY());
	}

	Matrixcd sigmaMinus() {
		/// Ref. [1] Eq. (19)
		return 0.5 * (sigmaX() + QM::im * sigmaY());
	}

	MLOcd creation(size_t p) {
		/// Ref. [1], Eq. (17), M = (\prod_(i<p) sigmaZ_i) * sigmaMinus_p
		MLOcd M(sigmaPlus(), p);
		for (size_t i = 0; i < p; ++i) {
			M.push_back(sigmaZ(), i);
		}
		return M;
	}

	MLOcd annihilation(size_t p) {
		/// Ref. [1], Eq. (18), M = (\prod_(i<p) sigmaZ_i) * sigmaMinus_p
		MLOcd M(sigmaMinus(), p);
		for (size_t i = 0; i < p; ++i) {
			M.push_back(sigmaZ(), i);
		}
		return M;
	}

	MLOcd twoIndexOperator(size_t p, size_t q, double eps) {
		size_t max = p;
		if (q > max) max = q;

		MLOcd A; /// A = (a_p^+) (a_q^+) (a_r) (a_s)
		for (size_t i = 0; i < max; ++i) {
			Matrixcd op = identityMatrixcd(2);

			if (i < q) {
				op = sigmaZ() * op;
			} else if (i == q) {
				op = sigmaMinus() * op;
			}

			if (i < p) {
				op = sigmaZ() * op;
			} else if (i == p) {
				op = sigmaPlus() * op;
			}
			if (residual(op, identityMatrixcd(2)) > eps) {
				A.push_back(op, i);
			}
		}
		return A;
	}

	MLOcd fourIndexOperator(size_t p, size_t q, size_t r, size_t s, double eps) {
		size_t max = p;
		if (q > max) max = q;
		if (r > max) max = r;
		if (s > max) max = s;

		MLOcd A; /// A = (a_p^+) (a_q^+) (a_r) (a_s)
		/// Rational: swipe over each
		for (size_t i = 0; i < max; ++i) {
			Matrixcd op = identityMatrixcd(2);
			if (i < s) {
				op = sigmaZ() * op;
			} else if (i == s) {
				op = sigmaMinus() * op;
			}

			if (i < r) {
				op = sigmaZ() * op;
			} else if (i == r) {
				op = sigmaMinus() * op;
			}

			if (i < q) {
				op = sigmaZ() * op;
			} else if (i == q) {
				op = sigmaPlus() * op;
			}

			if (i < p) {
				op = sigmaZ() * op;
			} else if (i == p) {
				op = sigmaPlus() * op;
			}

			if (residual(op, identityMatrixcd(2)) > eps) {
				A.push_back(op, i);
			}
		}
		return A;
	}

	SOPcd electronicHamiltonian(const Matrixd& Hpq, const Tensord& Hpqrs) {
		double eps = 1e-12;
		SOPcd H;
		for (size_t q = 0; q < Hpq.dim2(); ++q) {
			for (size_t p = 0; p < Hpq.dim1(); ++p) {
				H.push_back(twoIndexOperator(p, q, eps), Hpq(p, q));
			}
		}

		const TensorShape& shape = Hpqrs.shape();
		for (size_t I = 0; I < shape.totalDimension(); ++I) {
			auto idx = indexMapping(I, shape); /// I -> (p, q, r, s)
			size_t p = idx[0];
			size_t q = idx[1];
			size_t r = idx[2];
			size_t s = idx[3];
			H.push_back(fourIndexOperator(p, q, r, s, eps), Hpqrs(I));
		}
		return H;
	}
}

