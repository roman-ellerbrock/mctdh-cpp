//
// Created by Roman Ellerbrock on 3/8/20.
//

#ifndef CMFINTEGRATOR_H
#define CMFINTEGRATOR_H
#include "Core/tHamiltonianRepresentation.h"
#include "Core/tLayerInterface.h"
#include "Core/IntegratorVariables.h"
#include "Util/BS_integrator.h"

/**
 * \defgroup Integrator
 * \ingroup mctdh
 * \brief This module holds all classes related to mctdh-Integrators.
 *
 * */

/**
 * \class CMFIntegrator
 * \ingroup Integrator
 * \brief The constant mean-field integrator for solving the
 * mctdh-EOM.
 *
 * In the CMF-Integrator the mctdh-Matrices are assumed to be constant for
 * a given time to propagate a Wavefunction.
 *
 * */

template <typename T>
using bs_integrator = BS_integrator<tLayerInterface<T>&, Tensor<T>, T>;

template <typename T>
class tCMFIntegrator {
public:
	tCMFIntegrator(const SOP<T>& H,
		const Tree& tree, T phase);

	~tCMFIntegrator() = default;

	void Integrate(IntegratorVariables& ivars, ostream& os = cout);

	double Error(const TensorTree<T>& Psi, const TensorTree<T>& Chi, const MatrixTree<T>& rho, const Tree& tree) const;

	void Output(double time, const TensorTree<T>& Psi,
		const TensorTree<T>& Psistart, const SOP<T>& H,
		const Tree& tree, ostream& os);

private:

	void CMFstep(TensorTree<T>& Psi, double time, double timeend, double accuracy_leaf, const Tree& tree);

	tHamiltonianRepresentation<T> matrices_;

	vector<bs_integrator<T>> bs_integrators_;
	vector<double> dt_bs_;
	vector<tLayerInterface<T>> interfaces_;

	double tconst_;
	double max_increase_;
};


#endif //CMFINTEGRATOR_H
