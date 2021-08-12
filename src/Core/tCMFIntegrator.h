//
// Created by Roman Ellerbrock on 3/8/20.
//

#ifndef TCMFINTEGRATOR_H
#define TCMFINTEGRATOR_H
#include "Core/tHamiltonianRepresentation.h"
#include "Core/tLayerInterface.h"
#include "Core/tIntegratorVariables.h"
#include "Util/BS_integrator.h"
#include "Util/RungeKutta4.h"

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
using tbs_integrator = BS_integrator<tLayerInterface<T>&, Tensor<T>, T>;
template <typename T>
using trk_integrator = RungeKutta4::RK_integrator<tLayerInterface<T>&, Tensor<T>, T>;

template <typename T>
class tCMFIntegrator {
public:
	tCMFIntegrator(const SOP<T>& H,
		const Tree& tree, T phase);

	~tCMFIntegrator() = default;

	void Integrate(tIntegratorVariables<T>& ivars, ostream& os = cout);

	double Error(const TensorTree<T>& Psi, const TensorTree<T>& Chi, const MatrixTree<T>& rho, const Tree& tree) const;

	void Output(double time, const TensorTree<T>& Psi,
		const TensorTree<T>& Psistart, const SOP<T>& H,
		const Tree& tree, ostream& os);

private:

	void CMFstep(TensorTree<T>& Psi, double time, double timeend, double accuracy_leaf, const Tree& tree);

	tHamiltonianRepresentation<T> matrices_;
	vector<tbs_integrator<T>> bs_integrators_;
	vector<trk_integrator<T>> rk_integrators_;
	vector<double> dt_bs_;
	vector<tLayerInterface<T>> interfaces_;

	double tconst_;
	double max_increase_;
};


#endif //TCMFINTEGRATOR_H
