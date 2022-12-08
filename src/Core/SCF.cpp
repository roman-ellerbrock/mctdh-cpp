//
// Created by Roman Ellerbrock on 11/14/22.
//

#include "SCF.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "Core/HamiltonianRepresentation.h"
#include "TreeClasses/TreeIO.h"

void addNodes(vector<const Node *>& sweep, const Node *p) {
	sweep.push_back(p);
	if (p->isBottomlayer()) {
		return;
	} else {
		for (int k = 0; k < p->nChildren(); ++k) {
			const Node *child = &(p->child(k));
			addNodes(sweep, child);
			sweep.push_back(p);
		}
	}
}

vector<const Node *> scf_sweep(const Tree& tree) {
	const Node *p = &tree.topNode();

	vector<const Node *> sweep;
	addNodes(sweep, p);
	sweep.push_back(nullptr);
	return sweep;
}

Tensorcd apply(const Tensorcd& A,
	const SparseMatrixTreescd& hMats, const SparseMatrixTreescd& hCons,
	const SparseMatrixTreecd& hCorr, const SparseMatrixTreecd& hConCorr,
	const MatrixTreecd *rho, const Hamiltonian& H, const Node& node) {

	Tensorcd HA(A.shape());

	for (size_t l = 0; l < hMats.size(); ++l) {
		Tensorcd hA(A.shape());

		if (node.isBottomlayer()) {
			/// only apply part if active hole
			if (hMats[l].isActive(node) && hCons[l].isActive(node)) {
				const MLOcd& M = H[l];
				const Leaf& leaf = node.getLeaf();
				hA = M.apply(A, leaf);
				hA = matrixTensor(hCons[l][node], hA, node.parentIdx());
				HA += H.coeff(l) * hA;
			}
		} else {
			if (nActives(hMats[l], H[l], hCons[l], node) > 1) {
				TreeFunctions::apply(hA, hMats[l], &hCons[l], rho, A,
					hCons[l].sparseTree(), node, -1);
				HA += H.coeff(l) * hA;
			}
		}
	}

	if (!node.isBottomlayer()) {
		/// add correction
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			HA += matrixTensor(hCorr[child], A, child.childIdx());
		}
	}
	if (!node.isToplayer()) {
		HA += matrixTensor(hConCorr[node], A, node.parentIdx());
	}
	return HA;
}

Tensorcd apply(const Tensorcd& A, const HamiltonianRepresentation& hrep,
	const Hamiltonian& H, const Node& node) {
//	const MatrixTreecd *rho = &hrep.rho_;
	const MatrixTreecd *rho = nullptr;
	const SparseMatrixTreecd& hCorr = hrep.hCorr_;
	const SparseMatrixTreecd& hConCorr = hrep.hConCorr_;
	return apply(A, hrep.hMats_, hrep.hContractions_, hCorr, hConCorr, rho, H, node);
}

complex<double> fullContraction(const Tensorcd& A, const Tensorcd& B) {
	auto rho = A.dotProduct(B);
	return rho.trace();
}

double normalize(Tensorcd& A) {
	double norm = real(fullContraction(A, A));
	norm = sqrt(norm);
	A /= norm;
	return norm;
}

void testKrylovSpace(const KrylovSpace& space,
	const HamiltonianRepresentation& hrep,
	const Hamiltonian& H, const Node& node) {

	size_t dim = space.space_.size();
	Matrixcd S(dim, dim);

	for (size_t l = 0; l < dim; ++l) {
		for (size_t m = 0; m < dim; ++m) {
			S(l, m) = fullContraction(space.space_[l], space.space_[m]);
		}
	}
	cout << "|1 - S| = " << residual(identityMatrixcd(dim), S) << endl;

	Matrixcd hm(dim, dim);
	for (size_t l = 0; l < dim; ++l) {
		auto HPsi = apply(space.space_[l], hrep, H, node);
		for (size_t m = 0; m < dim; ++m) {
			hm(l, m) = fullContraction(HPsi, space.space_[m]);
		}
	}
	cout << "<<H>> = \n";
	(hm * QM::cm).print();
}

KrylovSpace solveKrylovSpace(Tensorcd Psi, const HamiltonianRepresentation& hrep,
	const Hamiltonian& H, const Node& node, size_t krylov_size) {

	vector<Tensorcd> krylovSpace(krylov_size, Tensorcd(Psi.shape()));

	/// normalize wavepacket
	auto norm = normalize(Psi);
	if (abs(norm - 1.) > 1e-5) {
		cerr << "Norm != 1.\n";
		cerr << "norm: " << norm << ", renormalizing." << endl;
//		exit(1);
	}
	krylovSpace[0] = Psi;

	vector<double> beta;
	vector<double> alpha;

	/// Do the first step separately, since no beta exists
	Tensorcd HPsi = apply(Psi, hrep, H, node);

	/// Calculate first alpha
	/// Note: in general <Psi_i|H|Psi_j> is a matrix if multiple states are stored in Psi
	complex<double> alphamat = fullContraction(Psi, HPsi);

	/// We use only a single state, so we read the only element of the 1x1 matrix
	alpha.push_back(real(alphamat));

	/// First step is simpler, because beta = 0
	Psi = HPsi - alpha.back() * Psi;

	/// remaining steps
	for (size_t l = 1; l < krylov_size; ++l) {
		beta.push_back(normalize(Psi));
		krylovSpace[l] = Psi;

		HPsi = apply(Psi, hrep, H, node);
		alphamat = fullContraction(Psi, HPsi);
		alpha.push_back(real(alphamat));

		const Tensorcd& v = krylovSpace[l];
		const Tensorcd& vlast = krylovSpace[l - 1];
		Psi = HPsi - alpha.back() * v - beta.back() * vlast;
	}

	/// Build tri-diagonal hamiltonian matrix
	Matrixcd Hrep(krylov_size, krylov_size);
	for (size_t l = 0; l < krylov_size; ++l) {
		Hrep(l, l) = alpha[l];
		if (l == krylov_size - 1) { break; }
		Hrep(l + 1, l) = beta[l];
		Hrep(l, l + 1) = beta[l];
	}
//	cout << "<H>:\n";
//	(Hrep * QM::cm).print();

	/// Diagonalize Hamiltonian matrix and save transformation matrix U and eigenvalues E.
	auto spectrum = diagonalize(Hrep);
//	cout << "E_krylov: ";
	(spectrum.second * QM::cm).print();
	return KrylovSpace(krylovSpace, spectrum);
}

void imaginaryTimePropagation(Tensorcd& Psi, const KrylovSpace& kry, double beta) {
	const Matrixcd& U = kry.spectrum_.first;
	const Vectord& E = kry.spectrum_.second;

	size_t krylov_size = kry.space_.size();
	Vectorcd c(krylov_size);
	for (size_t l = 0; l < krylov_size; ++l) {
		c(l) += U(l, 0) * conj(U(0, 0));
		for (size_t k = 0; k < krylov_size; ++k) {
			complex<double> exp_k = exp(-beta * (E(k) - E(0)));
			c(l) += U(l, k) * exp_k * conj(U(0, k));
		}
	}

	Psi.zero();
	for (size_t l = 0; l < krylov_size; ++l) {
		const Tensorcd& psi_l = kry.space_[l];
		Psi += c(l) * psi_l;
	}
	normalize(Psi);
}

void diagonalizeKrylov(Tensorcd& Psi, const KrylovSpace& kry) {
	const Matrixcd& U = kry.spectrum_.first;
	const Vectord& E = kry.spectrum_.second;

	size_t krylov_size = kry.space_.size();
	Psi.zero();
	for (size_t l = 0; l < krylov_size; ++l) {
		const Tensorcd& psi_l = kry.space_[l];
		Psi += U(l, 0) * psi_l;
	}
}

size_t adjacentIndex(const Node& from, const Node *to) {
	/// returns index pointing from node *from* to node *to*
	/// caution when touching the logic/order of conditions here!
	if (to == nullptr) {
		return from.parentIdx();
	}
	for (size_t k = 0; k < from.nChildren(); ++k) {
		if (from.child(k).address() == to->address()) {
			return k;
		}
	}
	if (from.parent().address() == to->address()) {
		return from.parentIdx();
	}
	cerr << "Node *to* is not adjacent to node *from*.\n";
	exit(1);
	return 0;
}

void outputSCF(Tensorcd Psi, const HamiltonianRepresentation& hrep,
	const Hamiltonian& H, const Node& node) {

	auto HPsi = apply(Psi, hrep, H, node);
	cout << setprecision(12);
	cout << "<H> = " << real(fullContraction(Psi, HPsi)) * QM::cm << endl;
}

void scf(SCF_parameters& par) {
	cout << "=====    SCF    =====\n";
	cout << "nIter = " << par.nIter << endl;
	cout << "nKrylov = " << par.nKrylov << endl;
	double time = 0;
	size_t krylov_size = par.nKrylov;
	double beta = par.beta;
	const Tree& tree = *par.tree;
	const Hamiltonian& H = *par.h;
	TensorTreecd& Psi = *par.psi;

	auto sweeper = scf_sweep(tree);
	double eps = 1. / ((double) (sweeper.size() - 1)); /// just for printing output

	HamiltonianRepresentation hrep(H, tree, tree, false);
//	hrep.initializeDense(H, tree); /// create this time with dense matrices
//	hrep.build(H, Psi, tree, time);
	cout << setprecision(12);
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			const Node& next = node.parent();
			hrep.buildSCF(H, Psi, node, next, time);
		}
	}

	for (size_t it = 0; it < par.nIter; ++it) {
//		cout << "Iteration " << it << endl;

		for (size_t l = 0; l < sweeper.size() - 1; ++l) {
			if (sweeper[l] == nullptr) { break; }

			/// Build Krylov Space
			const Node& node = *sweeper[l];
			size_t tot = node.shape().totalDimension();
			size_t ksize = krylov_size > tot ? tot : krylov_size;
			cout << it + l * eps << " ";
			KrylovSpace krylov = solveKrylovSpace(Psi[node], hrep, H, node, ksize);
			diagonalizeKrylov(Psi[node], krylov);

			/// isometrize Psi towards next node
			const Node *next_ptr = sweeper[l + 1];
			size_t outIdx = adjacentIndex(node, next_ptr);
			Tensorcd PsiW = Psi[node];
			Psi[node] = qr(Psi[node], outIdx);

			/// overlap with <I | A> and multiply into adjacent node
			if (next_ptr == nullptr) { break; }
			const Node& next = *next_ptr;
			auto r = contraction(Psi[node], PsiW, outIdx);
			size_t inIdx = adjacentIndex(next, &node);
			Psi[next] = matrixTensor(r, Psi[next], inIdx);

			/// rebuild matrix elements pointing towards next node
			hrep.buildSCF(H, Psi, node, next, time);
		}
	}
//	TreeIO::output(Psi, tree);
}
