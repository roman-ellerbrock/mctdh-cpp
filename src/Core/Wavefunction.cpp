//
// Created by Roman Ellerbrock on 1/13/22.
//
#include "Wavefunction.h"

void occupySingles(Tensorcd& A, size_t& n, const Node& node) {
	const TensorShape& shape = A.shape();

	// only do singles, if all fit in!
	size_t req = 1;
	for (size_t i = 0; i < shape.order() - 1; ++i) {
		req += shape[i] - 1;
	}
	if (req > shape.lastDimension()) { return; }

	// first excited states
	for (size_t i = 0; i < shape.order() - 1; ++i) {
		for (size_t m = 1; m < shape[i]; ++m) {
			vector<size_t> idx(shape.order(), 0);
			idx[i] = m;
			idx[shape.lastIdx()] = n++;
			A(idx) = 1.;
		}
	}
}

void occupyRandom(Tensorcd& A, size_t& n) {
	const TensorShape& shape = A.shape();
	// occupy with random values
	mt19937 gen(time(nullptr));
	uniform_real_distribution<double> dist(-1., 1.);
	for (; n < shape.lastDimension(); ++n) {
		for (size_t i = 0; i < shape.lastBefore(); ++i) {
			A(i, n) = dist(gen);
		}
		gramSchmidt(A);
	}
}

void occupyCIS(Tensorcd& A, const Node& node) {
	A.zero();
	// gs
	size_t n = 1;
	A(0) = 1.;

	occupySingles(A, n, node);
	occupyRandom(A, n);


/*	Tensorcd spf({2, 2});
	//	spf[0] = 0.844866; // 13 Fragments
	//	spf[1] = -0.534978;
	spf[0] = 0.844866;
	spf[1] = -0.534978;
	spf[2] = -spf[1];
	spf[3] = spf[0];
	gramSchmidt(spf);
	if (node.isBottomlayer())  { A = spf; }*/

	gramSchmidt(A);
}

void occupyCIS(TensorTreecd& Psi, const Tree& tree) {

	for( const Node& node : tree) {
		if (!node.isBottomlayer()) {
			occupyCIS(Psi[node], node);
		}
	}
}


