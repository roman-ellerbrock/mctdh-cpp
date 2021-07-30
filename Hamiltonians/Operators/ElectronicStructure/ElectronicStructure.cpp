//
// Created by Roman Ellerbrock on 7/27/21.
//

#include "ElectronicStructure.h"

using namespace JordanWigner;

TwoIndex readTwoIndexIntegral(const string& filename) {
	ifstream file(filename);

	TwoIndex hpq;
	while (file.peek() != EOF) {
		size_t p, q = 0;
		double h = 0.;
		file >> p >> q >> h;
//		cout << p << " " << q << " " << h << endl;
		hpq.push_back({p, q, h});
	}
	return hpq;
}

FourIndex readFourIndexIntegral(const string& filename) {
	ifstream file(filename);

	FourIndex hpqrs;
	while (file.peek() != EOF) {
		size_t p, q, r, s = 0;
		double h = 0.;
		file >> p >> q >> r >> s >> h;
//		cout << p << " " << q << " " << r << " " << s << " " << h << endl;
		hpqrs.push_back({p, q, r, s, h});
	}
	return hpqrs;
}

Matrixd convertTwoIndex(const TwoIndex& h) {
	size_t dim = 1;
	for (auto pq : h) {
		size_t p = get<0>(pq);
		if (p >= dim) { dim = p + 1; }
	}
	Matrixd ht(dim, dim);
	for (auto pq : h) {
		size_t p = get<0>(pq);
		size_t q = get<1>(pq);
		double val = get<2>(pq);
		ht(p, q) = val;
	}
	return ht;
}

Tensord convertFourIndex(const FourIndex& h) {
	size_t dim = 1;
	for (auto pqrs : h) {
		size_t p = get<0>(pqrs);
		if (p >= dim) { dim = p + 1; }
	}
	Tensord ht({dim, dim, dim, dim});
	for (auto pqrs : h) {
		size_t p = get<0>(pqrs);
		size_t q = get<1>(pqrs);
		size_t r = get<2>(pqrs);
		size_t s = get<3>(pqrs);
		double val = get<4>(pqrs);
		ht({p, q, r, s}) = val;
	}
	return ht;
}

SOPcd electronicStructure(const string& twoindex, const string& fourindex) {
	auto hpq_sparse = readTwoIndexIntegral(twoindex);
	Matrixd hpq = convertTwoIndex(hpq_sparse);
	cout << "hpq:\n";
	hpq.print();

	auto hpqrs_sparse = readFourIndexIntegral(fourindex);
	Tensord hpqrs = convertFourIndex(hpqrs_sparse);
	cout << "hpqrs:\n";
	hpqrs.print();

	return JordanWigner::electronicHamiltonian(hpq, hpqrs);
}
