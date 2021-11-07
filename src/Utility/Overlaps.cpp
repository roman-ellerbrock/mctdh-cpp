//
// Created by Roman Ellerbrock on 10/21/21.
//

#include "Overlaps.h"
#include "TreeClasses/MatrixTreeFunctions.h"

vector<TensorTreecd> readTensorTrees(const string& file1) {
	ifstream is(file1);
	vector<TensorTreecd> psi;
	int i = 0;
	while (!is.eof()) {
		psi.emplace_back(TensorTreecd(is));
		cout << i++ << endl;
		if (is.peek() == EOF) { break; }
	}
	return psi;
}

void wavefunctionOverlap(const string& file1, const string& file2,
	const Tree& tree, OverlapType type) {
	auto bra = readTensorTrees(file1);
	auto ket = readTensorTrees(file2);
	cout << "Number of wavefunctions found for bra:" << bra.size() << endl;
	cout << "Number of wavefunctions found for ket:" << ket.size() << endl;
	wavefunctionOverlap(bra, ket, tree, type);
}

void wavefunctionOverlap(const vector<TensorTreecd>& Psi,
	const vector<TensorTreecd>& Chi, const Tree& tree, OverlapType type) {
	const Node& top = tree.topNode();
	switch (type) {
		case full:
			for (size_t i = 0; i < Psi.size(); ++i) {
				for (size_t j = 0; j < Chi.size(); ++j) {
					MatrixTreecd S(tree);
					for (const Node& node : tree) {
						size_t n = Psi[i][node].shape().lastDimension();
						size_t m = Chi[j][node].shape().lastDimension();
						S[node] = Matrixcd(n, m);
					}
					TreeFunctions::dotProduct(S, Psi[i], Chi[j], tree);
					auto& s = S[top];
					cout << i << "\t" << j << "\t" << pow(abs(s(0, 0)), 2) << "\t";
					s.print();
				}
			}
			break;
		case diagonal:
			cout << fixed << setprecision(8);
			for (size_t i = 0; i < Psi.size(); ++i) {
				MatrixTreecd S(tree);
				for (const Node& node : tree) {
					size_t n = Psi[i][node].shape().lastDimension();
					size_t m = Chi[i][node].shape().lastDimension();
					S[node] = Matrixcd(n, m);
				}
				TreeFunctions::dotProduct(S, Psi[i], Chi[i], tree);
				auto& s = S[top];
				cout << i << "\t" << pow(abs(s(0, 0)), 2) << "\t" << real(s(0,0)) << "\n";
			}
			break;
	}
}

