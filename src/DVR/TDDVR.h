//
// Created by Roman Ellerbrock on 3/9/20.
//

#ifndef TDDVR_H
#define TDDVR_H
#include "DVR/TreeGrids.h"
#include "Core/Wavefunction.h"
#include "DVR/XMatrixTrees.h"
#include "TreeClasses/MatrixTree.h"

class TDDVR {
public:
	explicit TDDVR(const Tree& tree) : Xs_(tree), rho_(tree), trafo_(tree), hole_trafo_(tree), grids_(tree){
		auto xs = Xsop(tree);
		for (const auto& x : xs) {
			hole_grids_.emplace_back(VectorTreed(x, tree));
		}
	}
	~TDDVR() = default;

	void Calculate(const Wavefunction& Psi, const Tree& tree);

	void print(const Tree& tree) const {
		cout << "TDDVR: " << endl;
		cout << "Grids:" << endl;
		for (const Node& node : tree) {
			size_t dim = trafo_[node].Dim1();
			node.info();
			for (size_t i = 0; i < dim; ++i) {
				for (const VectorTreed& grid : grids_) {
					if (grid.Active(node)) {
						const Vectord& g = grid[node];
						cout << g(i) << "\t";
					}
				}
				cout << endl;
			}
		}
		cout << "Hole grids:" << endl;
		for (const Node& node : tree) {
			size_t dim = trafo_[node].Dim1();
			node.info();
			for (size_t i = 0; i < dim; ++i) {
				for (const VectorTreed& grid : hole_grids_) {
					if (grid.Active(node)) {
						const Vectord& g = grid[node];
						cout << g(i) << "\t";
					}
				}
				cout << endl;
			}
		}
	}

	TreeGrids grids_;
	MatrixTreecd trafo_;

	TreeGrids hole_grids_;
	MatrixTreecd hole_trafo_;
private:
	XMatrixTrees Xs_;
	MatrixTreecd rho_;

};


#endif //TDDVR_H
