//
// Created by Roman Ellerbrock on 5/3/22.
//

#include "PortfolioOptimization.h"
#include "TreeOperators/SumOfProductsOperator.h"


Tensord readAssets(const vector<string>& names, size_t m, size_t Delta) {
	size_t N = names.size();

	Tensord assets({m + Delta, N});
	size_t j = 0;
	for (const string& name : names) {
		ifstream is(name);
		for (size_t i = 0; i < m + Delta; ++i) {
			is >> assets(i, j);
		}
		j++;
	}
	return assets;
}

Tensord returns(const Tensord& A, size_t M, size_t Delta) {
	const TensorShape& shape = A.shape_;
	size_t N = shape.lastDimension(); // n_assets
	Tensord mu({M, N});
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < M; ++m) {
			mu(m, n) = (A(m + Delta, n) - A(m, n)) / A(m, n);
		}
	}
	return mu;
}

Tensord average(const Tensord& A, size_t M, size_t Delta) {
	const TensorShape& shape = A.shape_;
	size_t N = shape.lastDimension(); // n_assets
	Tensord e({M, N});
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < M; ++m) {
			for (size_t d = 0; d < Delta; ++d) {
				e(m, n) += A(m + d, n);
			}
		}
	}
	e /= (double) Delta;
	return e;
}

Tensord sigma(Tensord A, size_t M, size_t Delta) {
	auto e = average(A, M, Delta);

	size_t N = A.shape_.lastDimension(); // n_assets

	TensorShape shape({M, N, N});
	Tensord sigma(shape);
	for (size_t n = 0; n < N; ++n) {
		for (size_t l = 0; l < N; ++l) {
			for (size_t m = 0; m < M; ++m) {
				for (size_t d = 0; d < Delta; ++d) {
					sigma({m, l, n}) += (A({m + d, l}) - e({m, l})) * (A({m + d, n}) - e({m, n}));
				}
			}
		}
	}
	sigma /= (double) Delta;
	return sigma;
}

Tensord timeSeries(const Tensord& A) {
	size_t m = A.shape_.lastBefore() - 1;
	size_t N = A.shape_.lastDimension();
	TensorShape shape({m, N});
	Tensord dA(shape);
	for (size_t j = 0; j < N; ++j) {
		for (size_t i = 0; i < m; ++i) {
			dA(i, j) = (A(i, j) - A(i + 1, j)) / A(i, j);
//			dA(i, j) = A(i, j) - A(i + 1, j);
		}
	}
	return dA;
}

SOPcd meanVarianceAnalysis() {

	cout << "Mean Variance Analysis\n";
	vector<string> names = {"BTC-USD.csv", "ETH-USD.csv", "ADA-USD.csv"};
	size_t N = names.size(); // number of assets
	size_t m = 91; // number of time steps
	size_t Delta = 1;

	/// read
	auto A = readAssets(names, m, Delta);
	auto dA = timeSeries(A);
	cout << "time series of assets:\n";
	dA.print();

	auto sig = sigma(dA, m, Delta);
	cout << "cov:\n";
	sig.print();
	getchar();

	/// How do you obtain mu, and sigma? (forecast return and covariance)
	/// For proof of concepts: 3 month until now on a daily basis

	return {SOPcd()};
}
