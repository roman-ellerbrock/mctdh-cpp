//
// Created by Roman Ellerbrock on 5/3/22.
//

#include "PortfolioOptimization.h"
#include "TreeOperators/SumOfProductsOperator.h"


Tensord readAssets(const vector<string>& names, size_t m, size_t Delta) {
	size_t N = names.size();

	Tensord assets({m, N});
	size_t j = 0;
	for (const string& name : names) {
		ifstream is(name);
		for (size_t i = 0; i < m; ++i) {
			is >> assets(i, j);
		}
		j++;
	}
	return assets;
}

/*Tensord returns(const Tensord& A, size_t M, size_t Delta) {
	const TensorShape& shape = A.shape_;
	size_t N = shape.lastDimension(); // n_assets
	Tensord mu({M-Delta, N});
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < M; ++m) {
			mu(m, n) = (A(m + Delta, n) - A(m, n)) / A(m, n);
		}
	}
	return mu;
}*/

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

Tensord timeSeries(const Tensord& A, size_t Delta) {
	size_t m = A.shape_.lastBefore() - 1;
	size_t N = A.shape_.lastDimension();
	TensorShape shape({m - Delta, N});
	Tensord dA(shape);
	for (size_t j = 0; j < N; ++j) {
		for (size_t i = 0; i < m - Delta; ++i) {
			dA(i, j) = A(i + Delta, j) / A(i, j);
		}
	}

	/*	for (size_t j = 0; j < N; ++j) {
			for (size_t i = 0; i < o; ++i) {
				for (size_t k = 0; k < Delta; ++k) {
					size_t l = i * Delta;
					dA(i, j) = A(l + Delta, j) / A(l, j);
				}
			}
		}*/
	return dA;
}

template<typename T>
Tensor<T> log(Tensor<T> A) {
	for (size_t i = 0; i < A.shape_.totalDimension(); ++i) {
		A(i) = log(A(i));
	}
	return A;
}

template<typename T>
Tensor<T> exp(Tensor<T> A) {
	for (size_t i = 0; i < A.shape_.totalDimension(); ++i) {
		A(i) = exp(A(i));
	}
	return A;
}

Tensord returns(const Tensord& A, size_t Delta) {
	size_t m = A.shape_.lastBefore();
	size_t N = A.shape_.lastDimension();
	TensorShape shape({m - Delta, N});
	Tensord dA(shape);
	for (size_t j = 0; j < N; ++j) {
		for (size_t i = 0; i < m - Delta; ++i) {
			dA(i, j) = A(i + Delta, j) / A(i, j);
		}
	}
	return dA;
}

Tensord covariance(const Tensord& mu, const Tensord& avg, size_t Delta) {
	size_t m = mu.shape_.lastBefore();
	size_t N = mu.shape_.lastDimension();
	size_t o = m / Delta;
	TensorShape shape({o, N, N});
	Tensord cov(shape);
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			for (size_t l = 0; l < o; ++l) {
				for (size_t k = 0; k < Delta; ++k) {
					size_t idx = l * Delta + k;
					cov({l, i, j}) += (mu(idx, i) - avg(l, i)) * (mu(idx, j) - avg(l, j));
				}
				cov({l, i, j}) /= (double) Delta;
				cov({l, i, j}) = sqrt(cov({l, i, j}));
			}
		}
	}
	return cov;
}

Tensord avg(const Tensord& A, size_t Delta) {
	size_t m = A.shape_.lastBefore();
	size_t N = A.shape_.lastDimension();
	size_t o = m / Delta;
	TensorShape shape({o, N});
	Tensord mu(shape);
	for (size_t i = 0; i < N; ++i) {
		for (size_t l = 0; l < o; ++l) {
			for (size_t k = 0; k < Delta; ++k) {
				size_t idx = l * Delta + k;
				mu(l, i) += A(idx, i);
			}
		}
	}
	mu /= (double) Delta;
	return mu;
}

SOPcd meanVarianceAnalysis() {

	cout << "Mean Variance Analysis\n";
//	vector<string> names = {"BTC-USD.csv", "ETH-USD.csv", "SOL-USD.csv"};
	vector<string> names = {"BTC-USD.csv"};
	size_t N = names.size(); // number of assets
	size_t m = 180; // number of time steps
	size_t Delta = 30;

	/// read
	auto A = readAssets(names, m, Delta);

	auto dA = returns(A, Delta);
	auto mu_avg = avg(dA, Delta);
	dA.print();
	mu_avg.print();
/*	for (size_t j = 0; j < dA.shape_.lastBefore(); ++j) {
		for (size_t i = 0; i < dA.shape_.lastDimension(); ++i) {
			cout << dA(j, i) << "\t";
		}
		cout << "\n";
	}*/

	auto sig = covariance(dA, mu_avg, Delta);
	cout << "sigma:\n";
	sig.print();

/*	auto sig = sigma(dA, m, Delta);
	cout << "cov:\n";
	sig.print();
	getchar();*/
	getchar();

	/// How do you obtain mu, and sigma? (forecast return and covariance)
	/// For proof of concepts: 3 month until now on a daily basis

	return {SOPcd()};
}
