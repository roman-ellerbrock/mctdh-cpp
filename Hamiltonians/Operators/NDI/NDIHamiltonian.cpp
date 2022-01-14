//
// Created by Roman Ellerbrock on 12/12/21.
//

#include "NDIHamiltonian.h"

namespace Operator {
	std::string fix(const std::string& in) {
		return std::regex_replace(
			in,
			std::regex("(?:([a-zA-Z])([0-9]))|(?:([0-9])([a-zA-Z]))"),
			"\\1\\3 \\2\\4",
			std::regex_constants::format_sed
		);
	}

	void toPauli(MLOcd& M, string op, size_t t) {
		cout << op << t;
		if (op == "X") {
			M.push_back(sigma_x(), t);
		} else if (op == "Z") {
			M.push_back(sigma_z(), t);
		} else {
			cerr << "unknown operator for pauli conversion.\n";
			exit(1);
		}
	}

	void readLineMPO(SOPcd& S, string line, size_t nprod) {
		std::replace(line.begin(), line.end(), '*', ' ');
		stringstream ss(line);

		double coeff = 0;
		ss >> coeff;
		cout << coeff;

		line = line.erase(0, line.find(' '));
		line = fix(line);
		ss = stringstream(line);
		MLOcd M;
		for (size_t i = 0; i < nprod; ++i) {
			string op;
			ss >> op;
			size_t t;
			ss >> t;
			cout << " * ";
			toPauli(M, op, t);
		}
		for (size_t i = 0; i < nprod; ++i) {
//			if (M.mode(i) >= 13) {
			if (M.mode(i) >= 52) {
				cout << endl;
				return;
			}
		}
		cout << " in" << endl;
		S.push_back(M, coeff);
	}

	SOPcd ndi() {
		string filename = "pauli_52q.txt";
		ifstream file(filename);
		size_t n_single = 104;
		size_t n_double = 208;
		size_t nline = 1 + n_single + n_double;
		SOPcd S;
		string line;

		/// read & ignore energy shift
		getline(file, line);

		for (size_t i = 0; i < n_single; ++i) {
			getline(file, line);
			readLineMPO(S, line, 1);
		}
		for (size_t i = 0; i < n_double; ++i) {
			getline(file, line);
			readLineMPO(S, line, 2);
		}
		cout << "S.size() = " << S.size() << endl;
		return S;
	}
}
