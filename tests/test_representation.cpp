//
// Created by hoppe on 21.01.21.
//

#include "UnitTest++/UnitTest++.h"
#include "Parser/yaml_parser.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "Util/QMConstants.h"
#include "Core/Eigenstates.h"
#include <iomanip>

SUITE(OperatorRepresentation){

    double epsilon = 1.0e-16;

    class ORParser{
    public:
        ORParser() = default;
        explicit ORParser(const int type){
            string a;
            if(type == 0){
                // multi-layer
                a = "../tests/representation_tests/ML.yaml";
            }else if(type == 2){
                // single-layer
                a = "../tests/representation_tests/SL.yaml";
            }

            state_ = parser::run(a);
            mt19937 gen(0);
            state_.wavefunctions_["Psi"] = Wavefunction(gen, state_.tree_);

        }
        mctdh_state state_;
    };

    TEST(SL_HEXPECTATION){

        auto SORtest = ORParser(2);
        auto& Hamiltonian = *SORtest.state_.hamiltonian_;
        auto Psi = SORtest.state_.wavefunctions_["Psi"];
        auto tree = SORtest.state_.tree_;

        Wavefunction dPsi(tree);
        HamiltonianRepresentation HRep(Hamiltonian, tree);
        HRep.build(Hamiltonian, Psi, tree, 0.);

        auto eval = Expectation(HRep, Psi, Hamiltonian, tree);

        complex<double> compare[4],a[4];
        a[0] = eval(0,0);
        a[1] = eval(0,1);
        a[2] = eval(1,0);
        a[3] = eval(1,1);
        compare[0] = complex<double>(0.05468025458072921,0);
        compare[1] = complex<double>(-1.043262566935352e-05,0);
        compare[2] = complex<double>(-1.043262566935291e-05,0);
        compare[3] = complex<double>(0.1668074591154756,0);

        double diff = 0.0;
        for(size_t i = 0; i < 4; i++){
            diff += std::pow(abs(compare[0]-a[0]),2.0);
        }

        CHECK_CLOSE(0.0,diff,epsilon);

    }

    TEST(ML_HEXPECTATION){

        auto SORtest = ORParser(0);
        auto& Hamiltonian = *SORtest.state_.hamiltonian_;
        auto Psi = SORtest.state_.wavefunctions_["Psi"];
        auto tree = SORtest.state_.tree_;

        HamiltonianRepresentation HRep(Hamiltonian, tree);
        HRep.build(Hamiltonian, Psi, tree, 0.);

        auto eval = Expectation(HRep, Psi, Hamiltonian, tree);

        complex<double> compare[4],a[4];
        a[0] = eval(0,0);
        a[1] = eval(0,1);
        a[2] = eval(1,0);
        a[3] = eval(1,1);
        compare[0] = complex<double>(0.05468025458072921,0);
        compare[1] = complex<double>(8.516403493066244e-06);
        compare[2] = complex<double>(8.516403493068843e-06,0);
        compare[3] = complex<double>(0.139836057134736,0);

        double diff = 0.0;
        for(size_t i = 0; i < 4; i++){
            diff += std::pow(abs(compare[0]-a[0]),2.0);
        }

        CHECK_CLOSE(0.0,diff,epsilon);

    }
}