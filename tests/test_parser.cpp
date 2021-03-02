//
// Created by hoppe on 26.01.21.
//

#include "UnitTest++/UnitTest++.h"
#include "Parser/yaml_parser.h"
#include "TreeClasses/MatrixTreeFunctions.h"

SUITE(Parser){
    double epsilon = 1e-16;

    // Write |Psi1_> to file, read into |Psi2_>
    // then |diff> = |Psi_1> - |Psi_2>
    // if <diff|diff> > epsilon, then reading and/or writing was not successful
    TEST(Parser_SL_ReadWrite){

        class UTParser {
        public:
            UTParser() {
                string yaml_filename("../tests/parser_tests/sl_read_write.yaml");
                state_ = parser::run(yaml_filename);
                Psi1_ = state_.wavefunctions_["Psi1"];
                Psi2_ = state_.wavefunctions_["Psi2"];
            }
            mctdh_state state_;
            Wavefunction Psi1_, Psi2_;
        };

        auto testsys = UTParser();

        // difference of the read and written file
        auto diff = testsys.Psi1_ - testsys.Psi2_;

        // dot product is tested by qutree, <x|x> is real and >= 0 for all |x> in the hilbert space
        // therefore abs(<diff|diff>) = <diff|diff>
        const Tree& tree = testsys.state_.tree_;
        auto dotp = TreeFunctions::DotProduct(diff,diff,tree);

        double res = 0;
        for(size_t i = 0; i < tree.TopNode().shape().size(); i++){
            res += abs(dotp[tree.TopNode()](i,i));
        }

        CHECK_CLOSE(res,0,epsilon);
    }

    TEST(Parser_ML_ReadWrite){

        class UTParser {
        public:
            UTParser() {
                string yaml_filename("../tests/parser_tests/ml_read_write.yaml");
                state_ = parser::run(yaml_filename);
                Psi1_ = state_.wavefunctions_["Psi1"];
                Psi2_ = state_.wavefunctions_["Psi2"];
            }
            mctdh_state state_;
            Wavefunction Psi1_, Psi2_;
        };

        auto testsys = UTParser();

        // difference of the read and written file
        auto diff = testsys.Psi1_ - testsys.Psi2_;

        // dot product is tested by qutree, <x|x> is real and >= 0 for all |x> in the hilbert space
        // therefore abs(<diff|diff>) = <diff|diff>
        const Tree& tree = testsys.state_.tree_;
        auto dotp = TreeFunctions::DotProduct(diff,diff,tree);

        double res = 0;
        for(size_t i = 0; i < tree.TopNode().shape().size(); i++){
            res += abs(dotp[tree.TopNode()](i,i));
        }

        CHECK_CLOSE(res,0,epsilon);
    }

}

