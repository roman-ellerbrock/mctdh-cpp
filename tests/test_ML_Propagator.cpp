//
// Created by hoppe on 21.01.21.
//

#include "UnitTest++/UnitTest++.h"
#include "Parser/yaml_parser.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "Util/QMConstants.h"
#include "Core/Eigenstates.h"

#include <iomanip>

SUITE(MultiLayer){

    double epsilon = 1.0e-14;

    class SOPParser {
    public:
        SOPParser() {
            string yaml_filename("../tests/ml_tests/SOPPot.yaml");
            state_ = parser::run(yaml_filename);

            // take care that tested WF is always the same
            mt19937 gen(0);
            state_.wavefunctions_["Psi"] = Wavefunction(gen, state_.tree_);
        }
        mctdh_state state_;
    };

    TEST(ML_SOP_Imagtime){

        auto SOPtest = SOPParser();

        auto state = SOPtest.state_;
        auto Psi_0 = state.wavefunctions_["Psi"];

        // prepare integration data
        auto t_end = 10.0;
        auto t = 0.0;
        auto out = 10.0;
        auto dt = 0.1;
        auto cmf = 1e-8;
        auto bs = 1e-8;
        auto file_in = "in.dat";
        auto file_out = "out.dat";
        IntegratorVariables ivar(t, t_end, dt, out, cmf, bs,
                                 state.wavefunctions_["Psi"],
                                 *state.hamiltonian_,
                                 state.tree_,
                                 state.cdvrtree_,
                                 file_out,
                                 file_in,
                                 false);

        // TODO: it would be nice if Eigenstates would have a flag where it does not output anything (or an empty stream to pass to)
        // TODO: all its stuff. maybe an ivar variable?

        // TODO: hack: surpress trsahy output
        cout.setstate(ios_base::badbit);
        Eigenstates(ivar);
        cout.clear();

        // check autocorrelation of propagated state with start state
        // least squares difference of dot products
        auto Psi_1 = *ivar.psi;

        auto dotp = TreeFunctions::DotProduct(Psi_0,Psi_1,state.tree_);

        complex<double> a[4],results[4];
        auto overlap = dotp[state.tree_.TopNode()];

        a[0] = overlap(0,0);
        a[1] = overlap(0,1);
        a[2] = overlap(1,0);
        a[3] = overlap(1,1);

        // t_end = 10.0
        results[0] = complex<double>(0.99923124502831084381426762774936,0);
        results[1] = complex<double>(0.025921664708893770584507976195709,0);
        results[2] = complex<double>(-0.025049014898462007572277343570022,0);
        results[3] = complex<double>(0.96384759569531430756228473910596,0);

        double res = 0.0;
        for(int i = 0; i < 4; i++){
            res +=  abs(results[i] - a[i]);
        }
        res = sqrt(res);

        CHECK_CLOSE(0.0,res,epsilon);

    }

    TEST(ML_SOP_Realtime){

        auto SOPtest = SOPParser();

        auto state = SOPtest.state_;
        auto Psi_0 = state.wavefunctions_["Psi"];

        // prepare integration data
        auto t_end = 10.0;
        auto t = 0.0;
        auto out = 10.0;
        auto dt = 0.1;
        auto cmf = 1e-8;
        auto bs = 1e-8;
        auto file_in = "in.dat";
        auto file_out = "out.dat";
        IntegratorVariables ivar(t, t_end, dt, out, cmf, bs,
                                 state.wavefunctions_["Psi"],
                                 *state.hamiltonian_,
                                 state.tree_,
                                 state.cdvrtree_,
                                 file_out,
                                 file_in,
                                 false);

        // TODO: it would be nice if Eigenstates would have a flag where it does not output anything (or an empty stream to pass to)
        // TODO: all its stuff. maybe an ivar variable?

        // TODO: hack: surpress trashy output
        cout.setstate(ios_base::badbit);
        CMFIntegrator propagator(*state.hamiltonian_,state.tree_,state.cdvrtree_,1.0);
        propagator.Integrate(ivar);
        cout.clear();

        // check autocorrelation of propagated state with start state
        // least squares difference of dot products
        auto Psi_1 = *ivar.psi;

        auto dotp = TreeFunctions::DotProduct(Psi_0,Psi_1,state.tree_);

        complex<double> a[4],results[4];
        auto overlap = dotp[state.tree_.TopNode()];

        a[0] = overlap(0,0);
        a[1] = overlap(0,1);
        a[2] = overlap(1,0);
        a[3] = overlap(1,1);

        // t_end = 10.0
        results[0] = complex<double>(0.933672393570491143677259060496,-0.35617727579872449084064101043623);
        results[1] = complex<double>(-0.0052598257300509321243708882320789,-0.0086540683605269667910686592904312);
        results[2] = complex<double>(-0.0064039889768262251573682775074303,-0.0075824365151104843388307763518696);
        results[3] = complex<double>(0.55576509419983610715831900961348,-0.77703882165397331061029717602651);

        double res = 0.0;
        for(int i = 0; i < 4; i++){
            res +=  abs(results[i] - a[i]);
        }
        res = sqrt(res);

        CHECK_CLOSE(0.0,res,epsilon);

    }



    class nonSOPParser {
    public:
        nonSOPParser() {
            string yaml_filename("../tests/ml_tests/nonSOPPot.yaml");
            state_ = parser::run(yaml_filename);

            // take care that tested WF is always the same
            mt19937 gen(0);
            state_.wavefunctions_["Psi"] = Wavefunction(gen, state_.tree_);

        }
        mctdh_state state_;
    };

    TEST(ML_nonSOP_Imagtime) {

        auto nonSOPtest = nonSOPParser();

        auto state = nonSOPtest.state_;
        auto Psi_0 = state.wavefunctions_["Psi"];

        // prepare integration data
        auto t_end = 10.0;
        auto t = 0.0;
        auto out = 10.0;
        auto dt = 0.1;
        auto cmf = 1e-8;
        auto bs = 1e-8;
        auto file_in = "in.dat";
        auto file_out = "out.dat";
        IntegratorVariables ivar(t, t_end, dt, out, cmf, bs,
                                 state.wavefunctions_["Psi"],
                                 *state.hamiltonian_,
                                 state.tree_,
                                 state.cdvrtree_,
                                 file_out,
                                 file_in,
                                 false);

        // TODO: it would be nice if Eigenstates would have a flag where it does not output anything (or an empty stream to pass to)
        // TODO: all its stuff. maybe an ivar variable?

        // TODO: hack: surpress trsahy output
        cout.setstate(ios_base::badbit);
        Eigenstates(ivar);
        cout.clear();

        // check autocorrelation of propagated state with start state
        // least squares difference of dot products
        auto Psi_1 = *ivar.psi;

        auto dotp = TreeFunctions::DotProduct(Psi_0,Psi_1,state.tree_);

        complex<double> a[4],results[4];
        auto overlap = dotp[state.tree_.TopNode()];

        a[0] = overlap(0,0);
        a[1] = overlap(0,1);
        a[2] = overlap(1,0);
        a[3] = overlap(1,1);

        // t_end = 10.0
        results[0] = complex<double>(0.97907832931776717089178418973461,-6.6801110465371377112252487923959e-12);
        results[1] = complex<double>(0.046515520681890153154025568937868,2.2288332582213727655871335087958e-10);
        results[2] = complex<double>(-0.044671884858588035982496222686677,-2.5614961364989316161967116637173e-11);
        results[3] = complex<double>(0.80479785281357729687101709714625,4.5505039019537538243650751719656e-09);

        double res = 0.0;
        for(int i = 0; i < 4; i++){
            res +=  abs(results[i] - a[i]);
        }
        res = sqrt(res);

        CHECK_CLOSE(0.0,res,epsilon);

    }

    TEST(ML_nonSOP_Realtime) {

        auto nonSOPtest = nonSOPParser();

        auto state = nonSOPtest.state_;
        auto Psi_0 = state.wavefunctions_["Psi"];

        // prepare integration data
        auto t_end = 10.0;
        auto t = 0.0;
        auto out = 10.0;
        auto dt = 0.1;
        auto cmf = 1e-8;
        auto bs = 1e-8;
        auto file_in = "in.dat";
        auto file_out = "out.dat";
        IntegratorVariables ivar(t, t_end, dt, out, cmf, bs,
                                 state.wavefunctions_["Psi"],
                                 *state.hamiltonian_,
                                 state.tree_,
                                 state.cdvrtree_,
                                 file_out,
                                 file_in,
                                 false);

        // TODO: it would be nice if Eigenstates would have a flag where it does not output anything (or an empty stream to pass to)
        // TODO: all its stuff. maybe an ivar variable?

        // TODO: hack: surpress trashy output
        cout.setstate(ios_base::badbit);
        CMFIntegrator propagator(*state.hamiltonian_,state.tree_,state.cdvrtree_,1.0);
        propagator.Integrate(ivar);
        cout.clear();

        // check autocorrelation of propagated state with start state
        // least squares difference of dot products
        auto Psi_1 = *ivar.psi;

        auto dotp = TreeFunctions::DotProduct(Psi_0,Psi_1,state.tree_);

        complex<double> a[4],results[4];
        auto overlap = dotp[state.tree_.TopNode()];

        a[0] = overlap(0,0);
        a[1] = overlap(0,1);
        a[2] = overlap(1,0);
        a[3] = overlap(1,1);

        // t_end = 10.0
        results[0] = complex<double>(0.80068279971679834350339888260351,-0.513668689748916440507287006767);
        results[1] = complex<double>(-0.029098130507602695815627669162495,-0.013162063064506460824398814679626);
        results[2] = complex<double>(-0.028504119139095940999029110685115,-0.010467511031466964974145916755788);
        results[3] = complex<double>(-0.068591382377383575130735948732763,-0.58789395915707576278919077594765);

        double res = 0.0;
        for(int i = 0; i < 4; i++){
            res +=  abs(results[i] - a[i]);
        }
        res = sqrt(res);

        CHECK_CLOSE(0.0,res,epsilon);

    }

}