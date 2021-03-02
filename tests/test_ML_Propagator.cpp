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

    double epsilon = 1.0e-8;

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
        results[1] = complex<double>(0.025921664708893683848334177355355,0);
        results[2] = complex<double>(-0.025049014898461941652785256451352,0);
        results[3] = complex<double>(0.96384759569531419653998227659031,0);

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
        results[0] = complex<double>(0.93367239357048759096358025999507,-0.35617727579872293652840653521707);
        results[1] = complex<double>(-0.0052598257300510362077794468405045,-0.0086540683605268054617853934473715);
        results[2] = complex<double>(-0.0064039889768260837774049853976521,-0.0075824365151104947471716322127122);
        results[3] = complex<double>(0.55576509419983544102450423451955,-0.77703882165397275549878486344824);

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
        results[0] = complex<double>(0.97907832931863114644954748655437,-6.9596758669730746045240590651538e-12);
        results[1] = complex<double>(0.04651552065395099522238098188609,2.1073603670097236746304525516633e-10);
        results[2] = complex<double>(-0.044671884846911209809849907514945,-2.2819198696866986788614857983512e-11);
        results[3] = complex<double>(0.80479785277064386228573766857153,4.266770312452388349724356281711e-09);

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

        // t_end = 5.0
        results[0] = complex<double>(0.80068279971538403039232889568666,-0.5136686897492928061126349348342);
        results[1] = complex<double>(-0.029098130507253391896504979285965,-0.013162063064847351334663017041748);
        results[2] = complex<double>(-0.028504119139673936983436774994516,-0.010467511030803920701881537524969);
        results[3] = complex<double>(-0.06859138237627232514981301392254,-0.58789395915735298547843967753579);

        double res = 0.0;
        for(int i = 0; i < 4; i++){
            res +=  abs(results[i] - a[i]);
        }
        res = sqrt(res);

        CHECK_CLOSE(0.0,res,epsilon);

    }

}