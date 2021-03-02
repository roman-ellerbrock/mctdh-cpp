//
// Created by hoppe on 21.01.21.
//

#include "UnitTest++/UnitTest++.h"
#include "Parser/yaml_parser.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "Util/QMConstants.h"
#include "Core/Eigenstates.h"
#include <iomanip>

SUITE(SingleLayer){

    double epsilon = 1.0e-14;

    class SOPParser {
    public:
        SOPParser() {
            string yaml_filename("../tests/sl_tests/SOPPot.yaml");
            state_ = parser::run(yaml_filename);

            // take care that tested WF is always the same
            mt19937 gen(0);
            state_.wavefunctions_["Psi"] = Wavefunction(gen, state_.tree_);
        }
        mctdh_state state_;
    };

    TEST(SL_SOP_Imagtime){

        auto SOPtest = SOPParser();

        auto state = SOPtest.state_;
        auto Psi_0 = state.wavefunctions_["Psi"];

        // prepare integration data
        auto t_end = 20.0;
        auto t = 0.0;
        auto out = 20.0;
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

        // TODO: hack: suppress trashy output
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

        // t_end = 20.0
        results[0] = complex<double>(0.99782673944451849568082479891018,0.0);
        results[1] = complex<double>(-0.002920592262570642378555207940849,0.0);
        results[2] = complex<double>(-0.0042536967281211453609235384476506,0.0);
        results[3] = complex<double>(-0.82282252322424753732121871507843,0.0);

        double res = 0.0;
        for(int i = 0; i < 4; i++){
            res +=  abs(results[i] - a[i]);
        }

        CHECK_CLOSE(0.0,res,epsilon);

    }

    TEST(SL_SOP_Realtime){

        auto SOPtest = SOPParser();

        auto state = SOPtest.state_;
        auto Psi_0 = state.wavefunctions_["Psi"];

        // prepare integration data
        auto t_end = 20.0;
        auto t = 0.0;
        auto out = 20.0;
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

        // TODO: hack: suppress trashy output
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

        // t_end = 20.0
        results[0] = complex<double>(0.74337671881260825035298012153362,-0.66266382020244685158161246363306);
        results[1] = complex<double>(-0.0060760574929585153131084140909479,-0.0017908026181746264805100254946524);
        results[2] = complex<double>(-0.0062270631920600203607074618616934,-0.0017557967862555569783095110736326);
        results[3] = complex<double>(-0.44981486881960081047893140748783,-0.6552581740893791417690295020293);

        double res = 0.0;
        for(int i = 0; i < 4; i++){
            res +=  abs(results[i] - a[i]);
        }

        CHECK_CLOSE(0.0,res,epsilon);

    }



    class nonSOPParser {
    public:
        nonSOPParser() {
            string yaml_filename("../tests/sl_tests/nonSOPPot.yaml");
            state_ = parser::run(yaml_filename);

            // take care that tested WF is always the same
            mt19937 gen(0);
            state_.wavefunctions_["Psi"] = Wavefunction(gen, state_.tree_);

        }
        mctdh_state state_;
    };

    TEST(SL_nonSOP_Imagtime) {

        auto nonSOPtest = nonSOPParser();

        auto state = nonSOPtest.state_;
        auto Psi_0 = state.wavefunctions_["Psi"];

        // prepare integration data
        auto t_end = 20.0;
        auto t = 0.0;
        auto out = 20.0;
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

        // TODO: hack: suppress trashy output
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

        // t_end = 20.0
        results[0] = complex<double>(0.96310633347673291204671386367409,9.8879973006993237837463860291818e-13);
        results[1] = complex<double>(0.027050246357936080093820052638875,4.2725299763902857501994016384918e-11);
        results[2] = complex<double>(-0.02090054004021815389768867987641,6.4981501135331774613230088195984e-12);
        results[3] = complex<double>(0.57172440004303592964163271972211,1.1185265778611335809036006334297e-09);

        double res = 0.0;
        for(int i = 0; i < 4; i++){
            res +=  abs(results[i] - a[i]);
        }
        res = sqrt(res);

        CHECK_CLOSE(0.0,res,epsilon);

    }

    TEST(SL_nonSOP_Realtime) {

        auto nonSOPtest = nonSOPParser();

        auto state = nonSOPtest.state_;
        auto Psi_0 = state.wavefunctions_["Psi"];

        // prepare integration data
        auto t_end = 20.0;
        auto t = 0.0;
        auto out = 20.0;
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

        // TODO: hack: suppress trashy output
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

        // t_end = 20.0
        results[0] = complex<double>(0.42250058101621512651746570554678,-0.77470742466056263797469227938564);
        results[1] = complex<double>(-0.024461582851163003737848811169897,0.0054377417914377342011977845004367);
        results[2] = complex<double>(-0.024221779958592447312692996774786,0.0061315012641119295139224831814317);
        results[3] = complex<double>(-0.15663395927903642212974943959125,0.037685101329096380717853520536664);

        double res = 0.0;
        for(int i = 0; i < 4; i++){
            res +=  abs(results[i] - a[i]);
        }
        res = sqrt(res);

        CHECK_CLOSE(0.0,res,epsilon);

    }

}