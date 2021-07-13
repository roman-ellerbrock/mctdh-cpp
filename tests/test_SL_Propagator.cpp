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

    double epsilon = 1.0e-8;

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

/*    TEST(SL_SOP_Imagtime){

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

        auto dotp = TreeFunctions::dotProduct(Psi_0,Psi_1,state.tree_);

        complex<double> a[4],results[4];
        auto overlap = dotp[state.tree_.topNode()];
        overlap.print();

        a[0] = overlap(0,0);
        a[1] = overlap(0,1);
        a[2] = overlap(1,0);
        a[3] = overlap(1,1);

        // t_end = 20.0
        results[0] = complex<double>(0.99782673944451860670312726142583,0.0);
        results[1] = complex<double>(-0.00292059226257011935942720093351,0.0);
        results[2] = complex<double>(-0.0042536967281209736232994167437482,0.0);
        results[3] = complex<double>(-0.82282252322424764834352117759408,0.0);

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

        auto dotp = TreeFunctions::dotProduct(Psi_0,Psi_1,state.tree_);

        complex<double> a[4],results[4];
        auto overlap = dotp[state.tree_.topNode()];
        a[0] = overlap(0,0);
        a[1] = overlap(0,1);
        a[2] = overlap(1,0);
        a[3] = overlap(1,1);

        // t_end = 20.0
        results[0] = complex<double>(0.74337671881260669604074564631446,-0.66266382020244574135858783847652);
        results[1] = complex<double>(-0.006076057492958681846562107864429,-0.0017908026181745130729627835108886);
        results[2] = complex<double>(-0.0062270631920599587780240646850416,-0.0017557967862557599409562003600627);
        results[3] = complex<double>(-0.44981486881960586199369345195009,-0.65525817408938291652731322756154);

        double res = 0.0;
        for(int i = 0; i < 4; i++){
            res +=  abs(results[i] - a[i]);
        }

        CHECK_CLOSE(0.0,res,epsilon);

    }

*/

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

/*    TEST(SL_nonSOP_Imagtime) {

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

        auto dotp = TreeFunctions::dotProduct(Psi_0,Psi_1,state.tree_);

        complex<double> a[4],results[4];
        auto overlap = dotp[state.tree_.topNode()];
        a[0] = overlap(0,0);
        a[1] = overlap(0,1);
        a[2] = overlap(1,0);
        a[3] = overlap(1,1);

        // t_end = 20.0
        results[0] = complex<double>(0.96310633346587604108890445786528,1.7404825275731834579784864262881e-13);
        results[1] = complex<double>(0.027050246709701773056844231746254,2.3527849348515721354041609327738e-10);
        results[2] = complex<double>(-0.020900540253599622642077804357541,-4.188288698007809908154716473135e-11);
        results[3] = complex<double>(0.57172440002811508730218292839709,6.9269861378086840470709865708217e-09);

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

        auto dotp = TreeFunctions::dotProduct(Psi_0,Psi_1,state.tree_);

        complex<double> a[4],results[4];
        auto overlap = dotp[state.tree_.topNode()];
        a[0] = overlap(0,0);
        a[1] = overlap(0,1);
        a[2] = overlap(1,0);
        a[3] = overlap(1,1);

        // t_end = 20.0
        results[0] = complex<double>(0.42250058101832194124014563385572,-0.77470742465834352419307151649264);
        results[1] = complex<double>(-0.024461582854740742537513753518397,0.0054377417926397188935316862057334);
        results[2] = complex<double>(-0.024221779958220276268709980627136,0.0061315012642777742823962938700788);
        results[3] = complex<double>(-0.15663395927497214543677728215698,0.037685101334036776032920812440352);

        double res = 0.0;
        for(int i = 0; i < 4; i++){
            res +=  abs(results[i] - a[i]);
        }
        res = sqrt(res);

        CHECK_CLOSE(0.0,res,epsilon);

    }
    */

}