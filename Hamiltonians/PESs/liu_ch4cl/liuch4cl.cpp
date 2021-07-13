//
// Created by hoppe on 3/18/21.
//

#include "liuch4cl.h"

extern "C" {
void ch4clpot_(double* v, double* q, double* m);
void potini_();
}


liuch4cl::liuch4cl(vector<double> massvec, vector<int> coupling) {
    potini_();
    massvec_ = massvec;
    coupling_ = coupling;
}

double liuch4cl::evaluate(const Vectord &Xv, size_t part) const {

    double* v = new double;
    double* q = new double[12];
    double* m = new double[6];

    m[0] = massvec_[1];

    m[1] = massvec_[0];
    m[2] = massvec_[0];
    m[3] = massvec_[0];

    m[4] = massvec_[2];
    m[5] = massvec_[3];

    for(int i = 0; i < 12; i++){
        q[i] = Xv[i];
    }

    ch4clpot_(v,q,m);

    double vres = *v;

    delete[] q;
    delete[] m;
    delete v;

    return vres;
}
