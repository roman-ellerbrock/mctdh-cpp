//
// Created by Roman Ellerbrock on 3/13/20.
//

#ifndef CDVR_H
#define CDVR_H
#include "DVR/TDDVR.h"
#include "DVR/Potential.h"
#include "DVR/CDVRNodeTensor.h"
#include "DVR/CDVREdgeTensor.h"

class CDVR {
public:
	CDVR(const Wavefunction& Psi, const Potential& V, const Tree& tree, size_t part = 0);
	~CDVR() = default;

	void Update(const Wavefunction& Psi, const Potential& V, const Tree& tree, size_t part = 0);

	TDDVR tddvr_;

private:
	TensorTreecd dvr_;
	MatrixTreed cdvr_;

	CDVRNodeTensor node_tensor_;
	CDVREdgeTensor edge_tensor_;
};

#endif //CDVR_H
