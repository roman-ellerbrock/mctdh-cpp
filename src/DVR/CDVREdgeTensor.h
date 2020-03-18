//
// Created by Roman Ellerbrock on 3/18/20.
//

#ifndef CDVREDGETENSOR_H
#define CDVREDGETENSOR_H
#include "TreeClasses/NodeAttribute.h"
#include "TreeShape/Tree.h"

class CDVREdgeTensor : public NodeAttribute<Tensorcd>{
public:
	CDVREdgeTensor() = default;
	~CDVREdgeTensor() = default;

	void Initialize(const Tree& tree);

};


#endif //CDVREDGETENSOR_H
