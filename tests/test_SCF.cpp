//
// Created by Hannes Hoppe on 23.06.22.
//

#include <gtest/gtest.h>
#include "TreeShape/LinearizedSCFNodes.h"
#include "TreeShape/TreeFactory.h"
#include "Core/Wavefunction.h"
#include "Core/SCF.h"

TEST(SCF, ShiftDensity){

    Tree tree = TreeFactory::balancedTree(8, 4, 11);

    mt19937 rng(187);
    Wavefunction test(rng, tree);

    LinearizedSCFNodes nodelist;
    nodelist.linearize(tree);

    for(const auto& i : nodelist.getAddresses()){
        std::cout << i << "\t";
    }
    std::cout << std::endl;

    SCF::SCF_ShiftTopNode(tree, test, nodelist.getAddresses()[0], nodelist.getAddresses()[1]);

}