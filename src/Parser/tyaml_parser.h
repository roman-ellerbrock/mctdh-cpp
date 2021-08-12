//
// Created by Roman Ellerbrock on 2/27/20.
//

#ifndef YAML_PARSER_H
#define YAML_PARSER_H
#include "Parser/tmctdh_state.h"
#include "Core/tIntegratorVariables.h"
#include "yaml-cpp/yaml.h"

namespace tparser {

	template<typename T>
	tmctdh_state<T> run(const string& yaml_filename);

	Tree create_tree(const YAML::Node& node);

}

#endif //YAML_PARSER_H
