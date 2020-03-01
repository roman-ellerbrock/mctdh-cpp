//
// Created by Roman Ellerbrock on 2/27/20.
//

#ifndef YAML_PARSER_H
#define YAML_PARSER_H
#include "Parser/mctdh_state.h"

namespace parser {

	mctdh_state read(const string& yaml_filename);

	mctdh_state run(const string& yaml_filename);
}

#endif //YAML_PARSER_H
