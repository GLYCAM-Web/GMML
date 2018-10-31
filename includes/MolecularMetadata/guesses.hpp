#ifndef GMML_GUESS_META_HPP
#define GMML_GUESS_META_HPP

#include <string>
#include <iostream>
#include "../MolecularModeling/atom.hpp"

namespace Assembly
{
  bool guessIfC_CDoubleBond(MolecularModeling::Atom* carbon1, MolecularModeling::Atom* carbon2);
}

#endif //GMML_GUESS_META_HPP
