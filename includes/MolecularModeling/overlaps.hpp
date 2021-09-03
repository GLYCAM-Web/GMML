#ifndef OVERLAPS_HPP
#define OVERLAPS_HPP

//#include "atom.hpp"
//#include "../common.hpp"
#include <vector>

// Wondering about forward declaring versus #including atom.hpp
namespace MolecularModeling
{
class Atom;
typedef std::vector<MolecularModeling::Atom*> AtomVector;
}

namespace gmml
{
    double CalculateAtomicOverlaps(MolecularModeling::AtomVector atomsA, MolecularModeling::AtomVector atomsB, bool print = false);
    double CalculateAtomicOverlapsBetweenNonBondedAtoms(MolecularModeling::AtomVector atomsA, MolecularModeling::AtomVector atomsB);
    double CalculateAtomicOverlaps(MolecularModeling::Atom *atomA, MolecularModeling::Atom *atomB, double radiusA = 0.0, double radiusB = 0.0);
}


#endif // OVERLAPS_HPP
