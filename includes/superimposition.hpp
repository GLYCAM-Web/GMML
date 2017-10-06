#ifndef SUPERIMPOSITION_HPP
#define SUPERIMPOSITION_HPP

#include "Eigen_Algebra_Template_Library/Geometry"
#include "GeometryTopology/coordinate.hpp"
#include "MolecularModeling/assembly.hpp"
#include "MolecularModeling/atom.hpp"

//*******************************************

typedef std::vector<MolecularModeling::Atom*> AtomVector;
typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;
typedef std::vector<MolecularModeling::Assembly*> AssemblyVector;

//*******************************************

using namespace MolecularModeling;

namespace gmml
{
    void GenerateMatrixFromAssembyCoordinates(Assembly *assembly, Eigen::Matrix3Xd *matrix);

    void ReplaceAssemblyCoordinatesFromMatrix(Assembly *assembly, Eigen::Matrix3Xd *matrix);
    //Atom Vector version, may be removed
    void GenerateMatrixFromAtomVectorCoordinates(AtomVector *atoms, Eigen::Matrix3Xd *matrix);

    void ReplaceAtomVectorCoordinatesFromMatrix(AtomVector *atoms, Eigen::Matrix3Xd *matrix);

    Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out);

    void Superimpose(AtomVector moving, AtomVector target);

    void Superimpose(AtomVector moving, AtomVector target, AtomVector alsoMoving);

    void Superimpose(Assembly *moving, Assembly *target);

    void Superimpose(Assembly *moving, Assembly *target, Assembly *alsoMoving);

    void Superimpose(Assembly *moving, Assembly *target, AssemblyVector *alsoMoving);
    // A function to test Find3DAffineTransform()
    // void TestFind3DAffineTransform();
}


#endif // SUPERIMPOSITION_HPP
