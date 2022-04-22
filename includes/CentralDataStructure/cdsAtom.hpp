#ifndef INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP

#include <string>
#include <iostream>
#include <vector>

#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
//#include "includes/CentralDataStructure/cdsCoordinate.hpp"
#include "includes/GeometryTopology/coordinate.hpp"
using GeometryTopology::Coordinate;
namespace cds
{
class cdsAtom : public glygraph::Node<cdsAtom>
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTORS                   //
    cdsAtom() {}
    cdsAtom(const std::string& name, const Coordinate& coord);
    //////////////////////////////////////////////////////////
    //                       ACCESSORS                      //
    //////////////////////////////////////////////////////////
    Coordinate* getCoordinate();
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    void setCoordinate(const Coordinate& c);
    void addCoordinate(const Coordinate& c);
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
private:
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    std::vector<std::unique_ptr<Coordinate>> coordinates_;     /*!< Position of the atom >*/
};
}
#endif // ATOM_HPP
