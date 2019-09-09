#ifndef ATOM_HPP
#define ATOM_HPP

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "../GeometryTopology/coordinate.hpp"
#include "moleculardynamicatom.hpp"
#include "quantommechanicatom.hpp"
#include "dockingatom.hpp"
#include "oligosaccharidedetectionatom.hpp"

namespace MolecularModeling
{
    class Residue;
    class AtomNode;
    class Atom;
    typedef std::vector<MolecularModeling::Atom*> AtomVector; // Do this here so it can be used within the class, and is known
    class Atom : public MolecularDynamicAtom, public QuantomMechanicAtom, public DockingAtom, public OligoSaccharideDetectionAtom
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
            * List of Atoms
            */
            //typedef std::vector<Atom*> AtomVector;
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTORS                   //
            //////////////////////////////////////////////////////////
            /*! \fn
            * Default constructor
            */
            Atom();
            /*! \fn
            * Constructors
            */
            Atom(MolecularModeling::Residue* residue, std::string name, GeometryTopology::Coordinate::CoordinateVector coordinates);
            Atom(MolecularModeling::Residue* residue, std::string name, GeometryTopology::Coordinate coordinate);
            /*! \fn
            * Copy constructor(*)
            */
            Atom(const Atom* atom); // TODO This is redundant. The Copy Constructor below can handle this.
            /*! \fn
            * Copy constructor(&)
            */
            Atom(const Atom& atom);
            //////////////////////////////////////////////////////////
            //                       DESTRUCTOR                     //
            //////////////////////////////////////////////////////////
            // ~Atom();
            //////////////////////////////////////////////////////////
            //                       ACCESSORS                      //
            //////////////////////////////////////////////////////////
            /** \addtogroup Molecular_Data_Structure
            * @{
            */
            /*! \fn
            * An accessor function in order to access to the residue
            * @return residue_ attribute of the current object of this class
            */
            MolecularModeling::Residue* GetResidue() const;
            /*! \fn
            * An accessor function in order to access to the name
            * @return name_ attribute of the current object of this class
            */
            std::string GetName() const;
            /*! \fn
            * An accessor function in order to access to the coordinates
            * @return coordinates_ attribute of the current object of this class
            */
            GeometryTopology::Coordinate::CoordinateVector GetCoordinates() const;
            /*! \fn
            * An accessor function in order to access to the first coordinate
            * @return coordinates_.at(0) attribute of the current object of this class
            */
            GeometryTopology::Coordinate* GetCoordinate();
            /*! \fn
            * An accessor function in order to access to the chemical_type
            * @return chemical_type_ attribute of the current object of this class
            */
            std::string GetChemicalType() const;
            /*! \fn
            * An accessor function in order to access to the description
            * @return description_ attribute of the current object of this class
            */
            std::string GetDescription() const;
            /*! \fn
            * An accessor function in order to access to the element symbol
            * @return element_symbol_ attribute of the current object of this class
            */
            std::string GetElementSymbol() const;
            /*! \fn
            * An accessor function in order to access to the node
            * @return node_ attribute of the current object of this class
            */
            MolecularModeling::AtomNode* GetNode() const;
            /*! \fn
            * An accessor function in order to access to the id
            * @return id_ attribute of the current object of this class
            */
            std::string GetId() const;
            /*! \fn
            * An accessor function in order to access to the is_ring_ attribute of the current object
            * @return is_ring_ attribute of the current object of this class
            */
            bool GetIsRing() const;
            /*! \fn
            * An accessor function in order to access to the index
            * @return index_ attribute of the current object of this class
            */
            unsigned long long GetIndex() const;
            /*! \fn                                     //Added by ayush on 12/11/17 for molecules in assembly for atom type like O,H
            * An accessor function in order to access to the atom type
            * @return atom_type_ attribute of the current object of this class
            */
            std::string GetAtomType() const;
            /*! \fn                                     //Added by Dave on 03/23/18 for adding B-Factor into ontology
            * An accessor function in order to access to the atom b factor
            * @return b_factor_ attribute of the current object of this class
            */
            float GetBFactor() const;
            /** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /** \addtogroup Manipulators
            * @{
            */
            /*! \fn
            * A mutator function in order to set the residue of the current object
            * Set the residue_ attribute of the current atom
            * @param residue The residue attribute of the current object
            */
            void SetResidue(MolecularModeling::Residue* residue);
            /*! \fn
            * A mutator function in order to set the name of the current object
            * Set the name_ attribute of the current atom
            * @param name The name attribute of the current object
            */
            void SetName(std::string name);
            /*! \fn
            * A mutator function in order to set the coordinates of the current object
            * Set the coordinates_ attribute of the current atom
            * @param coordinates The coordinates attribute of the current object
            */
            void SetCoordinates(GeometryTopology::Coordinate::CoordinateVector coordinates);
            /*! \fn
            * A function in order to add the coordinate to the current object
            * Set the coordinates_ attribute of the current atom
            * @param coordinate The coordinate of the current object
            */
            void AddCoordinate(GeometryTopology::Coordinate* coordinate);
            /*! \fn
            * A mutator function in order to set the chemical type of the current object
            * Set the chemical_type_ attribute of the current atom
            * @param chemical_type The chemical type attribute of the current object
            */
            void SetChemicalType(std::string chemical_type);
            /*! \fn
            * A mutator function in order to set the description of the current object
            * Set the description_ attribute of the current atom
            * @param description The description attribute of the current object
            */
            void SetDescription(std::string description);
            /*! \fn
            * A mutator function in order to set the element symbol of the current object
            * Set the element_symbol_ attribute of the current atom
            * @param element_symbol The element symbol attribute of the current object
            */
            void SetElementSymbol(std::string element_symbol);
            /*! \fn
            * A mutator function in order to set the node of the current object
            * Set the node_ attribute of the current atom
            * @param node The node attribute of the current object
            */
            void SetNode(MolecularModeling::AtomNode* node);
            /*! \fn
            * A mutator function in order to set the id of the current object
            * Set the id_ attribute of the current atom
            * @param id The identification of the current object
            */
            void SetId(std::string id);
            /*! \fn
            * A mutator function in order to set the is_ring_ attribute of the current object
            * Set the is_ring_ attribute of the current atom
            * @param is_ring The boolean value representing if the current atom object is in aring or not
            */
            void SetIsRing(bool is_ring);
            //Added by ayush on 11/12/17 for molecules in assembly like atom type as O,C,H
            /*! \fn
            * A mutator function in order to set the atom type of the current object
            * Set the atom_type_ attribute of the current atom
            * @param atom_type The std::string stypes attribute of the current object
            */
            void SetAtomType(std::string atome_type);
            /*! \fn
            * A function to generate the index for an Atom.
            */
            unsigned long long generateAtomIndex();
            //Added by Dave on 03/23/18 for adding B Factor to ontology
            /*! \fn
            * A mutator function in order to set the b factor of the current object
            * Set the b_factor_ attribute of the current atom
            * @param b_factor The b factor attribute of the current object
            */
            void SetBFactor(float b_factor);
            /** @}*/
            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////
            /*! \fn
            * A function to find connected Atoms to this Atom.
            * @param visitedAtoms The AtomVector used to find connected atoms.
            */
            void FindConnectedAtoms(AtomVector& visitedStoms);
            /*! \fn
            * A function to aget the distance between this Atom and another.
            * @param otherAtoms The Atom used to get the distance from this Atom.
            */
            double GetDistanceToAtom(Atom* otherAtom);
            /*! \fn
            * A function to get the distance between this Atom's Coordinates
            * and another set of Coordinates.
            * @param coordinate The GeometryTopology::Coordinate used to get the distance from this Atom's Coordinates.
            */
            double GetDistanceToCoordinate(GeometryTopology::Coordinate* coordinate);


            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
            * A function to print out the atom contents in a structural format
            * Print out the information in a defined structure
            * @param out An output stream, the print result will be written in the given output stream
            */
            void Print(std::ostream& = std::cout); // @TODO DT - See the TODO below for operator<<
            //////////////////////////////////////////////////////////
            //                   OVERLOADED OPERATORS               //
            //////////////////////////////////////////////////////////
            // @TODO DT - Not sure if these are needed, but I am going to leave them
            // 				here for potential furthur discussion.
            //			The potential benefit would be to allow someone to make deep copies of an
            //				Atom object to an already initialized Atom object.
            // void operator=(const Atom&);
            // void operator=(const Atom*);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            MolecularModeling::Residue* residue_;                 /*!< A pointer back to the residue that this atom belongs to >*/
            std::string name_ = "";                 /*!< Name of the atom >*/
            GeometryTopology::Coordinate::CoordinateVector coordinates_;     /*!< Position of the atom >*/
            std::string chemical_type_;        /*!< A descriptor to determines the chemical type of the atom >*/
            std::string description_;          /*!< Short description of the atom >*/
            std::string element_symbol_;       /*!< Element symbol of the atom >*/
            MolecularModeling::AtomNode* node_;                   /*!< A Pointer to a node of the graph structure that indicates this atom >*/
            std::string id_;                   /*!< An identifier for an atom which is generated based on the type of the input file from which the structure has to be built
                    			Mostly it is like "residue_name:atom_name" >*/
            bool is_ring_;                     /*!< A boolean value which represents if an atom is involved in a sugar ring or not. This attribute is set during the Sugar ID process >*/
            unsigned long long index_;         /*!< A unqiue index for each atom in an assembly >*/
            std::string atom_type_;                /*!< List the atom type in an assembly >*/      //Added by ayush on 13/11/17 for molecules in assembly to set the atom type as an attribute like O,H
            float b_factor_;                  /*!< Gives the B Factor for the atom >*/            //Added by Dave on 03/23/18 for adding B Factor to ontology

            //////////////////////////////////////////////////////////
            //                   HELPER FUNCTIONS                   //
            //////////////////////////////////////////////////////////
            void Copy(const Atom* atom);
            void SetAttributes( MolecularModeling::Residue* residue, std::string name, GeometryTopology::Coordinate::CoordinateVector coordinates,
                                std::string chemical_type, std::string description, std::string element_symbol,
                                AtomNode* atomnode, std::string id, bool is_ring, std::string atom_type);
    };
    // @TODO DT - Get this working. For some reason it causes a Seg Fault.
    //			Ideally this function would be where the printing code is and the
    //			class Print function would call this operator.
    //			This would allow someone to do: std::cout << AtomObject << std::endl;
    // std::ostream& operator<<( std::ostream& out, const Atom& atom );
}
#endif // ATOM_HPP
