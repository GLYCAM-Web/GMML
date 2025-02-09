#ifndef OFFFILEATOM_HPP
#define OFFFILEATOM_HPP

#include <string>
#include <vector>
#include <iostream>
#include "../../GeometryTopology/coordinate.hpp"

namespace OffFileSpace
{
    class OffFileAtom
    {
      public:
        //////////////////////////////////////////////////////////
        //                       Constructor                    //
        //////////////////////////////////////////////////////////
        /*! \fn
         * Default constructor
         */
        OffFileAtom();

        //////////////////////////////////////////////////////////
        //                           ACCESSOR                   //
        //////////////////////////////////////////////////////////

        /*! \fn
         * An accessor function in order to access to atom name of the current object
         * @return  name_ attribute of the current atom
         */
        std::string GetName();
        /*! \fn
         * An accessor function in order to access to atom type of the current object
         * @return type_ attribute of the current atom
         */
        std::string GetType();
        /*! \fn
         * An accessor function in order to access to residue index of the current object
         * @return residue_index_ attribute of the current atom
         */
        int GetResidueIndex();
        /*! \fn
         * An accessor function in order to access to atom index of the current object
         * @return atom_index_ attribute of the current atom
         */
        int GetAtomIndex();
        /*! \fn
         * An accessor function in order to access to atomic number of the current object
         * @return atomic_numebr_ attribute of the current atom
         */
        int GetAtomicNumber();
        /*! \fn
         * An accessor function in order to access to charge of the current object
         * @return charge_ attribute of the current atom
         */
        double GetCharge();
        /*! \fn
         * An accessor function in order to access to coordinate of the current object
         * @return coordinate_ attribute of the current atom
         */
        GeometryTopology::Coordinate GetCoordinate();
        /*! \fn
         * An accessor function in order to access to the list of atom indices which are bonded to the current object
         * @return bonded_atom_indices attribute of the current atom
         */
        std::vector<int> GetBondedAtomsIndices();
        /*! \fn
         * An accessor function in order to access to atom order of the current object
         * @return atom_order_ attribute of the current atom
         */
        int GetAtomOrder();

        //////////////////////////////////////////////////////////
        //                           MUTATOR                    //
        //////////////////////////////////////////////////////////

        /*! \fn
         * A mutator function in order to set the atom name of the current object
         * Set the name_ attribute of the current atom
         * @param name Atom name
         */
        void SetName(const std::string name);

        /*! \fn
         * A mutator function in order to set the atom type of the current object
         * Set the type_ attribute of the current atom
         * @param type Atom type
         */
        void SetType(const std::string type);

        /*! \fn
         * A mutator function in order to set the residue index of the current object
         * Set the residue_index_ attribute of the current atom
         * @param residue_index The index of the residue that the current object belongs to
         */
        void SetResidueIndex(int residue_index);
        /*! \fn
         * A mutator function in order to set the atom index of the current object
         * Set the atom_index__ attribute of the current atom
         * @param atom_index The index of the current object in the belonging residue
         */
        void SetAtomIndex(int atom_index);
        /*! \fn
         * A mutator function in order to set the atomic number of the current object
         * Set the atomic_number_ attribute of the current atom
         * @param atomic_number Atomic number of the current object
         */
        void SetAtomicNumber(int atomic_number);
        /*! \fn
         * A mutator function in order to set the coordinate of the current object
         * Set the coordinate_ attribute of the current atom
         * @param coordinate The position of the current object in the cartesian coordinate
         */
        void SetCoordinate(GeometryTopology::Coordinate* coordinate);
        /*! \fn
         * A mutator function in order to set the list of bonded atom indices of the current object
         * Set the bonded_atom_indices_ attribute of the current atom
         * @param bonded_atom_indices a list of indices of atoms which are bonded to the current object
         */
        void SetBondedAtomsIndices(const std::vector<int> bonded_atoms_indices);
        /*! \fn
         * A mutator function in order to add a new atom index into the list of bonded atom indices of the current
         * object Add a new index to the bonded_atom_indices_ attribute of the current atom
         * @param index Index of the new atom which is connected to the current object
         */
        void AddBondedAtomIndex(int index);
        /*! \fn
         * A mutator function in order to set the atom order of the current object
         * Set the atom_order_ attribute of the current atom
         * @param atom_order_ written order of the current object in a library file
         */
        void SetAtomOrder(int atom_order);

        /*! \fn
         * A mutator function in order to set the atom charge of the current object
         * Set the charge_ attribute of the current atom
         * @param charge charge of the atom in the off file
         */
        void SetAtomCharge(double charge);

        //////////////////////////////////////////////////////////
        //                     DISPLAY FUNCTIONS                //
        //////////////////////////////////////////////////////////
        /*! \fn
         * A function to print out the current object information in a structural format
         * Print out the information in a defined structure
         * @param out An output stream, the print result will be written in the given output stream
         */
        void Print(std::ostream& out = std::cerr);

      private:
        //////////////////////////////////////////////////////////
        //                         ATTRIBUTES                   //
        //////////////////////////////////////////////////////////
        std::string name_ = ""; /*!< Atom name: Set by the 1st column of the atom section of a off file */
        std::string type_ = ""; /*!< Atom type; Set by the 2nd column of the atom section of a off file */
        int residue_index_ =
            (-1); /*!< Residue index that the atom belongs to; Set by the 4th column of the atom section of a off file*/
        int atom_index_ = (-1); /*!< Index of the atom in the belonging residue; Set by the 6th column of the atom
                                   section of a off file */
        int atomic_number_ =
            (0); /*!< Atomic number of the atom; Set by the 7th column of the atom section of a off file */
        double charge_ = (0.0); /*!< Charge of the atom; Set by the 8th column of the atom section of a off file */
        GeometryTopology::Coordinate coordinate_; /*!< Coordinate of the atom; Set by the corresponding line in the
                                                     positions section of a off file */
        std::vector<int>
            bonded_atoms_indices_; /*!< List of atom indices that are bonded to the atom;
                               Set by the corresponding atom index in the connectivity section of a off file */
        int atom_order_ = (0); /*!< Order of atoms of the residue in the atom section; Set by a line counter iterates on
                                  lines of the atom section of a off file */
    };
} // namespace OffFileSpace

#endif // OFFFILEATOM_HPP
