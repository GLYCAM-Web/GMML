// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBMODEL_HPP
#define PDBMODEL_HPP

#include <string>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbModelResidueSet;

    class PdbModel
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////


            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbModel();
            /*! \fn
              * Constructor with required parameters
              * @param model_block
              */
            PdbModel(std::stringstream& model_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the model serial number in a pdb model
              * @return model_serial_number_ attribute of the current object of this class
              */
            int GetModelSerialNumber();
            /*! \fn
              * An accessor function in order to access to the model residue set in a pdb model
              * @return model_residue_set_ attribute of the current object of this class
              */
            PdbModelResidueSet* GetModelResidueSet();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the model serial number of the current object
              * Set the model_serial_number_ attribute of the current pdb model
              * @param model_serial_number The model serial number attribute of the current object
              */
            void SetModelSerialNumber(int model_serial_number);
            /*! \fn
              * A mutator function in order to set the model residue set of the current object
              * Set the model_residue_set_ attribute of the current pdb model
              * @param model_residue_set The model residue set attribute of the current object
              */
            void SetModelResidueSet(PdbModelResidueSet* model_residue_set);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the model contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int model_serial_number_;               /*!< Serial number of a model >*/
            PdbModelResidueSet* model_residue_set_; /*!< Residue sets involving in a model >*/

    };
}

#endif // PDBMODEL_HPP
