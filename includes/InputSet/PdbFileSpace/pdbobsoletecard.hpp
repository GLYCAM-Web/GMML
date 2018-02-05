// Created by: Dave Montgomery

#ifndef PDBOBSOLETECARD_HPP
#define PDBOBSOLETECARD_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbObsoleteCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbObsoleteCard();
            /*! \fn
              * Constructor with required parameters
              * @param replacement_date Replacement Date
              * @param replacement_entries Replacement Entries
              */
            PdbObsoleteCard(const std::string& replacement_date, const std::vector<std::string>& replacement_entries);
            /*! \fn
              * Constructor with required parameters
              * @param specification_block
              */
            PdbObsoleteCard(std::stringstream& obsolete_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */

            /*! \fn
              * An accessor function in order to access to the replacement date in a pdb obsolete
              * @return replacement_date_ attribute of the current object of this class
              */
            std::string GetReplacementDate();
            /*! \fn
              * An accessor function in order to access to the replacement entries in a pdb obsolete
              * @return replacement_entries_ attribute of the current object of this class
              */
            std::vector<std::string> GetReplacementEntries();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //

           //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the replacement date of the current object
              * Set the replacement_date_ attribute of the current obsolete
              * @param replacement_date The replacement date of the current object
              */
            void SetReplacementDate(const std::string replacement_date);
            /*! \fn
              * A mutator function in order to set the replacement entries of the current object
              * Set the replacement_entries_ attribute of the current obsolete
              * @param replacement_entries The replacement entries of the current object
              */
            void SetReplacementEntries(const std::vector<std::string> replacement_entries);

/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the obsolete contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string replacement_date_;             /*!< Replacement Date >*/
            std::vector<std::string> replacement_entries_;          /*!< Replacement Entries >*/
    };
}
#endif // PDBOBSOLETECARD_HPP
