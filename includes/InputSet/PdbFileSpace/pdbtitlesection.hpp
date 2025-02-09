// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBTITLESECTION_HPP
#define PDBTITLESECTION_HPP

#include <string>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbTitleSection
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        /*! \fn
         * Default constructor
         */
        PdbTitleSection();
        /*! \fn
         * Constructor with required parameters
         * @param record_name Name for a title card record which appears in the first column of each line in a pdb file
         * @param title Title of a pdb file
         */
        PdbTitleSection(const std::string& record_name, const std::string& title);
        /*! \fn
         * Constructor with required parameters
         * @param stream_block
         */
        PdbTitleSection(std::stringstream& stream_block);

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        /** \addtogroup Molecular_Data_Structure
         * @{
         */
        /*! \fn
         * An accessor function in order to access to the record name in a title card
         * @return record_name_ attribute of the current object of this class
         */
        std::string GetRecordName();
        /*! \fn
         * An accessor function in order to access to the title in a title card
         * @return title_ attribute of the current object of this class
         */
        std::string GetTitle();
        /** @}*/
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        /** \addtogroup Manipulators
         * @{
         */
        /*! \fn
         * A mutator function in order to set the record name of the current object
         * Set the record_name_ attribute of the current title card
         * @param record_name The record name of the current object
         */
        void SetRecordName(const std::string record_name);
        /*! \fn
         * A mutator function in order to set the title of the current object
         * Set the title_ attribute of the current title card
         * @param title The title of the current object
         */
        void SetTitle(const std::string title);
        /** @}*/
        //////////////////////////////////////////////////////////
        //                        FUNCTIONS                     //
        //////////////////////////////////////////////////////////

        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        /*! \fn
         * A function to print out the pdb title card contents in a structural format
         * Print out the information in a defined structure
         * @param out An output stream, the print result will be written in the given output stream
         */
        void Print(std::ostream& out = std::cerr);

      private:
        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        std::string record_name_; /*!< Record name which appears in the first column of each line in a pdb file >*/
        std::string title_;       /*!< Title that appears in TITLE record of a pdb file >*/
    };
} // namespace PdbFileSpace

#endif // PDBTITLESECTION_HPP
