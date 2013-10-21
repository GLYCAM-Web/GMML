#ifndef PARAMETERFILE_HPP
#define PARAMETERFILE_HPP

#include <string>
#include <vector>
#include <map>


///////////////////// FORWARD DECLARATION //////////////////
namespace ParameterFileSpace
{
    class ParameterFileAtom;
    class ParameterFileBond;
    class ParameterFileAngle;
    class ParameterFileDihedral;
    class ParameterFile
    {
        public:
            //////////////////////////////////////////////////////////
            //                     TYPE DEFINITION                  //
            //////////////////////////////////////////////////////////
            /*! \def
              * A mapping between an atom type and its atom object
              */
            typedef std::map<std::string, ParameterFileAtom*> AtomTypeMap;
            /*! \def
              * A mapping between a connection of two atom types and its bond object
              */
            typedef std::map<std::vector<std::string>, ParameterFileBond*> BondMap;
            /*! /def
              * A mapping between a connection of three atom types and its angle object
              */
            typedef std::map<std::vector<std::string>, ParameterFileAngle*> AngleMap;
            /*! \def
              * A mapping between a connection of four atom types and its dihedral object
              */
            typedef std::map<std::vector<std::string>, ParameterFileDihedral*> DihedralMap;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Constructor
              * @param param_file An existing library file path to be read
              */
            ParameterFile(const std::string& param_file);

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to parameter file path of the current object
              * @return path_ attribute of the current object of this class
              */
            const std::string& GetFilePath() const;
            /*! \fn
              * An accessor function in order to access to the title of the current object
              * @return title_ attribute of the current object of this class
              */
            const std::string& GetTitle() const;
            /*! \fn
              * An accessor function in order to access to the map of atom types (as string) to their atom objects of the current object
              * @return atom_types_ attribute of the current object of this class
              */
            const AtomTypeMap& GetAtomTypes() const ;
            /*! \fn
              * An accessor function in order to access to the map of bonds (dual string of atom types) to their bond objects of the current object
              * @return bonds_ attribute of the current object of this class
              */
            const BondMap& GetBonds() const;
            /*! \fn
              * An accessor function in order to access to the map of angles (triple string of atom types) of the current object
              * @return angles_ attribute of the current object of this class
              */
            const AngleMap& GetAngles() const;
            /*! \fn
              * An accessor function in order to access to the map of dihedrals (fourple string of atom types) of the current object
              * @return path_ attribute of the current object of this class
              */
            const DihedralMap& GetDihedrals() const;

            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a parameter file
              */
            void Read(std::ifstream& in_file);
            /*! \fn
              * A function that parses a line of atom type section of the current object
              * Process the atom type lines of the parameter file
              * @param line A line of the atom type section of the current parameter file
              */
            void ProcessAtomType(const std::string& line);
            /*! \fn
              * A function to paras a line of hydrophicil atoms section of the current object
              * Process the hydrophilic atom type line of the parameter file
              * @param line A line of the hydrophilic atom types section of the current parameter file
              */
            void ProcessHydrophilicAtomType(const std::string& line);
            /*! \fn
              * A function that parses a line of bond section of the current object
              * Process the bond lines of the parameter file
              * @param line A line of the bond section of the current parameter file
              */
            void ProcessBond(const std::string& line);
            /*! \fn
              * A function to parse a line corresponding to the angle section of the current object
              * Process the angle lines of the parameter file
              * @param line A line of the angle section of the current parameter file
              */
            void ProcessAngle(const std::string& line);
            /*! \fn
              * A function that parses a line of dihedral section of the current object
              * Process the dihedral lines of the parameter file
              * @param line A line of dihedral section of the current parameter file
              * @param line_number The line number of the current read line from the parameter file
              * @param in_file A stream of the current parameter file
              */
            void ProcessDihedral(std::string& line, int& line_number, std::ifstream& in_file);
            /*! \fn
              * A function that parses a line of improper dihedral section of the current object
              * Process the improper dihedral lines of the parameter file
              * @param line A line of improper dihedral section of the current parameter file
              * @param line_number The line number of the current read line from the parameter file
              * @param in_file A stream of the current parameter file
              */
            void ProcessImproperDihedral(std::string& line, int& line_number, std::ifstream& in_file);
            /*! \fn
              * A function that parses a line of hydrogen bond section of the current object
              * Process the hydrogen bond lines of the parameter file
              * @param line A line of hydrogen bond section of the current parameter file
              */
            void ProcessHydrogenBond(const std::string& line);
            /*! \fn
              * A function that parses a line of equivalent symbols section of the current object
              * Process the equivalent symbols lines of the parameter file
              * @param line A line of equivalent symbols section of the current parameter file
              */
            void ProcessEquivalentSymbols(const std::string& line);
            /*! \fn
              * A function that parses a line of potential parameter section of the current object
              * Process the potential parameter lines of the parameter file located at the end of the file
              * @param line A line of potential parameter section of the current parameter file
              */
            void ProcessPotentialParameter(const std::string& line);
            /*! \fn
              * A function that parses a description of a line in dihedral section of the current object to extract double value of the given key
              * Process the description of a line in dihedral section of the parameter file
              * @param dscr The description part of a line in dihedral section of the current parameter file
              * @param key A key string which a value is assigned to it in the first parameter
              * @return A double value of the given key string
              */
            double ProcessDoubleDihedralDescription(const std::string& dscr, const std::string& key);

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            std::string path_;           /*!< Path of the given parameter file */
            std::string title_;          /*!< Title of the parameter file */
            AtomTypeMap atom_types_;     /*!< A collection of mapping between atom type and its attributes*/
            BondMap bonds_;              /*!< A collection of mapping between bond (double atom types) and its attributes*/
            AngleMap angles_;            /*!< A collection of mapping between angle (tripple atom types) and its attributes*/
            DihedralMap dihedrals_;      /*!< A collection of mapping between dihedral (quad atom types) and its attributes*/

    };
}

#endif // PARAMETERFILE_HPP
