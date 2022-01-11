#include <fstream>
#include <iostream>
#include "../../../includes/InputSet/CifFileSpace/ciffile.hpp"
//#include "../../../includes/utils.hpp"
#include "../../../includes/InputSet/CifFileSpace/ciffileprocessingexception.hpp"
#include "includes/CodeUtils/logging.hpp"

using CifFileSpace::CifFile;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CifFile::CifFile() {}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
std::string CifFile::GetPath()
{
    return path_;
}

std::string CifFile::GetResidueName()
{
    return residue_name_;
}

CifFile::CifFileAtomVector CifFile::GetAtoms()
{
    return atoms_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void CifFile::SetPath(const std::string path)
{
    path_ = path;
}

void CifFile::SetResidueName(const std::string residue_name)
{
    residue_name_ = residue_name;
}

void CifFile::SetAtoms(CifFile::CifFileAtomVector atoms)
{
    atoms_ = atoms;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////
///FRU C1  C1  C 0 1 N N N -6.763 20.037 85.842 0.791  1.055  -1.776 C1  FRU 1
void CifFile::Read(std::ifstream& in_file)
{
    std::string line;
    // int line_number = 0;

    /// Unable to read file
    if (!getline(in_file, line))
    {
        throw CifFileProcessingException("Error reading file");
    }

    /// Skip first lines at the begining of the file
    while(line[0] == '#' || line[0] == '_' || line.find("loop_") != std::string::npos || line.find("data_") != std::string::npos)
    {
        getline(in_file, line);
    }


    while((line[0] != '#' || line[0] != '_' || line.find("loop_") == std::string::npos || line.find("data_") == std::string::npos))
    {
        getline(in_file, line);
//        atoms_.push_back(new CifFileAtom(line));
//        std::cout << line << std::endl;

    }
}
//     int listing_index = 1;
//     getline(in_file, line);
//     while(line[0] != '!')
//     {
//         try
//         {
//             /// Process index section
//             RemoveQuotes(line);
//             RemoveSpaces(line);
//             residues_[line] = new LibraryFileResidue(line, listing_index);
//             listing_index++;
//             getline(in_file,line);      /// Read the next line
//         } catch(...)
//         {
//             throw LibraryFileProcessingException(__LINE__, "Error processing index section");
//         }
//     }
// }

// ///READER

// bool PdbqtFile::ParseCards(ifstream &in_stream)
// {
//     std::string line;

//     /// Unable to read file
//     if (!getline(in_stream, line))
//     {
//         gmml::log(__LINE__, __FILE__,  gmml::ERR,"Wrong input file format");
//         std::cout << "Wrong input file format" << std::endl;
//         throw PdbqtFileProcessingException("Error reading file");
//     }
//     std::string record_name = Split(line, " ").at(0);
//     if(record_name.compare("MODEL") == 0)
//     {
//         if(!ParseModelCard(in_stream, line))
//             return false;
//     }
//     else
//     {
//         gmml::log(__LINE__, __FILE__,  gmml::ERR,"Wrong input file format");
//         std::cout << "Wrong input file format" << std::endl;
//         return false;
//     }
//     return true;
// }
 



//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void CifFile::Print(std::ostream &out)
{
    out << "";
}
