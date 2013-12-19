// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbhelixresidue.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHelixResidue::PdbHelixResidue() : residue_name_(""), residue_chain_id_(' '), residue_sequence_number_(kNotSet), residue_insertion_code_(' ' ) {}
PdbHelixResidue::PdbHelixResidue(const string &residue_name, char residue_chain_id, int residue_sequence_number, char residue_insertion_code)
    : residue_name_(residue_name), residue_chain_id_(residue_chain_id), residue_sequence_number_(residue_sequence_number),
      residue_insertion_code_(residue_insertion_code) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHelixResidue::GetResidueName()
{
    return residue_name_;
}

char PdbHelixResidue::GetResidueChainId()
{
    return residue_chain_id_;
}

int PdbHelixResidue::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}

char PdbHelixResidue::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHelixResidue::SetResidueName(const string residue_name)
{
    residue_name_ = residue_name;
}

void PdbHelixResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}

void PdbHelixResidue::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}

void PdbHelixResidue::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

