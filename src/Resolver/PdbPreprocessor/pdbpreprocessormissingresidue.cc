
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorMissingResidue::PdbPreprocessorMissingResidue() {}

PdbPreprocessorMissingResidue::PdbPreprocessorMissingResidue(char chain_id, int start_index, int end_index, int index_before_gap, int index_after_gap, char starting_residue_insertion_code, char ending_residue_insertion_code) :
    residue_chain_id_(chain_id), starting_residue_sequence_number_(start_index), ending_residue_sequence_number_(end_index),
    residue_before_gap_(index_before_gap), residue_after_gap_(index_after_gap),
    starting_residue_insertion_code_(starting_residue_insertion_code), ending_residue_insertion_code_(ending_residue_insertion_code)
{
    selected_n_termination_ = (PossibleNChainTermination)1;
    selected_c_termination_ = (PossibleCChainTermination) 2;
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbPreprocessorMissingResidue::GetResidueChainId()
{
    return residue_chain_id_;
}
int PdbPreprocessorMissingResidue::GetStartingResidueSequenceNumber()
{
    return starting_residue_sequence_number_;
}
int PdbPreprocessorMissingResidue::GetEndingResidueSequenceNumber()
{
    return ending_residue_sequence_number_;
}
int PdbPreprocessorMissingResidue::GetResidueBeforeGap()
{
    return residue_before_gap_;
}
int PdbPreprocessorMissingResidue::GetResidueAfterGap()
{
    return residue_after_gap_;
}
PossibleNChainTermination PdbPreprocessorMissingResidue::GetSelectedNTermination()
{
    return selected_n_termination_;
}
PossibleCChainTermination PdbPreprocessorMissingResidue::GetSelectedCTermination()
{
    return selected_c_termination_;
}
char PdbPreprocessorMissingResidue::GetStartingResidueInsertionCode()
{
    return starting_residue_insertion_code_;
}
char PdbPreprocessorMissingResidue::GetEndingResidueInsertionCode()
{
    return ending_residue_insertion_code_;
}

string PdbPreprocessorMissingResidue::GetStringFormatOfSelectedNTermination()
{
    switch (selected_n_termination_)
    {
        case 1:
            return "COCH3";
        case 2:
            return "NH3+";
    }
}

string PdbPreprocessorMissingResidue::GetStringFormatOfSelectedCTermination()
{
    switch(selected_c_termination_)
    {
        case 1:
            return "NH2";
        case 2:
            return "NHCH3";
        case 3:
            return "CO2-";
    }
}

string PdbPreprocessorMissingResidue::GetStringFormatOfNTermination(PossibleNChainTermination n_termination)
{
    switch (n_termination)
    {
        case 1:
            return "COCH3";
        case 2:
            return "NH3+";
    }
}

string PdbPreprocessorMissingResidue::GetStringFormatOfCTermination(PossibleCChainTermination c_termination)
{
    switch(c_termination)
    {
        case 1:
            return "NH2";
        case 2:
            return "NHCH3";
        case 3:
            return "CO2-";
    }
}

vector<string> PdbPreprocessorMissingResidue::GetAllPossibleNChainTerminationAsString()
{
    vector<string> all_possible_n_chian_termination_as_string;
    for(int possible_n_chian_termination = COCH3; possible_n_chian_termination != NH3; possible_n_chian_termination++)
        all_possible_n_chian_termination_as_string.push_back(GetStringFormatOfNTermination((PossibleNChainTermination)possible_n_chian_termination));
    return all_possible_n_chian_termination_as_string;
}

vector<string> PdbPreprocessorMissingResidue::GetAllPossibleCChainTerminationAsString()
{
    vector<string> all_possible_c_chian_termination_as_string;
    for(int possible_c_chian_termination = NH2; possible_c_chian_termination != CO2; possible_c_chian_termination++)
        all_possible_c_chian_termination_as_string.push_back(GetStringFormatOfCTermination((PossibleCChainTermination)possible_c_chian_termination));
    return all_possible_c_chian_termination_as_string;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorMissingResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorMissingResidue::SetStartingResidueSequenceNumber(int starting_residue_sequence_number)
{
    starting_residue_sequence_number_ = starting_residue_sequence_number;
}
void PdbPreprocessorMissingResidue::SetEndingResidueSequenceNumber(int ending_residue_sequence_number)
{
    ending_residue_sequence_number_ = ending_residue_sequence_number;
}
void PdbPreprocessorMissingResidue::SetResidueBeforeGap(int residue_before_gap)
{
    residue_before_gap_ = residue_before_gap;
}
void PdbPreprocessorMissingResidue::SetResidueAfterGap(int residue_after_gap)
{
    residue_after_gap_ = residue_after_gap;
}
void PdbPreprocessorMissingResidue::SetSelectedNTermination(PossibleNChainTermination selected_n_termination)
{
    selected_n_termination_ = selected_n_termination;
}
void PdbPreprocessorMissingResidue::SetSelectedCTermination(PossibleCChainTermination selected_c_termination)
{
    selected_c_termination_ = selected_c_termination;
}
void PdbPreprocessorMissingResidue::SetStartingResidueInsertionCode(char starting_residue_insertion_code)
{
    starting_residue_insertion_code_ = starting_residue_insertion_code;
}
void PdbPreprocessorMissingResidue::SetEndingResidueInsertionCode(char ending_residue_insertion_code)
{
    ending_residue_insertion_code_ = ending_residue_insertion_code;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorMissingResidue::Print(ostream &out)
{
    out << "Chain: " << starting_residue_sequence_number_ << "-" << ending_residue_sequence_number_ << "(" << residue_chain_id_ << ")"
           << ", Residue insertion code before gap: " << starting_residue_insertion_code_
              << ", Residue insertion code after gap: " << ending_residue_insertion_code_
        << ", Residue before gap: " << residue_before_gap_
        << ", Selected c termination: " << GetStringFormatOfSelectedCTermination()
        << ", Residue after gap: " << residue_after_gap_
        << ", Selected n termination: " << GetStringFormatOfSelectedNTermination()           
        << endl;
}







