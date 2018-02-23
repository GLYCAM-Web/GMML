#include "../../../includes/InputSet/PdbFileSpace/pdbcispeptidecard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbFileSpace;


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCISPeptideCard::PdbCISPeptideCard() {}
PdbCISPeptideCard::PdbCISPeptideCard(string &line)
{
    if (!Trim(line).empty())
    {
      record_name_ = line.substr(0, 6);
      serial_number_ = atoi(line.substr(7, 3).c_str());
      pep_1_ = line.substr(11, 3);
      chain_id_1_ = line.substr(15, 1);
      seq_num_1_ = atoi(line.substr(17, 4).c_str());
      i_code_1_ = line.substr(21, 1);
      pep_2_ = line.substr(25, 3);
      chain_id_2_ = line.substr(29, 1);
      seq_num_2_ = atoi(line.substr(31, 4).c_str());
      i_code_2_ = line.substr(35, 1);
      mod_num_ = atoi(line.substr(43, 3).c_str());
      measure_ = atof(line.substr(53, 6).c_str());
    }
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

string PdbCISPeptideCard::GetRecordName()
{
    return record_name_;
}

int PdbCISPeptideCard::GetSerialNumber()
{
    return serial_number_;
}

string PdbCISPeptideCard::GetPeptide1ResidueName()
{
    return pep_1_;
}
string PdbCISPeptideCard::GetPeptide1ChainId()
{
    return chain_id_1_;
}

int PdbCISPeptideCard::GetPeptide1SequenceNumber()
{
    return seq_num_1_;
}
string PdbCISPeptideCard::GetPeptide1InsertionCode()
{
    return i_code_1_;
}

string PdbCISPeptideCard::GetPeptide2ResidueName()
{
    return pep_2_;
}
string PdbCISPeptideCard::GetPeptide2ChainId()
{
    return chain_id_2_;
}

int PdbCISPeptideCard::GetPeptide2SequenceNumber()
{
    return seq_num_2_;
}
string PdbCISPeptideCard::GetPeptide2InsertionCode()
{
    return i_code_2_;
}
int PdbCISPeptideCard::GetModelNumber()
{
    return mod_num_;
}

float PdbCISPeptideCard::GetMeasure()
{
    return measure_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbCISPeptideCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbCISPeptideCard::SetSerialNumber(int serial_number)
{
     serial_number_ = serial_number;
}

void PdbCISPeptideCard::SetPeptide1ResidueName(const string pep_1)
{
     pep_1_ = pep_1;
}

void PdbCISPeptideCard::SetPeptide1ChainId(const string chain_id_1)
{
     chain_id_1_ = chain_id_1;
}

void PdbCISPeptideCard::SetPeptide1SequenceNumber(int seq_num_1)
{
     seq_num_1_ = seq_num_1;
}

void PdbCISPeptideCard::SetPeptide1InsertionCode(const string i_code_1)
{
     i_code_1_ = i_code_1;
}

void PdbCISPeptideCard::SetPeptide2ResidueName(const string pep_2)
{
     pep_2_ = pep_2;
}

void PdbCISPeptideCard::SetPeptide2ChainId(const string chain_id_2)
{
     chain_id_2_ = chain_id_2;
}

void PdbCISPeptideCard::SetPeptide2SequenceNumber(int seq_num_2)
{
     seq_num_2_ = seq_num_2;
}

void PdbCISPeptideCard::SetPeptide2InsertionCode(const string i_code_2)
{
     i_code_2_ = i_code_2;
}

void PdbCISPeptideCard::SetModelNumber(int mod_num)
{
     mod_num_ = mod_num;
}

void PdbCISPeptideCard::SetMeasure(float measure)
{
     measure_ = measure;
}


//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbCISPeptideCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_;
    out << "Serial Number: " << serial_number_;
    out << "Peptide 1 Residue Name: " << pep_1_;
    out << "Peptide 1 Chain ID: " << chain_id_1_;
    out << "Peptide 1 Sequence Number: " << seq_num_1_;
    out << "Peptide 1 Insertion Code: " << i_code_1_;
    out << "Peptide 2 Residue Name: " << pep_2_;
    out << "Peptide 2 Chain ID: " << chain_id_2_;
    out << "Peptide 2 Sequence Number: " << seq_num_2_;
    out << "Peptide 2 Insertion Code: " << i_code_2_;
    out << "Model Number: " << mod_num_;
    out << "Angle Measurement: " << measure_;
    out << endl;
}
