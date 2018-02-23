#include "../../../includes/InputSet/PdbFileSpace/pdbmastercard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbMasterCard::PdbMasterCard() {}

PdbMasterCard::PdbMasterCard(stringstream& stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        record_name_ = line.substr(0,6);
        num_remark_ = atoi(line.substr(10,5).c_str());
        num_het_ = atoi(line.substr(20,5).c_str());
        num_helix_ = atoi(line.substr(25,5).c_str());
        num_sheet_ = atoi(line.substr(30,5).c_str());
        num_site_ = atoi(line.substr(40,5).c_str());
        num_x_form_ = atoi(line.substr(45,5).c_str());
        num_coord_ = atoi(line.substr(50,5).c_str());
        num_ter_ = atoi(line.substr(55,5).c_str());
        num_connect_ = atoi(line.substr(60,5).c_str());
        num_seq_ = atoi(line.substr(65,5).c_str());
        getline(stream_block, line);
        temp = line;
    }
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

string PdbMasterCard::GetRecordName()
{
    return record_name_;
}

int PdbMasterCard::GetNumRemark()
{
    return num_remark_;
}

int PdbMasterCard::GetNumHet()
{
    return num_het_;
}

int PdbMasterCard::GetNumHelix()
{
    return num_helix_;
}

int PdbMasterCard::GetNumSheet()
{
    return num_sheet_;
}

int PdbMasterCard::GetNumSite()
{
    return num_site_;
}

int PdbMasterCard::GetNumXForm()
{
    return num_x_form_;
}

int PdbMasterCard::GetNumCoord()
{
    return num_coord_;
}

int PdbMasterCard::GetNumTer()
{
    return num_ter_;
}

int PdbMasterCard::GetNumConnect()
{
    return num_connect_;
}

int PdbMasterCard::GetNumSeq()
{
    return num_seq_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbMasterCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbMasterCard::SetNumRemark(int num_remark)
{
    num_remark_ = num_remark;
}

void PdbMasterCard::SetNumHet(int num_het)
{
    num_het_ = num_het;
}

void PdbMasterCard::SetNumHelix(int num_helix)
{
    num_helix_ = num_helix;
}

void PdbMasterCard::SetNumSheet(int num_sheet)
{
    num_sheet_ = num_sheet;
}

void PdbMasterCard::SetNumSite(int num_site)
{
    num_site_ = num_site;
}

void PdbMasterCard::SetNumXForm(int num_x_form)
{
    num_x_form_ = num_x_form;
}

void PdbMasterCard::SetNumCoord(int num_coord)
{
    num_coord_ = num_coord;
}

void PdbMasterCard::SetNumTer(int num_ter)
{
    num_ter_ = num_ter;
}

void PdbMasterCard::SetNumConnect(int num_connect)
{
    num_connect_ = num_connect;
}

void PdbMasterCard::SetNumSeq(int num_seq)
{
    num_seq_ = num_seq;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbMasterCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_;
    out << "REMARK Cards: " << num_remark_;
    out << "HET Cards: " << num_het_;
    out << "HELIX Cards: " << num_helix_;
    out << "SHEET Cards: " << num_sheet_;
    out << "SITE Cards: " << num_site_;
    out << "ORIGX+SCALE+MTRIX Cards: " << num_x_form_;
    out << "ATOM+HETATM Cards: " << num_coord_;
    out << "TER Cards: " << num_ter_;
    out << "CONECT Cards: " << num_connect_;
    out << "SEQRES Cards: " << num_seq_;
    out << endl;
}
