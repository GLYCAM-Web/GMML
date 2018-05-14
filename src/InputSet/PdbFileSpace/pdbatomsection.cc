
#include "../../../includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbAtomSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbAtomSection::PdbAtomSection() : record_name_("ATOM") {}

PdbAtomSection::PdbAtomSection(std::stringstream &stream_block, std::string index)
{
    atom_cards_ = PdbAtomMap();
    ordered_atom_cards_ = PdbAtomCardOrderVector();
    std::string line;
    bool is_record_name_set = false;
//    cout << stream_block.str() << std::endl;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_record_name_set){
//            record_name_ = line.substr(0,6);
            record_name_ = "ATOM";
            gmml::Trim(record_name_);
            is_record_name_set=true;
        }

        PdbFileSpace::PdbAtomCard* atom_card = new PdbFileSpace::PdbAtomCard(line);
        //int ch = 65 + ConvertString<int>(Split(index, "_")[1]);
        atom_card->SetAtomCardIndexInResidueSet(index);
        //atom->SetAtomChainId((char)ch);
        atom_cards_[atom_card->GetAtomSerialNumber()] = atom_card;
        ordered_atom_cards_.push_back(atom_card);

        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbAtomSection::GetRecordName()
{
    return record_name_;
}
PdbAtomSection::PdbAtomMap PdbAtomSection::GetAtomCards()
{
    return atom_cards_;
}
PdbAtomSection::PdbAtomCardOrderVector PdbAtomSection::GetOrderedAtomCards()
{
    return ordered_atom_cards_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbAtomSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbAtomSection::SetAtomCards(PdbAtomMap atom_cards)
{
    atom_cards_.clear();
    for(PdbAtomMap::iterator it = atom_cards.begin(); it != atom_cards.end(); it++)
    {
        PdbFileSpace::PdbAtomCard* atom_card = (*it).second;
        int serial_number = (*it).first;
        atom_cards_[serial_number] = atom_card;
    }
}
void PdbAtomSection::SetOrderedAtomCards(PdbAtomCardOrderVector ordered_atom_cards)
{
    ordered_atom_cards_.clear();
    for(PdbAtomCardOrderVector::iterator it = ordered_atom_cards.begin(); it != ordered_atom_cards.end(); it++)
    {
        PdbFileSpace::PdbAtomCard* atom_card = (*it);
        ordered_atom_cards_.push_back(atom_card);
    }
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbAtomSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "_________________ Atoms _______________" << std::endl;
    for(PdbAtomSection::PdbAtomMap::iterator it = atom_cards_.begin(); it != atom_cards_.end(); it++)
    {
        out << "Atom Serial Number: ";
        if((it)->first != gmml::iNotSet)
            out << (it)->first << std::endl;
        else
            out << " " << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
    out << std::endl;
}
