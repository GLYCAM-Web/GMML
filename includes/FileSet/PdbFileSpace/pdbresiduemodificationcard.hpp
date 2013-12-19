// Author: Alireza Khatamian

#ifndef PDBRESIDUEMODIFICATIONCARD_HPP
#define PDBRESIDUEMODIFICATIONCARD_HPP

#include <map>
#include <string>

namespace PdbFileSpace
{
    class PdbResidueModification;
    class PdbResidueModificationCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, PdbResidueModification*> ResidueModificationMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbResidueModificationCard();
            PdbResidueModificationCard(const std::string& record_name);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            ResidueModificationMap GetResidueModifications();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetRecordName(const std::string record_name);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////


        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            ResidueModificationMap residue_modifications_;
    };
}

#endif // PDBRESIDUEMODIFICATIONCARD_HPP
