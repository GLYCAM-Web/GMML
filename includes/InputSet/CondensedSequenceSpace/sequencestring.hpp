#ifndef GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCESPACE_SEQUENCESTRING_HPP
#define GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCESPACE_SEQUENCESTRING_HPP

#include <string>

namespace CondensedSequenceSpace
{
    class SequenceString
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            SequenceString (std::string sequence);
            SequenceString (const SequenceString& another_sequence);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetOriginalSequence () const;
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetOriginalSequence (std::string sequence);
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////
            
            //////////////////////////////////////////////////////////
            //                       OPERATORS                      //
            //////////////////////////////////////////////////////////
            void operator= (const SequenceString& another_sequence);
        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string original_sequence;
    }; // end SequenceString
} // end CondensedSequenceSpace

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
std::ostream& operator<<(std::ostream&, const CondensedSequenceSpace::SequenceString&);

#endif // GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCESPACE_SEQUENCESTRING_HPP
