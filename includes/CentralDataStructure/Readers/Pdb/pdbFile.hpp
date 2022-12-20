#ifndef INCLUDES_INPUTSET_PDBFILE_PDBFILE_HPP
#define INCLUDES_INPUTSET_PDBFILE_PDBFILE_HPP
// ToDo split preprocessor into separate function.
// ToDo Warning about gaps being 1 residue.
// ToDo make more direct queries here instead of giving out HeaderRecord etc.
// ToDo ACE/NME between residues with same number but an insertion code.
#include "includes/CentralDataStructure/ensemble.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/authorRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/databaseReferenceRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/headerRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/journalRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/remarkRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/titleRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"

#include <string>
#include <vector>

namespace pdb
{
const int iPdbLineLength = 80;
class PdbModel;
class PdbFile : public cds::Ensemble
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    PdbFile();
    PdbFile(const std::string &pdbFilePath);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    inline std::string GetInputFilePath() const {return inFilePath_;}
    // These should be private and whatever info they give out should be directly queryable here.
    inline const HeaderRecord& GetHeaderRecord() const {return headerRecord_;}
    inline const TitleRecord& GetTitleRecord() const {return titleRecord_;}
    inline const AuthorRecord& GetAuthorRecord() const {return authorRecord_;}
    inline const JournalRecord& GetJournalRecord() const {return journalRecord_;}
    inline const RemarkRecord& GetRemarkRecord() const {return remarkRecord_;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    std::string GetUniprotIDs() const;
    const float& GetResolution() const;
    const float& GetBFactor() const;
    pdb::PreprocessorInformation PreProcess(PreprocessorOptions options);
    //////////////////////////////////////////////////////////
    //                        DISPLAY                       //
    //////////////////////////////////////////////////////////
    void Write(const std::string outName) const;
    void Write(std::ostream& stream) const;
private:
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    void ParseInFileStream(std::ifstream& pdbFileStream);
    std::stringstream ExtractHeterogenousRecordSection(std::ifstream &pdbFileStream, std::string &line, const std::vector<std::string> recordNames);
    std::stringstream ExtractHomogenousRecordSection(std::ifstream &pdbFileStream, std::string &line, std::string previousName);
    inline const std::vector<DatabaseReference>& GetDatabaseReferences() const {return databaseReferences_;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                        ATTRIBUTES                    //
    //////////////////////////////////////////////////////////
    std::string inFilePath_;
    HeaderRecord headerRecord_;
    TitleRecord titleRecord_;
    AuthorRecord authorRecord_;
    JournalRecord journalRecord_;
    RemarkRecord remarkRecord_;
    std::vector<DatabaseReference> databaseReferences_;
};
}
#endif /* INCLUDES_INPUTSET_PDBFILE_PDBFILE_HPP */
