#ifndef ASSEMBLY_HPP
#define ASSEMBLY_HPP

#include <string>
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm> // added for std::erase remove

#include "../GeometryTopology/coordinate.hpp"
#include "../GeometryTopology/plane.hpp"
#include "../GeometryTopology/ResidueLinkages/residue_linkage.hpp"
#include "../common.hpp"
#include "../Glycan/chemicalcode.hpp"
#include "../Glycan/sugarname.hpp"
#include "../Glycan/monosaccharide.hpp"
#include "../Glycan/ontologyvocabulary.hpp"
#include "../InputSet/PdbFileSpace/pdbfile.hpp"
#include "../InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../InputSet/TopologyFileSpace/topologyfile.hpp"
#include "../InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../InputSet/PdbFileSpace/pdbremarksection.hpp"
#include "../InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../InputSet/PdbFileSpace/pdbmodelsection.hpp"
#include "../InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../Glycan/oligosaccharide.hpp"
#include "../Glycan/monosaccharide.hpp"
#include "../Glycan/note.hpp"
#include "../InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../InputSet/PdbFileSpace/inputfile.hpp"
#include "molecule.hpp"
#include "Graph/graph.hpp"
// #include "residue.hpp"
// #include "atom.hpp"
// #include "atomnode.hpp"
#include "residuenode.hpp"

// has been moved to residue_linkage.hpp
// class Residue_linkage;
// typedef std::vector<Residue_linkage> ResidueLinkageVector;

namespace MolecularModeling
{
    class Assembly;
    typedef std::vector<MolecularModeling::Assembly*> AssemblyVector;

    class Assembly
    {
      public:
        //////////////////////////////////////////////////////////
        //                    TYPE DEFINITION                   //
        //////////////////////////////////////////////////////////
        //       typedef std::vector<Assembly*> AssemblyVector;
        //       typedef std::vector<Residue*> ResidueVector;
        //            typedef std::vector<Atom*> AtomVector;
        typedef std::map<std::string, gmml::GraphSearchNodeStatus> AtomStatusMap;
        typedef std::map<std::string, Atom*> AtomIdAtomMap;
        //            typedef std::vector<AtomVector> AtomVectorVector;
        typedef std::map<std::string, AtomVector> CycleMap;
        typedef std::map<std::string, std::map<std::string, std::vector<std::string>>> SelectPatternMap;
        typedef std::map<std::string, ResidueVector> HierarchicalContainmentMap;
        //            typedef std::map<Residue*, ResidueVector> ResidueAttachmentMap;
        typedef std::vector<Glycan::Oligosaccharide*> OligosaccharideVector;
        typedef std::map<int, int> AssemblytoPdbSequenceNumberMap;
        typedef std::map<int, int> AssemblytoPdbSerialNumberMap;
        typedef std::map<std::string, std::string> DerivativeModificationMap;
        typedef std::vector<std::vector<std::string>> AttachedGlycanStructuresVector;
        typedef std::vector<Glycan::Note*> NoteVector;
        //   typedef std::vector<AtomNode*>AtomNodeVector; //Added by ayush on 04/11/18 for TopologyFix in assembly
        //  typedef std::vector<ResidueNode*>ResidueNodeVector; //Added by ayush on 11/16/17 for identifying
        //  residuenodes in assembly
        //       typedef std::vector<MolecularModeling::Molecule*> MoleculeVector; //Added by ayush on 11/12/17 for
        //       molecules in assembly
        typedef Assembly TemplateAssembly; // typedef for marking a template assembly, which contains all necessary
                                           // template residues extracted from 3D template library.

        typedef struct
        {
            // pair.first: descriptor type
            // pair.second: descriptor name
            std::vector<std::pair<std::string, std::string>> linear_descriptors;

            // pair.first: residue index
            // pair.second: residue identifier
            std::vector<std::pair<std::string, std::string>> indices;

            // vector[0]: residue identifier 1
            // vector[1]: atom name 1
            // vector[2]: residue identifier 2
            // vector[3]: atom name 2
            std::vector<std::vector<std::string>> residue_links;

            // not sure what kinds of error/warning messages gmml will provide
            std::vector<std::string> error_warning_messages;
        } gmml_api_output;

        struct parsedString
        {
            std::string label;
            GraphDS::Node* node = NULL;
            GraphDS::Edge* edge = NULL;

            parsedString(std::string Label, GraphDS::Node* Node = NULL)
            { // Node is NULL by default for branch brackets
                label = Label;
                node  = Node;
            }

            parsedString(std::string Label, GraphDS::Edge* Edge)
            {
                label = Label;
                edge  = Edge;
            }
        };

        // This structure is used for the MatchPdbAtoms2Glycam function for better recording and error handling. Yao
        // 2022-03-26
        struct Pdb2glycamMatchingTracker;

        struct Pdb2glycamMatchingFailInfo
        {
            Pdb2glycamMatchingTracker* pdb2glycam_match_session_tracker = NULL;
            int iteration_length                                        = 0;
            MolecularModeling::Atom* failed_atom                        = NULL;
            std::string failure_notice;

            Pdb2glycamMatchingFailInfo(Pdb2glycamMatchingTracker* tracker, int depth, MolecularModeling::Atom* atom,
                                       std::string notice)
            {
                this->iteration_length                 = depth;
                this->failed_atom                      = atom;
                this->failure_notice                   = notice;
                this->pdb2glycam_match_session_tracker = tracker;

                if (this->pdb2glycam_match_session_tracker->largest_iteration_length < this->iteration_length)
                {
                    this->pdb2glycam_match_session_tracker->largest_iteration_length = this->iteration_length;
                }
            }
        };

        struct Pdb2glycamMatchingTracker
        {
            std::vector<std::map<MolecularModeling::Atom*, MolecularModeling::Atom*>>
                all_isomorphisms; // Records all possible matches
            // Records the 1st atom that cannot be matched in the last failed match attempt. Can be used for debugging
            // if match failed. Should be ignored if at least one match can be found. std::multimap<int,
            // MolecularModeling::Atom*> iteration_length_first_mismatched_atom_map;
            std::vector<Pdb2glycamMatchingFailInfo*> failures;
            int largest_iteration_length = 0;

            Pdb2glycamMatchingTracker()
            {}
        };

        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        /*! \fn
         * Default constructor
         */
        Assembly();
        /*! \fn
         * Constructor to build a structure from the given set of input files
         * @param file_paths Set of file paths that are required to build a structure
         *         Most of the time it has just one element (topology, prep, lib, pdb) but at some point it needs more
         * than one file (topology+coordinate)
         * @param type Type of the input which is selected from InputFileType enumerator
         */
        Assembly(std::string file_path, gmml::InputFileType type);
        /*! \fn
         * Constructor to build a structure from the given set of input files
         * @param file_path (topology, prep, lib, pdb)
         * @param type Type of the input which is selected from InputFileType enumerator
         */
        Assembly(
            std::vector<std::string> file_paths,
            gmml::InputFileType type); // redundant in parts due to above (new) constructor, but kept for compatibility.
        /*! \fn
         * Constructor to build a structure from multiple file types, a general version of the previous one
         * @param file_paths Set of set of file paths that are required to build a structure
         * @param types Set of inoput file types of the inputs which are selected from InputFileType enumerator
         */
        Assembly(std::vector<std::vector<std::string>> file_paths, std::vector<gmml::InputFileType> types);
        Assembly(Assembly* assembly);

        // An assembly built from a stringstream of atoms in PDB format
        Assembly(std::stringstream& atomStream);

        // An assembly created by residues of another assembly
        Assembly(std::vector<MolecularModeling::Residue*> residueVector);

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        /** \addtogroup Molecular_Data_Structure
         * @{
         */
        /*! \fn
         * An accessor function in order to access to the name
         * @return name_ attribute of the current object of this class
         */
        std::string GetName();
        /*! \fn
         * An accessor function in order to access to the assemblies
         * @return assemblies_ attribute of the current object of this class
         */
        AssemblyVector GetAssemblies();
        /*! \fn
         * An accessor function in order to access to the residues
         * @return residues_ attribute of the current object of this class
         */
        ResidueVector GetResidues();
        /*! \fn
         * An accessor function in order to access to the chemical type
         * @return chemical_type_ attribute of the current object of this class
         */
        std::string GetChemicalType();
        /*! \fn
         * An accessor function in order to access to the sequence number
         * @return sequence_number_ attribute of the current object of this class
         */
        int GetSequenceNumber();
        std::string GetId();
        /*! \fn
         * An accessor function in order to access to the center of description
         * @return description_ attribute of the current object of this class
         */
        std::string GetDescription();
        /*! \fn
         * An accessor function in order to access to the source file
         * @return source_file_ attribute of the current object of this class
         */
        std::string GetSourceFile();
        /*! \fn
         * An accessor function in order to access to the source file type
         * @return source_file_type_ attribute of the current object of this class
         */
        gmml::InputFileType GetSourceFileType();
        /*! \fn
         * An accessor function in order to access to the model index
         * @return model_index_ attribute of the current object of this class
         */
        int GetModelIndex();
        /*! \fn
         * A functions that extracts all atoms of an assembly
         * @return Vector of all atoms in the current object of assembly
         */
        AtomVector GetAllAtomsOfAssembly();
        /*! \fn
         * A functions that extracts all atoms of an assembly that are within protein residues and are backbone atoms
         * @return Vector of all atoms in the current object of assembly that within protein residues and are backbone
         * atoms
         */
        AtomVector GetAllAtomsOfAssemblyWithinProteinBackbone();
        /*! \fn
         * A functions that extracts all atoms of an assembly that are within protein residues but are not backbone
         * atoms
         * @return Vector of all atoms in the current object of assembly that are within protein residues but are not
         * backbone atoms
         */
        AtomVector GetAllAtomsOfAssemblyWithinProteinSidechain();
        /*! \fn
         * A functions that extracts all atoms of an assembly that are within protein residues
         * @return Vector of all atoms in the current object of assembly that are within protein residues
         */
        AtomVector GetAllAtomsOfAssemblyWithinProteinResidues();
        /*! \fn
         * A functions that extracts all atoms of an assembly that are not within protein residues
         * @return Vector of all atoms in the current object of assembly that are not within protein residues
         */
        AtomVector GetAllAtomsOfAssemblyNotWithinProteinResidues();
        /*! \fn
         * A functions that extracts all atoms of an assembly except atoms of water residues
         * @return Vector of all atoms in the current object of assembly except atoms of water residues
         */
        AtomVector GetAllAtomsOfAssemblyExceptProteinWaterResiduesAtoms();
        /*! \fn
         * A functions that extracts all residues of an assembly
         * @return Vector of all residues in the current object of assembly
         */
        ResidueVector GetAllResiduesOfAssembly();
        /*! \fn
         * A functions that extracts all protein residues of an assembly
         * @return Vector of all protein residues in the current object of assembly
         */
        ResidueVector GetAllProteinResiduesOfAssembly();
        /*! \fn
         * A functions that extracts all nucleic acid residues of an assembly
         * @return Vector of all nucleic acid residues in the current object of assembly
         */
        ResidueVector GetAllNucleicAcidResiduesOfAssembly();
        /*! \fn
         * A function to return all coordinates of all atoms in all residues and assemblies of an assembly
         * @return List of all coordinates of all atoms in all residues and assemblies of an assembly
         */
        GeometryTopology::CoordinateVector GetAllCoordinates();
        /*! \fn
         * A function to extract all the coordinates of all the cycle atoms of the monosaccharide.
         * @param mono The Monosaccharide object
         * @return coordinates The CoordinateVector with all the Coordinates
         */
        GeometryTopology::CoordinateVector GetCycleAtomCoordinates(Glycan::Monosaccharide* mono);
        /*! \fn
         * A function to return all issues/notes within an assembly
         * @return List of all notes of an assembly
         */
        NoteVector GetNotes();

        /*! \fn                                                       //Added by ayush on 11/16/17 for residuenodes in
         * assembly
         * * A functions that extracts all residuesnodes in an assembly
         * @return Vector of all residuesnode in the current object of assembly
         */
        MolecularModeling::ResidueNodeVector GetAllResidueNodesOfAssembly();

        /*! \fn                                                       //Added by ayush on 11/12/17 for molecules in
         * assembly
         * * A functions that extracts all molecules of an assembly
         * @return Vector of all molecules in the current object of assembly
         */
        MolecularModeling::MoleculeVector GetMolecules();

        /** @}*/
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        /** \addtogroup Manipulators
         * @{
         */
        /*! \fn
         * A mutator function in order to create a MolecularModeling::Residue
         * This line repeats the same thing above, but in a different way
         * @param PrepFileSpace::PrepFileResidue* and a Residue::Type that is defaulted to type "Undefined"
         */
        Residue& CreateResidue(PrepFileSpace::PrepFileResidue*,
                               Abstract::ResidueType residueType = Abstract::ResidueType::Undefined);
        /*! \fn
         * A mutator function in order to set the name of the current object
         * Set the name_ attribute of the current assembly
         * @param name The name attribute of the current object
         */
        void SetName(std::string name);
        /*! \fn
         * A mutator function in order to set the assemblies of the current object
         * Set the assemblies_ attribute of the current assembly
         * @param assemblies The assemblies attribute of the current object
         */
        void SetAssemblies(AssemblyVector assemblies);
        /*! \fn
         * A function in order to add the assembly to the current object
         * Set the assemblies_ attribute of the current assembly
         * @param assembly The assembly of the current object
         */
        void AddAssembly(Assembly* assembly);
        /*! \fn
         * A function in order to remove the assembly from the current object
         * Deletes from the assemblies_ attribute of this assembly
         * @param assembly The assembly of the current object
         */
        void RemoveAssembly(Assembly* assembly);
        void UpdateIds(std::string new_id);
        /*! \fn
         * A mutator function in order to set the residues of the current object
         * Set the residues_ attribute of the current assembly
         * @param residues The residues attribute of the current object
         */
        void SetResidues(ResidueVector residues);
        /*! \fn
         * A function in order to add the residue to the current object
         * Set the residues_ attribute of the current assembly
         * @param residue The residue of the current object
         */
        void AddResidue(Residue* residue);
        /*! \fn
         * A function in order to insert the residue to the current object
         * Set the residues_ attribute of the current assembly
         * @param residue The residue of the current object
         */
        void InsertResidue(Residue* insert_infront_of, Residue* to_be_inserted); // Added by Yao @ 06-25-2018
        /*! \fn
         * A function in order to remove a residue from the current object
         * Set the residues_ attribute of the current assembly
         * @param residue The residue to be removed.
         */
        void RemoveResidue(Residue* residue);
        /*! \fn
         * A mutator function in order to set the chemical type of the current object
         * Set the chemical_type_ attribute of the current assembly
         * @param chemical_type The chemical_type attribute of the current object
         */
        void SetChemicalType(std::string chemical_type);
        /*! \fn
         * A mutator function in order to set the sequence number of the current object
         * Set the sequence_number_ attribute of the current assembly
         * @param sequence_number The sequence number attribute of the current object
         */
        void SetSequenceNumber(int sequence_number);
        void SetId(std::string id);
        /*! \fn
         * A mutator function in order to set the description of the current object
         * Set the description_ attribute of the current assembly
         * @param description The description attribute of the current object
         */
        void SetDescription(std::string description);
        /*! \fn
         * A mutator function in order to set the source file of the current object
         * Set the source_file_ attribute of the current assembly
         * @param source_file The source file attribute of the current object
         */
        void SetSourceFile(std::string source_file);
        /*! \fn
         * A mutator function in order to set the source file type of the current object
         * Set the source_file_type_ attribute of the current assembly
         * @param source_file_type The source file type attribute of the current object
         */
        void SetSourceFileType(gmml::InputFileType source_file_type);
        /*! \fn
         * A mutator function in order to set the selected model index of the current object
         * Set the model_index_ attribute of the current assembly
         * @param model_index The target model index attribute of the current object
         */
        void SetModelIndex(int model_index);
        /*! \fn
         * A mutator function in order to set the list of notes of the current object
         * Set the notes_ attribute of the current assembly
         * @param notes The notes attribute of the current object
         */
        void SetNotes(NoteVector notes);
        /*! \fn
         * A function in order to add a note instance to the current object
         * Set the notes_ attribute of the current assembly
         * @param note The note instance of the current object
         */
        void AddNote(Glycan::Note* note);

        /*! \fn                                                         //Added by ayush on 11/16/17 for residuenodes in
         * assembly A function in order to set residuenodes to the current object of Assembly Set the residuenodes_
         * attribute of the current assembly
         * @param residuenodes The residuenodes of the current object
         */
        void SetResidueNodes(MolecularModeling::ResidueNodeVector residuenodes);

        /*! \fn                                                         //Added by ayush on 11/12/17 for molecules in
         * assembly A function in order to add the molecules to the current object Set the molecules_ attribute of the
         * current assembly
         * @param molecules The molecules of the current object
         */
        void SetMolecules(MolecularModeling::MoleculeVector molecules);

        void MergeAssembly(Assembly* other);
        /** @}*/
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        bool CheckCondensedSequenceSanity(
            std::string sequence,
            CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree& prep_residues);
        /** @addtogroup Molecular_Data_Structure_Builders
         * @{
         */
        void BuildAssemblyFromCondensedSequence(std::string sequence, std::string prep_file, std::string parameter_file,
                                                bool structure = false);

        /*! \fn
         * A function in order to build an assembly from a condensed sequence, using a prep file as source data
         * Convert a condensed sequence into assembly residues/atoms. Set ResidueNode, AtomNode. Set default
         * phi/psi/omega angles. Crude adjust torsion values to resolve clashes.
         * @param condensed_sequence The condensed sequence string, for example DManpa1-4DGlcpb1-OH
         * @param prep_file A prep file object.
         */
        void BuildAssemblyFromCondensedSequence(
            std::string condensed_sequence,
            PrepFileSpace::PrepFile* prep_file); // Created by Yao 06/25/2018, replace old version above
        /*! \fn
                  * A function in order to build an assembly from a condensed sequence, using a prep file as source data
                  * A wrapper function for python to call. I(Yao) don't know how to initiate a prep file object in
           python. So let python generate a string, convert into prep file object in the wrapper. Then the wrapper calls
           the above version.
                  * @param condensed_sequence The condensed sequence string, for example DManpa1-4DGlcpb1-OH
              * @param prep_file_path A prep file string to a prep file path.
                  */
        void BuildAssemblyFromCondensedSequence(
            std::string condensed_sequence,
            std::string prep_file_path); // Created by Yao 07/02/2018, a temporary wrapepr for Python to call.

        /*! \fn
         * A function in order to obtain all the glycam06 residue names in a condensed sequence, and turn each glycam06
         * residue into an assembly residue This function is called by the new version of
         * "BuildAssemblyFromCondensedSequence" function. It's part of the build process.
         * @param glycam06_residue_tree A vector of glycam06_residue*
         * @param template_assembly A template assembly(see below) contains only one of each type of residue that
         * appears in condensed sequence.
         */
        std::map<int, std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*>>
        ConvertCondensedSequence2AssemblyResidues(
            CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree& glycam06_residue_tree,
            TemplateAssembly* template_assembly);
        /*! \fn
              * A function that sets residue node and atom node connections for assembly residues based on the
           parent/child head/tail atom relationship of the corresponding condensed sequence residue.
              * This function is called by the new version of "BuildAssemblyFromCondensedSequence" function. It's part
           of the build process.
              * @param condensed_sequence_assembly_residue_map This map is the output of
           "ConvertCondensedSequence2AssemblyResidues". It's a map between a glycam06 residue and its corresponding
            assembly residue
              */
        void SetGlycam06ResidueBonding(
            std::map<int,
                     std::pair<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*, MolecularModeling::Residue*>>&
                condensed_sequence_assembly_residue_map);
        /*! \fn
              * A function that recursively go through an oligosaccharide tree structure from the reducing end(e.g.
           ROH,OME etc). It utilizes head-tail relationship, set phi/psi/omega torsions along with the recursion.
              * This function is called by the new version of "BuildAssemblyFromCondensedSequence" function. It's part
           of the build process.
              * @param parent_residue The current residue along the recursion process. It's named parent_residue,
           because this function next goes through one of its child residues, which becomes the next parent residue.
           Termination condition: the current parent_residue has no child, or tail atom equals head atom(dead end).
              */
        void RecursivelySetAngleGeometry(MolecularModeling::Residue* parent_residue);
        void RecursivelySetDihedralAngleGeometry(MolecularModeling::Residue* parent_residue);
        void SetDihedralAngleGeometryWithMetadata();
        void FigureOutResidueLinkagesInGlycan(Residue* from_this_residue1, Residue* to_this_residue2,
                                              ResidueLinkageVector* residue_linkages);
        /*! \fn
              * A function that recursively go through an oligosaccharide tree structure from the reducing end, tagging
           dihedrals without performing any rotation. This function has the same code organization as
           RecursivelySetGeometry. RecursivelySetGeometry sets torsion without tagging, but this function tags without
           performing any rotation.
              * This function is called by the new version of "BuildAssemblyFromCondensedSequence" function. It's part
           of the build process.
              * @param parent_residue The current residue along the recursion process. It's named parent_residue,
           because this function next goes through one of its child residues,
              * @param index_dihedral_map A map between an integer and a pair of AtomVector* and string. Integer is a
           bond index, starts from 0, starts numbering from reducing end. AtomVector* is pointer to an AtomVector
           containing the four atoms of a dihedral. String is the dihedral type, either phi, psi or omega.
              * @param linkage_index An integer that equals the current bond index, starting from the reducing end.
              */
        void RecursivelyTagDihedrals(MolecularModeling::Residue* parent_residue,
                                     std::multimap<int, std::pair<AtomVector*, std::string>>& index_dihedral_map,
                                     int& linkage_index);
        /*! \fn
         * A function that find residues that clash with other residues in an assembly.
         * This function is called by the new version of "BuildAssemblyFromCondensedSequence" function, after
         * RecursivelySetGeometry. It follows up setting geometry, and initiate the attempt  		     to crudely
         * resolve clashes resulting from default geometry setttings
         */
        ResidueVector FindClashingResidues();
        /*! \fn
         * A function that finds all pathways, which starts from a clashing residue, and ends at a branching
         * point(residue with mutiple tail atoms). This function is called by the new version of
         * "BuildAssemblyFromCondensedSequence" function, after FindClashingResidues.
         * @param all_clashing_residues A ResidueVector containing all the clashing residues. It is the output of
         * FindClashingResidues
         */
        std::vector<ResidueVector> FindPathToCommonAncestors(ResidueVector& all_clashing_residues);
        /*! \fn
         * A function that resolves clashes in assembly atoms.
         * This is an orchestrator function that calls FindAllOmegaTorsionsInPathway and FindBestSetOfTorsions
         * @param fused_clashing_paths A map between a common ancestor residue(branching residue), and the residue
         * vector containg all residues in a clashing pathway.
         * @param index_dihedral_map See above RecursivelyTagDihedrals @param documentation
         */
        void ResolveClashes(std::vector<ResidueVector>& fused_clashing_paths,
                            std::multimap<int, std::pair<AtomVector*, std::string>>& index_dihedral_map);

        /*! \fn
         * A function that finds all omega torsions in  pathway
         * The output of this function is a set of omega dihedrals. They will be rotated subsequently to resolve
         * clases(limited grid search)
         * @param pathway A residue vector containing a clashing path. The code finds dihedrals within this path.
         * @param dihedral_map See above.
         */
        std::vector<AtomVector*>
        FindAllOmegaTorsionsInPathway(ResidueVector& pathway,
                                      std::multimap<int, std::pair<AtomVector*, std::string>>& dihedral_map);
        /*! \fn
         * A function that returns the set of coordinates with the lowest clashes that the algorithm can generate.
         * This function is called by ResolveClashes
         * @param available_dihedrals a vector of AtomVector*. Each AtomVector contains four atoms for a dihdral.
         */
        GeometryTopology::CoordinateVector FindBestSetOfTorsions(std::vector<AtomVector*>& available_dihedrals);
        /*! \fn
              * A function that takes in all dihedrals, and their corresponding rotation values. Exhaustively generates
           individual sets of rotation values.
              * This function is called in the "FindBestSetOfTorsions" function
              * @param all_dihedral_rotation_values A vector of pair, each pair is an atomvector containing four atoms,
           and a vector of double containing allowable rotation values e.g pair<C4-O4-C1-C2, -5,0,5>; pair <C5-C6-O6-C1,
           -5,0,5>
              * @param current_dihedral_index The index of the current dihedral in the vector
           "all_dihedral_rotation_values". To initiate recursion from the nth dihedral, supply n-1 as argument
              * @param container_for_combinations A vector of individual rotation combination. Each individual
           combination is a vector of pair<AtomVector*, double>. Each pair represents a
              * particular dihedral and its corresponding rotation value.
              * angle_index_per_dihedral A temporary container holding the rotation value at each dihedral. To initiate
           a new recursion, supply an empty vector.
              */
        void GenerateAllTorsionCombinations(
            std::vector<std::pair<AtomVector*, std::vector<double>>>& all_dihedral_rotation_values,
            unsigned int current_dihedral_index,
            std::vector<std::vector<std::pair<AtomVector*, double>>>& container_for_combinations,
            std::vector<double>& angle_index_per_dihedral);

        /*! \fn
              * A function that "grafts" a child residue onto a parent residue. It puts the head atoms of child a
           certain distance away from the tail atom of parent.In this process, the code keeps parent residue intact, but
           moves the child residue accordingly.
              * This function is called within "RecursiveSetGeometry". It is part of the geometry setting process.
              * @param parent_tail_atom A bond involves the tail atom of parent, and the head atom of child. This is the
           Atom* of parent tail atom
              * @param child_head_atom A bond involves the tail atom of parent, and the head atom of child. This is the
           Atom* of child head atom
              */
        void SetResidueResidueBondDistance(MolecularModeling::Atom* parent_tail_atom,
                                           MolecularModeling::Atom* child_head_atom);
        /** @}*/
        /*! \fn
              * This function is the backend code that handles the website build rotamer feature
              * @param working_assembly The assembly from whose atoms rotamers will be populated, typically "this"
           assembly.
              * @param rotamers_glycosidic_angles_info A vector of "linkages". A linkage is a bond between
           monosaccharide residues in an oligosaccharide. Each linkage contains phi,psi, and sometime omega torsions.
           Detailed documentation check typedef.
              * @index_dihedral_map See above.
              */
        void GenerateRotamersForCondensedSequence(
            Assembly* working_assembly,
            CondensedSequenceSpace::CondensedSequence::CondensedSequenceRotamersAndGlycosidicAnglesInfo
                rotamers_glycosidic_angles_info,
            std::multimap<int, std::pair<AtomVector*, std::string>>& index_dihedral_map);
        void Grafting(MolecularModeling::Atom* parent_tail_atom, MolecularModeling::Atom* child_head_atom,
                      double head_tail_bond_length); // Added by Yao 08/20/2019.
        void AttachResidues(Residue* residue, Residue* parent_residue, int branch_index, std::string parameter_file);
        void RemoveHydrogenAtAttachedPosition(Residue* residue, int branch_index);
        void RemoveAllHydrogenAtoms();
        void RemoveAllProteinResidues();
        void RemoveAllNucleicAcidResidues();
        void SetDerivativeAngle(Residue* residue, Residue* parent_residue, int branch_index);
        void AdjustCharge(Residue* residue, Residue* parent_residue, int branch_index);
        void SetAttachedResidueBond(Residue* residue, Residue* parent_residue, int branch_index,
                                    std::string parameter_file);
        void SetAttachedResidueAngle(Residue* residue, Residue* parent_residue, int branch_index,
                                     std::string parameter_file);
        void SetAttachedResidueTorsion(Residue* residue, Residue* parent_residue, int branch_index);
        void SetPhiTorsion(Residue* residue, Residue* parent_residue, int branch_index, double torsion);
        void SetPsiTorsion(Residue* residue, Residue* parent_residue, int branch_index, double torsion,
                           bool crystallographic_definition = true);
        void SetOmegaTorsion(Residue* residue, Residue* parent_residue, int branch_index, double torsion, int type = 6);
        void SetOmegaDerivativeTorsion(Residue* residue, Residue* parent_residue, int branch_index, double torsion);
        double GetDihedral(MolecularModeling::Atom* atom1, MolecularModeling::Atom* atom2,
                           MolecularModeling::Atom* atom3, MolecularModeling::Atom* atom4); // Yao 09/07/2019
        void SetDihedral(Atom* atom1, Atom* atom2, Atom* atom3, Atom* atom4, double torsion);
        void SetAngle(Atom* atom1, Atom* atom2, Atom* atom3, double angle);
        /*! \fn
         * A function to build a structure from a single pdb file
         * Imports data from pdb file data structure into central data structure
         * @param pdb_file_path Path to a pdb file
         */
        /** \addtogroup Manipulators
         * @{
         */
        void BuildAssemblyFromPdbFile(std::string pdb_file_path,
                                      std::vector<std::string> amino_lib_files  = std::vector<std::string>(),
                                      std::vector<std::string> glycam_lib_files = std::vector<std::string>(),
                                      std::vector<std::string> other_lib_files  = std::vector<std::string>(),
                                      std::vector<std::string> prep_files       = std::vector<std::string>(),
                                      std::string parameter_file                = "");
        /*! \fn
         * A function to build a structure from a single pdb file
         * Imports data from pdb file data structure into central data structure
         * @param pdb_file Pdb file object
         */
        void BuildAssemblyFromPdbFile(PdbFileSpace::PdbFile* pdb_file,
                                      std::vector<std::string> amino_lib_files  = std::vector<std::string>(),
                                      std::vector<std::string> glycam_lib_files = std::vector<std::string>(),
                                      std::vector<std::string> other_lib_files  = std::vector<std::string>(),
                                      std::vector<std::string> prep_files       = std::vector<std::string>(),
                                      std::string parameter_file                = "");

        /*! \fn
         * A function to build a structure from a stringstream of atoms
         * Imports data from stringstream into central data structure
         * @param atomStream stringstream of atom cards from a pdb file
         */
        void BuildAssemblyFromAtomStream(std::stringstream& atomStream);
        /*! \fn
         * A function to build a structure from a single pdbqt file
         * Imports data from pdbqt file data structure into central data structure
         * @param pdbqt_file_path Path to a pdbqt file
         */

        void BuildAssemblyFromPdbqtFile(std::string pdbqt_file_path, std::string parameter_file = "");
        /*! \fn
         * A function to build a structure from a single pdbqt file
         * Imports data from pdbqt file data structure into central data structure
         * @param pdbqt_file Pdbqt file object
         */
        void BuildAssemblyFromPdbqtFile(PdbqtFileSpace::PdbqtFile* pdbqt_file, std::string parameter_file = "");
        /*! \fn
         * A function to build a structure from a single topology file
         * Imports data from topology file data structure into central data structure
         * @param topology_file_path Path to a topology file
         */
        void BuildAssemblyFromTopologyFile(std::string topology_file_path, std::string parameter_file = "");
        /*! \fn
         * A function to build a structure from a single topology file
         * Imports data from topology file data structure into central data structure
         * @param topology_file Topology file object
         */
        void BuildAssemblyFromTopologyFile(TopologyFileSpace::TopologyFile* topology_file,
                                           std::string parameter_file = "");
        /*! \fn
         * A function to build a structure from a single library file
         * Imports data from library file data structure into central data structure
         * @param library_file_path Path to a library file
         */
        void BuildAssemblyFromLibraryFile(std::string library_file_path, std::string parameter_file = "");
        /*! \fn
         * A function to build a structure from a single library file
         * Imports data from library file data structure into central data structure
         * @param library_file Library file object
         */
        void BuildAssemblyFromLibraryFile(LibraryFileSpace::LibraryFile* library_file, std::string parameter_file = "");
        /*! \fn
         * A function to build a structure from a combination of a topology file and its corresponding coordinate file
         * Imports data from topology file data structure into central data structure and assign the atom coordinates
         * based on the coordinate file
         * @param topology_file_path Path to a topology file
         * @param coordinate_file_path Path to a coordinate file corresponding to the given topology file
         */
        void BuildAssemblyFromTopologyCoordinateFile(std::string topology_file_path, std::string coordinate_file_path,
                                                     std::string parameter_file = "");
        /*! \fn
         * A function to build a structure from a combination of a topology file and its corresponding coordinate file
         * Imports data from topology file data structure into central data structure and assign the atom coordinates
         * based on the coordinate file
         * @param topology_file Topology file object
         * @param coordinate_file Coordinate file object corresponding to the given topology file
         */
        void BuildAssemblyFromTopologyCoordinateFile(TopologyFileSpace::TopologyFile* topology_file,
                                                     CoordinateFileSpace::CoordinateFile* coordinate_file,
                                                     std::string parameter_file = "");

        /*! \fn
         * A funcion that builds a template assembly that contains all template residues extracted from prep file.
         * @param prep_file A pointer to prep file
         * @param query_residue_names the names of all residues needed
         */
        TemplateAssembly* BuildTemplateAssemblyFromPrepFile(
            CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree& glycam06_residue_tree,
            PrepFileSpace::PrepFile* prep_file);
        TemplateAssembly* BuildTemplateAssemblyFromPrepFile(std::vector<std::string>& query_residue_names,
                                                            PrepFileSpace::PrepFile* prep_file);
        /*! \fn
         * A function to build a structure from a single prep file
         * Imports data from prep file data structure into central data structure
         * @param prep_file_path Path to a prep file
         */
        void BuildAssemblyFromPrepFile(std::string prep_file_path, std::string parameter_file = "");
        /*! \fn
         * A function to build a structure from a single prep file
         * Imports data from prep file data structure into central data structure
         * @param prep_file Prep file object
         */
        void BuildAssemblyFromPrepFile(PrepFileSpace::PrepFile* prep_file, std::string parameter_file = "");
        //            PdbFileSpace::PdbFile* BuildPdbFileStructureFromCondensedSequence
        /*! \fn
         * A function to build a pdb file structure from the current assembly object
         * Exports data from assembly data structure into pdb file structure
         * @param link_card_direction An integer to define the direction in the link cards (-1: O -> C and 1: C -> O)
         */
        PdbFileSpace::PdbFile* BuildPdbFileStructureFromAssembly(int link_card_direction    = -1,
                                                                 int connect_card_existance = 1, int model_index = -1,
                                                                 bool useInputPDBResidueNumbers = true);
        /*! \fn
         * A function to build a pdbqt file structure from the current assembly object
         * Exports data from assembly data structure into pdbqt file structure
         */
        PdbqtFileSpace::PdbqtFile* BuildPdbqtFileStructureFromAssembly();

        void ExtractPdbModelSectionFromAssembly(PdbFileSpace::PdbModelResidueSet* residue_set, int& serial_number,
                                                int& sequence_number, int model_number, bool useInputPDBResidueNumbers,
                                                AssemblytoPdbSequenceNumberMap& assembly_to_pdb_sequence_number_map,
                                                AssemblytoPdbSerialNumberMap& assembly_to_pdb_serial_number_map);
        void ExtractPdbqtModelCardFromAssembly(PdbqtFileSpace::PdbqtModelResidueSet* residue_set, int& serial_number,
                                               int& sequence_number, int model_number);
        void ExtractPdbLinkSectionFromAssembly(PdbFileSpace::PdbLinkSection* link_card, int model_index,
                                               AssemblytoPdbSequenceNumberMap assembly_to_pdb_sequence_number,
                                               int link_card_direction = -1);
        void ExtractPdbConnectSectionFromAssembly(PdbFileSpace::PdbConnectSection* connect_card,
                                                  AssemblytoPdbSerialNumberMap assembly_to_pdb_serial_number);
        /*! \fn
         * A function to extract bonds from the current assembly object
         * @param inserted_bond_types Bond types that have been already detected in the assembly structure
         * @param assembly_atom A source atom in a bond in the assembly structure
         * @param neighbor Second atom in a bond which is a neighbor of assembly_atom
         * @param bonds All known bonds from parameter file
         * @param bond_type_counter A counter that indicates the number of bond types that have been already detected
         * and also determines the index associated with it
         * @param topology_file Output topology file structure that the detected bond types belong to
         */
        void ExtractTopologyBondTypesFromAssembly(std::vector<std::vector<std::string>>& inserted_bond_types,
                                                  Atom* assembly_atom, Atom* neighbor,
                                                  ParameterFileSpace::ParameterFile::BondMap& bonds,
                                                  int& bond_type_counter,
                                                  TopologyFileSpace::TopologyFile* topology_file);
        /*! \fn
         * A function to extract bond types from the current assembly object
         * @param inserted_bonds Bonds that have been already detected in the assembly structure
         * @param inserted_bond_types Bond types that have been already detected in the assembly structure
         * @param assembly_atom A source atom in a bond in the assembly structure
         * @param neighbor Second atom in a bond which is a neighbor of assembly_atom
         * @param topology_file Output topology file structure that the detected bond types belong to
         */
        void ExtractTopologyBondsFromAssembly(std::vector<std::vector<std::string>>& inserted_bonds,
                                              std::vector<std::vector<std::string>>& inserted_bond_tyoes,
                                              Atom* assembly_atom, Atom* neighbor,
                                              TopologyFileSpace::TopologyFile* topology_file);

        /*! \fn
         * A function to build a prep file structure from the current assembly object
         * Exports data from assembly data structure into prep file structure
         * @param parameter_file_path A path to a parameter file
         */
        PrepFileSpace::PrepFile* BuildPrepFileStructureFromAssembly(std::string parameter_file_path);
        /*! \fn
         * A function to build improper dihedral types of prep file structure from the current assembly object
         * Exports data from assembly data structure to generate prep improper dihedral types
         * @param assembly_atom A source atom in a dihedral in the assembly structure
         * @param inserted_improper_dihedrals_types Improper dihedral types that have been already detected in an
         * assembly structure
         * @param dihedrals All known dihedrals from parameter file
         */
        void ExtractPrepImproperDihedralTypesFromAssembly(Atom* assembly_atom,
                                                          std::vector<std::string>& inserted_improper_dihedral_types,
                                                          ParameterFileSpace::ParameterFile::DihedralMap& dihedrals);
        /*! \fn
         * A function to build dihedrals of prep file structure from the current assembly object
         * Exports data from assembly data structure to generate prep dihedrals
         * @param assembly_atom A source atom in a dihedral in the assembly structure
         * @param inserted_improper_dihedrals Improper dihedrals that have been already detected in an assembly
         * structure
         * @param inserted_improper_dihedral_types Improper dihedral types that have been already detected in an
         * assembly structure
         * @param dihedrals All known dihedrals from parameter file
         */
        void
        ExtractPrepImproperDihedralsFromAssembly(Atom* assembly_atom,
                                                 std::vector<std::vector<std::string>>& inserted_improper_dihedrals,
                                                 std::vector<std::string> inserted_improper_dihedral_types,
                                                 ParameterFileSpace::ParameterFile::DihedralMap& dihedrals);
        /*! \fn
         * A function to extract topological types of atoms of a residue in assembly object
         * @param assembly_atoms Atoms of a residue in the assembly structure
         * @return List of all topological types of atoms in a residue of an assembly
         */
        std::vector<gmml::TopologicalType>
        GetAllTopologicalTypesOfAtomsOfResidue(AtomVector assembly_atoms, PrepFileSpace::PrepFileResidue::Loop& loops,
                                               std::vector<int>& bond_index,
                                               int dummy_atoms = gmml::DEFAULT_DUMMY_ATOMS);
        /*! \fn
         * A function to build a topology file structure from the current assembly object
         * Exports data from assembly data structure into topology file structure
         */
        TopologyFileSpace::TopologyFile*
        BuildTopologyFileStructureFromAssembly(std::string parameter_file_path,
                                               std::string ion_parameter_file_path = "");
        /*! \fn
         * A function to build agnle types of topology file structure from the current assembly object
         * Exports data from assembly data structure to generate topology atom types
         * @param assembly_atom A source atom in an angle in the assembly structure
         * @param neighbor Second atom in an angle which is a neighbor of assembly_atom
         * @param neighbor_of_neighbor Third atom in an angle which is a neighbor of neighbor atom and it is not
         * identical to assmebly_atom
         * @param inserted_angle_types Angle types that have been already detected in an assembly structure
         * @param angle_type_counter A counter that indicates the number of angle types that have been already detected
         * and also determines the index associated with it
         * @param topology_file Output topology file structure that the detected angle types belong to
         * @param angles All known angles from parameter file
         */
        void ExtractTopologyAngleTypesFromAssembly(Atom* assembly_atom, Atom* neighbor, Atom* neighbor_of_neighbor,
                                                   std::vector<std::vector<std::string>>& inserted_angle_types,
                                                   int& angle_type_counter,
                                                   TopologyFileSpace::TopologyFile* topology_file,
                                                   ParameterFileSpace::ParameterFile::AngleMap& angles);
        /*! \fn
         * A function to build agnle types of topology file structure from the current assembly object
         * Exports data from assembly data structure to generate topology atom types
         * @param assembly_atom A source atom in an angle in the assembly structure
         * @param neighbor Second atom in an angle which is a neighbor of assembly_atom
         * @param neighbor_of_neighbor Third atom in an angle which is a neighbor of neighbor atom and it is not
         * identical to assmebly_atom
         * @param inserted_angles Angles that have been already detected in an assembly structure
         * @param inserted_angle_types Angle types that have been already detected in an assembly structure
         * @param angle_type_counter A counter that indicates the number of angle types that have been already detected
         * and also determines the index associated with it
         * @param topology_file Output topology file structure that the detected angle types belong to
         */
        void ExtractTopologyAnglesFromAssembly(Atom* assembly_atom, Atom* neighbor, Atom* neighbor_of_neighbor,
                                               std::vector<std::vector<std::string>>& inserted_angles,
                                               std::vector<std::vector<std::string>>& inserted_angle_types,
                                               TopologyFileSpace::TopologyFile* topology_file);
        /*! \fn
         * A function to build dihedral types of topology file structure from the current assembly object
         * Exports data from assembly data structure to generate topology dihedral types
         * @param assembly_atom A source atom in a dihedral in the assembly structure
         * @param neighbor Second atom in a dihedral which is a neighbor of assembly_atom
         * @param neighbor_of_neighbor Third atom in a dihedral which is a neighbor of neighbor atom and it is not
         * identical to assmebly_atom
         * @param neighbor_of_neighbor_of_neighbor Fourth atom in a dihedral which is a neighbor of neighbor of neighbor
         * atom and it is not identical to neighbor atom
         * @param inserted_dihedrals_types Dihedral types that have been already detected in an assembly structure
         * @param dihedral_type_counter A counter that indicates the number of dihedral types that have been already
         * detected and also determines the index associated with it
         * @param topology_file Output topology file structure that the detected angle types belong to
         * @param dihedrals All known dihedrals from parameter file
         */
        void ExtractTopologyDihedralTypesFromAssembly(Atom* assembly_atom, Atom* neighbor, Atom* neighbor_of_neighbor,
                                                      Atom* neighbor_of_neighbor_of_neighbor,
                                                      std::vector<std::string>& inserted_dihedral_types,
                                                      int& dihedral_type_counter,
                                                      TopologyFileSpace::TopologyFile* topology_file,
                                                      ParameterFileSpace::ParameterFile::DihedralMap& dihedrals);
        /*! \fn
         * A function to build dihedrals of topology file structure from the current assembly object
         * Exports data from assembly data structure to generate topology dihedrals
         * @param assembly_atom A source atom in a dihedral in the assembly structure
         * @param neighbor Second atom in a dihedral which is a neighbor of assembly_atom
         * @param neighbor_of_neighbor Third atom in a dihedral which is a neighbor of neighbor atom and it is not
         * identical to assmebly_atom
         * @param neighbor_of_neighbor_of_neighbor Fourth atom in a dihedral which is a neighbor of neighbor of neighbor
         * atom and it is not identical to neighbor atom
         * @param inserted_dihedrals Dihedrals that have been already detected in an assembly structure
         * @param inserted_dihedral_types Dihedral types that have been already detected in an assembly structure
         * @param dihedrals All known dihedrals from parameter file
         * @param topology_file Output topology file structure that the detected angle types belong to
         */
        void ExtractTopologyDihedralsFromAssembly(Atom* assembly_atom, Atom* neighbor, Atom* neighbor_of_neighbor,
                                                  Atom* neighbor_of_neighbor_of_neighbor,
                                                  std::vector<std::vector<std::string>>& inserted_dihedrals,
                                                  std::vector<std::string>& inserted_dihedral_types,
                                                  ParameterFileSpace::ParameterFile::DihedralMap& dihedrals,
                                                  TopologyFileSpace::TopologyFile* topology_file);

        /*! \fn
         * A function to build a coordinate file structure from the current assembly object
         * Exports data from assembly data structure into coordinate file structure
         */
        /** @}*/
        /** @addtogroup Molecular_Data_Structure_Builders
         * @{ */
        CoordinateFileSpace::CoordinateFile* BuildCoordinateFileStructureFromAssembly();
        /*! \fn
         * A function to build a library file structure from the current assembly object
         * Exports data from assembly data structure into library file structure
         */
        LibraryFileSpace::LibraryFile* BuildLibraryFileStructureFromAssembly();
        /*! \fn
         * A function to build a graph structure (bonding information) for the current object of central data structure
         * @param building_option A building option that can be selected from BuildingStructureOption enumerator
         * @param options List of additional options that can be defined by user
         * @param file_paths In the case that building structure based on the information of database files is selected
         * this argument is a list of database file names
         */
        void BuildStructure(gmml::BuildingStructureOption building_option, std::vector<std::string> options,
                            std::vector<std::string> file_paths);
        /*! \fn
         * A function to build a graph structure for the current object of central data structure based on the distance
         * between the atoms of the structure
         * @param cutoff Threshold of closeness of the atoms to be considered as bonded
         * @param model_index In the case that the structure has multiple model (multiple coordinates for atoms, such as
         * pdb) this arguments indicates the desired model index
         */
        void BuildStructureByDistance(int number_of_threads = 1, double cutoff = gmml::maxCutOff, int model_index = 0);
        /*! \fn
         * A function to build a graph structure for the current object of central data structure based on the bonding
         * information provided in the original file
         */
        void BuildStructureByOriginalFileBondingInformation();
        /*! \fn
         * A function to build a graph structure for the current object of central data structure based on the bonding
         * information provided in the original file in the case that the original file is a pdb file
         */
        void BuildStructureByPDBFileInformation();
        /*! \fn
         * A function to build a graph structure for the current object of central data structure based on the bonding
         * information provided in the original file in the case that the original file is a topology file
         */
        void BuildStructureByTOPFileInformation();
        /*! \fn
         * A function to build a graph structure for the current object of central data structure based on the bonding
         * information provided in the original file in the case that the original file is a lib file
         */
        void BuildStructureByLIBFileInformation();
        /*! \fn
         * A function to build a graph structure for the current object of central data structure based on the bonding
         * information provided in the original file in the case that the original file is a prep file
         */
        void BuildStructureByPrepFileInformation();
        /*! \fn
         * A function to build a graph structure for the current object of central data structure based on the bonding
         * information provided in the given database files
         * @param types List of types of the database files
         * @param file_paths List of the database file paths
         */
        void BuildStructureByDatabaseFilesBondingInformation(std::vector<gmml::InputFileType> types,
                                                             std::vector<std::string> file_paths);
        /** @}*/
        /** \addtogroup Data_Sets
         * @{
         */
        /*! \fn
         * A function that counts the number of atoms in all assemblies and residues of the assembly
         * @return counter Number of atoms in all assemblies and residues in the current object of assembly
         */
        int CountNumberOfAtoms();
        /*! \fn
         * A function that counts the number of atoms types in all assemblies and residues of the assembly
         * @return counter Number of atoms in all assemblies and residues in the current object of assembly
         */
        int CountNumberOfAtomTypes();
        /*! \fn
         * A function that counts the number of residues in all assemblies and residues of the assembly
         * @return type_list Number of residues in all assemblies and residues in the current object of assembly
         */
        int CountNumberOfResidues();
        /*! \fn
         * A function that counts the number of bonds including hydrogen in all assemblies and residues of the assembly
         * @param parameter_file_path The path of the parameter file
         * @return counter Number of bonds including hydrogen in all assemblies and residues in the current object of
         * assembly
         */
        int CountNumberOfBondsIncludingHydrogen(std::string parameter_file_path);
        /*! \fn
         * A function that counts the number of bonds excluding hydrogen in all assemblies and residues of the assembly
         * @param parameter_file_path The path of the parameter file
         * @return counter Number of bonds excluding hydrogen in all assemblies and residues in the current object of
         * assembly
         */
        int CountNumberOfBondsExcludingHydrogen(std::string parameter_file_path);
        /*! \fn
         * A function that counts the number of bonds in all assemblies and residues of the assembly
         * @return counter Number of bonds in all assemblies and residues in the current object of assembly
         */
        int CountNumberOfBonds();
        /*! \fn
         * A function that counts the number of bond types in all assemblies and residues of the assembly
         * @param parameter_file_path The path of the parameter file
         * @return type_list Number of bond types in all assemblies and residues in the current object of assembly
         */
        int CountNumberOfBondTypes(std::string parameter_file_path);
        /*! \fn
         * A function that counts the number of angles including hydrogen in all assemblies and residues of the assembly
         * @param parameter_file_path The path of the parameter file
         * @return counter Number of angles including hydrogen in all assemblies and residues in the current object of
         * assembly
         */
        int CountNumberOfAnglesIncludingHydrogen(std::string parameter_file_path);
        /*! \fn
         * A function that counts the number of angles excluding hydrogen in all assemblies and residues of the assembly
         * @param parameter_file_path The path of the parameter file
         * @return counter Number of angles excluding hydrogen in all assemblies and residues in the current object of
         * assembly
         */
        int CountNumberOfAnglesExcludingHydrogen(std::string parameter_file_path);
        /*! \fn
         * A function that counts the number of angles in all assemblies and residues of the assembly
         * @return counter Number of angles in all assemblies and residues in the current object of assembly
         */
        int CountNumberOfAngles();
        /*! \fn
         * A function that counts the number of angle types in all assemblies and residues of the assembly
         * @param parameter_file_path The path of the parameter file
         * @return type_list Number of angle types in all assemblies and residues in the current object of assembly
         */
        int CountNumberOfAngleTypes(std::string parameter_file_path);
        /*! \fn
         * A function that counts the number of dihedrals includig hydrogen in all assemblies and residues of the
         * assembly
         * @param parameter_file_path The path of the parameter file
         * @return type_list Number of dihedrals including hydrogen in all assemblies and residues in the current object
         * of assembly
         */
        int CountNumberOfDihedralsIncludingHydrogen(std::string parameter_file_path);
        /*! \fn
         * A function that counts the number of dihedrals excluding hydrogen in all assemblies and residues of the
         * assembly
         * @param parameter_file_path The path of the parameter file
         * @return counter Number of dihedrals excluding hydrogen in all assemblies and residues in the current object
         * of assembly
         */
        int CountNumberOfDihedralsExcludingHydrogen(std::string parameter_file_path);
        /*! \fn
         * A function that counts the number of dihedrals in all assemblies and residues of the assembly
         * @param parameter_file_path The path of the parameter file
         * @return counter Number of dihedrals in all assemblies and residues in the current object of assembly
         */
        int CountNumberOfDihedrals(std::string parameter_file_path);
        std::vector<std::vector<std::string>> CreateAllAtomTypePermutationsforDihedralType(std::string atom_type1,
                                                                                           std::string atom_type2,
                                                                                           std::string atom_type3,
                                                                                           std::string atom_type4);
        std::vector<std::vector<std::string>>
        CreateAllAtomTypePermutationsforImproperDihedralType(std::string atom_type1, std::string atom_type2,
                                                             std::string atom_type3, std::string atom_type4);
        /*! \fn
         * A function that counts the number of dihedral types in all assemblies and residues of the assembly
         * @param parameter_file_path The path of the parameter file
         * @return type_list Number of dihedral types in all assemblies and residues in the current object of assembly
         */
        int CountNumberOfDihedralTypes(std::string parameter_file_path);
        /*! \fn
         * A function that returns all atoms with at least three neighbors in all assemblies and residues of the
         * assembly
         * @return atoms_with_at_least_three_neighbors A list of atoms with at least three neighbors in all assemblies
         * and residues in the current object of assembly
         */
        AtomVector GetAllAtomsOfAssemblyWithAtLeastThreeNeighbors();
        /*! \fn
         * A function that counts the number of excluded atoms in all assemblies and residues of the assembly
         * @return counter Number of excluded atoms in all assemblies and residues in the current object of assembly
         */
        int CountNumberOfExcludedAtoms();
        /*! \fn
         * A function that counts the number of atoms in largest residue in all assemblies and residues of the assembly
         * @return max Number of atoms in largest residue in all assemblies and residues in the current object of
         * assembly
         */
        int CountMaxNumberOfAtomsInLargestResidue();

        /*! \fn
         * A function to select a set of atoms in an assembly by using a string pattern
         * @param pattern Pattern of selection inclusing assembly number, residue name and/or atom name
         * Sample pattern: "1.*:#520,MAN,GAL@#3740-3750,^C:ALA@O1;1.1.2,1.1.3:NAG@O$"
         * @return A list of atoms which satisfy the condition specified in the pattern
         */
        AtomVector Select(std::string pattern);
        SelectPatternMap ParsePatternString(std::string pattern);
        void GetHierarchicalMapOfAssembly(HierarchicalContainmentMap& hierarchical_map, std::stringstream& index);
        /** @}*/
        void ClearAssembly();

        //            void CycleDetection();
        //            std::vector<std::vector<std::string> > CreateAllCyclePermutations(std::string id1, std::string
        //            id2, std::string id3, std::string id4, std::string id5, std::string id6);

        //            CycleMap DetectCyclesByDFS(std::string cycle_size = "5|6");

        /*! \fn
         * A function in order to access to the list of all residue from lib files
         * @param lib_files The list of paths to library files
         * @return all_residues_
         */
        LibraryFileSpace::LibraryFile::ResidueMap
        GetAllResiduesFromMultipleLibFilesMap(std::vector<std::string> lib_files);
        /*! \fn
         * A function in order to access to the list of all residue from prep files
         * @param prep_files The list of paths to prep files
         * @return all_residues_
         */
        PrepFileSpace::PrepFile::ResidueMap GetAllResiduesFromMultiplePrepFilesMap(std::vector<std::string> prep_files);

        /*! \fn
         * A function in order to extract all residue names existing in the given lib files
         * @param lib_files The list of paths to library files
         * @return all_residue_names
         */
        gmml::ResidueNameMap GetAllResidueNamesFromMultipleLibFilesMap(std::vector<std::string> lib_files);
        /** \addtogroup Manipulators
         * @{
         */

        void ObtainAllOligasaccharideResidues(Glycan::Oligosaccharide* oligo, ResidueVector& oligo_residues);
        void PutAglyconeInNewResidueAndRearrangeGlycanResidues(
            std::vector<Glycan::Oligosaccharide*> oligosaccharides,
            std::map<Glycan::Oligosaccharide*, ResidueVector>& oligo_residue_map);

        gmml::GlycamResidueNamingMap
        ExtractResidueGlycamNamingMap(OligosaccharideVector oligosaccharides,
                                      std::map<Glycan::Oligosaccharide*, std::vector<std::string>>& oligo_id_map,
                                      std::map<Glycan::Oligosaccharide*, ResidueVector>& oligo_residue_map);
        void ExtractOligosaccharideNamingMap(
            gmml::GlycamResidueNamingMap& pdb_glycam_map, Glycan::Oligosaccharide* oligosaccharide,
            CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree
                condensed_sequence_glycam06_residue_tree,
            int& index, std::map<Glycan::Oligosaccharide*, std::vector<std::string>>& oligo_id_map,
            Glycan::Oligosaccharide* root_oligo, std::map<Glycan::Oligosaccharide*, ResidueVector>& oligo_residue_map);
        void UpdateResidueName2GlycamName(gmml::GlycamResidueNamingMap residue_glycam_map, std::string prep_file);
        void TestUpdateResidueName2GlycamName(gmml::GlycamResidueNamingMap residue_glycam_map, std::string prep_file);
        void MatchPdbAtoms2Glycam(std::map<Glycan::Oligosaccharide*, ResidueVector>& oligo_residue_map,
                                  std::string prep_file,
                                  std::map<Glycan::Oligosaccharide*, Pdb2glycamMatchingTracker*>& match_tracker);
        // Avoid map ptr as key. Yao 20240826
        /*int RecursiveMoleculeSubgraphMatching(Atom* target_atom, AtomVector& target_atoms, Atom* template_atoms,
                                              std::map<Atom*, std::string>& target_atom_label_map,
                                              std::map<Atom*, std::string>& template_atom_label_map,
                                              std::vector<std::map<Atom*, Atom*>>& target_template_vertex_match,
                                              std::vector<std::map<Atom*, Atom*>>& template_target_vertex_match,
                                              std::vector<Atom*>& insertion_order, Pdb2glycamMatchingTracker* tracker,
                                              int depth);*/
        int RecursiveMoleculeSubgraphMatching(Atom* target_atom, AtomVector& target_atoms, int target_atom_index,
                                              Atom* template_atom, AtomVector& template_atoms, int template_atom_index,
                                              std::vector<std::string>& target_atom_label_map,
                                              std::vector<std::string>& template_atom_label_map,
                                              std::vector<std::map<Atom*, Atom*>>& target_template_vertex_match,
                                              std::vector<std::map<Atom*, Atom*>>& template_target_vertex_match,
                                              std::vector<Atom*>& insertion_order, Pdb2glycamMatchingTracker* tracker,
                                              int depth);
        bool CheckAndAcceptMatches(MolecularModeling::AtomVector& target_atoms,
                                   std::vector<std::map<Atom*, Atom*>>& target_template_vertex_match,
                                   std::vector<std::map<Atom*, Atom*>>& template_target_vertex_match,
                                   Pdb2glycamMatchingTracker* tracker);
        /*bool AllAtomEdgesMatch(Atom* target_atom, Atom* template_atom,
                               std::map<Atom*, std::string>& target_atom_label_map,
                               std::map<Atom*, std::string>& template_atom_label_map);*/
        bool AllAtomEdgesMatch(Atom* target_atom, AtomVector& target_atoms, Atom* template_atom,
                               AtomVector& template_atoms, std::vector<std::string>& target_atom_label_map,
                               std::vector<std::string>& template_atom_label_map);
        /*bool AtomVertexMatch(Atom* target_atom, Atom* template_atom,
                             std::map<Atom*, std::string>& target_atom_label_map,
                             std::map<Atom*, std::string>& template_atom_label_map);*/
        bool AtomVertexMatch(int target_atom_index, int template_atom_index,
                             std::vector<std::string>& target_atom_label_map,
                             std::vector<std::string>& template_atom_label_map);
        void Reorder06TreeTempDeleteLater(
            CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree& unordered_tree,
            CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree& ordered_tree,
            int current_tree_index, std::vector<int>& visited);
        bool IfMatchScenarioAlreadyExists(std::map<Atom*, Atom*>& target_template_vertex_match,
                                          Pdb2glycamMatchingTracker* tracker);
        void RemoveDownstreamMatches(Atom* target_atom,
                                     std::vector<std::map<Atom*, Atom*>>& target_template_vertex_match,
                                     std::vector<std::map<Atom*, Atom*>>& template_target_vertex_match,
                                     std::vector<Atom*> target_insertion_order);
        bool IfVertexAlreadyMatched(Atom* vertex_atom, std::vector<std::map<Atom*, Atom*>>& match_map);

        /// Pattern mathing
        bool PatternMatching(Residue* residue, ResidueVector query_residues, gmml::GlycamAtomNameMap& pdb_glycam_map,
                             gmml::GlycamAtomNameMap& glycam_atom_map);
        void CreateLabelGraph(Residue* residue, Residue* query_residue);
        void PruneLabelGraphByNeighboringLabels(Residue* query_residue);
        void CreatePrunedMatchingGraph(Residue* residue, ResidueVector query_residues);

        /*! \fn
         * A wrapper function for the next function. I do this to make the vector of monosaccharides external. I want to
         * do some manipulations to them.
         */
        OligosaccharideVector ExtractSugars(std::vector<std::string> amino_lib_files, bool glyporbity_report = false,
                                            bool populate_ontology = false, bool individualOntologies = false,
                                            std::string CCD_Path = " ");
        /*! \fn
         * A function in order to extract all the saccharide structures
         * @param amino_lib_files The list of paths to amino library files, used for identifying terminal residues
         * @param gyprobity_report A flag to prompt information for glyprobity report
         * @param populate_ontology A flag to prompt ontology population
         * @return oligosaccharides A list of extarcted oligosaccharide structures
         */
        OligosaccharideVector ExtractSugars(std::vector<std::string> amino_lib_files,
                                            std::vector<Glycan::Monosaccharide*>& monos, bool glyporbity_report = false,
                                            bool populate_ontology = false, bool individualOntologies = false,
                                            std::string CCD_Path = " ");

        std::vector<gmml_api_output> PDBExtractSugars(std::vector<std::string> amino_lib_files, std::string CCD_Path);

        /*! \fn
         * A funstion in order to initiate population of turtle formatted triples (subject-predicate-object) for
         * creating the GMMO ontology
         * @param outfile The output stream of the ontology file
         * @param oligosaccharides All the oligosaccharide structures that have been identified by the ExtractSugars
         * function
         */
        void PopulateOntology(std::ofstream& out_file, OligosaccharideVector oligosaccharides);
        /*! \fn
         * A function in order to populate the Note class of the ontology which is designed for storing the
         * issues/notes/problems identified in a PDB file
         * @param pdb_stream The output stream of PDB triples to be added to the main output stream
         * @param note_stream The output stream of Note triples to be added to the main output stream
         * @param pdb_uri The URI for the PDB instance to be used in the ontology. e.g http://gmmo.uga.edu/#3H32
         * @param notes The list of notes that have created by ExtractSugars function
         * @param id_prefix The specific prefix for the URIs related to aspecific PDB. e.g 3H32_
         * @param note_id The numeric id to be used for URI of a note
         */
        void PopulateNotes(std::stringstream& pdb_stream, std::stringstream& note_stream, std::string pdb_uri,
                           NoteVector notes, std::string id_prefix, int note_id);
        /*! \fn
         * A function in order to populate the Oligosaccharide class of the ontology
         * @param pdb_stream The output stream of PDB triples to be added to the main output stream
         * @param oligo_stream The output stream of Oligosaccharide triples to be added to the main output stream
         * @param mono_stream The output stream of Monosaccharide triples to be added to the main output stream
         * @param linkage_stream The output stream of Linkage class triples to be added to the main output stream
         * @param pdb_uri The URI for the PDB instance to be used in the ontology. e.g http://gmmo.uga.edu/#3H32
         * @param id_prefix The specific prefix for the URIs related to aspecific PDB. e.g 3H32_
         * @param link_id The numeric id to be used for URI of a link
         * @param oligosaccharides All the oligosaccharide structures that have been identified by the ExtractSugars
         * function
         * @param side_or_ring_atoms The list of side atoms and ring atoms of a monosaccharide
         * @param visited_oligos The list of oligos (monos) that have been visited and processed by traversing the tree
         * like structure of main oligosaccharide. each oligosaccharide has a core of monosaccharide. The collection of
         * linked oligosaccharides forms the main oligosaccharide structure.
         * @param mono_to_short_name_map The map containing short-names of monosaccharides (of current oligo)
         * @param oligo_to_res_uri_map The map containing URIs of monosaccharides (of current oligo)
         * @param root_oligo_id The id of the main oligosaccharide's core mono.
         */
        void PopulateOligosaccharide(std::stringstream& pdb_stream, std::stringstream& oligo_stream,
                                     std::stringstream& oligo_sequence_stream, std::stringstream& residue_stream,
                                     std::stringstream& mono_stream, std::stringstream& linkage_stream,
                                     std::string pdb_uri, std::string id_prefix, int& link_id,
                                     OligosaccharideVector oligos, std::vector<std::string>& side_or_ring_atoms,
                                     std::vector<int>& visited_oligos,
                                     std::map<std::string, std::string>& mono_to_short_name_map,
                                     std::map<std::string, std::string>& oligo_to_res_uri_map, int& root_oligo_id);
        /*! \fn
         * A function in order to populate the Linkage class of the ontology
         * @param linkage_stream The output stream of Linkage triples to be added to the main output stream
         * @param oligo_uri The URI for the Oligosaccharide instance to be used in the ontology. e.g
         * http://gmmo.uga.edu/#3H32_oligo1
         * @param parent_res_uri
         * @param child_res_uri
         * @param linkNum
         * @param thisLinkage
         * @param thisMono
         * @param thisMonoNeighbor
         * each oligosaccharide has a core of monosaccharide. The collection of linked oligosaccharides forms the main
         * oligosaccharide structure.
         */
        void PopulateLinkage(std::stringstream& linkage_stream, std::stringstream& oligo_stream, std::string oligo_uri,
                             std::string parent_mono_uri, std::string child_mono_uri, int& linkNum,
                             Glycan::GlycosidicLinkage* thisLinkage, Glycan::Monosaccharide* thisMono,
                             Glycan::Monosaccharide* thisMonoNeighbor);
        /*! \fn
         * A function in order to populate the Linkage class of the ontology
         * @param oligo_sequence_stream The output stream of SequenceLinkage triples to be added to the main output
         * stream
         * @param oligo An Assembly Oligosaccharide structure to be used to create linkage instances
         * @param oligo_uri The URI for the Oligosaccharide instance to be used in the ontology. e.g
         * http://gmmo.uga.edu/#3H32_oligo1
         * @param id_prefix The specific prefix for the URIs related to a specific PDB. e.g 3H32_
         * @param visited_oligos The list of oligos (monos) that have been visited and processed by traversing the tree
         * like structure of main oligosaccharide. each oligosaccharide has a core of monosaccharide. The collection of
         * linked oligosaccharides forms the main oligosaccharide structure.
         * @param mono_to_short_name_map The map containing short-names of monosaccharides (of current oligo)
         * @param oligo_to_res_uri_map The map containing URIs of monosaccharides (of current oligo)
         * @param root_oligo_id The id of the main oligosaccharide's core mono
         */
        void PopulateSequenceLinkage(std::stringstream& oligo_sequence_stream, Glycan::Oligosaccharide* oligo,
                                     std::string oligo_uri, std::string id_prefix, std::vector<int>& visited_oligos,
                                     std::map<std::string, std::string>& mono_to_short_name_map,
                                     std::map<std::string, std::string>& oligo_to_res_uri_map, int& root_oligo_id);

        /*! \fn
         * A function in order to populate the derivatives of the oligosaccharide sequence
         * @param oligo_sequence_stream The output stream of SequenceLinkage triples to be added to the main output
         * stream
         * @param mono_short_name The short-name of current monosaccharide (of current oligo)
         * @param oligo_uri The URI for the Oligosaccharide instance to be used in the ontology. e.g
         * http://gmmo.uga.edu/#3H32_oligo1
         * @param res_uri The URI for the current monosaccharide to be checked for derivatives
         */
        void CheckDerivativesAndPopulate(std::stringstream& oligo_sequence_stream, std::stringstream& residue_stream,
                                         std::string oligo_uri, std::string res_uri, Glycan::Monosaccharide* mono);

        /*! \fn
         * A function in order to check if the monosaccharide has derivates
         * @param mono_short_name The short-name of current monosaccharide (of current oligo)
         */
        bool hasDerivative(std::string mono_short_name);

        /*! \fn
         * A function in order to get the derivates of monosaccharide
         * @param mono_short_name The short-name of current monosaccharide (of current oligo)
         * @param derivatives The vector of derivatives of current monosaccharide (of current oligo)
         */
        void getDerivatives(std::string& mono_short_name, std::vector<std::string>& derivatives);

        /*! \fn
         * A function in order to extract the index of the carbon atom of of a monosaccharides that is linked to another
         * monosaccharide. 2 for DNeupNAca in DNeupNAca2-3DGalp
         * @param linkage_carbon_id The Assembly atom identifier of the carbon atom involved in a linkage
         * @return c_index The index of the carbon atom based on the naming conventions.
         */
        int ExtractLinkageCarbonIndex(Glycan::Oligosaccharide* oligo, std::string linkage_carbon_id);
        /*! \fn
         * A function in order to populate the Monosaccharide class of the ontology
         * @param oligo_stream The output stream of Oligosaccharide triples to be added to the main output stream
         * @param oligo_uri The URI for the Oligosaccharide instance to be used in the ontology. e.g
         * http://gmmo.uga.edu/#3H32_oligo1
         * @param id_prefix The specific prefix for the URIs related to aspecific PDB. e.g 3H32_
         * @param mono An Assembly Monosaccharide object to be used for populating the Monosaccharide triples
         * @param side_or_ring_atoms The list of side atoms and ring atoms of a monosaccharide
         */
        void PopulateMonosaccharide(std::stringstream& mono_stream, std::stringstream& oligo_stream,
                                    std::stringstream& pdb_stream, std::string oligo_uri, std::string id_prefix,
                                    Glycan::Monosaccharide* mono, std::vector<std::string>& side_or_ring_atoms,
                                    std::string pdb_uri, int numR);
        /*! \fn
         * A function in order to populate the RingAtom class of the ontology
         * @param ring_atom_stream The output stream of RingAtom triples to be added to the main output stream
         * @param id_prefix The specific prefix for the URIs related to aspecific PDB. e.g 3H32_
         * @param ring_uri The URI for the RingAtom instance to be used in the ontology. e.g
         * http://gmmo.uga.edu/#5BO9_C2_4752_SIA_A_409_n_n_1
         * @param ring_resource The ring atom resource in the ontology. e.g. 5BO9_C2_4752_SIA_A_409_n_n_1
         * @param ring_index The ring index of the atom. can be any of 1, 2, 3, 4, 5
         * @param mono An Assembly Monosaccharide object
         * @param ring_atom An Assembly Atom object to be used for populating the RingAtom triples
         * @param side_or_ring_atoms The list of side atoms and ring atoms of a monosaccharide
         */
        void PopulateRingAtom(std::stringstream& ring_atom_stream, std::string id_prefix, std::string ring_uri,
                              std::string ring_resource, unsigned int ring_index, Atom* ring_atom,
                              Glycan::Monosaccharide* mono, std::vector<std::string>& side_or_ring_atoms);
        /*! \fn
         * A function in order to populate the SideAtom class of the ontology
         * @param side_atom_stream The output stream of SideAtom triples to be added to the main output stream
         * @param id_prefix The specific prefix for the URIs related to aspecific PDB. e.g 3H32_
         * @param side_uri The URI for the SideAtom instance to be used in the ontology. e.g
         * http://gmmo.uga.edu/#5BO9_O3_4778_GAL_A_410_n_n_1
         * @param side_resource The ring atom resource in the ontology. e.g. 5BO9_O3_4778_GAL_A_410_n_n_1
         * @param side_index The index of the side atom. can be any of 1, 2, 3, 4, 5, 6, 7, 8 which will be stores as
         * any of -1, 1, 2, 3, 4, 5, +1, +2, +3 in the ontology ring atom with index 1 can have side atoms with indeces
         * -1 and 1. ring atom with index 4 in a furanose or 5 in a pyranose can have side atoms with indeces +1, +2 and
         * +3
         * @param side_atom An Assembly Atom object to be used for populating the SideAtom triples
         * @param mono An Assembly Monosaccharide object
         * @param side_or_ring_atoms The list of side atoms and ring atoms of a monosaccharide
         */
        void PopulateSideAtom(std::stringstream& side_atom_stream, std::string id_prefix, std::string side_uri,
                              std::string side_resource, int ring_index, int side_index, Atom* side_atom,
                              Glycan::Monosaccharide* mono, std::vector<std::string>& side_or_ring_atoms);
        /*! \fn
         * A function in order to populate the SugarName class of the ontology
         * @param mono_stream The output stream of SugarName triples to be added to the main output stream
         * @param mono_uri The URI for the Monosaccharide instance to be used in the ontology. e.g
         * http://gmmo.uga.edu/#5BO9_mono1
         * @param id_prefix The specific prefix for the URIs related to aspecific PDB. e.g 3H32_
         * @param mono_id The numeric id to be used as part of the URI of a SugarName
         * @param sugar_name An Assembly SugarName object to be used for populating SugarName class of the ontology
         */
        void PopulateSugarName(std::stringstream& mono_stream, std::string mono_uri, std::string id_prefix, int mono_id,
                               Glycan::SugarName sugar_name);
        /*! \fn
         * A function in order to populate the Residue class of the ontology
         * @param pdb_stream The output stream of PDB triples to be added to the main output stream
         * @param residue_stream The output stream of Residue triples to be added to the main output stream
         * @param pdb_uri The URI for the PDB instance to be used in the ontology. e.g http://gmmo.uga.edu/#3H32
         * @param id_prefix The specific prefix for the URIs related to aspecific PDB. e.g 3H32_
         * @param residues A list of Assembly residue objects to be used for populating the Residue class of the
         * ontology
         * @param sugar_name An Assembly SugarName object to be used for populating SugarName class of the ontology
         * @param side_or_ring_atoms The list of side atoms and ring atoms of a monosaccharide
         */
        void PopulateResidue(std::stringstream& pdb_stream, std::stringstream& residue_stream, std::string pdb_uri,
                             std::string id_prefix, ResidueVector residues,
                             std::vector<std::string> side_or_ring_atoms);
        /*! \fn
         * A function in order to populate the Atom class of the ontology. This function populates the atoms that
         * haven't been populated by PopulateRingAtom and PopulateSideAtom functions
         * @param atom_stream The output stream of Atom triples to be added to the main output stream
         * @param atom_uri The URI for the Atom instance to be used in the ontology.
         * @param atom_resource The atom resource in the ontology.
         * @param id_prefix The specific prefix for the URIs related to aspecific PDB. e.g 3H32_
         * @param atom An Assembly Atom object to be used for populating the Atom triples
         */
        void PopulateAtom(std::stringstream& atom_stream, std::string atom_uri, std::string atom_resource,
                          std::string id_prefix, Atom* atom);
        /*! \fn
         * A function in order to create a title for the specific PDB file triples
         * @param pdb_resource The PDB resource which the title is created for
         * @param pdb_stream The output stream of PDB triples to be added to the main output stream
         */
        void CreateTitle(std::string pdb_resource, std::stringstream& pdb_stream);
        void CreateSubtitle(std::string pdb_resource, std::stringstream& pdb_stream);
        // /*! \fn
        // * A function in order to create a turtle formatted triple (subject predicate object) and appending it to the
        // output file stream
        // * @param s The subject part of the triple
        // * @param p The predicate part of the triple
        // * @param o The object part of the triple
        // * @param stream The output stream which is going to be written in the ontology turtle file
        // */
        // void AddTriple(std::string s, std::string p, std::string o, std::stringstream& stream);
        // /*! \fn
        // * A function in order to create a turtle formatted triple (subject predicate object=literal value) and
        // appending it to the output file stream
        // * @param s The subject part of the triple
        // * @param p The predicate part of the triple
        // * @param o The object part of the triple, the object is not a resource in this case. It can be a literal
        // value e.g. string , int
        // * @param stream The output stream which is going to be written in the ontology turtle file
        // */
        // void AddLiteral(std::string s, std::string p, std::string o, std::stringstream& stream);
        // /*! \fn
        // * A function in order to create a turtle formatted triple (subject predicate object=literal value) and
        // appending it to the output file stream
        // * @param s The subject part of the triple
        // * @param p The predicate part of the triple
        // * @param o The object part of the triple, the object is not a resource in this case. It can be a literal
        // value e.g. string , int
        // * @param stream The output stream which is going to be written in the ontology turtle file
        // */
        // void AddDecimal(std::string s, std::string p, float o, std::stringstream& stream);
        /*! \fn
         * A function in order to create the an ontology resource based on the given resource type
         * @param resource The resource type e.g. PDB, Residue, Atom
         * @param number The numeric id that shoud be included in the URI resource
         * @param id_prefix The specific prefix for the URIs related to aspecific PDB. e.g 3H32_
         * @param id The resource that is going to be created e.g. O3_4778_GAL_A_410_n_n_1
         * @return uri_resource The resource that has been created by the function
         */
        std::string CreateURIResource(gmml::URIType resource, int number, std::string id_prefix, std::string id);
        /*! \fn
         * A function in order to create the unique URI for a resource of the ontology
         * @param uri_resource The resource that has been created by the CreateURIResource function
         * @return uri The unique URI that has been created by the function
         */
        std::string CreateURI(std::string uri_resource);

        /*! \fn
         * A function in order to formulate and execute a cURL command based on the given SPARQL query
         * @param output_file_type The format of the result from query execution. e.g. csv, json, xml
         * @param query The SPARQL query that is going to be used in the cURL command
         */
        void FormulateCURL(std::string output_file_type, std::string query);
        /*! \fn
         * A function in order to extract information from ontology based on the name of the sugar
         * @param stereo_name The stereochemistry name of the sugar e.g. b-D-mannopyranose
         * @param stereo_condensed_name The condensed version of stereochemistry name of the sugar e.g. DManpb
         * @param name The complete name of the sugar e.g. 2-sulfo-b-D-mannopyranose
         * @param condensed_name The condensed version of the complete name of the sugar e.g. DManp[2s]b
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        std::string FormulateCURLGF(std::string output_file_type, std::string query, std::string url);
        /*! \fn
         * A function used in GlyFinder project in order to extract information from ontology based on the name of the
         * sugar
         * @param stereo_name The stereochemistry name of the sugar e.g. b-D-mannopyranose
         * @param stereo_condensed_name The condensed version of stereochemistry name of the sugar e.g. DManpb
         * @param name The complete name of the sugar e.g. 2-sulfo-b-D-mannopyranose
         * @param condensed_name The condensed version of the complete name of the sugar e.g. DManp[2s]b
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        void ExtractOntologyInfoByNameOfGlycan(std::string stereo_name, std::string stereo_condensed_name,
                                               std::string name, std::string condensed_name,
                                               std::string output_file_type = "csv");
        /*! \fn
         * A function in order to extract information from ontology based on the different parts of the sugar name
         * @param isomer The isomer(D/L) part of the monosaccharide name
         * @param ring_type The ring type(p/f) part of the monosaccharide name
         * @param name The complete name of the sugar e.g. 2-sulfo-b-D-mannopyranose
         * @param configuration The configuration(a/b/x) part of the monosaccharide name
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        void ExtractOntologyInfoByNamePartsOfGlycan(std::string isomer, std::string ring_type,
                                                    std::string configuration, std::string output_file_type = "csv");
        /*! \fn
         * A function in order to extract information from ontology based on a specific PDB ID
         * @param pdb_id The unique PDB identifier
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        void ExtractOntologyInfoByPDBID(std::string pdb_id, std::string output_file_type = "csv");
        /*! \fn
         * A function in order to extract information from ontology based on a specific chemical code
         * inspired http://glycam.org/docs/gmml/2016/03/31/glycode-internal-monosaccharide-representation
         * @param chemical_code The chemical code structure of the sugar e.g. _2^3^4P_a^+1
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        std::string ExtractOntologyInfoByPDBIDGF(std::string pdb_id, std::string url,
                                                 std::string output_file_type = "csv");
        /*! \fn
         * A function used in GlyFinder project in order to extract information from ontology based on a specific
         * chemical code inspired http://glycam.org/docs/gmml/2016/03/31/glycode-internal-monosaccharide-representation
         * @param chemical_code The chemical code structure of the sugar e.g. _2^3^4P_a^+1
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        void ExtractOntologyInfoByStringChemicalCode(std::string chemical_code, std::string output_file_type = "csv");
        /*! \fn
         * A function in order to extract information from ontology based on a specific oligosaccharide name
         * @param oligo_name The name of the oligosaccharide. e.g. DNeupNAca2-3DGalpb1-4DGlcpNAcb1-ROH
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        void ExtractOntologyInfoByOligosaccharideNameSequence(std::string oligo_name,
                                                              std::string output_file_type = "csv");
        /*! \fn
         * A function in order to extract information from ontology based on a given specific pattern
         * @param oligo_name_pattern The oligosaccharide pattern that is going to be used in the query
         * e.g. DGlcpNAcb1-4DGlc*, *b1-4L*, *GlcpNAcb1-4DGlcpNAcb, DGlcpNAcb*4DGlcpNAca, *DGlcpNAcb1-4DGlc*,
         * DGlcpNAcb*DGlc*, *DManpa1-6[DManpa1-2DManpa1-3]D*
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        void ExtractOntologyInfoByOligosaccharideNameSequenceByRegex(std::string oligo_name_pattern,
                                                                     std::string output_file_type = "csv");
        /*! \fn
         * A function in order to extract information from ontology based on the orientations of the side atoms of a
         * monosaccharide structure
         * @param ring_type The ring type(p/f) part of the monosaccharide name
         * @param anomeric_orientation The orientation of the side oxygen attached to anomeric carbon of the ring
         * @param minus_one_orientation The orientation of the side carbon attached to anomeric carbon of the ring
         * @param index_two_orientation The orientation of the side oxygen attached to second carbon of the ring
         * @param index_three_orientation The orientation of the side oxygen attached to third carbon of the ring
         * @param index_four_orientation The orientation of the side oxygen attached to fourth carbon of the ring
         * @param plus_one_orientation The orientation of the side carbon attached to last carbon of the ring
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        std::string ExtractOntologyInfoByOligosaccharideNameSequenceGF(std::string oligo_name,
                                                                       std::string output_file_type = "csv");
        /*! \fn
         * A function used in GlyFinder project in order to extract information from ontology based on a given specific
         * pattern
         * @param oligo_name_pattern The oligosaccharide pattern that is going to be used in the query
         * e.g. DGlcpNAcb1-4DGlc*, *b1-4L*, *GlcpNAcb1-4DGlcpNAcb, DGlcpNAcb*4DGlcpNAca, *DGlcpNAcb1-4DGlc*,
         * DGlcpNAcb*DGlc*, *DManpa1-6[DManpa1-2DManpa1-3]D*
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        std::string ExtractOntologyInfoByOligosaccharideNameSequenceByRegexGF(std::string oligo_name_pattern,
                                                                              std::string url,
                                                                              std::string output_file_type = "csv");
        /*! \fn
         * A function used in GlyFinder project in order to extract information from ontology based on the orientations
         * of the side atoms of a monosaccharide structure
         * @param ring_type The ring type(p/f) part of the monosaccharide name
         * @param anomeric_orientation The orientation of the side oxygen attached to anomeric carbon of the ring
         * @param minus_one_orientation The orientation of the side carbon attached to anomeric carbon of the ring
         * @param index_two_orientation The orientation of the side oxygen attached to second carbon of the ring
         * @param index_three_orientation The orientation of the side oxygen attached to third carbon of the ring
         * @param index_four_orientation The orientation of the side oxygen attached to fourth carbon of the ring
         * @param plus_one_orientation The orientation of the side carbon attached to last carbon of the ring
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        void ExtractOntologyInfoByGlycanStructure(std::string ring_type, std::string anomeric_orientation,
                                                  std::string minus_one_orientation, std::string index_two_orientation,
                                                  std::string index_three_orientation,
                                                  std::string index_four_orientation = "",
                                                  std::string plus_one_orientation   = "",
                                                  std::string output_file_type       = "csv");
        /*! \fn
         * A function in order to extract information from ontology based on a specific modification/derivative pattern
         * @param ring_type The ring type(p/f) part of the monosaccharide name
         * @param derivative_modification_map A mapping between the monosaccharide's atom position/index and the
         * derivative/modification attached to it e.g. derivative_modification_map['2'] = 'xC-N-C=OCH3'
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        void ExtractOntologyInfoByDerivativeModificationMap(std::string ring_type,
                                                            DerivativeModificationMap derivative_modification_map,
                                                            std::string output_file_type = "csv");
        /*! \fn
         * A function in order to extract information from ontology based on the sugar structures
         * @param attached_structures The orientation of side atoms of attache monosaccharide
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        void ExtractOntologyInfoByAttachedGlycanStructures(AttachedGlycanStructuresVector attached_structures,
                                                           std::string output_file_type = "csv");
        /*! \fn
         * A function in order to extract information from ontology based on identified note(s) of PDB file(s)
         * @param pdb_file The PDB ID that the notes have been extracted from
         * @param note_type The type of the note. e.g. error, warning, comment
         * @param note_category The category of the note. e.g. anomeric, residue name, monosaccharide
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        void ExtractOntologyInfoByNote(std::string pdb_file, std::string note_type, std::string note_category,
                                       std::string output_file_type = "csv");
        /*! \fn
         * A function in order to extract information from ontology based on the given SPARQL query
         * @param query The custom SPARQL query
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        void ExtractOntologyInfoByCustomQuery(std::string query, std::string output_file_type = "csv");
        /*! \fn
         * A function in order to extract information from ontology based on user input
         * @param searchType The type of search submitted from GlyFinder
         * @param searchTerm The search term submitted from GlyFinder
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        std::string QueryOntology(std::string searchType, std::string searchTerm, float resolution_min,
                                  float resolution_max, float b_factor_min, float b_factor_max,
                                  float oligo_b_factor_min, float oligo_b_factor_max, int isError = -1,
                                  int isWarning = -1, int isComment = -1, int isLigand = -1, int isGlycomimetic = -1,
                                  int isNucleotide = -1, std::string aglycon = "", std::string count = "False",
                                  int page = 1, int resultsPerPage = 25, std::string sortBy = "", std::string url = "",
                                  std::string output_file_type = "csv");

        std::string MoreQuery(std::string pdb_id, std::string oligo_sequence, std::string oligo, std::string url,
                              std::string output_file_type);

        std::string ontologyDownload(std::string searchType, std::string searchTerm, float resolution_min,
                                     float resolution_max, float b_factor_min, float b_factor_max,
                                     float oligo_b_factor_min, float oligo_b_factor_max, int isError, int isWarning,
                                     int isComment, int isLigand, int isGlycomimetic, int isNucleotide,
                                     std::string aglycon, std::string count, int page, int resultsPerPage,
                                     std::string sortBy, std::string url, std::string output_file_type);

        std::string ontologyPDBDownload(std::string searchType, std::string searchTerm, float resolution_min,
                                        float resolution_max, float b_factor_min, float b_factor_max,
                                        float oligo_b_factor_min, float oligo_b_factor_max, int isError, int isWarning,
                                        int isComment, int isLigand, int isGlycomimetic, int isNucleotide,
                                        std::string aglycon, std::string count, int page, int resultsPerPage,
                                        std::string sortBy, std::string url, std::string output_file_type);

        void ConvertGraphToQuery(std::stringstream& queryStream, GraphDS::Graph queryGraph);
        GraphDS::Graph CreateQueryStringGraph(std::string queryString);
        void ConnectNodes(int start, int end, std::vector<parsedString>& parsedVector, GraphDS::Graph& graph);

        /*! \fn
         * A function in order to extract necessary atom coordinates from ontology to calculate phi/psi/omega torsion
         * angles
         * @param disaccharide_pattern The disaccharide pattern that is going to be searched in ontology
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        void ExtractAtomCoordinatesForTorsionAnglesFromOntologySlow(
            std::string disaccharide_pattern /*, std::string output_file_type = "csv"*/);
        /*! \fn
         * A function in order to extract necessary atom coordinates from ontology to calculate phi/psi/omega torsion
         * angles
         * @param disaccharide_pattern The disaccharide pattern that is going to be searched in ontology
         * @param output_file_type The format of the result to expect from query execution. e.g. csv, json, xml
         */
        void ExtractAtomCoordinatesForTorsionAnglesFromOntologyFast(
            std::string disaccharide_pattern /*, std::string output_file_type = "csv"*/);
        /*! \fn
         * A function in order to calculate torsion angles based on the result of the query for extracting cooridnates
         * from ontology
         */
        void ExtractTorsionAnglesFromSlowQueryResult();
        /*! \fn
         * A function in order to calculate torsion angles based on the result of the query for extracting cooridnates
         * from ontology
         */
        void ExtractTorsionAnglesFromFastQueryResult();
        /*! \fn
         * A function in order to extract atom coordinates from ontology based on the given arguments, calculate the
         * bond lenghts, mean and standard deviation
         * @param atom_name1 The name of the first atom
         * @param atom_name1 The name of the second atom
         * @param is_atom2_ring A boolean value indicating if the second atom is a side (exocyclic) atom
         * @param mono_name The name of the monosaccharide which contains the given atoms
         * @return statistics A list of calculated statistics. Mean and standard deviation
         */
        std::vector<double> CalculateBondlengthsStatisticsBasedOnOntologyInfo(std::string atom_name1,
                                                                              std::string atom_name2,
                                                                              bool is_atom2_ring,
                                                                              std::string mono_name);
        /*! \fn
         * A function in order to extract atom coordinates from ontology based on the given arguments, calculate the
         * bond angles, mean and standard deviation
         * @param atom_name1 The name of the first atom
         * @param atom_name1 The name of the second atom
         * @param atom_name3 The name of the third atom
         * @param is_atom2_ring A boolean value indicating if the second atom is a side (exocyclic) atom
         * @param mono_name The name of the monosaccharide which contains the given atoms
         * @return statistics A list of calculated statistics. Mean and standard deviation
         */
        std::vector<double> CalculateBondAnglesStatisticsBasedOnOntologyInfo(std::string atom_name1,
                                                                             std::string atom_name2,
                                                                             std::string atom_name3, bool is_atom3_ring,
                                                                             std::string mono_name);
        /*! \fn
         * A function in order to extract torsion angles from a PDB file for a given disaccharide pattern
         * @param amino_lib_files The list of paths to amino library files to process PDB file
         * @param disaccharide The disaccharide pattern that is going to be used to extract torsio angles
         */
        void ExtractTorsionAnglesFromPDB(std::vector<std::string> amino_lib_files, std::string disaccharide);
        /*! \fn
         * A function in order to check if a parent and child oligosaccharide matches the given monosaccharide names and
         * linkage indeces if the values matches, the function calculates the phi/psi angle(s)
         * @param oligo An Assembly Oligosaccharide object that is going to be checked for matching the given values
         * @param phi_angle The phi angle that is going to be calculated by this function
         * @param phi_angle The phi angle that is going to be calculated by this function
         * @param first_mono The first monosaccharide in the disaccharide pattern which is the child monosaccharide in
         * the tree-like structure of the main oligosaccharide e.g. DGalpb in DNeupNAca2-3DGalpb
         * @param second_mono The second monosaccharide in the disaccharide pattern which is the parent monosaccharide
         * in the tree-like structure of the main oligosaccharide e.g. DNeupNAca in DNeupNAca2-3DGalpb
         * @return true if disaccharide pattern was found and false otherwise.
         */
        bool MatchDisaccharide(std::queue<Glycan::Oligosaccharide*> oligo, double& phi_angle, double& psi_angle,
                               std::string first_mono, char mono1_carbon_index, std::string second_mono,
                               char mono2_carbon_index);
        /*! \fn
         * A function in order to calculate bond angles based on the coordinate objects
         * @param atom1_crd The geometric coordinate of the first atom of the bond angle
         * @param atom2_crd The geometric coordinate of the second atom of the bond angle
         * @param atom3_crd The geometric coordinate of the third atom of the bond angle
         * @return The calculated bond angle (radian)
         */
        double CalculateBondAngleByCoordinates(GeometryTopology::Coordinate* atom1_crd,
                                               GeometryTopology::Coordinate* atom2_crd,
                                               GeometryTopology::Coordinate* atom3_crd);
        /*! \fn
         * A function in order to calculate bond angles based on the coordinate objects
         * @param atom1 The first atom of the bond angle
         * @param atom2 The second atom of the bond angle
         * @param atom3 The third atom of the bond angle
         * @return The calculated bond angle (radian)
         */
        double CalculateBondAngleByAtoms(Atom* atom1, Atom* atom2, Atom* atom3);
        /*! \fn
         * A function in order to calculate torsion angles based on the coordinate objects
         * @param atom1_crd The geometric coordinate of the first atom of the torsion angle
         * @param atom2_crd The geometric coordinate of the second atom of the torsion angle
         * @param atom3_crd The geometric coordinate of the third atom of the torsion angle
         * @param atom4_crd The geometric coordinate of the fourth atom of the torsion angle
         * @return current_dihedral The calculated torsion angle (radian)
         */
        //  double CalculateTorsionAngleByCoordinates(GeometryTopology::Coordinate* atom1_crd,
        //  GeometryTopology::Coordinate* atom2_crd,
        //                                            GeometryTopology::Coordinate* atom3_crd,
        //                                            GeometryTopology::Coordinate* atom4_crd);
        /*! \fn
         * A function in order to calculate torsion angles based on the Assembly atom objects
         * @param atom1 The first atom of the torsion angle
         * @param atom2 The second atom of the torsion angle
         * @param atom3 The third atom of the torsion angle
         * @param atom4 The fourth atom of the torsion angle
         * @return current_dihedral The calculated torsion angle (radian)
         */
        // OG Feb 2021. DO NOT USE CalculateTorsionAngleByAtoms. Only used by guesses.cc for now until bug sorted. USE
        // GeometryTopology::CalculateDihedralAngle
        double CalculateTorsionAngleByAtoms(Atom* atom1, Atom* atom2, Atom* atom3, Atom* atom4);
        /*! \fn
         * A function in order to calculate a torsion angle matrix based on a file that contains torison angles
         * calculated from PDB file(s) the matrix shows the number of angles that falls in a certain scope, the scope is
         * calculated based on the given range
         * @param torsion_file The files that contains all the calculated torsion angles for PDB file(s)
         * @return low_range The lowest number of the range
         * @return high_range The highest number of the range
         */
        void CalculateTorsionStatistics(std::string torsion_file = "", int low_range = -180, int high_range = 180);
        /*! \fn
         * A function in order to calculate the average B Factor for an Oligosaccharide
         * @param oligo The oligosaccharide
         * @param monos The monosaccharide vector that makes up the oligosaccharide
         */
        void CalculateOligosaccharideBFactor(Glycan::Oligosaccharide* oligo,
                                             std::vector<Glycan::Monosaccharide*> monos);
        double CalculateOmegaAngle(Glycan::Oligosaccharide* parent_oligo, std::string parent_atom_id,
                                   std::string glycosidic_atom_id);
        double CalculatePsiAngle(Glycan::Oligosaccharide* child_oligo, std::string parent_atom_id,
                                 std::string child_atom_id, std::string glycosidic_atom_id);
        double CalculatePhiAngle(Glycan::Oligosaccharide* child_oligo, std::string parent_atom_id,
                                 std::string child_atom_id, std::string glycosidic_atom_id);
        bool guessIfC_CDoubleBond(MolecularModeling::Atom* carbon1, MolecularModeling::Atom* carbon2);
        std::pair<double, double> guessBondLengthByAtomType(MolecularModeling::Atom* atom1,
                                                            MolecularModeling::Atom* atom2);
        //    std::pair<double,double> guessBondLengthByAtomType(MolecularModeling::Atom* atom1,
        //    MolecularModeling::Atom* atom2);
        void GetAuthorNaming(std::vector<std::string> amino_lib_files, Glycan::Monosaccharide* mono,
                             std::string CCD_Path);
        bool checkIfNucleotide(Glycan::Monosaccharide* mono);

        /*! \fn
         * A function in order to extract and print out all saccharides ring atoms information
         */
        /** @}*/
        /** \addtogroup Verifiers_and_Issue_Resolvers
         * @{
         */
        // void ExtractRingAtomsInformation();
        /*! \fn
         * A function in order to detect cycles in the molecular graph using the exhaustive ring perception algorithm
         * The algorithm is derived from http://pubs.acs.org/doi/pdf/10.1021/ci960322f
         * @return cycles A map between the string version of atoms of cycles and the list of cycle atom objects
         */
        Assembly::CycleMap DetectCyclesByExhaustiveRingPerception();
        /*! \fn
         * A function in order to prune the graph (recursively removing nodes with zero or 1 neighbors)
         * @param all_atoms The list of atoms of the graph which is going to be updated by the function
         */
        void PruneGraph(AtomVector& all_atoms);
        /*! \fn
         * A function in order to convert the graph into a path graph (creating list of edges between the nodes and a
         * list of labels for those edges )
         * @param path_graph_edges The list of edges between the nodes to be filled by the function
         * @param path_graph_labels The list of edge labels to be filled by the function
         * @param atoms The list of atoms of the graph
         */
        void ConvertIntoPathGraph(std::vector<std::string>& path_graph_edges,
                                  std::vector<std::string>& path_graph_labels, AtomVector atoms);
        /*! \fn
                  * A function in order to reduce the path graph such that if there is a path/walk a-b-c in the graph
           converting it to a-c and creating a new label for the new edge and checking if the new edge makes a cycle
                  * @param path_graph_edges The list of edges between the nodes
                  * @param path_graph_labels The list of edge labels
                  * @param reduced_path_graph_edges The list of edges of the (reduced) path graph
                  * @param reduced_path_graph_labels The list of edge labels of the (reduced) path graph
                  * @param common_atom The atom that needs to be checked if it is involved in a walk (a path like a-b-c)
                  * @param cycles The list of cycles to be filled by the function
                  */
        void ReducePathGraph(std::vector<std::string> path_graph_edges, std::vector<std::string> path_graph_labels,
                             std::vector<std::string>& reduced_path_graph_edges,
                             std::vector<std::string>& reduced_path_graph_labels, std::string common_atom,
                             std::vector<std::string>& cycles);

        /*! \fn
         * A function in order to detect cycles in the molecular graph using depth first search algorithm
         * @return cycles A map between the string version of atoms of cycles and the list of cycle atom objects
         */
        CycleMap DetectCyclesByDFS();
        /*! \fn
         * A function of depth first search algorithm in order to traverse the graph
         * @param path_graph_edges The list of edges between the nodes to be filled by the function
         * @param path_graph_labels The list of edge labels to be filled by the function
         * @param atoms The list of atoms of the graph
         */
        void DFSVisit(AtomVector atoms, AtomStatusMap& atom_status_map, AtomIdAtomMap& atom_parent_map, Atom* atom,
                      int& counter, AtomIdAtomMap& dest_srd_map);
        /*! \fn
         * A recursive function in order to back track a path from current atom to a source atom to return the atom
         * objects using the information extracted by the DFS algorithm
         * @param src_id The identifier of the source atom
         * @param current_atom The current atom that the function has been called on
         * @param atom_parent_map A map between atom identifier and its parent atom object
         * @param cycle The list of atom objects involved in a cycle
         * @param cycle_stream The cycle path that has been back traversed so far
         */
        void ReturnCycleAtoms(std::string src_id, Atom* current_atom, AtomIdAtomMap& atom_parent_map, AtomVector& cycle,
                              std::stringstream& cycle_stream);
        /*! \fn
         * A function in order to discard rings/cycles that have double bonds
         * @param cycles A map between the string version of atoms of cycles and the list of cycle atom objects
         */
        CycleMap FilterCyclesWithDoubleBonds(CycleMap& cycles);
        /*! \fn
         * A function in order to discard rings/cycles that are only made from carbons atoms
         * @param cycles A map between the string version of atoms of cycles and the list of cycle atom objects
         */
        void FilterAllCarbonCycles(CycleMap& cycles);
        /*! \fn
         * A function in order to discard cycles/rings that are sharing an edge (fused cycles)
         * @param cycles A map between the string version of atoms of cycles and the list of cycle atom objects
         */
        void RemoveFusedCycles(CycleMap& cycles);
        /*! \fn
         * A function in order to detect the anomeric carbon of the ring (the carbon which has two oxygon neighbors)
         * @param anomeric_carbons_note The Note for the anomeric carbon to be filled by the function
         * @param anomeric_carbons_status The detection status of the anomeric carbon to be filled by the function
         * @param cycle The list of cycle atoms
         * @param cycle_atom_str The string version of atom identifiers of the cycle
         * @return o_neighbor The anomeric carbon which is the neighbor of oxygen
         */
        // Moved to Monosaccharide class
        //  Atom* FindAnomericCarbon(Glycan::Note* anomeric_carbons_note, std::vector<std::string>&
        //  anomeric_carbons_status, AtomVector cycle, std::string cycle_atoms_str);
        /*! \fn
         * A function in order to sort atom objects of the cycle starting from the anomeric carbon of the ring (ring
         * oxygen will be last atom)
         * @param cycle The list of cycle atoms
         * @param anomeric_atom The anomeric carbon of the ring
         * @return sorted_cycle The sorted list of cycle atom objects
         */
        AtomVector SortCycle(AtomVector cycle, Atom* anomeric_atom, std::stringstream& sorted_cycle_stream);
        /*! \fn
         * A function in order to calculate geometry outliers for glyprobity report (e.g. bond lengths, bond angles,
         * torsion angles)
         * @param mono The monosaccharide object which is processed by this function to calculate its outliers
         */
        void CalculateGlyprobityGeometryOutliers(Glycan::Monosaccharide* mono);
        /** @}*/
        /*! \fn
         * A function in order to extract the relative oriantation of the side atoms that are attached to ring atoms
         * against the ring
         * @param mono The monosaccharide object
         * @param sorted_cycle_stream The sorted atom of the cycle so far (to be filled with the fuction)
         * @param cycle_atom_str The string version of atom identifiers of the cycle
         * @return orientations The list of side atoms orinetations
         */
        // Moved to monosaccharide class
        //  std::vector<std::string> GetSideGroupOrientations(Glycan::Monosaccharide* mono, std::string
        //  cycle_atoms_str);

        /*! \fn
         *  A function to start complete side group atoms detection from all side chains detected.
         * @param mono A vector of all monosaccharides detected in input file.
         */
        // Moved to monosaccharide class
        //  void InitiateDetectionOfCompleteSideGroupAtoms (std::vector<Glycan::Monosaccharide*> monos);
        /*! \fn
         * A function in order to check if plus one side atom belongs to the current monosaccharide ,called within
         * GetSideGroupOrientations.
         * @param SideAtomArm,a reference to the vector element housing this side chain atoms(i.e.
         * Monosaccharide.side_atoms_).
         * @param working_atom The atom whose node neighbors are currently being checked for attached anomeric carbons.
         */
        // Moved to monosaccharide class
        //  bool CheckIfPlusOneSideAtomBelongsToCurrentMonosaccharide(AtomVector& SideAtomArm, AtomVector & cycle_atoms,
        //  Atom* working_atom);
        /*! \fn
         * A function in order to obtain complete side group atoms in a monosaccharide,involving recursive calls. Called
         * within GetSideGroupOrientations.
         * @param SideAtomArm, a reference to the vector element housing this side chain atoms(i.e.
         * Monosaccharide.side_atoms_).
         * @param working_atom, the atom whose node neighbors are currently being checked for new side group atoms.
         * @param cycle_atoms, AtomVector of the ring atoms of the current monosaccharide, used to prevent recursion
         * from proceeding towards the ring.
         * @param visited_atoms, AtomVector of atoms that have served as working atoms, used to prevent recursion from
         * going backwards.
         */
        // Moved to monosaccharide class
        //  void SetCompleteSideGroupAtoms(AtomVector& SideAtomArm, Atom* working_atom, AtomVector & cycle_atoms,
        //  AtomVector & visited_atoms);
        /*! \fn
         * A function to make all atoms of a monosaccharide a new residue, replacing the corresponding old one in input
         * file. This is to solve the problem wher one residue contains more than one sugar.
         * @param monos a vector of Monosaccharide*, containing all identified monosaccharides in the input files.
         */
        void UpdateMonosaccharides2Residues(std::vector<Glycan::Monosaccharide*>& monos);
        /*! \fn
         * A function in order to check if the current side chain in the previous function is a terminal side chain.
         * @param starting_atom A pointer to the atom where recursive checking initiates.
         * @param current_cycle_and_visited_atoms The cycle whose side chains are being extracted.Recursion shall not
         * proceed towards it.
         */
        // Moved to monosaccharide class
        //  void CheckIfSideChainIsTerminal(Atom* starting_atom, AtomVector & current_cycle_and_visited_atoms, bool &
        //  is_terminal);

        /** \addtogroup Manipulators
         * @{
         */
        /*! \fn
         * A function in order to build a chemical code structure for the monosaccharide based on the side atom
         * orientations
         * @param orientations The list of side atoms orinetations
         * @return code The chemical code of the monosaccharide structure
         */
        // moved to monosaccharide class
        //  Glycan::ChemicalCode* BuildChemicalCode(std::vector<std::string> orientations);
        /*! \fn
         * A function in order to extract additional side atoms of the saccharide (+2 and +3 positions)
         * @param mono The monosaccharide object
         * @return plus_sides The list of additional side atoms of the monosaccharide
         */
        // Moved to monosaccharide class
        //  AtomVector ExtractAdditionalSideAtoms(Glycan::Monosaccharide* mono);
        /*! \fn
         * A function in order to extract the probable derivatives that are attached to the side atoms of the ring and
         * sets the derivative_map attribute of the monosaccharide Entry examples for the derivative map: [-1,
         * derivative pattern] [a, derivative pattern] [2, derivative pattern] ... [+1, derivative pattern] [+2,
         * derivative pattern]
         * @param mono The monosaccharide object
         */
        // Moved to monosaccharide class
        //  void ExtractDerivatives(Glycan::Monosaccharide* mono);
        /*! \fn
         * A function in order to generate a complete name for the monosaccharide structure based on its derivatives
         * @param mono The monosaccharide object
         */
        // Moved to monosaccharide class
        //  void GenerateCompleteSugarName(Glycan::Monosaccharide* mono);
        /*! \fn
         * A function in order to add the modification info to the name of the sugar based on the first group of rules
         * (No Bracket -> 2(r:6&!-1), Warning position -> a, Error Position -> 5(r6),4(r5) )
         * @param key The index of the atom which is a key in the derivative map
         * @param pattern The pattern of the modification which is based on the value in the derivative map
         * @param mono The monosaccharide object
         * @param long_name_pattern The long name pattern of the modification pattern which should be added to the
         * sugar's long name
         * @param cond_name_pattern The short name pattern of the modification pattern which should be added to the
         * sugar's short name
         * @param tail The stream that will be added to the end of monosaccharide name and will be updated by this
         * function
         * @param head The stream that will be added to the beginning of monosaccharide name and will be updated by this
         * function
         * @param minus_one A boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
         * @param in_bracket The boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
         */
        void AddModificationRuleOneInfo(std::string key, std::string pattern, Glycan::Monosaccharide* mono,
                                        std::string long_name_pattern, std::string cond_name_pattern,
                                        std::stringstream& head, std::stringstream& tail, bool minus_one,
                                        std::stringstream& in_bracket);
        /*! \fn
         * A function in order to add the modification info to the name of the sugar based on the first group of rules
         * (No Bracket -> 2(r:6&!-1), Warning position -> a, Error Position -> 5(r6),4(r5) )
         * @param key The index of the atom which is a key in the derivative map
         * @param pattern The pattern of the modification which is based on the value in the derivative map
         * @param mono The monosaccharide object
         * @param long_name_pattern_at_minus_one The long name modification pattern for non-ring carbon attached to
         * anomeric carbon which should be added to the sugar's long name
         * @param long_name_pattern_at_minus_one The long name modification pattern for non-ring carbon attached to last
         * ring carbon which should be added to the sugar's long name
         * @param cond_name_pattern The short name pattern of the modification pattern which should be added to the
         * sugar's short name
         * @param tail The stream that will be added to the end of monosaccharide name and will be updated by this
         * function
         * @param minus_one A boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
         * @param in_bracket The boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
         */
        void AddModificationRuleTwoInfo(std::string key, std::string pattern, Glycan::Monosaccharide* mono,
                                        std::string long_name_pattern_at_minus_one,
                                        std::string long_name_pattern_at_plus_one, std::string cond_name_pattern,
                                        std::stringstream& tail, bool minus_one, std::stringstream& in_bracket);

        /*! \fn
         * A function in order to add unknown derivative info to the name of the sugar
         * @param key The index of the atom which is a key in the derivative map
         * @param pattern The pattern of the modification which is based on the value in the derivative map
         * @param mono The monosaccharide object
         * @param long_name_pattern The long name pattern of the modification pattern which should be added to the
         * sugar's long name
         * @param cond_name_pattern The short name pattern of the modification pattern which should be added to the
         * sugar's short name
         * @param head The stream that will be added to the beginning of monosaccharide name and will be updated by this
         * function
         * @param minus_one A boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
         * @param in_bracket The boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
         */
        void AddUnknownDerivativeRuleInfo(std::string key, std::string pattern, Glycan::Monosaccharide* mono,
                                          std::string long_name_pattern, std::string cond_name_pattern,
                                          std::stringstream& head, bool minus_one, std::stringstream& in_bracket);
        /*! \fn
         * A function in order to add the derivative info to the name of the sugar based on the first group of rules
         * (Error Position -> 5(r6),4(r5) )
         * @param key The index of the atom which is a key in the derivative map
         * @param pattern The pattern of the modification which is based on the value in the derivative map
         * @param mono The monosaccharide object
         * @param long_name_pattern The long name pattern of the modification pattern which should be added to the
         * sugar's long name
         * @param cond_name_pattern The short name pattern of the modification pattern which should be added to the
         * sugar's short name
         * @param head The stream that will be added to the beginning of monosaccharide name and will be updated by this
         * function
         * @param minus_one A boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
         * @param in_bracket The boolean value which indicates whether anomeric carbon is attached to a non-ring carbon
         */
        void AddDerivativeRuleInfo(std::string key, std::string pattern, Glycan::Monosaccharide* mono,
                                   std::string long_name_pattern, std::string cond_name_pattern,
                                   std::stringstream& head, bool minus_one, std::stringstream& in_bracket);
        /*! \fn
         * A function in order to update the chemical code structure of a complex monosaccharide (monosaccharide with
         * side atoms at position +2 and +3)
         * @param mono The monosaccharide object
         */
        // Moved to monosaccharide class
        //  void UpdateComplexSugarChemicalCode(Glycan::Monosaccharide* mono);
        /*! /fn
         * A function used to update the PDB code of a complex monosaccharide after derivatives have been found.
         * @param mono The monosaccharide object
         */
        // Moved to monosaccharide class
        //  void UpdatePdbCode( Glycan::Monosaccharide * mono );
        /*! \fn
         * A function in order to extract oligosacchride structure based on the linkages between monosaccharides
         * @param monos The list of extracted monosaccharide object
         * @param dataset_residue_names List of known terminal residue names
         * @param number_of_covalent_links Number of sugars that covalenly bond to a protein
         * @param number_of_probable_non_covalent_complex Number of sugars that might be non-covalantly interacting with
         * a protein
         * @return oligosacchrides The list of extracted oligosacchrides
         */
        void createOligosaccharideGraphs(std::vector<Glycan::Monosaccharide*> detected_monos,
                                         gmml::ResidueNameMap dataset_residue_names, int& number_of_covalent_links,
                                         int& number_of_probable_non_covalent_complexes);
        std::vector<Glycan::Oligosaccharide*>
        createOligosaccharides(std::vector<Glycan::Monosaccharide*> detected_monos);

        OligosaccharideVector ExtractOligosaccharides(std::vector<Glycan::Monosaccharide*> monos,
                                                      gmml::ResidueNameMap dataset_residue_names,
                                                      int& number_of_covalent_links,
                                                      int& number_of_probable_non_covalent_complexes);
        /** @}*/
        /** \addtogroup Verifiers_and_Issue_Resolvers
         * @{
         */
        /*! \fn
         * A function in order to check if the target atom is attached to OME terminal
         * @param target_atom The atom which will be checked for a terminal
         * @return pattern The discovered pattern attached to the target atom
         */
        std::string CheckOMETerminal(Atom* target, AtomVector& terminal_atoms);
        /*! \fn
         * A function in order to check if the target atom is attached to ROH terminal
         * @param target_atom The atom which will be checked for a terminal
         * @return pattern The discovered pattern attached to the target atom
         */
        std::string CheckROHTerminal(Atom* target, AtomVector& terminal_atoms);
        /*! \fn
         * A function in order to check if the target atom is attached to TBT terminal
         * @param target_atom The atom which will be checked for a terminal
         * @return pattern The discovered pattern attached to the target atom
         */
        std::string CheckTBTTerminal(Atom* target, AtomVector& terminal_atoms);
        /*! \fn
         * A function in order to check if the target atom is attached to TBT, ROH or OME patterns when the residue name
         * is different from the terminal residue names
         * @param target_atom The atom which will be checked for a terminal
         * @return pattern The discovered pattern attached to the target atom
         */
        std::string CheckTerminals(Atom* target, AtomVector& terminal_atoms);
        void CheckLinkageNote(Glycan::Monosaccharide* mono1, Glycan::Monosaccharide* mono2, std::string linkage,
                              std::vector<std::string>& checked_linkages);

        /*! \fn
         * A function in order to calculate the orientation(S/R) of additional side atoms
         * @param prev_atom The previous neighbor atom of the current atom
         * @param target_atom The atom which the orientation of its oxygen neighbor is going to be calculated
         * @param next_atom The next neighbor atom of the current atom
         * @return orientation The relative orientation of the additional side atom's oxygen neighbor
         */
        std::string CalculateRSOrientations(Atom* prev_atom, Atom* target, Atom* next_atom);
        /** @}*/
        /** \addtogroup Molecular_Data_Structure_Builders
         * @{
         */
        /*! \fn
                  * A recursive function in order to build the tree like structure of the oligosacchrides (a directed
           graph made from the molecular graph structure of monosaccharides. each monosaccharide is a node and their
           linkages are the edges)
                  * @param key The root monosaccharide(node) of each oligosacchride
                  * @param value The monosaccharide attached to the key monosacharide
                  * @param oligo The oligosaccharide object to be constructed by the function
                  * @param visited_monos The monosaccharide objects that have been visited by previous calls to the
           function
                  * @param monos_table A mapping between a monosaccharide object and the list of monosaccharides
           attached to it
                  * @param monos_table_linkages A mapping between a monosaccharide object and the list of atoms that are
           involved in the linkages with the monosaccharides attached to the current monosaccharide object
                  * @param visited_linkages The list of linkages(atoms involved in the linkages between monosaccharides)
           that have been visited so far by calls to the function
                  */
        void BuildOligosaccharideTreeStructure(
            Glycan::Monosaccharide* key, std::vector<Glycan::Monosaccharide*> val, Glycan::Oligosaccharide* oligo,
            std::vector<int>& visited_monos,
            std::map<Glycan::Monosaccharide*, std::vector<Glycan::Monosaccharide*>> monos_table,
            std::map<Glycan::Monosaccharide*, std::vector<std::string>> monos_table_linkages,
            std::vector<std::string>& visited_linkages);
        /** @}*/
        /** \addtogroup Verifiers_and_Issue_Resolvers
         * @{
         */
        /*! \fn
         * A function in order to check if the target atom is attached to a derivative with the pattern xCH-N
         * @param target_atom The atom which will be checked for a derivative
         * @param cycle_atom_str The string version of atom identifiers of the cycle
         * @param pattern_atoms The list of atoms that involved in the pattern
         * @return pattern The discovered pattern of the attached derivative
         */
        std::string CheckxC_N(Atom* target, std::string cycle_atoms_str /*, AtomVector& pattern_atoms*/);
        /*! \fn
         * A function in order to check if the target atom is attached to a derivative with the pattern xC-O-C=OCH3 or
         * xC-N-C=OCH3
         * @param target_atom The atom which will be checked for a derivative
         * @param cycle_atom_str The string version of atom identifiers of the cycle
         * @param NxO The option for checking the pattern with oxygen or nitrogen
         * @param pattern_atoms The list of atoms that involved in the pattern
         * @return pattern The discovered pattern of the attached derivative
         */
        std::string CheckxC_NxO_CO_C(Atom* target, std::string cycle_atoms_str, char NxO, AtomVector& pattern_atoms);
        /*! \fn
         * A function in order to check if the target atom is attached to a derivative with the pattern xC-N-C=OCH2OH or
         * xC-O-C=OCH2OH
         * @param target_atom The atom which will be checked for a derivative
         * @param cycle_atom_str The string version of atom identifiers of the cycle
         * @param NxO The option for checking the pattern with oxygen or nitrogen
         * @param pattern_atoms The list of atoms that involved in the pattern
         * @return pattern The discovered pattern of the attached derivative
         */
        std::string CheckxC_NxO_CO_CO(Atom* target, std::string cycle_atoms_str, char NxO, AtomVector& pattern_atoms);
        /*! \fn
         * A function in order to check if the target atom is attached to a derivative with the pattern xC-N-SO3 or
         * xC-O-SO3
         * @param target_atom The atom which will be checked for a derivative
         * @param cycle_atom_str The string version of atom identifiers of the cycle
         * @param NxO The option for checking the pattern with oxygen or nitrogen
         * @param pattern_atoms The list of atoms that involved in the pattern
         * @return pattern The discovered pattern of the attached derivative
         */
        std::string CheckxC_NxO_SO3(Atom* target, std::string cycle_atoms_str, char NxO, AtomVector& pattern_atoms);
        /*! \fn
         * A function in order to check if the target atom is attached to a derivative with the pattern xC-N-PO3 or
         * xC-PO3
         * @param target_atom The atom which will be checked for a derivative
         * @param cycle_atom_str The string version of atom identifiers of the cycle
         * @param NxO The option for checking the pattern with oxygen or nitrogen
         * @param pattern_atoms The list of atoms that involved in the pattern
         * @return pattern The discovered pattern of the attached derivative
         */
        std::string CheckxC_NxO_PO3(Atom* target, std::string cycle_atoms_str, char NxO, AtomVector& pattern_atoms);
        /*! \fn
         * A function in order to check if the target atom is attached to a derivative with the pattern xC-N-CH3 or
         * xC-O-CH3
         * @param target_atom The atom which will be checked for a derivative
         * @param cycle_atom_str The string version of atom identifiers of the cycle
         * @param NxO The option for checking the pattern with oxygen or nitrogen
         * @param pattern_atoms The list of atoms that involved in the pattern
         * @return pattern The discovered pattern of the attached derivative
         */
        std::string CheckxC_NxO_C(Atom* target, std::string cycle_atoms_str, char NxO, AtomVector& pattern_atoms);
        /*! \fn
         * A function in order to check if the target atom is attached to a derivative with the pattern C-(O,O),
         * C-(O.OH), rC-(O.O) and rC-(O.OH)
         * @param target_atom The atom which will be checked for a derivative
         * @param cycle_atom_str The string version of atom identifiers of the cycle
         * @param NxO The option for checking the pattern with oxygen or nitrogen
         * @param pattern_atoms The list of atoms that involved in the pattern
         * @return pattern The discovered pattern of the attached derivative
         */
        std::string CheckxCOO(Atom* target, std::string cycle_atoms_str /*, AtomVector& pattern_atoms*/);
        void AddIon(std::string ion_name, std::string lib_file, std::string parameter_file, int ion_count = 0);
        void AddSolvent(double extension, double closeness, Assembly* solvent_component_assembly, std::string lib_file);
        void SplitSolvent(Assembly* solvent, Assembly* solute);
        void SplitIons(Assembly* assembly, ResidueVector ions);

        double GetTotalCharge();
        double GetRadius();
        double GetTotalMass();
        void GetCenterOfMass(GeometryTopology::Coordinate* center_of_mass);
        void GetCenterOfGeometry(GeometryTopology::Coordinate* center_of_geometry);
        void GetBoundary(GeometryTopology::Coordinate* lower_left_back_corner,
                         GeometryTopology::Coordinate* upper_right_front_corner);
        /*! \fn
         * A function in order to calculate the surface area of overlap between atoms of two assemblies
         * @param assemblyB is the second assembly. Overlaps between atoms of the same assembly are not counted
         * @return Total overlap between assemblies, relative to the surface area of a buried C atom.
         */
        double CalculateAtomicOverlaps(Assembly* assemblyB);

        /*! \fn                                                                             //Added by ayush on 11/16/17
         * for residuenodes in assembly A function to add a residuenode to the ResidueNodeVector in an assembly
         * @param residue a pointer to the residue object of class Residue
         */
        void AddResidueNode(MolecularModeling::ResidueNode* residuenode);

        /*! \fn                                                                             //Added by ayush on 11/17/17
         * for residuenodes in assembly A function to generate residuenodes in an assembly
         * @return A residuenode vector of all the residues nodes present in the assembly
         */
        MolecularModeling::ResidueNodeVector GenerateResidueNodesInAssembly();

        /*! \fn                                                                             //Added by ayush on 11/12/17
         * for molecules in assembly A function to add a molecule to the MoleculeVector in an assembly
         * @param molecule a pointer to the molecule object of class Molecule
         */
        void AddMolecule(MolecularModeling::Molecule* molecule);

        /*! \fn                                                                             //Added by ayush on 12/7/17
         * for molecules in assembly A function to generate molecules in an assembly
         * @return A molecule vector of all the molecules present in the assembly
         */
        void GenerateMoleculesInAssembly();

        /*! \fn                                                                             //Added by ayush on 12/7/17
         * for molecules in assembly A DFS function to traverse the residuenodes the assembly
         */
        void GenerateMoleculesDFSUtil(MolecularModeling::ResidueNode* residuenode);

        double CalculateAtomicOverlaps(AtomVector assemblyBAtoms);
        AtomVector GetAllAtomsOfAssemblyWithinXAngstromOf(GeometryTopology::Coordinate* coordinate, double distance);

        /*! \fn                                                                              //Added by ayush on
         * 04/11/18 for TopologyFix in assembly A function that returns list of atoms bonded to each other by start and
         * direction in the Assembly. eg:(start)Atom1->(direction)Atom2
         * @param start_atom The starting point in the assembly list of Atoms.
         * @param direction_atom The direction for traversing the bonding among atom list of Assembly
         * @param ignore_list The list of atoms which is ignored during traversal
         * @return bonded_atoms_bystartdirection_ A list of assembly atoms bonded to each other based on start point and
         * direction
         */
        AtomVector GetAllBondedAtomsByStartDirection(Atom* start_atom, Atom* direction_atom, AtomVector ignore_list);

        /*! \fn                                                                          //Added by ayush on 04/11/18
         * for TopologyFix in assembly A function to check if an atom exists in Assebly AtomList
         * @param An atom to check
         * @start_atom_neighbors Neighbors of the start atom in considertion
         * @ignore_list Ignore list of atoms provided by the user
         * @return True/False based on existence
         */
        bool CheckIfAtomExistInAssembly(Atom* atom);

        /*! \fn                                                                              //Added by ayush on
         * 04/11/18 for TopologyFix in assembly A function that performs the Depth First Search traversal to find the
         * bonded atoms based on start and direction atoms.
         * @param atom The current atom under consideration.
         * @param start_atom_neighbors The list of neighbors of the current atom.
         * @param ignore_list The list of atoms which is ignored during traversal
         */
        void BondedAtomsByStartDirectionDFSUtil(Atom* atom, AtomVector start_atom_neighbors, AtomVector ignore_list);

        /*! \fn                                                                          //Added by ayush on 04/16/18
         * for TopologyFix in assembly A function that returns the cooridnate vector corresponding to the atom vector
         * based on the index specified.
         * @param atomList Vector of atoms (Atoms)
         * @param CoordinateIndex The index of the coordinate set that should be extracted
         * @return Vector of pointers to the coordinates for the vector of atoms, in the same order as the vector of
         * atoms.
         */
        GeometryTopology::CoordinateVector GetCoordinatesFromAtomVector(AtomVector atomList, int CoordinateIndex);

        /*! \fn                                                                          //Added by ayush on 06/16/18
         * for OffFile A function that create an OFF file from Assembly
         * @param CoordinateIndex The index of the coordinate set
         */
        void CreateOffFileFromAssembly(std::string file_name, int CoordinateIndex);
        void SerializeResidueNumbers();
        void SetChargesAndAtomTypes();
        void EnsureIntegralCharge();
        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        /*! \fn
         * A function to print out the assembly contents in a structural format
         * Print out the information in a defined structure
         * @param out An output stream, the print result will be written in the given output stream
         */
        void Print(std::ostream& out = std::cerr);

        void PrettyPrintHet(std::ostream& out = std::cerr);
        void PrintHetResidues(std::ostream& out = std::cerr);
        void PrintHetAtoms(std::ostream& out = std::cerr);

        void WriteHetResidues(std::string file_name);
        void WriteHetAtoms(std::string file_name);

      private:
        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        std::string name_;          /*!< Assembly name >*/
        AssemblyVector assemblies_; /*!< List of assemblies that may involved in a big structure as a new assembly >*/
        ResidueVector residues_;    /*!< List of residues involved in the current object of assembly >*/
        std::string chemical_type_; /*!< A descriptor for the chemical type of the current object of assembly >*/
        int sequence_number_;       /*!< An integer number to indicates sequence of importing the current assembly >*/
        std::string id_;
        std::string description_;              /*!< Short description for the current assembly >*/
        std::string source_file_;              /*!< File name that the current assembly has been built upon >*/
        gmml::InputFileType source_file_type_; /*!< Type of the file that the current assembly has been built upon >*/
        int model_index_;  /*!< In case that there are more than one models for an assembly, this attribute indicated
                              which model is the target model >*/
        NoteVector notes_; /*!< A list of note instances from the Note struct in Glycan name space which is used for
                              representing the potential issues within a structure >*/
        MolecularModeling::ResidueNodeVector residuenodes_;
        /*!< List of residuenodes present in the current object of assembly >*/ // Added by ayush on 11/16/17 for
                                                                                // residuenodes in assembly
        MolecularModeling::MoleculeVector molecules_;
        /*!< List of molecules present in the current object of assembly >*/ // Added by ayush on 11/12/17 for
                                                                             // molecules in assembly
        AtomVector bonded_atoms_bystartdirection_;
        /*!< List of atoms bonded based on start point an ddirection in Assembly>*/ // Added by ayush on 04/11/18
                                                                                    // for Bonded Atoms based on
                                                                                    // start-direction in assembly
        PdbFileSpace::InputFile* input_file_; /*!< A pointer back to the Input object >*/
    };

    struct DistanceCalculationThreadArgument
    {
        int thread_index;
        int number_of_threads;
        int model_index;
        double cutoff;
        Assembly* a;

        DistanceCalculationThreadArgument()
        {
            thread_index      = 0;
            number_of_threads = 1;
            model_index       = 0;
            cutoff            = gmml::maxCutOff;
            a                 = NULL;
        }

        DistanceCalculationThreadArgument(int ti, int tn, int mi, double c, Assembly* assembly)
        {
            thread_index      = ti;
            number_of_threads = tn;
            model_index       = mi;
            cutoff            = c;
            a                 = assembly;
        }
    };

    struct DistanceCalculationByMatrixThreadArgument
    {
        int thread_index;
        int model_index;
        double cutoff;
        std::vector<Atom*>* first_chunk;
        std::vector<Atom*>* second_chunk;

        DistanceCalculationByMatrixThreadArgument()
        {
            thread_index = 0;
            model_index  = 0;
            cutoff       = gmml::maxCutOff;
            first_chunk  = new std::vector<Atom*>();
            second_chunk = new std::vector<Atom*>();
        }

        DistanceCalculationByMatrixThreadArgument(int ti, int mi, double c, std::vector<Atom*>* fc,
                                                  std::vector<Atom*>* sc)
        {
            thread_index = ti;
            model_index  = mi;
            cutoff       = c;
            first_chunk  = fc;
            second_chunk = sc;
        }
    };

    struct BacktrackingElements
    {
      public:
        BacktrackingElements(gmml::GlycamAtomNameMap pdb_glycam_map, gmml::GlycamAtomNameMap glycam_atom_map,
                             AtomVector atoms, int index, std::queue<Atom*> to_visit = std::queue<Atom*>())
        {
            pdb_glycam_map_  = pdb_glycam_map;
            glycam_atom_map_ = glycam_atom_map;
            atoms_           = atoms;
            index_           = index;
            to_visit_        = to_visit;
        }

        gmml::GlycamAtomNameMap pdb_glycam_map_;
        gmml::GlycamAtomNameMap glycam_atom_map_;
        AtomVector atoms_;
        std::queue<Atom*> to_visit_;
        int index_;
    };
} // namespace MolecularModeling

#endif // ASSEMBLY_HPP
