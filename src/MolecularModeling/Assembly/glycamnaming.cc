#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <limits>

#include "../../../includes/MolecularModeling/assembly.hpp"
#include "../../../includes/MolecularModeling/residue.hpp"
#include "../../../includes/MolecularModeling/atom.hpp"
#include "../../../includes/MolecularModeling/atomnode.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyfile.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologybond.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../../includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbtitlesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbconnectsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinksection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/GeometryTopology/grid.hpp"
#include "../../../includes/GeometryTopology/cell.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <unistd.h>
#include <errno.h>
#include <string.h>

using MolecularModeling::Assembly;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void Assembly::ObtainAllOligasaccharideResidues(Glycan::Oligosaccharide* oligo, ResidueVector& oligo_residues)
{
    MolecularModeling::Residue* sugar_residue = oligo->root_->cycle_atoms_[0]->GetResidue();
    oligo_residues.push_back(sugar_residue);

    for (unsigned int i = 0; i < oligo->child_oligos_.size(); i++)
    {
        Glycan::Oligosaccharide* child_oligo = oligo->child_oligos_[i];
        ObtainAllOligasaccharideResidues(child_oligo, oligo_residues);
    }

    return;
}

void Assembly::PutAglyconeInNewResidueAndRearrangeGlycanResidues(
    std::vector<Glycan::Oligosaccharide*> oligosaccharides,
    std::map<Glycan::Oligosaccharide*, ResidueVector>& oligo_residue_map)
{
    for (std::vector<Glycan::Oligosaccharide*>::iterator it = oligosaccharides.begin(); it != oligosaccharides.end();
         it++)
    {
        Glycan::Oligosaccharide* oligo      = *it;
        MolecularModeling::Atom* anomeric_o = NULL;
        MolecularModeling::Atom* anomeric_c = NULL;

        if (oligo->root_->cycle_atoms_.at(0) != NULL)
        {
            anomeric_c = oligo->root_->cycle_atoms_.at(0);
        }

        // Anomeric group has now been removed from sidegroup if it shouldn't be part of oligosaccharide (e.g. NLN
        // nitrogen). So now, find anomeric oxygen from node neighbors of anomeric carbon.
        AtomVector anomeric_c_neighbors = anomeric_c->GetNode()->GetNodeNeighbors();
        for (unsigned int i = 0; i < anomeric_c_neighbors.size(); i++)
        {
            MolecularModeling::Atom* this_neighbor = anomeric_c_neighbors[i];
            std::string element                    = this_neighbor->GetElementSymbol();
            if (!this_neighbor->GetIsCycle() && (element == "O" || element == "S" || element == "N"))
            {
                anomeric_o = this_neighbor;
            }
        }

        if (anomeric_o == NULL || anomeric_c == NULL)
        {
            std::cerr << "Either anomeric carbon or oxygen is not detected. Check your sugar structure. Aborting\n";
            throw;
        }

        AtomVector terminal_atoms = AtomVector();
        std::string terminal = CheckTerminals(anomeric_o, terminal_atoms);

        if(terminal.compare("") == 0)
        {
            std::cerr << "Unknown aglycone detected.Aborting" << std::endl;
            throw;
        }

        std::string oligo_name = oligo->IUPAC_name_;
        if (oligo->terminal_.compare("") == 0)
        {
            std::cerr << "Error, no aglycone detected.Aborting\n";
            throw;
        }

        CondensedSequenceSpace::CondensedSequence* condensed_sequence =
            new CondensedSequenceSpace::CondensedSequence(oligo_name);
        // Gets the three letter code of all carbohydrates involved in current oligosaccharide
        CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree
            condensed_sequence_glycam06_residue_tree = condensed_sequence->GetCondensedSequenceGlycam06ResidueTree();

        ResidueVector oligo_residues;
        ObtainAllOligasaccharideResidues(oligo, oligo_residues);

        // Determine the location of the previous residue from this oligosaccharide.
        bool no_preceding_residue                     = false;
        MolecularModeling::Residue* preceding_residue = NULL;
        ResidueVector AllResiduesInAssembly           = this->GetResidues();

        for (unsigned int i = 0; i < AllResiduesInAssembly.size(); i++)
        {
            MolecularModeling::Residue* this_residue = AllResiduesInAssembly[i];

            // For the 1st residue that is part of this oligo, get the previous residue.
            if (std::find(oligo_residues.begin(), oligo_residues.end(), this_residue) != oligo_residues.end())
            {
                if (i == 0)
                {
                    no_preceding_residue = true;
                }
                else
                {
                    preceding_residue = AllResiduesInAssembly[i - 1];
                }
                break;
            }
        }

        // Now remove all oligo residues. They will be reordered.
        for (unsigned int i = 0; i < oligo_residues.size(); i++)
        {
            this->RemoveResidue(oligo_residues[i]);
        }

        std::string aglycone_name                         = condensed_sequence_glycam06_residue_tree.at(0)->GetName();
        MolecularModeling::Residue* residue_with_aglycone = oligo_residues[0]; //Strange that the 1st residue is not always the aglycone
        /*MolecularModeling::Residue* residue_with_aglycone = NULL;
		for (unsigned int i = 0; i < oligo_residues.size(); i++)
		{
			MolecularModeling::Residue* oligo_res = oligo_residues[i];
			if (oligo_res->GetName() == aglycone_name){
				residue_with_aglycone = oligo_res;
				break;
			}
		}

		if (residue_with_aglycone == NULL){
			std::cout << "Error, aglycone residue with name " << aglycone_name << " cannot be detected in this oligsaccharide." << std::endl;
			std::exit(1);
		}*/

        MolecularModeling::Residue* new_terminal_residue  = NULL;

        MolecularModeling::Residue* new_terminal_residue = NULL;

        if (!(anomeric_o->GetResidue()
                  ->CheckIfProtein())) // If this terminal atom is not part of protein,for example,NLN, then it should
                                       // be in a new glycam residue, for example, ROH.
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "First sugar residue is: " + residue_with_aglycone->GetId());
            new_terminal_residue = new Residue(this,aglycone_name);
			new_terminal_residue->SetName(aglycone_name);

            oligo_residue_map[oligo].push_back(new_terminal_residue);

            // If there is no preceding residue, simply add the new terminal residue, which will be the 1st residue
            if (no_preceding_residue)
            {
                this->AddResidue(new_terminal_residue);
            }
            else
            {
                ResidueVector::iterator preceding_residue_it =
                std::find(this->residues_.begin(), this->residues_.end(), preceding_residue);
                ResidueVector::iterator insertion_location = preceding_residue_it + 1;
                this->residues_.insert(insertion_location, new_terminal_residue);
            }

            // For example, the ROH coming from terminal BGC, should go in front of BGC.
            // OG std::string new_terminal_residue_id = terminal_atom->GetId(); //This new residue takes the atom id as
            // residue id. // OG: huh? OG new_terminal_residue->SetId(new_terminal_residue_id); OG Ok let's see if doing
            // something more "sensible" breaks something else:
            int newResidueNumberInt = std::stoi(new_terminal_residue->GetNumber()) -
                                      1; // OG: You sometimes need to be -1, but I don't think that functionality
            // currently exists, and this thing is so bad I ain't putting any more time into it.
            new_terminal_residue->SetId(new_terminal_residue->CreateID(
                aglycone_name, residue_with_aglycone->GetChainID(), std::to_string(newResidueNumberInt)));

            gmml::log(__LINE__, __FILE__, gmml::INF, "new_terminal_residue_id is: " + new_terminal_residue->GetId());

            // Get all the atoms that should go to this new terminal residue. Find all connected atoms of the terminal
            // atom that's not on the sugar side of the old residue

            for (AtomVector::iterator atom_it = terminal_atoms.begin(); atom_it != terminal_atoms.end(); atom_it++)
            {
                MolecularModeling::Atom* new_atom = *atom_it;
                // std::cout << "New terminal residue " << new_terminal_residue->GetName() << " adding atom " <<
                // new_atom->GetName() << std::endl;
                new_terminal_residue->AddAtom(new_atom);
                new_atom->SetResidue(new_terminal_residue);
                residue_with_aglycone->RemoveAtom(new_atom, false);

                gmml::log(__LINE__, __FILE__, gmml::INF,
                          "New terminal residue id is: " + new_terminal_residue->GetId());
                gmml::log(__LINE__, __FILE__, gmml::INF, "Adding atom with id: " + new_atom->GetId());
            }
        }

        if (new_terminal_residue != NULL)
        {
            ResidueVector::iterator new_terminal_residue_location =
                std::find(this->residues_.begin(), this->residues_.end(), new_terminal_residue);
            ResidueVector::iterator insertion_location = new_terminal_residue_location + 1;
            this->residues_.insert(insertion_location, oligo_residues.begin(), oligo_residues.end());
        }
        else
        {

            ResidueVector::iterator insertion_location = this->residues_.begin();
            if (no_preceding_residue)
            {
                insertion_location = this->residues_.begin();
            }
            else
            {
                ResidueVector::iterator preceding_residue_location =
                    std::find(this->residues_.begin(), this->residues_.end(), preceding_residue);
                insertion_location = preceding_residue_location + 1;
            }
            this->residues_.insert(insertion_location, oligo_residues.begin(), oligo_residues.end());
        }
    }
    return;
}

void Assembly::Reorder06TreeTempDeleteLater(CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree& unordered_tree, CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree& ordered_tree, int current_tree_index, std::vector<int>& visited){
	CondensedSequenceSpace::CondensedSequenceGlycam06Residue* g06_residue_current = unordered_tree[current_tree_index];	
	ordered_tree.push_back(g06_residue_current);
	visited.push_back(current_tree_index);
	std::map<int, int> sorted_children;

    for (unsigned int i = 0; i < unordered_tree.size(); i++){
    	CondensedSequenceSpace::CondensedSequenceGlycam06Residue* g06_residue = unordered_tree[i];
        int parent_index = g06_residue->GetParentId();

        if (parent_index == current_tree_index){
        	std::string parent_oxygen_name = g06_residue->GetParentOxygen();
            int po_index = (parent_oxygen_name.length() <=1) ? 0 : std::stoi(parent_oxygen_name.substr(1));
            sorted_children[po_index] = i;
            //This assumes that atom name consists of only 1 character as element and the rest are numbers. If not this will break.
        }
    }

	for (std::map<int, int>::iterator mapit = sorted_children.begin(); mapit != sorted_children.end(); mapit++){
		int next_child_index = mapit->second;
		if (std::find(visited.begin(), visited.end(), next_child_index) != visited.end()) continue;
		this->Reorder06TreeTempDeleteLater(unordered_tree, ordered_tree, next_child_index, visited);
	}
}

gmml::GlycamResidueNamingMap
Assembly::ExtractResidueGlycamNamingMap(std::vector<Glycan::Oligosaccharide*> oligosaccharides,
                                        std::map<Glycan::Oligosaccharide*, std::vector<std::string>>& oligo_id_map,
                                        std::map<Glycan::Oligosaccharide*, ResidueVector>& oligo_residue_map)
{
    // TODO: Done
    // Update the map
    // Key: std::string -> residue_id or atom_id
    // Value: std::vector<std::string> -> list of possible glycam naming
    gmml::GlycamResidueNamingMap pdb_glycam_residue_map = gmml::GlycamResidueNamingMap();
    // Iterates on all oligosaccharides and update the naming map
    for (std::vector<Glycan::Oligosaccharide*>::iterator it = oligosaccharides.begin(); it != oligosaccharides.end();
         it++)
    {
        int index                      = 0;
        Glycan::Oligosaccharide* oligo = *it;
        std::string oligo_name         = oligo->IUPAC_name_;
        // std::cout << "Iupac name: " << oligo_name << std::endl;
        // In case that there is no terminal attached to the reducing end adds a temporary terminal residue to make the
        // sequence parser able to parse the sequence
        if (oligo->terminal_.compare("") == 0)
        {
            std::cerr << "Error, no aglycone detected.Aborting\n";
            throw;
        }

        CondensedSequenceSpace::CondensedSequence* condensed_sequence =
            new CondensedSequenceSpace::CondensedSequence(oligo_name);
        // Gets the three letter code of all carbohydrates involved in current oligosaccharide
        CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree
            condensed_sequence_glycam06_residue_tree_unordered = condensed_sequence->GetCondensedSequenceGlycam06ResidueTree();

        // If for some reason there are no residues in the condensed sequence, which really shouldn't happen. Terminate
        // the program as well.
        if (condensed_sequence_glycam06_residue_tree_unordered.empty())
        {
            std::cerr << "Error, no glycam 06 residues detected" << std::endl;
            throw;
        }

		//Yao temporarily patch code 20240827. Reorder 06 tree based on bond index;
		CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree condensed_sequence_glycam06_residue_tree;
		//Add aglycone
		condensed_sequence_glycam06_residue_tree.push_back(condensed_sequence_glycam06_residue_tree_unordered[0]);
		std::vector<int> visited(1, 0);
		this->Reorder06TreeTempDeleteLater(condensed_sequence_glycam06_residue_tree_unordered, condensed_sequence_glycam06_residue_tree, 1, visited);

        // Detect anomeric carbon and oxygen to determine aglycone atoms, which in regular PDB file is put in the same
        // residue of the 1st sugar. For Glycam, this goes into a new residue, such as ROH.
        MolecularModeling::Atom* anomeric_o = NULL;
        MolecularModeling::Atom* anomeric_c = NULL;
        if (oligo->root_->cycle_atoms_.at(0) != NULL)
        {
            anomeric_c = oligo->root_->cycle_atoms_.at(0);
        }

        // Anomeric group has now been removed from sidegroup if it shouldn't be part of oligosaccharide (e.g. NLN
        // nitrogen). So now, find anomeric oxygen from node neighbors of anomeric carbon.
        AtomVector anomeric_c_neighbors = anomeric_c->GetNode()->GetNodeNeighbors();
        for (unsigned int i = 0; i < anomeric_c_neighbors.size(); i++)
        {
            MolecularModeling::Atom* this_neighbor = anomeric_c_neighbors[i];
            std::string element                    = this_neighbor->GetElementSymbol();
            if (!this_neighbor->GetIsCycle() && (element == "O" || element == "S" || element == "N"))
            {
                anomeric_o = this_neighbor;
            }
        }

        if (anomeric_o == NULL || anomeric_c == NULL)
        {
            std::cerr << "Either anomeric carbon or oxygen is not detected. Check your sugar structure. Aborting\n";
            throw;
        }

        AtomVector terminal_atoms = AtomVector();
        std::string terminal = CheckTerminals(anomeric_o, terminal_atoms);
		
		std::cout << "In extract map detected terminal name: " << terminal << std::endl;
		
        if(terminal.compare("") == 0)
        {
            std::cerr << "Unknown aglycone detected.Aborting" << std::endl;
            throw;
        }

        // TODO:
        // Add the residue mismatch into a structure for the Ontology usage

        // Separating terminal residue in glycam naming
        std::string terminal_atom_id    = anomeric_o->GetId();
        std::string terminal_residue_id = anomeric_o->GetResidue()->GetId();
        gmml::log(__LINE__, __FILE__, gmml::INF, "terminal_atom_id is: " + terminal_atom_id);
        gmml::log(__LINE__, __FILE__, gmml::INF, "terminal_residue_id is: " + terminal_residue_id);
        if (pdb_glycam_residue_map.find(terminal_atom_id) == pdb_glycam_residue_map.end())
        {
            pdb_glycam_residue_map[terminal_atom_id] == std::vector<std::string>();
        }

        std::string aglycone_name = condensed_sequence_glycam06_residue_tree.at(0)->GetName();
        pdb_glycam_residue_map[terminal_atom_id].push_back(aglycone_name);


        oligo_id_map[oligo].push_back(terminal_atom_id);

        if (pdb_glycam_residue_map.find(terminal_residue_id) == pdb_glycam_residue_map.end())
        {
            pdb_glycam_residue_map[terminal_residue_id] = std::vector<std::string>();
        }

        pdb_glycam_residue_map[terminal_residue_id].push_back(aglycone_name);
        oligo_id_map[oligo].push_back(terminal_residue_id);

        gmml::log(__LINE__, __FILE__, gmml::INF, "aglycone_name is: " + aglycone_name);

        index++;
        this->ExtractOligosaccharideNamingMap(pdb_glycam_residue_map, oligo, condensed_sequence_glycam06_residue_tree,
                                              index, oligo_id_map, oligo, oligo_residue_map);
    }
    return pdb_glycam_residue_map;
}

void Assembly::ExtractOligosaccharideNamingMap(
    gmml::GlycamResidueNamingMap& pdb_glycam_map, Glycan::Oligosaccharide* oligosaccharide,
    CondensedSequenceSpace::CondensedSequence::CondensedSequenceGlycam06ResidueTree
        condensed_sequence_glycam06_residue_tree,
    int& index, std::map<Glycan::Oligosaccharide*, std::vector<std::string>>& oligo_id_map,
    Glycan::Oligosaccharide* root_oligo, std::map<Glycan::Oligosaccharide*, ResidueVector>& oligo_residue_map)
{
    std::string name       = condensed_sequence_glycam06_residue_tree.at(index)->GetName();
    // TODO: Done
    // Update to hold all possible three letter names for the specific residue_id
    std::string residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
    if (pdb_glycam_map.find(residue_id) == pdb_glycam_map.end())
    {
        pdb_glycam_map[residue_id] = std::vector<std::string>();
    }
    pdb_glycam_map[residue_id].push_back(name);

    oligo_id_map[oligosaccharide].push_back(residue_id);
    ResidueVector residue_vector = oligo_residue_map[root_oligo];
    if (std::find(residue_vector.begin(), residue_vector.end(),
                  oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()) == residue_vector.end())
    {
        oligo_residue_map[root_oligo].push_back(oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue());
    }
    index++;
    // Separating SO3 and PO3 residues in glycam naming
    while (index < (int)condensed_sequence_glycam06_residue_tree.size() &&
           condensed_sequence_glycam06_residue_tree.at(index)->GetIsDerivative())
    {
        int parent_index = condensed_sequence_glycam06_residue_tree.at(index)->GetParentId();
        int carbon_index =
            gmml::ConvertString<int>(condensed_sequence_glycam06_residue_tree.at(index)->GetAnomericCarbon().substr(1));
        MolecularModeling::Atom* carbon_atom = NULL;
        if (std::find_if(oligosaccharide->root_->derivatives_map_.begin(),
                         oligosaccharide->root_->derivatives_map_.end(),
                         [](const std::pair<std::string, std::string>& element)
                         {
                             return element.first == "-1";
                         }) == oligosaccharide->root_->derivatives_map_.end())
        // if(oligosaccharide->root_->derivatives_map_.find("-1") == oligosaccharide->root_->derivatives_map_.end())
        {
            if ((carbon_index == 6 || carbon_index == 7 || carbon_index == 8) && carbon_atom == NULL)
            {
                if (oligosaccharide->root_->sugar_name_.ring_type_.compare("P") == 0)
                {
                    carbon_atom = oligosaccharide->root_->side_atoms_.at(oligosaccharide->root_->side_atoms_.size() - 1)
                                      .at(carbon_index - 6);
                }
            }
            if ((carbon_index == 5 || carbon_index == 6 || carbon_index == 7) && carbon_atom == NULL)
            {
                if (oligosaccharide->root_->sugar_name_.ring_type_.compare("F") == 0)
                {
                    carbon_atom = oligosaccharide->root_->side_atoms_.at(oligosaccharide->root_->side_atoms_.size() - 1)
                                      .at(carbon_index - 5);
                }
            }
            if (carbon_atom == NULL)
            {
                carbon_atom = oligosaccharide->root_->side_atoms_.at(carbon_index).at(0);
            }
        }
        else
        {
            if ((carbon_index == 7 || carbon_index == 8 || carbon_index == 9) && carbon_atom == NULL)
            {
                if (oligosaccharide->root_->sugar_name_.ring_type_.compare("P") == 0)
                {
                    carbon_atom = oligosaccharide->root_->side_atoms_.at(oligosaccharide->root_->side_atoms_.size() - 1)
                                      .at(carbon_index - 7);
                }
            }
            if ((carbon_index == 6 || carbon_index == 7 || carbon_index == 8) && carbon_atom == NULL)
            {
                if (oligosaccharide->root_->sugar_name_.ring_type_.compare("F") == 0)
                {
                    carbon_atom = oligosaccharide->root_->side_atoms_.at(oligosaccharide->root_->side_atoms_.size() - 1)
                                      .at(carbon_index - 6);
                }
            }
            if (carbon_atom == NULL)
            {
                carbon_atom = oligosaccharide->root_->side_atoms_.at(carbon_index - 1).at(0);
            }
        }
        if (carbon_atom != NULL)
        {
            AtomVector n_linkage_derivative_atoms = AtomVector();
            std::string n_linkage_derivative_string =
                CheckxC_NxO_SO3(carbon_atom, oligosaccharide->root_->cycle_atoms_str_, 'N', n_linkage_derivative_atoms);
            if (n_linkage_derivative_string.compare("xC-N-SO3") == 0)
            {
                for (AtomVector::iterator it = n_linkage_derivative_atoms.begin();
                     it != n_linkage_derivative_atoms.end(); it++)
                {
                    MolecularModeling::Atom* atom  = *it;
                    std::string derivative_atom_id = atom->GetId();
                    if (pdb_glycam_map.find(derivative_atom_id) == pdb_glycam_map.end())
                    {
                        pdb_glycam_map[derivative_atom_id] = std::vector<std::string>();
                    }
                    pdb_glycam_map[derivative_atom_id].push_back("SO3");
                    oligo_id_map[oligosaccharide].push_back(derivative_atom_id);
                }
                std::string new_name = condensed_sequence_glycam06_residue_tree.at(parent_index)->GetName();
                std::string derivative_residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
                if (pdb_glycam_map.find(derivative_residue_id) == pdb_glycam_map.end())
                {
                    pdb_glycam_map[derivative_residue_id] = std::vector<std::string>();
                }
                pdb_glycam_map[derivative_residue_id].push_back(gmml::ConvertT<int>(carbon_index) + new_name.substr(1));
                oligo_id_map[oligosaccharide].push_back(derivative_residue_id);
            }
            AtomVector o_linkage_derivative_atoms = AtomVector();
            std::string o_linkage_derivative_string =
                CheckxC_NxO_SO3(carbon_atom, oligosaccharide->root_->cycle_atoms_str_, 'O', o_linkage_derivative_atoms);
            if (o_linkage_derivative_string.compare("xC-O-SO3") == 0)
            {
                for (AtomVector::iterator it = o_linkage_derivative_atoms.begin();
                     it != o_linkage_derivative_atoms.end(); it++)
                {
                    MolecularModeling::Atom* atom  = *it;
                    std::string derivative_atom_id = atom->GetId();
                    if (pdb_glycam_map.find(derivative_atom_id) == pdb_glycam_map.end())
                    {
                        pdb_glycam_map[derivative_atom_id] = std::vector<std::string>();
                    }
                    pdb_glycam_map[derivative_atom_id].push_back("SO3");
                    oligo_id_map[oligosaccharide].push_back(derivative_atom_id);
                }
                std::string new_name = condensed_sequence_glycam06_residue_tree.at(parent_index)->GetName();
                std::string derivative_residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
                if (pdb_glycam_map.find(derivative_residue_id) == pdb_glycam_map.end())
                {
                    pdb_glycam_map[derivative_residue_id] = std::vector<std::string>();
                }
                pdb_glycam_map[derivative_residue_id].push_back(gmml::ConvertT<int>(carbon_index) + new_name.substr(1));
                oligo_id_map[oligosaccharide].push_back(derivative_residue_id);
            }
            AtomVector n_linkage_derivative_atoms_1   = AtomVector();
            std::string n_linkage_derivative_string_1 = CheckxC_NxO_PO3(
                carbon_atom, oligosaccharide->root_->cycle_atoms_str_, 'N', n_linkage_derivative_atoms_1);
            if (n_linkage_derivative_string_1.compare("xC-N-PO3") == 0)
            {
                for (AtomVector::iterator it = n_linkage_derivative_atoms_1.begin();
                     it != n_linkage_derivative_atoms_1.end(); it++)
                {
                    MolecularModeling::Atom* atom  = *it;
                    std::string derivative_atom_id = atom->GetId();
                    if (pdb_glycam_map.find(derivative_atom_id) == pdb_glycam_map.end())
                    {
                        pdb_glycam_map[derivative_atom_id] = std::vector<std::string>();
                    }
                    pdb_glycam_map[derivative_atom_id].push_back("PO3");
                    oligo_id_map[oligosaccharide].push_back(derivative_atom_id);
                }
                std::string new_name = condensed_sequence_glycam06_residue_tree.at(parent_index)->GetName();
                std::string derivative_residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
                if (pdb_glycam_map.find(derivative_residue_id) == pdb_glycam_map.end())
                {
                    pdb_glycam_map[derivative_residue_id] = std::vector<std::string>();
                }
                pdb_glycam_map[derivative_residue_id].push_back(gmml::ConvertT<int>(carbon_index) + new_name.substr(1));
                oligo_id_map[oligosaccharide].push_back(derivative_residue_id);
            }
            AtomVector o_linkage_derivative_atoms_1   = AtomVector();
            std::string o_linkage_derivative_string_1 = CheckxC_NxO_PO3(
                carbon_atom, oligosaccharide->root_->cycle_atoms_str_, 'O', o_linkage_derivative_atoms_1);
            if (o_linkage_derivative_string_1.compare("xC-O-PO3") == 0)
            {
                for (AtomVector::iterator it = o_linkage_derivative_atoms_1.begin();
                     it != o_linkage_derivative_atoms_1.end(); it++)
                {
                    MolecularModeling::Atom* atom  = *it;
                    std::string derivative_atom_id = atom->GetId();
                    if (pdb_glycam_map.find(derivative_atom_id) == pdb_glycam_map.end())
                    {
                        pdb_glycam_map[derivative_atom_id] = std::vector<std::string>();
                    }
                    pdb_glycam_map[derivative_atom_id].push_back("PO3");
                    oligo_id_map[oligosaccharide].push_back(derivative_atom_id);
                }
                std::string new_name = condensed_sequence_glycam06_residue_tree.at(parent_index)->GetName();
                std::string derivative_residue_id = oligosaccharide->root_->cycle_atoms_.at(0)->GetResidue()->GetId();
                if (pdb_glycam_map.find(derivative_residue_id) == pdb_glycam_map.end())
                {
                    pdb_glycam_map[derivative_residue_id] = std::vector<std::string>();
                }
                pdb_glycam_map[derivative_residue_id].push_back(gmml::ConvertT<int>(carbon_index) + new_name.substr(1));
                oligo_id_map[oligosaccharide].push_back(derivative_residue_id);
            }
        }
        index++;
    }

    // Recursively assign glycam naming to the monosaccharides of an oligosaccharide
    // The problem is that the order of child_oligo doesn't necessary match that of the glycam06 tree!
    for (unsigned int i = 0; i < oligosaccharide->child_oligos_.size(); i++)
    {
        this->ExtractOligosaccharideNamingMap(pdb_glycam_map, oligosaccharide->child_oligos_.at(i),
                                              condensed_sequence_glycam06_residue_tree, index, oligo_id_map, root_oligo,
                                              oligo_residue_map);
    }
}

void Assembly::TestUpdateResidueName2GlycamName(gmml::GlycamResidueNamingMap residue_glycam_map, std::string prep_file)
{
    for (AssemblyVector::iterator it = this->GetAssemblies().begin(); it != this->GetAssemblies().end(); it++)
    {
        (*it)->TestUpdateResidueName2GlycamName(residue_glycam_map, prep_file);
    }

    PrepFileSpace::PrepFile* prep                     = new PrepFileSpace::PrepFile(prep_file);
    PrepFileSpace::PrepFile::ResidueMap prep_residues = prep->GetResidues();
    ResidueVector residues                            = this->GetResidues();

    for (ResidueVector::iterator it2 = residues.begin(); it2 != residues.end(); it2++)
    {
        Residue* residue         = *it2;
        std::string residue_name = residue->GetName();
        std::string residue_id   = residue->GetId();
        if (residue_glycam_map.find(residue_id) != residue_glycam_map.end())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Residue found in residue_glycam_map, adding to terminal_residues, id is: " + residue_id);
            std::vector<std::string> glycam_names = residue_glycam_map[residue_id];
            std::string real_glycam_name          = "";
            if (glycam_names.size() == 1) // if residue_id and glycam names is one to one relationship
            {
                real_glycam_name = residue_glycam_map[residue_id][0];
            }
            else if (glycam_names.size() > 1) // if not one to one
            {
                // prevent residue from being named to its attached aglycones,whose names are also in the name vector
                std::vector<std::string> known_aglycon_names = {"ROH", "OME", "TBT"};
                for (std::vector<std::string>::iterator it3 = glycam_names.begin(); it3 != glycam_names.end(); it3++)
                {
                    // if this name is not a known aglycone
                    if (std::find(known_aglycon_names.begin(), known_aglycon_names.end(), *it3) ==
                        known_aglycon_names.end())
                    {
                        real_glycam_name = *it3;
                    }
                }
            }
            residue->SetName(real_glycam_name);

            // Right now cannot handle X configuration, exit when this happens
            // Now rename atoms using graph matching.
            // this->RenameAtoms(residue, prep);
        }
    }
}

void Assembly::MatchPdbAtoms2Glycam(std::map<Glycan::Oligosaccharide*, ResidueVector>& oligo_residue_map,
                                    std::string prep_file,
                                    std::map<Glycan::Oligosaccharide*, Pdb2glycamMatchingTracker*>& match_tracker)
{
    int local_debug = -1;
    if (local_debug > 0)
        //index++;
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Initializing prep file");
    }
    PrepFileSpace::PrepFile* prep = new PrepFileSpace::PrepFile(prep_file);
    for (std::map<Glycan::Oligosaccharide*, ResidueVector>::iterator mapit = oligo_residue_map.begin();
         mapit != oligo_residue_map.end(); mapit++)
    {
        if (local_debug > 0)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Oligosaccharide: " + std::to_string(std::distance(oligo_residue_map.begin(), mapit)));
        }
        Glycan::Oligosaccharide* oligo       = mapit->first;
        ResidueVector corresponding_residues = mapit->second;
        Pdb2glycamMatchingTracker* tracker   = new Pdb2glycamMatchingTracker();
        match_tracker[oligo]                 = tracker;

        std::string condensed_sequence = oligo->IUPAC_name_;
        MolecularModeling::Assembly template_assembly;
        if (local_debug > 0)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Condensed sequence is: " + condensed_sequence);
            gmml::log(__LINE__, __FILE__, gmml::INF, "Building template assembly");
        }
        template_assembly.BuildAssemblyFromCondensedSequence(condensed_sequence, prep);
        AtomVector template_atoms = template_assembly.GetAllAtomsOfAssembly();
        for (unsigned int i = 0; i < template_atoms.size(); i++)
        {}

        //PdbFileSpace::PdbFile *outputPdbFile = template_assembly.BuildPdbFileStructureFromAssembly();
        //outputPdbFile->Write("template.pdb");

        AtomVector target_atoms;
        for (ResidueVector::iterator resit = corresponding_residues.begin(); resit != corresponding_residues.end();
             resit++)
        {
            Residue* target_residue     = *resit;
            AtomVector atoms_in_residue = target_residue->GetAtoms();
            target_atoms.insert(target_atoms.end(), atoms_in_residue.begin(), atoms_in_residue.end());
        }

        std::vector<std::string> target_atom_label_map;
        for (AtomVector::iterator it = target_atoms.begin(); it != target_atoms.end(); it++)
        {
            Atom* atom = *it;
            std::string label;
            label += atom->GetResidue()->GetName();
            label += "_";

            if (atom->GetIsCycle())
            {
                label += "Ring_";
            }
            else
            {
                label += "Side_";
            }

            label                       += atom->GetElementSymbol();
            // label += atom->DetermineChirality();
            target_atom_label_map.push_back(label);
            //std::cout << "Target atom ptr " << atom << " " << atom->GetResidue()->GetName() << "-" <<  atom->GetName() << " has label " << label << std::endl;
        }

        std::vector<std::string> template_atom_label_map;

        for (AtomVector::iterator it = template_atoms.begin(); it != template_atoms.end(); it++)
        {
            Atom* atom = *it;
            std::string label;
            label += atom->GetResidue()->GetName();
            label += "_";

            if (atom->GetIsCycle())
            {
                label += "Ring_";
            }
            else
            {
                label += "Side_";
            }

            label                         += atom->GetElementSymbol();
            // label += atom->DetermineChirality();
            template_atom_label_map.push_back(label);
            //std::cout << "Template atom " << atom->GetResidue()->GetName() << "-" << atom->GetName() << " has label " << label << std::endl;
        }

        Atom* target_start_atom            = target_atoms[0];
        std::map<Atom*, Atom*> empty_match = std::map<Atom*, Atom*>();
        std::vector<std::map<Atom*, Atom*>> target_template_vertex_match;

        target_template_vertex_match.push_back(empty_match);
        std::vector<std::map<Atom*, Atom*>> template_target_vertex_match;
        template_target_vertex_match.push_back(empty_match);

        // std::vector<std::map<Atom*, Atom*> > all_isomorphisms;
        std::vector<Atom*> target_insertion_order;

        for (MolecularModeling::AtomVector::iterator template_it = template_atoms.begin();
             template_it != template_atoms.end(); template_it++)
        {
            Atom* template_atom = *template_it;
			int template_index = std::distance(template_atoms.begin(), template_it);
            int depth           = 0;
            if (local_debug > 0)
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, "Template atom is: " + template_atom->GetName());
                gmml::log(__LINE__, __FILE__, gmml::INF, "Depth is: " + std::to_string(depth));
                gmml::log(__LINE__, __FILE__, gmml::INF, "Target atom is: " + target_start_atom->GetName());
                gmml::log(__LINE__, __FILE__, gmml::INF, "About to run RecursiveMoleculeSubgraphMatching()");
            }
            this->RecursiveMoleculeSubgraphMatching(
                target_start_atom, target_atoms, 0, template_atom, template_atoms, template_index, target_atom_label_map, template_atom_label_map,
                target_template_vertex_match, template_target_vertex_match, target_insertion_order, tracker, depth);
        }
    }
    return;
}

bool Assembly::IfVertexAlreadyMatched(Atom* vertex_atom, std::vector<std::map<Atom*, Atom*>>& match_map)
{
    /*if (match_map.find(vertex_atom) != match_map.end()){
        return true;
    }
     */
    for (std::vector<std::map<Atom*, Atom*>>::iterator mapit = match_map.begin(); mapit != match_map.end(); mapit++)
    {
        std::map<Atom*, Atom*>& existing_match = *mapit;
        bool vertex_atom_found                 = false;
        for (std::map<Atom*, Atom*>::iterator mapit2 = existing_match.begin(); mapit2 != existing_match.end(); mapit2++)
        {
            if (mapit2->first == vertex_atom)
            {
                vertex_atom_found = true;
                break;
            }
        }
        if (vertex_atom_found)
        {
            return true;
        }
    }
    return false;
}

void Assembly::RemoveDownstreamMatches(Atom* target_atom,
                                       std::vector<std::map<Atom*, Atom*>>& target_template_vertex_match,
                                       std::vector<std::map<Atom*, Atom*>>& template_target_vertex_match,
                                       std::vector<Atom*> target_insertion_order)
{
    // std::cout << "Removing downstream matches for target atom " << target_atom->GetName() << std::endl;
    MolecularModeling::AtomVector::iterator target_match_position = target_insertion_order.end();
    for (MolecularModeling::AtomVector::iterator atom_it = target_insertion_order.begin();
         atom_it != target_insertion_order.end(); atom_it++)
    {
        if (*atom_it == target_atom)
        {
            target_match_position = atom_it;
            break;
        }
    }
    target_match_position++;

    if (target_match_position < target_insertion_order.end())
    {
        for (MolecularModeling::AtomVector::iterator atom_it = target_match_position;
             atom_it != target_insertion_order.end(); atom_it++)
        {
            Atom* downstream_target = *atom_it;
            // Remove downstream matches
            for (std::vector<std::map<Atom*, Atom*>>::iterator vec_it = target_template_vertex_match.begin();
                 vec_it != target_template_vertex_match.end(); vec_it++)
            {
                std::map<Atom*, Atom*>& current_match = *vec_it;
                if (current_match.find(downstream_target) != current_match.end())
                {
                    Atom* corresponding_template = current_match[downstream_target];
                    for (std::vector<std::map<Atom*, Atom*>>::iterator vec_it2 = template_target_vertex_match.begin();
                         vec_it2 != template_target_vertex_match.end(); vec_it2++)
                    {
                        std::map<Atom*, Atom*>& current_match2 = *vec_it2;
                        current_match2.erase(corresponding_template);
                    }
                }
            }
        }

        // Yao Mar 02 2022, this for loop seems to be a bug. Pdb2glycam on 1umi failed. Commentting this out 1umi
        // worked. Uncomment this later on if it's not a bug. Permanently delete this if later it's certain that this is
        // a bug. Remove the match itself
        /*for (MolecularModeling::AtomVector::iterator atom_it = target_match_position; atom_it !=
        target_insertion_order.end(); atom_it++){ Atom* downstream_target = *atom_it; for (std::vector<std::map<Atom*,
        Atom*> >::iterator vec_it = target_template_vertex_match.begin(); vec_it != target_template_vertex_match.end();
        vec_it++){ std::map<Atom*, Atom*>& current_match = *vec_it; current_match.erase(downstream_target);
            }
        }*/

        target_insertion_order.erase(target_match_position, target_insertion_order.end());
    }
}

bool Assembly::IfMatchScenarioAlreadyExists(std::map<Atom*, Atom*>& target_template_vertex_match,
                                            Pdb2glycamMatchingTracker* tracker)
{
    bool match_already_exists                             = false;
    std::vector<std::map<Atom*, Atom*>>& all_isomorphisms = tracker->all_isomorphisms;

    for (std::vector<std::map<Atom*, Atom*>>::iterator existing_isos_it = all_isomorphisms.begin();
         existing_isos_it != all_isomorphisms.end(); existing_isos_it++)
    {

        std::map<Atom*, Atom*>& current_existing_match         = *existing_isos_it;
        bool current_existing_match_identical_to_current_match = true;

        for (std::map<Atom*, Atom*>::iterator current_match_it = target_template_vertex_match.begin();
             current_match_it != target_template_vertex_match.end(); current_match_it++)
        {
            Atom* current_match_target        = current_match_it->first;
            Atom* current_match_template      = current_match_it->second;
            bool current_match_pair_duplicate = false;

            for (std::map<Atom*, Atom*>::iterator existing_match_it = current_existing_match.begin();
                 existing_match_it != current_existing_match.end(); existing_match_it++)
            {

                Atom* existing_match_target   = existing_match_it->first;
                Atom* existing_match_template = existing_match_it->second;

                if (current_match_target == existing_match_target && current_match_template == existing_match_template)
                {
                    // std::cout << "Duplicate target-teplate pair: " << current_match_target->GetResidue()->GetName()
                    // << "-" << current_match_target->GetName() << std::endl;
                    current_match_pair_duplicate = true;
                    break;
                }
            }

            if (!current_match_pair_duplicate)
            {
                current_existing_match_identical_to_current_match = false;
                break;
            }
        }

        if (current_existing_match_identical_to_current_match)
        {
            match_already_exists = true;
            break;
        }
    }

    return match_already_exists;
}

int Assembly::RecursiveMoleculeSubgraphMatching(Atom* target_atom, AtomVector& target_atoms, int target_atom_index, 
												Atom* template_atom, AtomVector& template_atoms, int template_atom_index,
                                                std::vector<std::string>& target_atom_label_map,
                                                std::vector<std::string>& template_atom_label_map,
                                                std::vector<std::map<Atom*, Atom*>>& target_template_vertex_match,
                                                std::vector<std::map<Atom*, Atom*>>& template_target_vertex_match,
                                                std::vector<Atom*>& target_insertion_order,
                                                Pdb2glycamMatchingTracker* pdb2glycam_matching_tracker,
                                                int previous_call_depth)
{
	//std::cout << "Trying to match target atom " << target_atom->GetResidue()->GetName() << "-" << target_atom->GetName() << " and " << template_atom->GetResidue()->GetName() << "-" << template_atom->GetName() << std::endl;
    // This function is intended to find multiple matches, but right now, only the 1st match is correct, the rest are
    // wrong.

    // Return value:
    // 0-Target vertex already matched.(Current not implemented, returns 5 when it's supposed to return 0)
    // 1-Template vertex already matched
    // 2-No equivalent template vertex exists(does not satisfy either vertex or edge equivalence)
    // 3-At least one equivalent template vertex exists, but a irreconcilable mismatches occur downstream.
    // 4-An isomorphism is found downstream
    // 5-A partial match is in progress, i.e. so far all vertex matches, but has not yet reached a complete isomorphism.
    // 1-3 are mismatch, 0,4-5 are successful match.

    int this_call_depth = previous_call_depth + 1;
    if (this->IfVertexAlreadyMatched(target_atom, target_template_vertex_match))
    {
        return 0; // Actually should return 0. Target vertex already matched by a previous pathway
    }

    if (this->IfVertexAlreadyMatched(template_atom, template_target_vertex_match))
    {
        return 1;
    }

    if (this->AtomVertexMatch(target_atom_index, template_atom_index, target_atom_label_map, template_atom_label_map))
    {
        target_template_vertex_match[0][target_atom]   = template_atom;
        template_target_vertex_match[0][template_atom] = target_atom;
        target_insertion_order.push_back(target_atom);
    }
    else
    {
        return 2; // No equivalent template vertex exists.
    }

    bool at_least_one_isomorphism_found_downstream = false;
    int downstream_matches_count                   = 0;
    std::map<int, bool> match_status_tracker;

    if (target_template_vertex_match[0].size() == target_atoms.size())
    {
        bool match_scenario_already_exists =
            this->IfMatchScenarioAlreadyExists(target_template_vertex_match[0], pdb2glycam_matching_tracker);
        if (!match_scenario_already_exists)
        {
            pdb2glycam_matching_tracker->all_isomorphisms.push_back(target_template_vertex_match[0]);
        }
		//std::cout << "An isomorphism match is found upon " << target_atom->GetResidue()->GetName() << "-" << target_atom->GetName() << std::endl;
        return 4; // An isomorphism match is found.
    }

    AtomVector target_atom_neighbors = target_atom->GetNode()->GetNodeNeighbors();
    // Filter for unmatched target neighbors
    AtomVector unmatched_target_neighbors;
    for (AtomVector::iterator target_it = target_atom_neighbors.begin(); target_it != target_atom_neighbors.end();
         target_it++)
    {
		MolecularModeling::Atom* neighbor = *target_it;
		if (std::find(target_atoms.begin(), target_atoms.end(), neighbor) == target_atoms.end()) continue;
        if (!this->IfVertexAlreadyMatched(*target_it, target_template_vertex_match))
        {
			//std::cout << "Found unmatched target neighbor " << (*target_it)->GetResidue()->GetName() << "-" << (*target_it)->GetName() << std::endl;
            unmatched_target_neighbors.push_back(*target_it);
        }
    }
    // If current target atom match, but there are no unmatched neighbor (i.e. all neighbors matched or dead end),
    // downstream match is deemed successful.
    if (unmatched_target_neighbors.empty())
    {
        return 5; // If an isomorphism isn't found above, and there is no downstream neighbors to match(dead end), it is
                  // a partial match.
    }

    AtomVector template_atom_neighbors = template_atom->GetNode()->GetNodeNeighbors();
    // Filter for unmatched template neighbors
    AtomVector unmatched_template_neighbors;
    for (AtomVector::iterator template_it = template_atom_neighbors.begin();
         template_it != template_atom_neighbors.end(); template_it++)
    {
        if (!this->IfVertexAlreadyMatched(*template_it, template_target_vertex_match))
        {
            unmatched_template_neighbors.push_back(*template_it);
        }
    }
    // If there are less unmatched template neighbor than are there unmatched targets, downstream matching is guaranteed
    // to fail. No need to actually calculate.
    if (unmatched_template_neighbors.size() < unmatched_target_neighbors.size())
    {
        std::stringstream mismatch_notice;
        mismatch_notice << "Less bonded template neighbors than target neighbors";
        pdb2glycam_matching_tracker->failures.emplace_back(new Pdb2glycamMatchingFailInfo(
            pdb2glycam_matching_tracker, this_call_depth, target_atom, mismatch_notice.str()));
        if (this_call_depth > pdb2glycam_matching_tracker->largest_iteration_length)
        {
            pdb2glycam_matching_tracker->largest_iteration_length = this_call_depth;
        }
        return 3; // Current target vertex match, but downstream mismatch is impossible.
        //return 5; // Current target vertex match, but downstream mismatch is impossible.
    }

    AtomVector matched_template_neighbors;

    for (unsigned int i = 0; i < unmatched_target_neighbors.size(); i++)
    {
		MolecularModeling::Atom* unmatched = unmatched_target_neighbors[i];
        AtomVector::iterator unmatched_it = std::find(target_atoms.begin(), target_atoms.end(), unmatched);
        int unmatched_target_index = std::distance(target_atoms.begin(), unmatched_it);

		for (unsigned int j = 0; j < unmatched_template_neighbors.size(); j++){
			MolecularModeling::Atom* unmatched_template = unmatched_template_neighbors[j];	
			AtomVector::iterator unmatched_template_it = std::find(template_atoms.begin(), template_atoms.end(), unmatched_template);
        	int unmatched_template_index = std::distance(template_atoms.begin(), unmatched_template_it);

            if (std::find(matched_template_neighbors.begin(), matched_template_neighbors.end(), unmatched_template) == matched_template_neighbors.end())
            {
				//std::cout << "Descending into " << next_target_atom->GetResidue()->GetName() << "-" << next_target_atom->GetName() << " corresponding template: " << template_atom_to_match->GetResidue()->GetName() << "-" << template_atom_to_match->GetName() << std::endl;

                int downstream_match_status = this->RecursiveMoleculeSubgraphMatching(
                    unmatched, target_atoms, unmatched_target_index, unmatched_template, template_atoms, unmatched_template_index, target_atom_label_map,
                    template_atom_label_map, target_template_vertex_match, template_target_vertex_match,
                    target_insertion_order, pdb2glycam_matching_tracker, this_call_depth);
					//std::cout << "Done " << unmatched->GetResidue()->GetName() << "-" << unmatched->GetName() << " returned " << downstream_match_status << std::endl;

                if (downstream_match_status == 4 || downstream_match_status == 5 || downstream_match_status == 0)
                {
                    match_status_tracker[unmatched_target_index] = true;

                    bool all_matches_true = true;
                    for (std::map<int, bool>::iterator mapit2 = match_status_tracker.begin();
                         mapit2 != match_status_tracker.end(); mapit2++)
                    {
                        if (!mapit2->second)
                        {
                            all_matches_true = false;
                            break;
                        }
                    }
                    if (all_matches_true)
                    {
                        downstream_matches_count++;
                    }

                    if (downstream_match_status == 4)
                    {
                        at_least_one_isomorphism_found_downstream = true;
                        this->RemoveDownstreamMatches(target_atom, target_template_vertex_match,
                                                      template_target_vertex_match, target_insertion_order);
                    }

					if (downstream_match_status == 4 || downstream_match_status == 5){
                    	matched_template_neighbors.push_back(unmatched_template);
					}
                }
            }

        } // for each unmatched template
    } // for1

    if (downstream_matches_count > 0)
    {
        if (at_least_one_isomorphism_found_downstream)
        {
            this->RemoveDownstreamMatches(target_atom, target_template_vertex_match, template_target_vertex_match,
                                          target_insertion_order);
            return 4; // Return upon successful isomorphism match
        }
        else
        {
            return 5;
        }
    }

    // If there is downstream mismatch, even if this target vertex itself matches, it does not support an isomorphism
    // match.
    if (this->IfVertexAlreadyMatched(target_atom, target_template_vertex_match))
    {
        this->RemoveDownstreamMatches(target_atom, target_template_vertex_match, template_target_vertex_match,
                                      target_insertion_order);
    }

    std::stringstream mismatch_notice;
    mismatch_notice << "This target itself has a match, but no downstream matches can be found.";
    pdb2glycam_matching_tracker->failures.emplace_back(new Pdb2glycamMatchingFailInfo(
        pdb2glycam_matching_tracker, this_call_depth, target_atom, mismatch_notice.str()));
    if (this_call_depth > pdb2glycam_matching_tracker->largest_iteration_length)
    {
        pdb2glycam_matching_tracker->largest_iteration_length = this_call_depth;
    }
    return 3; // Return upon downstream mismatch, but this target vertex itself matches.
	//return 5;
}

bool Assembly::CheckAndAcceptMatches(MolecularModeling::AtomVector& target_atoms,
                                     std::vector<std::map<Atom*, Atom*>>& target_template_vertex_match,
                                     std::vector<std::map<Atom*, Atom*>>& template_target_vertex_match,
                                     Pdb2glycamMatchingTracker* tracker)
{
    bool matches_accepted = false;
    for (std::vector<std::map<Atom*, Atom*>>::iterator vec_it = target_template_vertex_match.begin();
         vec_it != target_template_vertex_match.end(); vec_it++)
    {
        std::map<Atom*, Atom*>& pending_match = *vec_it;
        if (pending_match.size() == target_atoms.size())
        {
            bool match_scenario_already_exists = this->IfMatchScenarioAlreadyExists(pending_match, tracker);
            if (!match_scenario_already_exists)
            {
                // std::cout << "One isomorphism added" << std::endl;
                tracker->all_isomorphisms.push_back(pending_match);
            }
        }
    }
    // If matches are accepted, remove all but the first match in the pending match container
    if (target_template_vertex_match.size() > 1)
    {
        target_template_vertex_match.erase(target_template_vertex_match.begin() + 1,
                                           target_template_vertex_match.end());
    }
    if (template_target_vertex_match.size() > 1)
    {
        template_target_vertex_match.erase(template_target_vertex_match.begin() + 1,
                                           template_target_vertex_match.end());
    }
    return matches_accepted;
}

bool Assembly::AllAtomEdgesMatch(Atom* target_atom, AtomVector& target_atoms, Atom* template_atom, AtomVector& template_atoms,
                                 std::vector<std::string>& target_atom_label_map,
                                 std::vector<std::string>& template_atom_label_map)
{
    bool all_edges_match          = true;
    AtomVector all_target_edges   = target_atom->GetNode()->GetNodeNeighbors();
    AtomVector all_template_edges = template_atom->GetNode()->GetNodeNeighbors();

    for (AtomVector::iterator template_edge_it = all_template_edges.begin();
         template_edge_it != all_template_edges.end(); template_edge_it++)
    {
        // Atom* template_edge = *template_edge_it;
        // std::cout << "All template edges: " << template_edge->GetName() << "-" <<
        // template_atom_label_map[template_edge] << std::endl;
    }

    for (AtomVector::iterator target_edge_it = all_target_edges.begin(); target_edge_it != all_target_edges.end();
         target_edge_it++)
    {
        Atom* target_edge = *target_edge_it;
		AtomVector::iterator this_target_it = std::find(target_atoms.begin(), target_atoms.end(), target_edge);
		int this_target_index = std::distance(target_atoms.begin(), this_target_it);
        bool match_found  = false;
        AtomVector matched_edges;

        for (AtomVector::iterator template_edge_it = all_template_edges.begin();
             template_edge_it != all_template_edges.end(); template_edge_it++)
        {
            Atom* template_edge = *template_edge_it;
			AtomVector::iterator this_edge_it = std::find(template_atoms.begin(), template_atoms.end(), template_edge);
			int this_template_index = std::distance(template_atoms.begin(), this_edge_it);

            if (this->AtomVertexMatch(this_target_index, this_template_index, target_atom_label_map, template_atom_label_map) &&
                std::find(matched_edges.begin(), matched_edges.end(), template_edge) == matched_edges.end())
            {

                match_found = true;
                matched_edges.push_back(template_edge);
                break;
            }
        }

        if (!match_found)
        {
            // std::cout << "target edge " << target_edge->GetName() << "-" << target_atom_label_map[target_edge] <<
            // "does not match any template edges" << std::endl;
            all_edges_match = false;
            break;
        }
    }

    return all_edges_match;
}

bool Assembly::AtomVertexMatch(int target_atom_index, int template_atom_index,
                               std::vector<std::string>& target_atom_label_map,
                               std::vector<std::string>& template_atom_label_map)
{
	//std::cout << "Target vs template label: " << target_atom_label_map[target_atom_index] << " vs " << template_atom_label_map[template_atom_index] << std::endl;
    if (target_atom_label_map[target_atom_index] == template_atom_label_map[template_atom_index])
    {
        return true;
    }
    return false;
}

void Assembly::UpdateResidueName2GlycamName(gmml::GlycamResidueNamingMap residue_glycam_map, std::string prep_file)
{
    for (AssemblyVector::iterator it = this->GetAssemblies().begin(); it != this->GetAssemblies().end(); it++)
    {
        (*it)->UpdateResidueName2GlycamName(residue_glycam_map, prep_file);
    }

    PrepFileSpace::PrepFile* prep                     = new PrepFileSpace::PrepFile(prep_file);
    PrepFileSpace::PrepFile::ResidueMap prep_residues = prep->GetResidues();
    ResidueVector residues                            = this->GetResidues();

    ResidueVector updated_residues  = ResidueVector();
    int residue_sequence_number     = 0;
    ResidueVector terminal_residues = ResidueVector();
    for (ResidueVector::iterator it2 = residues.begin(); it2 != residues.end(); it2++)
    {
        Residue* residue         = *it2;
        std::string residue_name = residue->GetName();
        std::string residue_id   = residue->GetId();
        gmml::log(__LINE__, __FILE__, gmml::INF, "Looking for residue in residue_glycam_map, id is: " + residue_id);
        if (residue_glycam_map.find(residue_id) != residue_glycam_map.end())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Residue found in residue_glycam_map, adding to terminal_residues, id is: " + residue_id);
            terminal_residues.push_back(new Residue(this, ""));
            residue_sequence_number--;
            bool terminal                          = false;
            int residue_name_size                  = residue_name.size();
            // TODO: Done
            // Update to match the atoms of the residue with all the three-letter code prep residues to set the updated
            // residue name correspondingly Create a map between the current atom names and the corresponding prep atom
            // names when the match is found Add the residue mismatch into a structure for the Ontology usage
            std::vector<std::string> glycam_names  = residue_glycam_map[residue_id];
            gmml::GlycamAtomNameMap pdb_glycam_map = gmml::GlycamAtomNameMap();
            gmml::GlycamAtomNameMap glycam_pdb_map = gmml::GlycamAtomNameMap();
            std::string glycam_name                = *glycam_names.begin();
            ResidueVector query_residues           = ResidueVector();
            for (std::vector<std::string>::iterator name_it = glycam_names.begin(); name_it != glycam_names.end();
                 name_it++)
            {
                std::string glycam_residue_name = *name_it;
                // Right now cannot handle X configuraton, exit when this happens
                if (glycam_residue_name[glycam_residue_name.size() - 1] == 'X')
                {
                    //		    std::cout << "Unable to determine alpha/beta configuration." << std::endl;
                    //		    std::cout << "Aborting." << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                PrepFileSpace::PrepFile::ResidueMap customized_prep_residues = PrepFileSpace::PrepFile::ResidueMap();
                customized_prep_residues[glycam_residue_name]                = prep_residues[glycam_residue_name];
                PrepFileSpace::PrepFile* temp_prep_file                      = new PrepFileSpace::PrepFile();
                temp_prep_file->SetResidues(customized_prep_residues);
                temp_prep_file->Write(glycam_residue_name + ".prep");

                Assembly* prep_assembly = new Assembly();
                prep_assembly->BuildAssemblyFromPrepFile(temp_prep_file);
                prep_assembly->SetSourceFile(glycam_residue_name + ".prep");
                prep_assembly->BuildStructureByPrepFileInformation();
                remove((glycam_residue_name + ".prep").c_str());

                query_residues.push_back(prep_assembly->GetResidues().at(0));
            }
            // TODO:
            // The two residues atom names are matched but the geometry needs to be checked
            // if GeometryCheck returns true
            //             PatternMatching(residue, query_residues, pdb_glycam_map, glycam_pdb_map);
            AtomVector atoms                            = residue->GetAtoms();
            std::map<std::string, Residue*> residue_set = std::map<std::string, Residue*>();
            for (AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
            {
                MolecularModeling::Atom* atom = *it1;
                std::string atom_id           = atom->GetId();
                // Terminal residue naming

                if (residue_glycam_map.find(residue_id) != residue_glycam_map.end())
                {
                    terminal_residues.at(terminal_residues.size() - 1)->SetAssembly(this);
                    if (residue_glycam_map[residue_id].size() ==
                        1) // if residue_id and glycam names is one to one relationship
                    {
                        glycam_name = *residue_glycam_map[residue_id].begin();
                    }
                    else if (residue_glycam_map[residue_id].size() > 1) // if not one to one
                    {
                        // prevent residue from being named to its attached aglycones,whose names are also in the name
                        // vector
                        std::vector<std::string> known_aglycon_names = std::vector<std::string>();
                        known_aglycon_names.push_back("ROH");
                        known_aglycon_names.push_back("OME");
                        for (std::vector<std::string>::iterator it = residue_glycam_map[residue_id].begin();
                             it != residue_glycam_map[residue_id].end(); it++)
                        {
                            std::string possible_glycam_name = *it;
                            // if this name is not a known aglycone
                            if (std::find(known_aglycon_names.begin(), known_aglycon_names.end(),
                                          possible_glycam_name) == known_aglycon_names.end())
                            {
                                glycam_name = possible_glycam_name;
                            }
                        }
                    }
                    terminal_residues.at(terminal_residues.size() - 1)->SetName(glycam_name);
                    terminal_residues.at(terminal_residues.size() - 1)->AddAtom(atom);
                    std::vector<std::string> residue_id_tokens = gmml::Split(residue_id, "_");
                    std::string terminal_residue_id = residue_id_tokens.at(0) + "_" + residue_id_tokens.at(1) + "_" +
                                                      gmml::ConvertT<int>(residue_sequence_number) + "_" +
                                                      residue_id_tokens.at(3) + "_" + residue_id_tokens.at(4);
                    terminal_residues.at(terminal_residues.size() - 1)->SetId(terminal_residue_id);
                    atom->SetResidue(terminal_residues.at(terminal_residues.size() - 1));

                    if (!terminal)
                    {
                        updated_residues.push_back(terminal_residues.at(terminal_residues.size() - 1));
                        terminal = true;
                    }
                    int index = atom_id.find(residue_name);
                    if (index >= 0)
                    {
                        atom->SetId(atom_id.replace(index, residue_name_size, glycam_name));
                    }
                    std::vector<std::string> atom_id_tokens = gmml::Split(atom->GetId(), "_");
                    std::string atom_id_new                 = atom_id_tokens.at(0) + "_" + atom_id_tokens.at(1) + "_" +
                                              atom_id_tokens.at(2) + "_" + atom_id_tokens.at(3) + "_" +
                                              gmml::ConvertT<int>(residue_sequence_number) + "_" +
                                              atom_id_tokens.at(5) + "_" + atom_id_tokens.at(6);
                    atom->SetId(atom_id_new);
                }
                // Non-terminal residues glycam naming
                else
                {
                    // TODO:
                    // Update to match the atoms with the corresponding prep residue to change the atom names with
                    // respect to prep file Use the map to update the atom naming Add the atom name mismatch into a
                    // structure for the Ontology usage
                    std::string prep_atom_id;
                    if (pdb_glycam_map.find(atom_id) != pdb_glycam_map.end())
                    {
                        prep_atom_id = pdb_glycam_map[atom_id];
                    }
                    else
                    {
                        prep_atom_id = atom_id;
                    }
                    glycam_name = gmml::Split(prep_atom_id, "_")[2];
                    if (residue_set.find(glycam_name) == residue_set.end())
                    {
                        residue_set[glycam_name] = (new Residue(this, glycam_name));
                        residue_sequence_number--;
                        std::vector<std::string> res_id_tokens = gmml::Split(residue_id, "_");
                        std::string res_id                     = glycam_name + "_" + res_id_tokens.at(1) + "_" +
                                             gmml::ConvertT<int>(residue_sequence_number) + "_" + res_id_tokens.at(3) +
                                             "_" + res_id_tokens.at(4);
                        residue_set[glycam_name]->SetId(res_id);
                    }
                    std::string atom_name     = atom->GetName();
                    std::string new_atom_name = gmml::Split(prep_atom_id, "_")[0];
                    std::string new_atom_id   = atom_id;
                    int index                 = new_atom_id.find(residue_name);
                    if (index >= 0)
                    {
                        new_atom_id = new_atom_id.replace(index, residue_name_size, glycam_name);
                    }
                    index = new_atom_id.find(atom_name);
                    if (index >= 0)
                    {
                        new_atom_id = new_atom_id.replace(index, atom_name.size(), new_atom_name);
                    }
                    atom->SetId(new_atom_id);
                    atom->SetResidue(residue_set[glycam_name]);
                    residue_set[glycam_name]->AddAtom(atom);
                }
            }
            for (std::map<std::string, Residue*>::iterator it1 = residue_set.begin(); it1 != residue_set.end(); it1++)
            {
                updated_residues.push_back((*it1).second);
            }
        }
        else
        {
            updated_residues.push_back(residue);
        }
    }
    this->SetResidues(updated_residues);
}

/*! \todo  Give this function a return value.

Here is the error:

src/MolecularModeling/Assembly/glycamnaming.cc: In member function 'bool
MolecularModeling::Assembly::PatternMatching(MolecularModeling::Residue*, MolecularModeling::Assembly::ResidueVector,
gmml::GlycamAtomNameMap&, gmml::GlycamAtomNameMap&)': src/MolecularModeling/Assembly/glycamnaming.cc:568:1: warning: no
return statement in function returning non-void [-Wreturn-type]
 }
 ^
 */
bool Assembly::PatternMatching(Residue* residue, ResidueVector query_residues, gmml::GlycamAtomNameMap& pdb_glycam_map,
                               gmml::GlycamAtomNameMap& glycam_atom_map)
{
    CreatePrunedMatchingGraph(residue, query_residues);
    for (ResidueVector::iterator it = query_residues.begin(); it != query_residues.end(); it++)
    {
        Residue* query_residue                                 = *it;
        AtomVector lowest_degree_atoms                         = query_residue->GetAtomsWithLowestIntraDegree();
        std::stack<BacktrackingElements*> backtracking_stack_1 = std::stack<BacktrackingElements*>();
        for (int i = lowest_degree_atoms.size() - 1; i >= 0; i--)
        {
            backtracking_stack_1.push(
                new BacktrackingElements(pdb_glycam_map, glycam_atom_map, lowest_degree_atoms, i));
        }
        if (!backtracking_stack_1.empty())
        {
            // backtracking point
            BacktrackingElements* backtracking_point_1 = backtracking_stack_1.top();
            backtracking_stack_1.pop();
            pdb_glycam_map                         = backtracking_point_1->pdb_glycam_map_;
            glycam_atom_map                        = backtracking_point_1->glycam_atom_map_;
            lowest_degree_atoms                    = backtracking_point_1->atoms_;
            MolecularModeling::Atom* source_atom   = lowest_degree_atoms.at(backtracking_point_1->index_);
            AtomVector source_atom_intra_neighbors = source_atom->GetNode()->GetIntraNodeNeighbors();
            std::stack<BacktrackingElements*> backtracking_stack_2 = std::stack<BacktrackingElements*>();
            for (int i = source_atom_intra_neighbors.size() - 1; i >= 0; i--)
            {
                backtracking_stack_2.push(
                    new BacktrackingElements(pdb_glycam_map, glycam_atom_map, source_atom_intra_neighbors, i));
            }
            if (!backtracking_stack_2.empty())
            {
                // backtracking point
                BacktrackingElements* backtracking_point_2 = backtracking_stack_2.top();
                backtracking_stack_2.pop();
                pdb_glycam_map                        = backtracking_point_2->pdb_glycam_map_;
                glycam_atom_map                       = backtracking_point_2->glycam_atom_map_;
                source_atom_intra_neighbors           = backtracking_point_2->atoms_;
                MolecularModeling::Atom* mapped_atom  = source_atom_intra_neighbors.at(backtracking_point_2->index_);
                pdb_glycam_map[mapped_atom->GetId()]  = source_atom->GetId();
                glycam_atom_map[source_atom->GetId()] = mapped_atom->GetId();
                AtomVector source_atom_neighbors      = source_atom->GetNode()->GetNodeNeighbors();
                std::queue<MolecularModeling::Atom*> to_visit = std::queue<MolecularModeling::Atom*>();
                std::queue<MolecularModeling::Atom*> parents  = std::queue<MolecularModeling::Atom*>();
                for (unsigned int i = 0; i < source_atom_neighbors.size(); i++)
                {
                    if (glycam_atom_map.find(source_atom_neighbors.at(i)->GetId()) == glycam_atom_map.end())
                    {
                        to_visit.push(source_atom_neighbors.at(i));
                        parents.push(mapped_atom);
                    }
                }
                std::stack<BacktrackingElements*> backtracking_stack_3 = std::stack<BacktrackingElements*>();
                while (!to_visit.empty())
                {
                    MolecularModeling::Atom* atom   = to_visit.front();
                    MolecularModeling::Atom* parent = parents.front();
                    to_visit.pop();
                    parents.pop();
                    AtomVector atom_intra_neighbors                      = atom->GetNode()->GetIntraNodeNeighbors();
                    AtomVector atom_intra_neighbors_filtered_with_parent = AtomVector();
                    for (unsigned int i = 0; i < atom_intra_neighbors.size(); i++)
                    {
                        AtomVector neighbors = atom_intra_neighbors.at(i)->GetNode()->GetNodeNeighbors();
                        if (find(neighbors.begin(), neighbors.end(), parent) != neighbors.end() &&
                            pdb_glycam_map.find(atom_intra_neighbors.at(i)->GetId()) == pdb_glycam_map.end())
                        {
                            atom_intra_neighbors_filtered_with_parent.push_back(atom_intra_neighbors.at(i));
                        }
                    }
                    for (int i = atom_intra_neighbors_filtered_with_parent.size() - 1; i >= 0; i--)
                    {
                        backtracking_stack_3.push(new BacktrackingElements(
                            pdb_glycam_map, glycam_atom_map, atom_intra_neighbors_filtered_with_parent, i, to_visit));
                    }
                    if (!backtracking_stack_3.empty())
                    {
                        // backtracking point
                        BacktrackingElements* backtracking_point_3 = backtracking_stack_3.top();
                        backtracking_stack_3.pop();
                        pdb_glycam_map                            = backtracking_point_3->pdb_glycam_map_;
                        glycam_atom_map                           = backtracking_point_3->glycam_atom_map_;
                        atom_intra_neighbors_filtered_with_parent = backtracking_point_3->atoms_;
                        MolecularModeling::Atom* mapped =
                            atom_intra_neighbors_filtered_with_parent.at(backtracking_point_3->index_);
                        pdb_glycam_map[mapped->GetId()] = atom->GetId();
                        glycam_atom_map[atom->GetId()]  = mapped->GetId();
                        AtomVector atom_neighbors       = atom->GetNode()->GetNodeNeighbors();
                        for (unsigned int i = 0; i < atom_neighbors.size(); i++)
                        {
                            if (glycam_atom_map.find(atom_neighbors.at(i)->GetId()) == glycam_atom_map.end())
                            {
                                to_visit.push(atom_neighbors.at(i));
                                parents.push(mapped);
                            }
                        }
                    }
                }
                // Check if the mapping is complete
                // continue to the next residue
                // else
                // backtrack
            }
            // Check if the mapping is complete
            // continue to the next residue
            // else
            // backtrack
        }
        // Check if the mapping is complete
        // continue to the next residue
        // else
        // backtrack
    }
    // Check if the mapping is complete
    // return true
    // else
    // return false
}

void Assembly::CreateLabelGraph(Residue* residue, Residue* query_residue)
{
    AtomVector query_atoms = query_residue->GetAtoms();
    AtomVector atoms       = residue->GetAtoms();
    for (AtomVector::iterator it = query_atoms.begin(); it != query_atoms.end(); it++)
    {
        MolecularModeling::Atom* query_atom   = *it;
        AtomNode* query_atom_node             = query_atom->GetNode();
        AtomVector query_atom_intra_neighbors = AtomVector();
        if (query_atom_node == NULL)
        {
            query_atom_node = new AtomNode();
        }
        for (AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
        {
            MolecularModeling::Atom* atom = *it1;
            AtomNode* atom_node           = atom->GetNode();
            if (atom_node == NULL)
            {
                atom_node = new AtomNode();
            }
            if (query_atom_node->GetElementLabel().compare(atom_node->GetElementLabel()) == 0)
            {
                query_atom_intra_neighbors.push_back(atom);
            }
        }
        query_atom_node->SetIntraNodeNeighbors(query_atom_intra_neighbors);
    }
}

void Assembly::PruneLabelGraphByNeighboringLabels(Residue* query_residue)
{
    AtomVector query_atoms = query_residue->GetAtoms();
    for (AtomVector::iterator it = query_atoms.begin(); it != query_atoms.end(); it++)
    {
        MolecularModeling::Atom* query_atom = *it;
        AtomNode* query_atom_node           = query_atom->GetNode();
        if (query_atom_node != NULL)
        {
            std::string query_atom_neighborhood_label     = query_atom_node->CreateNeighboringLabel();
            AtomVector query_atom_intra_neighbors         = query_atom_node->GetIntraNodeNeighbors();
            AtomVector updated_query_atom_intra_neighbors = AtomVector();
            for (AtomVector::iterator it1 = query_atom_intra_neighbors.begin(); it1 != query_atom_intra_neighbors.end();
                 it1++)
            {
                MolecularModeling::Atom* query_atom_intra_neighbor = *it1;
                AtomNode* query_atom_intra_neighbor_node           = query_atom_intra_neighbor->GetNode();
                if (query_atom_intra_neighbor_node != NULL)
                {
                    std::string query_atom_intra_neighbor_neighborhood_label =
                        query_atom_intra_neighbor_node->CreateNeighboringLabel();
                    if (query_atom_neighborhood_label.compare(query_atom_intra_neighbor_neighborhood_label) == 0)
                    {
                        updated_query_atom_intra_neighbors.push_back(query_atom_intra_neighbor);
                    }
                }
            }
            query_atom_node->SetIntraNodeNeighbors(updated_query_atom_intra_neighbors);
        }
    }
}

void Assembly::CreatePrunedMatchingGraph(Residue* residue, ResidueVector query_residues)
{
    if (residue->GraphElementLabeling())
    {
        for (ResidueVector::iterator it = query_residues.begin(); it != query_residues.end(); it++)
        {
            Residue* query_residue = *it;
            if (query_residue->GraphElementLabeling())
            {
                this->CreateLabelGraph(residue, query_residue);
                this->PruneLabelGraphByNeighboringLabels(query_residue);
                //                query_residue->Print();
            }
        }
    }
}
