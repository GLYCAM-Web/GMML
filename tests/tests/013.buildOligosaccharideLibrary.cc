#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h> // mkdir
#include "includes/gmml.hpp"

int main(int argc, char** argv)
{
    if (argc != 5)
    {
        std::cerr << "Usage:   ./buildOligosaccharideLibrary <List IDs and sequences> <Char delimiter used in list> <Folder to put outputs> <Prepfile>\n";
        std::cerr << "Example: ./buildOligosaccharideLibrary smallTestLibrary.txt _ testOutputs/ ../dat/prep/GLYCAM_06j-1_GAGS.prep\n";
        std::cerr << "Don't use a delimiter that appears in glycan sequences or ids. Like - or , or [] etc\n";
        std::exit(EXIT_FAILURE); 
    }
    //std::cout << "argv [1] is" << argv[1] << ", " << "\n";
    //Convert command line inputs to legible variables
    std::ifstream infile(argv[1]);
    char delimiter = argv[2][0]; // The second [0] gets me the first element of the argv which is type char**
    std::string outputFolderName = argv[3];
    std::string prepFile = argv[4];
    struct stat info;
    if( stat(argv[3], &info) != 0)
    {
        std::cerr << "Folder " << outputFolderName << "/ does not exist and it ain't my job to make it.\n";
        std::exit(EXIT_FAILURE);  
    }
    std::string line;
    while (std::getline(infile, line))
    {
        StringVector splitLine = gmml::splitStringByDelimiter(line, delimiter);
        std::string inputSequence = splitLine.at(1);
        std::cout << "\n*********************\nBuilding " << inputSequence << "\n*********************\n";
        CondensedSequence::carbohydrateBuilder carbBuilder(inputSequence, prepFile);
        if (carbBuilder.IsStatusOk())
        {
            carbBuilder.Print();
            std::string inputGlycanID = splitLine.at(0);
            carbBuilder.GenerateSingle3DStructureSingleFile(outputFolderName, "PDB", inputGlycanID);
        }
        else
        {
            std::cerr << "Error thrown by gmml was: " << carbBuilder.GetStatusMessage() << std::endl;
        }
    }
}
