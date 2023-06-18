#include "../../includes/CodeUtils/codetests.hpp"

using CodeUtils::CodeTests;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
CodeTests::CodeTests()
{}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
StringVector CodeTests::ListCodeTests()
{
    StringVector testlist;
    testlist.push_back("ProduceSegmentationFault");
    return testlist;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////

void CodeTests::ProduceSegmentationFault()
{
    int* a;
    a = NULL;
    std::cerr << "This program has requested a segfault.  Complying now." << std::flush;
    std::cout << *a << std::endl;
    return;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
