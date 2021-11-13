#include <sstream>
#include <stdexcept>

#include "../../../includes/common.hpp"
//#include "../../../includes/utils.hpp"
#include "../../../includes/InputSet/CifFileSpace/ciffileprocessingexception.hpp"
#include "includes/CodeUtils/logging.hpp"

using CifFileSpace::CifFileProcessingException;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
CifFileProcessingException::CifFileProcessingException(const std::string &message)
    : line_number_(gmml::dNotSet), message_(message) {}

CifFileProcessingException::CifFileProcessingException(int line_number, const std::string &message)
    : line_number_(line_number), message_(message) {}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
/// Exception handler for Cif file exceptions
const char* CifFileProcessingException::what() const throw()
{
    what_ = "CifFile: " + message_;
    if (line_number_ != gmml::dNotSet)
    {
        std::stringstream ss;
        if(ss << line_number_)
        {
            what_ += " (line " + ss.str() + ")";
            gmml::log(__LINE__, __FILE__, gmml::ERR, what_.c_str());
            return what_.c_str();
        }
        else
        {
            gmml::log(__LINE__, __FILE__, gmml::ERR, "to_string: invalid conversion");
            throw std::invalid_argument("to_string: invalid conversion");       /// Invalid conversion from int to string
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::ERR, what_.c_str());
    return what_.c_str();
}

CifFileProcessingException::~CifFileProcessingException() throw() {}
