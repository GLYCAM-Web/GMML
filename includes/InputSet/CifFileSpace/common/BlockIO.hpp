/*
FILE:     BlockIO.hpp
*/
/*
VERSION:  7.105
*/
/*
DATE:     10/21/2013
*/
/*
  Comments and Questions to: sw-help@rcsb.rutgers.edu
*/
/*
COPYRIGHT 1999-2013 Rutgers - The State University of New Jersey

This SOFTWARE has been modified from 
the version originally obtained from RUTGERS.  
*/


#ifndef BLOCKIO_H
#define BLOCKIO_H

#include "rcsb_types.hpp"

const unsigned int WORDSPERBLOCK = 2048;
const unsigned int BLKSIZE = 8192;

class BlockIO
{
public:
    BlockIO();
    ~BlockIO();

    void AssociateBuffer(char** newBuffer);

    unsigned int ReadBlock(const int fd, const UInt32 blockNum);
    unsigned int WriteBlock(const int fd, const UInt32 blockNum);

private:
    UInt32 _buffer[WORDSPERBLOCK]; // A buffer for reading/writing blocks

};

#endif

