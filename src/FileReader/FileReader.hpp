// ----------------------------------------------------------------------------
/**
   File: FileReader.hpp

   Status:         Version 1.0
   Language: C++

   License: UPB licence

   Copyright (c) <2013> <University of Paderborn>
   Permission is hereby granted, free of charge, to any person
   obtaining a copy of this software and associated documentation
   files (the "Software"), to deal in the Software without restriction,
   including without limitation the rights to use, copy, modify and
   merge the Software, subject to the following conditions:

   1.) The Software is used for non-commercial research and
       education purposes.

   2.) The above copyright notice and this permission notice shall be
       included in all copies or substantial portions of the Software.

   3.) Publication, Distribution, Sublicensing, and/or Selling of
       copies or parts of the Software requires special agreements
       with the University of Paderborn and is in general not permitted.

   4.) Modifications or contributions to the software must be
       published under this license. The University of Paderborn
       is granted the non-exclusive right to publish modifications
       or contributions in future versions of the Software free of charge.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
   OTHER DEALINGS IN THE SOFTWARE.

   Persons using the Software are encouraged to notify the
   Department of Communications Engineering at the University of Paderborn
   about bugs. Please reference the Software in your publications
   if it was used for them.


   Author: Oliver Walter

   E-Mail: walter@nt.uni-paderborn.de

   Description: io class for lattice reading, processing and writing

   Limitations: -

   Change History:
   Date         Author       Description
   2014         Walter       Initial
*/
// ----------------------------------------------------------------------------
#ifndef _FILEREADER_HPP_
#define _FILEREADER_HPP_

#include <fst/vector-fst.h>
#include "StringToIntMapper.hpp"
#include "../ParameterParser/ParameterParser.hpp"
#include "FileData.hpp"

/* class to read input files */
class FileReader {
  // parameters
  const ParameterStruct Params;

  static const bool PARSE_REFERENCES = true;

  /* internal functions: input */

  // Lattice reading commands
  void ReadHTKLattices(FileData & data) const;
  void ReadOpenFSTLattices(FileData & data) const;

  void ReadSegmentList(
    std::size_t InputFileId,
    std::string line, int debug_
  ) const;

  void ReadTextFiles(
    FileData & data,
    const std::vector<std::string > &InputFiles,
    std::vector<LogVectorFst> *InputFsts,
    std::vector<std::string> *FileNames,
    bool ParseReferences=false
  ) const;

  // read the symbols file mapping strings to integers
  void ReadSymbols(FileData & data) const;

  void ReadInputFilesFromList(FileData & data) const;

  // read reference transcription from file ReferenceFile
  void ReadReferenceTranscription(FileData & data) const;

  void ReadInputArcInfos(FileData & data) const;

  /* Modifications */
  void PruneLattices(
    FileData & data,
    double PruningFactor
  ) const;

  void ApplyAcousticModelScalingFactor(FileData & data) const;

  void ApplyWordEndTransducer(FileData & data) const;

  void ApplySentEndTransducer(FileData & data) const;

  void CalculateLatticePhonemeErrorRate(FileData & data) const;

  /* internal functions: output */
  void WriteOpenFSTLattices(const FileData & data) const;

  void WriteInputArcInfos(const FileData & data) const;

  bool IsSilence(
    std::string phone
  ) const;

  std::string GetSubstrAfterSep(
    std::string inStr,
    char sep
  ) const;

  void ReadPronDict(FileData & data) const;

  StateId AddPhoneToFst(
    FileData & data,
    LogVectorFst& LatticeFst, StateId State,
    const std::string& phone
  ) const;

public:
  /* constructor */
  // set some predefined symbols in string to int map and read data
  FileReader(
    ParameterStruct Params
  );

  void ReadData(FileData & data) const;

};

#endif
