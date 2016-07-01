// ----------------------------------------------------------------------------
/**
   File: FileData.cpp

   Status:         Version 1.0
   Language: C++

   License: UPB licence

   Copyright (c) <2016> <University of Paderborn>
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


   Author: Thomas Glarner

   E-Mail: glarner@nt.uni-paderborn.de

   Description: io class for lattice reading, processing and writing

   Limitations: -

   Change History:
   Date         Author       Description
   2016         Glarner      Initial
*/
// ----------------------------------------------------------------------------
#ifndef _FILEDATA_HPP_
#define _FILEDATA_HPP_

#include <fst/vector-fst.h>
#include "definitions.hpp"
#include "StringToIntMapper.hpp"

class FileData{
  // a global string to int mapper which is updated with each reading process
  StringToIntMapper GlobalStringToInt;

  // Pronunciation dictionary
  PronDictType PronDict;

  // members for initialization fsts
  StringToIntMapper InitStringToInt;
  std::vector<LogVectorFst> InitFsts;
  // filename of read initialization files
  std::vector<std::string> InitFileNames;

  // membersf for input fsts
  StringToIntMapper InputStringToInt;
  std::vector<LogVectorFst> InputFsts;
  std::vector<std::string> InputFileNames;
  std::vector<ArcInfo> InputArcInfos; // ArcInfo members: {label, start, end}

  // members for reference fsts
  StringToIntMapper ReferenceStringToInt;
  std::vector<LogVectorFst> ReferenceFsts;
  std::vector<std::string> ReferenceFileNames;

public:
  /* Constructor */
  FileData();

  /* Copy Constructor */
  FileData(const FileData& lhs);

  /* interface */
  // get vector mapping integer to string (characters)
  const std::vector<std::string> &GetInputIntToStringVector() const;

  // integer to string mapping for the reference transcription (characters)
  const std::vector<std::string> &GetReferenceIntToStringVector() const;

  // integer to string mapping for the init transcription (characters)
  const std::vector<std::string> &GetInitIntToStringVector() const;

  size_t GetNumInputs() const;

  // Simple accessors
  const LogVectorFst &GetInputFst(size_t index) const { return InputFsts[index]; }
  const std::vector<LogVectorFst> &GetInputFsts() const { return InputFsts; }
  const LogVectorFst &GetReferenceFst(size_t index) const { return ReferenceFsts[index]; }
  const std::vector<LogVectorFst> &GetReferenceFsts() const { return ReferenceFsts; }
  const LogVectorFst &GetInitFst(size_t index) const { return InitFsts[index]; }
  const std::vector<LogVectorFst> &GetInitFsts() const { return InitFsts; }

  LogVectorFst &GetInputFst(size_t index) { return InputFsts[index]; }
  std::vector<LogVectorFst> &GetInputFsts() { return InputFsts; }
  LogVectorFst &GetReferenceFst(size_t index) { return ReferenceFsts[index]; }
  std::vector<LogVectorFst> &GetReferenceFsts() { return ReferenceFsts; }
  LogVectorFst &GetInitFst(size_t index) { return InitFsts[index]; }
  std::vector<LogVectorFst> &GetInitFsts() { return InitFsts; }

  const std::vector<std::string> &GetInputFileNames() const { return InputFileNames; }
  const std::vector<std::string> &GetReferenceFileNames() const { return ReferenceFileNames; }
  const std::vector<std::string> &GetInitFileNames() const { return InitFileNames; }
  const std::vector<ArcInfo> &GetInputArcInfos() const { return InputArcInfos; }

  std::vector<std::string> &GetInputFileNames() { return InputFileNames; }
  std::vector<std::string> &GetReferenceFileNames() { return ReferenceFileNames; }
  std::vector<std::string> &GetInitFileNames() { return InitFileNames; }
  std::vector<ArcInfo> &GetInputArcInfos() { return InputArcInfos; }

  const StringToIntMapper &GetGlobalStringToInt() const { return GlobalStringToInt; }
  const StringToIntMapper &GetInitStringToInt() const { return InitStringToInt; }
  const StringToIntMapper &GetInputStringToInt() const { return InputStringToInt; }
  const StringToIntMapper &GetReferenceStringToInt() const { return ReferenceStringToInt; }

  StringToIntMapper &GetGlobalStringToInt() { return GlobalStringToInt; }
  StringToIntMapper &GetInitStringToInt() { return InitStringToInt; }
  StringToIntMapper &GetInputStringToInt() { return InputStringToInt; }
  StringToIntMapper &GetReferenceStringToInt() { return ReferenceStringToInt; }

  PronDictType &GetPronDict() { return PronDict; }

};

#endif // _FILEDATA_HPP_
