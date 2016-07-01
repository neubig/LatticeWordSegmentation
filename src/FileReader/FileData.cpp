// ----------------------------------------------------------------------------
/**
   File: FileData.cpp
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
*/
// ----------------------------------------------------------------------------
#include "FileData.hpp"

FileData::FileData() {
  GlobalStringToInt.Insert(EPS_SYMBOL);
  GlobalStringToInt.Insert(PHI_SYMBOL);
  GlobalStringToInt.Insert(UNKBEGIN_SYMBOL);
  GlobalStringToInt.Insert(UNKEND_SYMBOL);
  GlobalStringToInt.Insert(SENTSTART_SYMBOL);
  GlobalStringToInt.Insert(SENTEND_SYMBOL);
}

FileData::FileData(const FileData& lhs):
  GlobalStringToInt(lhs.GlobalStringToInt),
  InitStringToInt(lhs.InitStringToInt),
  InitFsts(lhs.InitFsts),
  InitFileNames(lhs.InitFileNames),
  InputStringToInt(lhs.InputStringToInt),
  InputDiscreteSeqs(lhs.InputDiscreteSeqs),
  InputFsts(lhs.InputFsts),
  InputFileNames(lhs.InputFileNames),
  InputArcInfos(lhs.InputArcInfos),
  ReferenceStringToInt(lhs.ReferenceStringToInt),
  ReferenceFsts(lhs.ReferenceFsts),
  ReferenceFileNames(lhs.ReferenceFileNames)
{
}


const std::vector<std::string> &FileData::GetInputIntToStringVector() const
{
  return InputStringToInt.GetIntToStringVector();
}

const std::vector<std::string> &FileData::GetReferenceIntToStringVector() const
{
  return ReferenceStringToInt.GetIntToStringVector();
}


const std::vector<std::string> &FileData::GetInitIntToStringVector() const
{
  return InitStringToInt.GetIntToStringVector();
}

size_t FileData::GetNumInputs() const
{
  return InputFsts.size();
}
