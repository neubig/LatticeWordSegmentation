/*
* Copyright 2016, Takahiro Shinozaki
* 
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
* 
*     http://www.apache.org/licenses/LICENSE-2.0
* 
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#ifndef DISCRETE_HMM_HPP__
#define DISCRETE_HMM_HPP__

#include <iostream>
#include <fst/fstlib.h>
#include <fst/vector-fst.h>
#include <cmath>

#include "definitions.hpp"

class DiscreteHMM {

private:

  int codebookSize_; // Assume discrete HMM
  int numUnits_;
  int numUnitStates_;
  std::vector<string> unitNames_;

  // Assume left-to-right topology. 
  std::vector<std::vector<double>> stateObsPriorParam_;
  std::vector<std::vector<double>> stateObsCount_;
  std::vector<std::vector<double>> stateTransPriorParam_;
  std::vector<std::vector<double>> stateTransCount_;

  std::vector<std::vector<double>> stateObsProb_;
  std::vector<std::vector<double>> stateTransProb_;



public:

  DiscreteHMM(int codebookSize, int numUnits, int numUnitStates=3, double obsPriorParam=1.0, double transPriorParam=1.0);
    
  void sampleObsProbs();

  void sampleTransProbs();

  double stateProb(int stId, int obsCode);

  double transProb(int stId, int dest);

  double addStateObsCount(int stId, int obsCode, float count);

  double addStateTransCount(int stId, int dest, float count);

  int addSampleCounts(const std::vector<int> &alignment, const std::vector<int> &feaSeq);

  double removeStateObsCount(int stId, int obsCode, float count);

  double removeStateTransCount(int stId, int dest, float count);

  int removeSampleCounts(const std::vector<int> &alignment, const std::vector<int> &feaSeq);

  LogVectorFst frameStateWfst(const std::vector<int> &feaSeq);

  LogVectorFst monoPhoneWfst();

  // Returns a sequence of HMM state ID (i.e, input labels of sample)
  std::vector<int> parseSample(const fst::Fst<fst::LogArc> & sample);

  std::vector<double> dirichletSample(std::vector<double> &alpha, std::vector<double> *diffAlpha);

  double gammaSample(double a, double scale);
  double exponSample(double l);

};

#endif
