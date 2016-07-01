#include "DiscreteHMM.hpp"

using namespace fst;

DiscreteHMM::DiscreteHMM(int codebookSize, int numUnits, int numUnitStates, double obsPriorParam, double transPriorParam)
  : codebookSize_(codebookSize), numUnits_(numUnits), numUnitStates_(numUnitStates),
    stateObsPriorParam_(vector<vector<double>>(numUnits*numUnitStates, vector<double>(codebookSize, obsPriorParam))),
    stateObsCount_(vector<vector<double>>(numUnits*numUnitStates, vector<double>(codebookSize, 0.0))),
    stateTransPriorParam_(vector<vector<double>>(numUnits*numUnitStates, vector<double>(2, transPriorParam))),
  stateTransCount_(vector<vector<double>>(numUnits*numUnitStates, vector<double>(2, 0.0)))
{
  stateObsProb_.resize(numUnits_*numUnitStates_);
  stateTransProb_.resize(numUnits_*numUnitStates_);
  sampleObsProbs();
  sampleTransProbs();
}
    
void DiscreteHMM::sampleObsProbs() {
  for (unsigned i=0; i<stateObsProb_.size(); i++) {
    stateObsProb_[i] = dirichletSample(stateObsPriorParam_[i], &(stateObsCount_[i]));
  }
}

void DiscreteHMM::sampleTransProbs() {
  for (unsigned i=0; i<stateTransProb_.size(); i++) {
    stateTransProb_[i] = dirichletSample(stateTransPriorParam_[i], &(stateTransCount_[i]));
  }
}

double DiscreteHMM::stateProb(int stId, int obsCode) {
  return stateObsProb_[stId][obsCode];
}

double DiscreteHMM::transProb(int stId, int dest) { // dest: 0 (self-loop), 1 (transition to next state)
  return stateTransProb_[stId][dest];
}


double DiscreteHMM::addStateObsCount(int stId, int obsCode, float count=1.0) {
  return stateObsCount_[stId][obsCode] += count;
}

double DiscreteHMM::addStateTransCount(int stId, int dest, float count=1.0) {
  return stateTransCount_[stId][dest] += count;
}

int DiscreteHMM::addSampleCounts(const vector<int> &alignment, const vector<int> &feaSeq) {
  unsigned i;
  int stId, nextStId, obsCode;

  if (alignment.size() != feaSeq.size()) {
    cerr << "Error: Inconsistent alignment and feature length." << endl;
    exit(1);
  }
  for (i=0; i<alignment.size(); i++) {
    stId = alignment[i];
    obsCode = feaSeq[i];
    stateObsCount_[stId][obsCode] += 1.0;
  }
  for (i=0; i<alignment.size()-1; i++) {
    stId = alignment[i];
    nextStId = alignment[i+1];
    if ((nextStId-stId==0) || (nextStId-stId==1)) {
      stateTransCount_[stId][nextStId-stId] += 1.0;
    } else if (nextStId%numUnitStates_ == 0) {
      stateTransCount_[stId][1] += 1.0;
    } else {
      cerr << "Error: Unexpected HMM state sequence. Not a left-to-right HMM?" << endl;
      exit(1);
    }
  }
  stId = alignment[i];
  stateTransCount_[stId][1] += 1.0;
  return alignment.size();
}

double DiscreteHMM::removeStateObsCount(int stId, int obsCode, float count=1.0) {
  return stateObsCount_[stId][obsCode] -= count;
}

double DiscreteHMM::removeStateTransCount(int stId, int dest, float count=1.0) {
  return stateTransCount_[stId][dest] -= count;
}

int DiscreteHMM::removeSampleCounts(const vector<int> &alignment, const vector<int> &feaSeq) {
  unsigned i;
  int stId, nextStId, obsCode;

  if (alignment.size() != feaSeq.size()) {
    cerr << "Error: Inconsistent alignment and feature length." << endl;
    exit(1);
  }
  // state observation counts
  for (i=0; i<alignment.size(); i++) {
    stId = alignment[i];
    obsCode = feaSeq[i];
    stateObsCount_[stId][obsCode] -= 1.0;
  }
  // state transition counts
  for (i=0; i<alignment.size()-1; i++) {
    stId = alignment[i];
    nextStId = alignment[i+1];
    if ((nextStId-stId==0) || (nextStId-stId==1)) {
      stateTransCount_[stId][nextStId-stId] -= 1.0;
    } else if (nextStId%numUnitStates_ == 0) {
      stateTransCount_[stId][1] -= 1.0;
    } else {
      cerr << "Error: Unexpected HMM state sequence. Not a left-to-right HMM?" << endl;
      exit(1);
    }
  }
  if (0<alignment.size()) {
    stId = alignment[i];
    stateTransCount_[stId][1] -= 1.0;
  }
  return alignment.size();
}

LogVectorFst DiscreteHMM::frameStateWfst(const vector<int> &feaSeq) {
  LogVectorFst fst;
  unsigned frame;
  int hmmStId;
  double prob;

  fst.AddState();
  fst.SetStart(0);
  for (frame=0; frame < feaSeq.size(); frame++) {
    fst.AddState();
  }
  fst.SetFinal(frame,0);
  for (frame=0; frame < feaSeq.size(); frame++) {
    for (hmmStId=0; hmmStId<numUnits_*numUnitStates_; hmmStId++) {
      prob=stateProb(hmmStId, feaSeq[frame]);
      // shift hmmStId by 1 since 0 is used for eps
      fst.AddArc(frame,LogArc(hmmStId+1,hmmStId+1,LogWeight(-1*log(prob)),frame+1));
      // fst.AddArc(frame,LogArc(hmmStId+1,hmmStId+1,-1*log(prob),frame+1));
    }
  }
  return fst;
}

LogVectorFst DiscreteHMM::monoPhoneWfst() {
  LogVectorFst fst;
  int utId, hmmStId, lDiscreteHMMStId, wfstStId;
  double prob;

  fst.AddState();
  fst.SetStart(0);
  fst.SetFinal(0,0);
  wfstStId = 0;
  for (utId=0; utId<numUnits_; utId++) {
    fst.AddState();
    wfstStId++;
    hmmStId = utId*numUnitStates_;
    // shift utId (phone name) by 2 since 0 is used for eps and 1 is used for phi
    // shift hmmStId (HMM state ID) by 1 since 0 is used for eps
    fst.AddArc(0,LogArc(hmmStId+1,utId+2,LogWeight(0),wfstStId));
    for (lDiscreteHMMStId=0; lDiscreteHMMStId<numUnitStates_-1; lDiscreteHMMStId++) {
      fst.AddState();
      prob=transProb(hmmStId, 0); // self-loop
      fst.AddArc(wfstStId,LogArc(hmmStId+1,0,LogWeight(-1*log(prob)),wfstStId));
      prob=transProb(hmmStId, 1); // next state
      fst.AddArc(wfstStId,LogArc(hmmStId+1+1,0,LogWeight(-1*log(prob)),wfstStId+1)); 
      hmmStId++;
      wfstStId++;
    }
    prob=transProb(hmmStId, 1);  
    fst.AddArc(wfstStId,LogArc(0,0,LogWeight(-1*log(prob)),0)); // next state (return to home)
    prob=transProb(hmmStId, 0);  
    fst.AddArc(wfstStId,LogArc(hmmStId+1,0,LogWeight(-1*log(prob)),wfstStId)); // self-loop
  }
  return fst;
}

// Returns a sequence of HMM state ID (i.e, input labels of sample)
vector<int> DiscreteHMM::parseSample(const Fst<LogArc> & sample) { 
  int hmmStId;
  LogArc arc;

  vector<int> ret;
  Fst<LogArc>::StateId sid = sample.Start();
  // continue until there are no more left
  while(true) {
    ArcIterator< Fst<LogArc> > ai(sample,sid);
    if(ai.Done())
      break;
    arc = ai.Value();
    hmmStId = arc.ilabel - 1;
    if (0<=hmmStId) { // exclude eps
      ret.push_back(hmmStId);
    }
    sid = arc.nextstate;
  }
  return ret;
}

vector<double> DiscreteHMM::dirichletSample(vector<double> &alpha, vector<double> *diffAlpha=NULL) {
  unsigned i;
  double sum;
  vector<double> ret;

  sum = 0.0;
  ret.resize(alpha.size());
  for (i=0; i<ret.size(); i++) {
    ret[i] = gammaSample(alpha[i], 1.0);
    sum += ret[i];
  }
  for (i=0; i<ret.size(); i++) {
    ret[i] /= sum;
  }
  return ret;
}

double DiscreteHMM::gammaSample(double a, double scale) {
  double b, c, e, u, v, w, y, x, z;
  if(a > 1) { // Best's XG method
    b = a-1;
    c = 3*a-.75;
    bool accept = false;
    do {
      u = (double)rand()/RAND_MAX;
      v = (double)rand()/RAND_MAX;
      w = u*(1-u);
      y = sqrt(c/w)*(u-.5);
      x = b+y;
      if(x >= 0) {
	z = 64*w*w*w*v*v;
	accept = (z <= 1-2*y*y/x || log(z) <= 2*(b*log(x/b)-y));
      }
    } while (!accept);
  } else { // Johnk's method
    do {
      u = (double)rand()/RAND_MAX;
      v = (double)rand()/RAND_MAX;
      x = pow(u,1/a);
      y = pow(v,1/(1-a));
    } while (x+y > 1);
    e = exponSample(1);
    x = e*x/(x+y);
  }
  return x * scale;
}

double DiscreteHMM::exponSample(double l) {
  return -1*log(1-(double)rand()/RAND_MAX)/l;
}
