/**
 *  @copyright Copyright 2020 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file SignalTransformerTools.cpp
 */

#include "SignalTransformerTools.h"

using namespace std;

/**
 * Map Raw Signals from all SiPMs according to matrix and side they belong to.
 * Returns map< scin ID, side < signals >>.
 * Side A is the first element int he vector, Side B is the second one.
 */
const pair<map<int, vector<vector<JPetRawSignal>>>, map<int, vector<JPetRawSignal>>>
SignalTransformerTools::getRawSigMtxWLSMap(const JPetTimeWindow* timeWindow)
{
  map<int, vector<vector<JPetRawSignal>>> rawSigMtxMap;
  map<int, vector<JPetRawSignal>> rawSigWLSMap;

  if (!timeWindow) {
    WARNING("Pointer of Time Window object is not set, returning empty map");
    return make_pair(rawSigMtxMap, rawSigWLSMap);
  }

  const unsigned int nRawSigs = timeWindow->getNumberOfEvents();
  for (unsigned int i = 0; i < nRawSigs; i++) {
    auto rawSig = dynamic_cast<const JPetRawSignal&>(timeWindow->operator[](i));

    if(rawSig.getPM().getDesc()=="scin") {
      auto scinID = rawSig.getPM().getScin().getID();
      auto pmSide = rawSig.getPM().getSide();
      auto search = rawSigMtxMap.find(scinID);

      if (search == rawSigMtxMap.end()) {
        // There is no element with searched scin ID in this map, adding new one
        vector<JPetRawSignal> tmpVec;
        vector<vector<JPetRawSignal>> tmpVecVec;
        tmpVecVec.push_back(tmpVec);
        tmpVecVec.push_back(tmpVec);
        if(pmSide==JPetPM::SideA){
          tmpVecVec.at(0).push_back(rawSig);
        } else if(pmSide==JPetPM::SideB){
          tmpVecVec.at(1).push_back(rawSig);
        }
        rawSigMtxMap.insert(pair<int, vector<vector<JPetRawSignal>>>(scinID, tmpVecVec));
      } else {
        if(pmSide==JPetPM::SideA){
          search->second.at(0).push_back(rawSig);
        }else if(pmSide==JPetPM::SideB){
          search->second.at(1).push_back(rawSig);
        }
      }
    } else if(rawSig.getPM().getDesc()=="wls") {
      auto pmID = rawSig.getPM().getID();
      auto search2 = rawSigWLSMap.find(pmID);
      if (search2 == rawSigWLSMap.end()) {
        vector<JPetRawSignal> tmpVec;
        tmpVec.push_back(rawSig);
      } else {
        search2->second.push_back(rawSig);
      }
    }
  }
  return make_pair(rawSigMtxMap, rawSigWLSMap);
}

/**
 * Method iterates over all matrices in the detector with signals,
 * calling merging procedure for each
 */
vector<JPetMatrixSignal> SignalTransformerTools::mergeSignalsAllSiPMs(
   map<int, vector<vector<JPetRawSignal>>>& rawSigMtxMap, double mergingTime
) {
  vector<JPetMatrixSignal> allMtxSignals;
  // Iterating over whole map
  for (auto& rawSigScin : rawSigMtxMap) {
    for (auto& rawSigSide : rawSigScin.second){
      auto mtxSignals = mergeRawSignalsOnSide(
        rawSigSide, mergingTime
      );
      allMtxSignals.insert(allMtxSignals.end(), mtxSignals.begin(), mtxSignals.end());
    }
  }
  return allMtxSignals;
}

/**
 * Method iterates over all Raw Sigals on some SiPMs on the same matrix,
 * matching them into groups on max. 4 as a MatrixSignal
 */
vector<JPetMatrixSignal> SignalTransformerTools::mergeRawSignalsOnSide(
  vector<JPetRawSignal>& rawSigVec, double mergingTime
) {
  vector<JPetMatrixSignal> mtxSigVec;
  sortByTime(rawSigVec);

  while (rawSigVec.size() > 0) {
    // Create Matrix Signal and add first Raw Signal by default
    JPetMatrixSignal mtxSig;
    mtxSig.setPM(rawSigVec.at(0).getPM());
    if(!mtxSig.addRawSignal(rawSigVec.at(0))){
      ERROR("Problem with adding the first signal to new object.");
      break;
    }

    unsigned int nextIndex = 1;
    while(true){

      if(rawSigVec.size() <= nextIndex){
        // nothing left to check
        break;
      }

      // signal matching condidion
      if(fabs(getRawSigBaseTime(rawSigVec.at(nextIndex))
        -getRawSigBaseTime(rawSigVec.at(0))) < mergingTime
      ){
        // mathing signal found
        if(mtxSig.addRawSignal(rawSigVec.at(nextIndex))){
          // added succesfully
          rawSigVec.erase(rawSigVec.begin()+nextIndex);
        } else {
          // this mtx pos is already occupied, check the next one
          nextIndex++;
        }
      } else {
        // next signal is too far from reference one, this MtxSig is finished
        break;
      }
    }
    mtxSig.setTime(calculateAverageTime(mtxSig));
    rawSigVec.erase(rawSigVec.begin());
    mtxSigVec.push_back(mtxSig);
  }
  return mtxSigVec;
}

/**
 * Returning time of leading Signal Channel on the first threshold from Raw Signal
 */
double SignalTransformerTools::getRawSigBaseTime(JPetRawSignal& rawSig)
{
  return rawSig.getPoints(JPetSigCh::Leading, JPetRawSignal::ByThrValue).at(0).getTime();
}

/**
 * Calculating average time of Matrix Signal based on times of contained Raw Signals
 */
double SignalTransformerTools::calculateAverageTime(JPetMatrixSignal& mtxSig)
{
  double averageTime = 0.0;
  auto rawSignals = mtxSig.getRawSignals();
  for(auto rawSig : rawSignals){
    averageTime += getRawSigBaseTime(rawSig.second);
  }
  averageTime = averageTime/((double) rawSignals.size());
  return averageTime;
}

/**
 * Sorting method for Raw Signals, based on time of leading THR1 Signal Channel
 */
void SignalTransformerTools::sortByTime(vector<JPetRawSignal>& input)
{
  std::sort(
    input.begin(), input.end(), [] (JPetRawSignal rawSig1, JPetRawSignal rawSig2) {
      return getRawSigBaseTime(rawSig1) < getRawSigBaseTime(rawSig2);
    }
  );
}
