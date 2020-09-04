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
vector<JPetMatrixSignal> SignalTransformerTools::mergeScinSiPMSignals(
  map<int, vector<vector<JPetRawSignal>>>& rawSigMtxMap, double mergingTime,
  boost::property_tree::ptree& scinSync, JPetStatistics& stats, bool saveHistos
) {
  vector<JPetMatrixSignal> allMtxSignals;
  // Iterating over whole map
  for (auto& rawSigScin : rawSigMtxMap) {
    if(rawSigScin.second.size()==0) { continue; }
    if(rawSigScin.second.at(0).size()==0) { continue; }
    auto scinID = rawSigScin.second.at(0).at(0).getPM().getScin().getID();
    // Getting offsets for this scintillator -
    // if calibrations are empty then default vaule is 0.0
    double offset = scinSync.get("scin_offsets."+to_string(scinID), 0.0);

    for (auto& rawSigSide : rawSigScin.second){
      auto mtxSignals = mergeRawSignalsOnSide(rawSigSide, mergingTime, offset);
      allMtxSignals.insert(allMtxSignals.end(), mtxSignals.begin(), mtxSignals.end());
    }
  }

  for (auto& mtxSig : allMtxSignals) {

    if(saveHistos){
      if(mtxSig.getPM().getDesc()=="scin") {
        auto scinID = mtxSig.getPM().getScin().getID();
        // if(scinID<fMinScinID || scinID>fMaxScinID) { continue; }

        stats.getHisto1D("mtxsig_multi")->Fill(mtxSig.getRawSignals().size());
        if(mtxSig.getPM().getSide()==JPetPM::SideA){
          stats.getHisto1D("mtxsig_per_scin_sideA")->Fill(scinID);
        } else if(mtxSig.getPM().getSide()==JPetPM::SideB){
          stats.getHisto1D("mtxsig_per_scin_sideB")->Fill(scinID);
        }
        auto rawSigVec = mtxSig.getRawSignals();
        if(rawSigVec.size() == 4) {
          auto side = rawSigVec.at(1).getPM().getSide();
          auto pm1ID = rawSigVec.at(1).getPM().getID();
          auto pm2ID = rawSigVec.at(2).getPM().getID();
          auto pm3ID = rawSigVec.at(3).getPM().getID();
          auto pm4ID = rawSigVec.at(4).getPM().getID();
          auto t1 = SignalTransformerTools::getRawSigBaseTime(rawSigVec.at(1));
          auto t2 = SignalTransformerTools::getRawSigBaseTime(rawSigVec.at(2));
          auto t3 = SignalTransformerTools::getRawSigBaseTime(rawSigVec.at(3));
          auto t4 = SignalTransformerTools::getRawSigBaseTime(rawSigVec.at(4));
          stats.getHisto1D(Form("offset_sipm_%d", pm2ID))->Fill(t2-t1);
          stats.getHisto1D(Form("offset_sipm_%d", pm3ID))->Fill(t3-t1);
          stats.getHisto1D(Form("offset_sipm_%d", pm4ID))->Fill(t4-t1);
        }
      }
    }
  }
  return allMtxSignals;
}

/**
 * Method iterates over all WLSs and creates vestor of signals, that are assigned to it.
 * For each created vector, the merging method is called
 */
vector<JPetMatrixSignal> SignalTransformerTools::mergeWLSSiPMSignals(
  const map<int, JPetWLS*>& wlsMap, map<int, vector<JPetRawSignal>>& rawSigWLSMap,
  double mergingTime, JPetStatistics& stats, bool saveHistos
) {
  vector<JPetMatrixSignal> allWLSSignals;

  for(auto& wlsEle : wlsMap) {
    auto wlsID = wlsEle.second->getID();
    auto pmIDs = wlsEle.second->getPMIDs();

    // cout << "WLS ID: " << wlsID << " sipms: " << pmIDs.size() << " ids: ";
    // for(auto id : pmIDs) { cout << id; }
    // cout << endl;

    vector<JPetRawSignal> wlsSignals;
    for(auto& pmEle : rawSigWLSMap) {
      for(auto id : pmIDs) {
        if(pmEle.first == id) {
           wlsSignals.insert(wlsSignals.end(), pmEle.second.begin(), pmEle.second.end());
        }
      }
    }

    auto matchedWLSSig = SignalTransformerTools::mergeRawSignalsOnSide(wlsSignals, mergingTime, 0.0);
    allWLSSignals.insert(allWLSSignals.end(), matchedWLSSig.begin(), matchedWLSSig.end());

    if(saveHistos){
      for(auto mtxSig : matchedWLSSig) {
        stats.getHisto1D(Form("wls_%d_tot", wlsID))->Fill(mtxSig.getTOT());

        for(auto sigEle : mtxSig.getRawSignals()){
          auto pmID = sigEle.second.getPM().getID();
          auto leads = sigEle.second.getPoints(JPetSigCh::Leading, JPetRawSignal::ByThrValue);
          auto trails = sigEle.second.getPoints(JPetSigCh::Trailing, JPetRawSignal::ByThrValue);
          double tot = trails.at(0).getTime()-leads.at(0).getTime();
          stats.getHisto1D(Form("wls_%d_sipm_%d_tot", wlsID, pmID))->Fill(tot);
        }
      }
    }
  }

  return allWLSSignals;
}

/**
 * Method iterates over all Raw Sigals on some SiPMs on the same matrix,
 * matching them into groups on max. 4 as a MatrixSignal
 */
vector<JPetMatrixSignal> SignalTransformerTools::mergeRawSignalsOnSide(
  vector<JPetRawSignal>& rawSigVec, double mergingTime, double offset
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
    mtxSig.setTime(calculateAverageTime(mtxSig)-offset);
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
