/**
 *  @copyright Copyright 2020 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  @file HitFinderTools.cpp
 */

using namespace std;

#include "HitFinderTools.h"
#include <TMath.h>
#include <vector>
#include <cmath>
#include <map>

/**
 * Helper method for sotring signals in vector
 */
void HitFinderTools::sortByTime(vector<JPetMatrixSignal>& sigVec)
{
  sort(sigVec.begin(), sigVec.end(),
    [](const JPetMatrixSignal & sig1, const JPetMatrixSignal & sig2) {
      return sig1.getTime() < sig2.getTime();
    }
  );
 }

/**
 * Method distributing Signals according to Scintillator they belong to
 */
map<string, map<int, vector<JPetMatrixSignal>>> HitFinderTools::getMappedSignals(const JPetTimeWindow* timeWindow)
{
  map<string, map<int, vector<JPetMatrixSignal>>> signalSidesMap;

  if (!timeWindow) {
    WARNING("Pointer of Time Window object is not set, returning empty map");
    return signalSidesMap;
  }

  map<int, vector<JPetMatrixSignal>> scinSigMap;
  map<int, vector<JPetMatrixSignal>> wlsSigMap;
  vector<JPetMatrixSignal> wlsSigVec;

  const unsigned int nSignals = timeWindow->getNumberOfEvents();
  for (unsigned int i = 0; i < nSignals; i++) {
    auto mtxSig = dynamic_cast<const JPetMatrixSignal&>(timeWindow->operator[](i));

    if(mtxSig.getMatrix().getType()=="SideA" || mtxSig.getMatrix().getType()=="SideB") {
      auto scinID = mtxSig.getMatrix().getScin().getID();
      auto search = scinSigMap.find(scinID);
      if(search==scinSigMap.end()){
        vector<JPetMatrixSignal> tmp;
        tmp.push_back(mtxSig);
        scinSigMap[scinID] = tmp;
      } else {
        search->second.push_back(mtxSig);
      }
    } else if(mtxSig.getMatrix().getType()=="WLS") {
      wlsSigVec.push_back(mtxSig);
    }
  }
  wlsSigMap[0] = wlsSigVec;
  signalSidesMap["Scin"] = scinSigMap;
  signalSidesMap["WLS"] = wlsSigMap;
  return signalSidesMap;
}

/**
 * Loop over all Scins invoking matching procedure
 */
vector<JPetHit> HitFinderTools::matchAllSignals(
  map<string, map<int, vector<JPetMatrixSignal>>>& signalSidesMap,
  double minTimeDiffAB, double maxTimeDiffAB, JPetStatistics& stats, bool saveHistos
) {
  vector<JPetHit> allHits;

  // Sorting wls signals in time
  auto wlsSignals = signalSidesMap.at("WLS").at(0);

  for (auto& scinSigMap : signalSidesMap.at("Scin")) {
    // Matching A-B sides signals for same scin
    auto hits = matchSignals(
      scinSigMap.second, wlsSignals, minTimeDiffAB, maxTimeDiffAB, stats, saveHistos
    );
    allHits.insert(allHits.end(), hits.begin(), hits.end());
  }
  if (saveHistos) {
    stats.getHisto1D("hits_tslot")->Fill(allHits.size());
  }
  return allHits;
}

/**
 * Method matching signals on the same Scintillator
 */
vector<JPetHit> HitFinderTools::matchSignals(
  vector<JPetMatrixSignal>& scinSignals, vector<JPetMatrixSignal>& wlsSignals,
  double minTimeDiffAB, double maxTimeDiffAB, JPetStatistics& stats, bool saveHistos
) {
  vector<JPetHit> hits;

  sortByTime(scinSignals);
  sortByTime(wlsSignals);

  // Matching signals on sides
  for (unsigned int i = 0; i < scinSignals.size(); i++) {
    if(i>=scinSignals.size()) { break; }
    for (unsigned int j = i+1; j < scinSignals.size(); j++) {
      if(j>=scinSignals.size()) { break; }

      auto tDiff = scinSignals.at(j).getTime() - scinSignals.at(i).getTime();
      if(saveHistos) { stats.getHisto1D("ab_tdiff_all")->Fill(tDiff); }

      // Time condition
      if (tDiff > minTimeDiffAB && tDiff < maxTimeDiffAB) {

        // Different sides condition
        if (scinSignals.at(i).getMatrix().getType() == scinSignals.at(j).getMatrix().getType()) { continue; }

        // Found A-B signals in conincidence
        if(saveHistos) { stats.getHisto1D("ab_tdiff_acc")->Fill(tDiff); }

        // If layer 1 (black module) - create AB Hit
        auto layerID = scinSignals.at(i).getMatrix().getScin().getSlot().getLayer().getID();
        if(layerID == 1){
          hits.push_back(createHit(scinSignals.at(i), scinSignals.at(j)));
        }

        // Layer 2 or 4 - Red Module, searching for WLS signals
        if(layerID == 2 || layerID == 4) {
          auto hitTime = (scinSignals.at(i).getTime()+scinSignals.at(j).getTime())/2.0;
          auto wlsSigIndex = matchWLSSignal(wlsSignals, hitTime, minTimeDiffAB, maxTimeDiffAB, stats, saveHistos);
          // Singal found if index different than -1
          if(wlsSigIndex==-1){
            hits.push_back(createHit(scinSignals.at(i), scinSignals.at(j)));
          } else {
            hits.push_back(createHit(scinSignals.at(i), scinSignals.at(j), wlsSignals.at(wlsSigIndex)));
            wlsSignals.erase(wlsSignals.begin() + wlsSigIndex);
          }
        }
      } else {
        if(saveHistos) { stats.getHisto1D("ab_tdiff_rej")->Fill(tDiff); }
        i = j;
        break;
      }
    }
  }
  return hits;
}

/**
 * Checking times of WLS signals if they mach with hit time in conincidence window
 */
int HitFinderTools::matchWLSSignal(
  std::vector<JPetMatrixSignal>& wlsSignals, double hitTime,
  double minTimeDiffAB, double maxTimeDiffAB, JPetStatistics& stats, bool saveHistos
){
  for(unsigned int i = 0; i < wlsSignals.size(); i++){
    auto tDiff = fabs(hitTime-wlsSignals.at(i).getTime());
    if(saveHistos) { stats.getHisto1D("hit_wls_tdiff_all")->Fill(tDiff); }
    if(tDiff > minTimeDiffAB && tDiff < maxTimeDiffAB) {
      if(saveHistos) { stats.getHisto1D("hit_wls_tdiff_acc")->Fill(tDiff); }
      return i;
    } else {
      if(saveHistos) { stats.getHisto1D("hit_wls_tdiff_rej")->Fill(tDiff); }
    }
  }
  return -1;
}

/**
 * Method for Hit creation with A-B signals
 */
JPetHit HitFinderTools::createHit(
  const JPetMatrixSignal& signal1, const JPetMatrixSignal& signal2
) {
  JPetMatrixSignal signalA;
  JPetMatrixSignal signalB;

  if (signal1.getMatrix().getType() == "SideA") {
    signalA = signal1;
    signalB = signal2;
  } else if(signal1.getMatrix().getType() == "SideB") {
    signalA = signal2;
    signalB = signal1;
  }

  JPetHit hit;
  hit.setSignalA(signalA);
  hit.setSignalB(signalB);
  hit.setTime((signalA.getTime() + signalB.getTime()) / 2.0);
  hit.setQualityOfTime(-1.0);
  hit.setTimeDiff(signalB.getTime() - signalA.getTime());
  hit.setQualityOfTimeDiff(-1.0);
  // TOT is a sum of over all threshold in all signals on both sides
  // As an quality of energy we temporaily put multiplicity of signals (2-8)
  hit.setEnergy(signalA.getTOT() + signalB.getTOT());
  hit.setQualityOfEnergy(signalA.getRawSignals().size()+signalB.getRawSignals().size());
  hit.setPosX(signalA.getMatrix().getScin().getCenterX());
  hit.setPosY(signalA.getMatrix().getScin().getCenterY());
  hit.setPosZ(-99.0);

  hit.setScin(signalA.getMatrix().getScin());
  hit.setWLS(JPetWLS::getDummyResult());

  return hit;
}

/**
 * Method for Hit creation with A-B signals and WLS signal for z position estimation
 */
JPetHit HitFinderTools::createHit(
  const JPetMatrixSignal& signal1, const JPetMatrixSignal& signal2, const JPetMatrixSignal& signalWLS
) {
  JPetMatrixSignal signalA;
  JPetMatrixSignal signalB;

  if (signal1.getMatrix().getType() == "SideA") {
    signalA = signal1;
    signalB = signal2;
  } else if(signal1.getMatrix().getType() == "SideB") {
    signalA = signal2;
    signalB = signal1;
  }

  JPetHit hit;
  hit.setSignalA(signalA);
  hit.setSignalB(signalB);
  hit.setSignalWLS(signalWLS);

  hit.setTime((signalA.getTime() + signalB.getTime()) / 2.0);
  hit.setQualityOfTime(-1.0);
  hit.setTimeDiff(signalB.getTime() - signalA.getTime());
  hit.setQualityOfTimeDiff(-1.0);

  // TOT is a sum of over all threshold in all signals on both sides
  // As an quality of energy we temporaily put multiplicity of signals (2-8)
  hit.setEnergy(signalA.getTOT() + signalB.getTOT());
  hit.setQualityOfEnergy(signalA.getRawSignals().size()+signalB.getRawSignals().size());

  hit.setPosX(signalA.getMatrix().getScin().getCenterX());
  hit.setPosY(signalA.getMatrix().getScin().getCenterY());

  auto wlsTOT = signalWLS.getTOT();
  auto rawSignals = signalWLS.getRawSignals();
  double sumWeights = 0.0;
  double sumPositions = 0.0;
  for(auto rawSignal : rawSignals){
    auto partToT = rawSignal.second.getTOT();
    auto zPos = rawSignal.second.getPM().getPosition();
    sumPositions += zPos*partToT/wlsTOT;
    sumWeights += partToT/wlsTOT;
  }
  if(sumWeights!=0.0){
    hit.setPosZ(sumPositions/sumWeights);
  } else {
    hit.setPosZ(-99.0);
  }

  hit.setScin(signalA.getMatrix().getScin());
  hit.setWLS(signalWLS.getMatrix().getWLS());

  return hit;
}


/**
 * Method for Hit creation in case of reference detector.
 * Setting only some necessary fields.
 */
JPetHit HitFinderTools::createDummyHit(const JPetMatrixSignal& signal)
{
  JPetHit hit;
  JPetMatrixSignal dummy;
  hit.setSignalA(dummy);
  hit.setSignalB(signal);
  hit.setTime(signal.getTime());
  hit.setQualityOfTime(-1.0);
  hit.setTimeDiff(0.0);
  hit.setQualityOfTimeDiff(-1.0);
  hit.setEnergy(signal.getTOT());
  hit.setQualityOfEnergy(-1.0);
  hit.setPosX(-99.0);
  hit.setPosY(-99.0);
  hit.setPosZ(-99.0);
  hit.setScin(JPetScin::getDummyResult());
  hit.setWLS(JPetWLS::getDummyResult());
  return hit;
}

/**
* Calculation of the total TOT of the hit - Time over Threshold:
* the sum of the TOTs on both thresholds and on the both sides (A,B)
*/
double HitFinderTools::calculateTOT(JPetHit& hit)
{
  double tot = 0.0;

  auto rawSignalsA = dynamic_cast<const JPetMatrixSignal&>(hit.getSignalA()).getRawSignals();
  auto rawSignalsB = dynamic_cast<const JPetMatrixSignal&>(hit.getSignalB()).getRawSignals();

  for(auto rawSig: rawSignalsA){
    tot += rawSig.second.getTOT();
  }

  for(auto rawSig: rawSignalsB){
    tot += rawSig.second.getTOT();
  }

  return tot;
}
