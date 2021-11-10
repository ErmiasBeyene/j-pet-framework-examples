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
#include <cmath>
#include <map>
#include <vector>

/**
 * Helper method for sotring signals in vector
 */
void HitFinderTools::sortByTime(vector<JPetMatrixSignal>& sigVec)
{
  sort(sigVec.begin(), sigVec.end(), [](const JPetMatrixSignal& sig1, const JPetMatrixSignal& sig2) { return sig1.getTime() < sig2.getTime(); });
}

/**
 * Method distributing Signals according to Scintillator they belong to
 */
map<string, map<int, vector<JPetMatrixSignal>>> HitFinderTools::getMappedSignals(const JPetTimeWindow* timeWindow)
{
  map<string, map<int, vector<JPetMatrixSignal>>> signalSidesMap;

  if (!timeWindow)
  {
    WARNING("Pointer of Time Window object is not set, returning empty map");
    return signalSidesMap;
  }

  map<int, vector<JPetMatrixSignal>> scinSigMap;
  map<int, vector<JPetMatrixSignal>> wlsSigMap;
  vector<JPetMatrixSignal> wlsSigVec;

  const unsigned int nSignals = timeWindow->getNumberOfEvents();
  for (unsigned int i = 0; i < nSignals; i++)
  {
    auto mtxSig = dynamic_cast<const JPetMatrixSignal&>(timeWindow->operator[](i));

    if (mtxSig.getMatrix().getType() == "SideA" || mtxSig.getMatrix().getType() == "SideB")
    {
      auto scinID = 1;
      auto search = scinSigMap.find(scinID);
      if (search == scinSigMap.end())
      {
        vector<JPetMatrixSignal> tmp;
        tmp.push_back(mtxSig);
        scinSigMap[scinID] = tmp;
      }
      else
      {
        search->second.push_back(mtxSig);
      }
    }
    else if (mtxSig.getMatrix().getType() == "WLS")
    {
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
vector<JPetHit> HitFinderTools::matchAllSignals(map<string, map<int, vector<JPetMatrixSignal>>>& signalSidesMap, double minTimeDiffAB,
                                                double maxTimeDiffAB, JPetStatistics& stats, bool saveHistos)
{
  vector<JPetHit> allHits;

  // Creating hits from SiPMs attached to WLS
  for (auto signal : signalSidesMap.at("WLS").at(0))
  {
    auto wlsHit = createWLSHit(signal);
    allHits.push_back(wlsHit);
    if (saveHistos)
    {
      stats.getHisto1D("hit_pos_z_wls")->Fill(wlsHit.getPosZ());
    }
  }

  auto scinSignals = signalSidesMap.at("Scin");

  // Standard Side A-B signal matching
  for (auto& signals : scinSignals)
  {
    // Match signals for scintillators
    auto scinHits = matchSignals(signals.second, minTimeDiffAB, maxTimeDiffAB, stats, saveHistos);
    allHits.insert(allHits.end(), scinHits.begin(), scinHits.end());
  }

  // Sorting wls signals in time
  // auto wlsSignals = signalSidesMap.at("WLS").at(0);
  //
  // for (auto& scinSigMap : signalSidesMap.at("Scin")) {
  //   // Matching A-B sides signals for same scin
  //   auto hits = matchSignals(
  //     scinSigMap.second, wlsSignals, minTimeDiffAB, maxTimeDiffAB, stats, saveHistos
  //   );
  //   allHits.insert(allHits.end(), hits.begin(), hits.end());
  // }
  // if (saveHistos) {
  //   stats.getHisto1D("hits_tslot")->Fill(allHits.size());
  // }

  return allHits;
}

vector<JPetHit> HitFinderTools::matchSignals(vector<JPetMatrixSignal>& scinSignals, double minTimeDiffAB, double maxTimeDiffAB, JPetStatistics& stats,
                                             bool saveHistos)
{
  vector<JPetHit> scinHits;
  vector<JPetMatrixSignal> remainSignals;
  sortByTime(scinSignals);

  while (scinSignals.size() > 0)
  {
    auto mtxSig = scinSignals.at(0);
    if (scinSignals.size() == 1)
    {
      remainSignals.push_back(mtxSig);
      break;
    }

    for (unsigned int j = 1; j < scinSignals.size(); j++)
    {
      if (minTimeDiffAB < scinSignals.at(j).getTime() - mtxSig.getTime() && scinSignals.at(j).getTime() - mtxSig.getTime() < maxTimeDiffAB)
      {
        if (mtxSig.getMatrix().getType() != scinSignals.at(j).getMatrix().getType())
          if ((mtxSig.getMatrix().getType() == "SideA" && scinSignals.at(j).getMatrix().getType() == "SideB") ||
              (mtxSig.getMatrix().getType() == "SideB" && scinSignals.at(j).getMatrix().getType() == "SideA"))
          {
            auto hit = createScinHit(mtxSig, scinSignals.at(j));

            if (saveHistos)
            {
              stats.getHisto2D("hit_pos_XY")->Fill(hit.getPosX(), hit.getPosY());
              stats.getHisto1D("hit_pos_z")->Fill(hit.getPosZ());
              stats.getHisto1D("hit_tdiff")->Fill(hit.getTimeDiff());
              stats.getHisto2D("time_diff_per_scin")->Fill(hit.getTimeDiff(), hit.getScin().getID());
              stats.getHisto1D("hit_per_scin")->Fill(mtxSig.getMatrix().getScin().getID());
              stats.getHisto1D("hit_per_scin")->Fill(scinSignals.at(j).getMatrix().getScin().getID());
              stats.getHisto1D("hit_sig_multi")->Fill(hit.getQualityOfEnergy());
            }

            scinHits.push_back(hit);
            scinSignals.erase(scinSignals.begin() + j);
            scinSignals.erase(scinSignals.begin() + 0);
            break;
          }
          else
          {
            if (j == scinSignals.size() - 1)
            {
              remainSignals.push_back(mtxSig);
              scinSignals.erase(scinSignals.begin() + 0);
              break;
            }
            else
            {
              continue;
            }
          }
      }
      else
      {
        if (saveHistos && mtxSig.getMatrix().getType() != scinSignals.at(j).getMatrix().getType())
        {
          stats.getHisto1D("remain_signals_tdiff")->Fill(scinSignals.at(j).getTime() - mtxSig.getTime());
        }
        remainSignals.push_back(mtxSig);
        scinSignals.erase(scinSignals.begin() + 0);
        break;
      }
    }
  }
  if (remainSignals.size() > 0 && saveHistos)
  {
    stats.getHisto1D("remain_signals_scin")->Fill((double)(remainSignals.at(0).getMatrix().getScin().getID()), remainSignals.size());
  }
  return scinHits;
}

/**
 * Method for Hit creation with A-B signals
 */
JPetHit HitFinderTools::createScinHit(const JPetMatrixSignal& signal1, const JPetMatrixSignal& signal2)
{
  JPetMatrixSignal signalA;
  JPetMatrixSignal signalB;

  if (signal1.getMatrix().getType() == "SideA")
  {
    signalA = signal1;
    signalB = signal2;
  }
  else if (signal1.getMatrix().getType() == "SideB")
  {
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
  hit.setQualityOfEnergy(signalA.getRawSignals().size() + signalB.getRawSignals().size());
  hit.setPosX(signalA.getMatrix().getScin().getCenterX());
  hit.setPosY(signalA.getMatrix().getScin().getCenterY());
  // Hardcoded velocity 12 cm/ns = 0.012 cm/ps
  double velocity = 0.012;
  hit.setPosZ(velocity * hit.getTimeDiff() / 2.0);

  hit.setScin(signalA.getMatrix().getScin());
  hit.setWLS(JPetWLS::getDummyResult());

  return hit;
}

JPetHit HitFinderTools::createWLSHit(const JPetMatrixSignal& signalWLS)
{
  JPetHit hit;
  hit.setSignalWLS(signalWLS);
  hit.setTime(signalWLS.getTime());
  hit.setQualityOfTime(-1.0);
  hit.setTimeDiff(-1.0);
  hit.setQualityOfTimeDiff(-1.0);
  hit.setEnergy(signalWLS.getTOT());
  hit.setQualityOfEnergy(signalWLS.getRawSignals().size());
  hit.setPosX(signalWLS.getMatrix().getWLS().getCenterX());
  hit.setPosY(signalWLS.getMatrix().getWLS().getCenterY());
  hit.setPosZ(signalWLS.getMatrix().getWLS().getCenterZ());
  hit.setScin(JPetScin::getDummyResult());
  hit.setWLS(signalWLS.getMatrix().getWLS());
  return hit;
}

/**
 * Method for Hit creation with A-B signals and WLS signal for z position estimation
 */
JPetHit HitFinderTools::createHit(const JPetMatrixSignal& signal1, const JPetMatrixSignal& signal2, const JPetMatrixSignal& signalWLS)
{
  JPetMatrixSignal signalA;
  JPetMatrixSignal signalB;

  if (signal1.getMatrix().getType() == "SideA")
  {
    signalA = signal1;
    signalB = signal2;
  }
  else if (signal1.getMatrix().getType() == "SideB")
  {
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
  hit.setQualityOfEnergy(signalA.getRawSignals().size() + signalB.getRawSignals().size());

  hit.setPosX(signalA.getMatrix().getScin().getCenterX());
  hit.setPosY(signalA.getMatrix().getScin().getCenterY());

  auto wlsTOT = signalWLS.getTOT();
  auto rawSignals = signalWLS.getRawSignals();
  double sumWeights = 0.0;
  double sumPositions = 0.0;
  for (auto rawSignal : rawSignals)
  {
    auto partToT = rawSignal.second.getTOT();
    auto zPos = rawSignal.second.getPM().getPosition();
    sumPositions += zPos * partToT / wlsTOT;
    sumWeights += partToT / wlsTOT;
  }
  if (sumWeights != 0.0)
  {
    hit.setPosZ(sumPositions / sumWeights);
  }
  else
  {
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

  for (auto rawSig : rawSignalsA)
  {
    tot += rawSig.second.getTOT();
  }

  for (auto rawSig : rawSignalsB)
  {
    tot += rawSig.second.getTOT();
  }

  return tot;
}

/**
 * Method matching signals on the same Scintillator
 */
// vector<JPetHit> HitFinderTools::matchSignals(vector<JPetMatrixSignal>& scinSignals, vector<JPetMatrixSignal>& wlsSignals, double minTimeDiffAB,
//                                              double maxTimeDiffAB, JPetStatistics& stats, bool saveHistos)
// {
//   vector<JPetHit> hits;
//
//   sortByTime(scinSignals);
//   sortByTime(wlsSignals);
//
//   // Matching signals on sides of a scintillator
//   for (unsigned int i = 0; i < scinSignals.size(); i++)
//   {
//     if (i >= scinSignals.size())
//     {
//       break;
//     }
//     for (unsigned int j = i + 1; j < scinSignals.size(); j++)
//     {
//       if (j >= scinSignals.size())
//       {
//         break;
//       }
//
//       // Different sides condition
//       if (scinSignals.at(i).getMatrix().getType() == scinSignals.at(j).getMatrix().getType())
//       {
//         continue;
//       }
//
//       auto tDiff = scinSignals.at(j).getTime() - scinSignals.at(i).getTime();
//       if (saveHistos)
//       {
//         stats.getHisto1D("ab_tdiff_all")->Fill(tDiff);
//       }
//
//       // Time condition
//       if (tDiff > minTimeDiffAB && tDiff < maxTimeDiffAB)
//       {
//
//         // Found A-B signals in conincidence
//         if (saveHistos)
//         {
//           stats.getHisto1D("ab_tdiff_acc")->Fill(tDiff);
//         }
//
//         // If layer 1 (black module) - create AB Hit
//         auto layerID = scinSignals.at(i).getMatrix().getScin().getSlot().getLayer().getID();
//         if (layerID == 1)
//         {
//           hits.push_back(createHit(scinSignals.at(i), scinSignals.at(j)));
//         }
//
//         // Layer 2 or 4 - Red Module, searching for WLS signals
//         if (layerID == 2 || layerID == 4)
//         {
//           auto hitTime = (scinSignals.at(i).getTime() + scinSignals.at(j).getTime()) / 2.0;
//           auto wlsSigIndex = matchWLSSignal(wlsSignals, hitTime, minTimeDiffAB, maxTimeDiffAB, stats, saveHistos);
//           // Singal found if index different than -1
//           if (wlsSigIndex == -1)
//           {
//             hits.push_back(createHit(scinSignals.at(i), scinSignals.at(j)));
//           }
//           else
//           {
//             hits.push_back(createHit(scinSignals.at(i), scinSignals.at(j), wlsSignals.at(wlsSigIndex)));
//             wlsSignals.erase(wlsSignals.begin() + wlsSigIndex);
//           }
//         }
//         i = j;
//       }
//       else
//       {
//         if (saveHistos)
//         {
//           stats.getHisto1D("ab_tdiff_rej")->Fill(tDiff);
//         }
//         i = j - 1;
//       }
//       // i = j;
//       break;
//     }
//   }
//   return hits;
// }

/**
 * Checking times of WLS signals if they mach with hit time in conincidence window
 */
// int HitFinderTools::matchWLSSignal(vector<JPetMatrixSignal>& wlsSignals, double hitTime, double minTimeDiffAB, double maxTimeDiffAB,
//                                    JPetStatistics& stats, bool saveHistos)
// {
//   for (unsigned int i = 0; i < wlsSignals.size(); i++)
//   {
//     auto tDiff = fabs(hitTime - wlsSignals.at(i).getTime());
//     if (saveHistos)
//     {
//       stats.getHisto1D("hit_wls_tdiff_all")->Fill(tDiff);
//     }
//     if (tDiff >= minTimeDiffAB && tDiff <= maxTimeDiffAB)
//     {
//       if (saveHistos)
//       {
//         stats.getHisto1D("hit_wls_tdiff_acc")->Fill(tDiff);
//       }
//       return i;
//     }
//     else
//     {
//       if (saveHistos)
//       {
//         stats.getHisto1D("hit_wls_tdiff_rej")->Fill(tDiff);
//       }
//     }
//   }
//   return -1;
// }
