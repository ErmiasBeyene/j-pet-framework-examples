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
 *  @file HitFinderTools.h
 */

#ifndef HITFINDERTOOLS_H
#define HITFINDERTOOLS_H

#include <JPetHit/JPetHit.h>
#include <JPetMatrixSignal/JPetMatrixSignal.h>
#include <JPetStatistics/JPetStatistics.h>
#include <JPetTimeWindow/JPetTimeWindow.h>
#include <string>
#include <vector>

/**
 * @brief Tools set fot HitFinder module
 *
 * Tols include methods of signal mapping and matching,
 * helpers of sorting, radian check and methods for reference detecctor
 *
 */
class HitFinderTools
{
public:
  static void sortByTime(std::vector<JPetMatrixSignal>& signals);

  static std::map<std::string, std::map<int, std::vector<JPetMatrixSignal>>> getMappedSignals(const JPetTimeWindow* timeWindow);

  static std::vector<JPetHit> matchAllSignals(std::map<std::string, std::map<int, std::vector<JPetMatrixSignal>>>& signalSidesMap,
                                              double minTimeDiffAB, double maxTimeDiffAB, JPetStatistics& stats, bool saveHistos);

  static std::vector<JPetHit> matchSignals(std::vector<JPetMatrixSignal>& scinSignals, double minTimeDiffAB, double maxTimeDiffAB,
                                           JPetStatistics& stats, bool saveHistos);

  static JPetHit createScinHit(const JPetMatrixSignal& signal1, const JPetMatrixSignal& signal2);

  static JPetHit createWLSHit(const JPetMatrixSignal& signalWLS);

  static JPetHit createHit(const JPetMatrixSignal& signal1, const JPetMatrixSignal& signal2, const JPetMatrixSignal& signalWLS);

  static JPetHit createDummyHit(const JPetMatrixSignal& signal);
  static double calculateTOT(JPetHit& hit);

  /// old
  // static std::vector<JPetHit> matchSignals(
  //   std::vector<JPetMatrixSignal>& scinSignals, std::vector<JPetMatrixSignal>& wlsSignals,
  //   double minTimeDiffAB, double maxTimeDiffAB, JPetStatistics& stats, bool saveHistos
  // );
  // static int matchWLSSignal(std::vector<JPetMatrixSignal>& wlsSignals, double hitTime, double minTimeDiffAB, double maxTimeDiffAB,
  //                           JPetStatistics& stats, bool saveHistos);
};

#endif /* !HITFINDERTOOLS_H */
