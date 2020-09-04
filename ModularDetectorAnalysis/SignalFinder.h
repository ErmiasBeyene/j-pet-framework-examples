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
 *  @file SignalFinder.h
 */

#ifndef SIGNALFINDER_H
#define SIGNALFINDER_H

#include <JPetRawSignal/JPetRawSignal.h>
#include <JPetUserTask/JPetUserTask.h>
#include <vector>

class JPetWriter;

#ifdef __CINT__
#define override
#endif

/**
 * @brief User Task: method organizing Signal Channels to Raw Signals
 *
 * Task organizes Signal Channels from every JPetTimeWindow to Raw Signals
 * Parameters for time window values used in tools can be specified in user options,
 * default are provided.
 */
class SignalFinder: public JPetUserTask
{
public:
  SignalFinder(const char* name);
  virtual ~SignalFinder();
  virtual bool init() override;
  virtual bool exec() override;
  virtual bool terminate() override;

protected:
  void saveRawSignals(const std::vector<JPetRawSignal>& sigChVec);
  const std::string kLeadTrailMaxTimeParamKey = "SignalFinder_LeadTrailMaxTime_double";
  const std::string kSaveControlHistosParamKey = "Save_Control_Histograms_bool";
  const std::string kEdgeMaxTimeParamKey = "SignalFinder_EdgeMaxTime_double";
  const std::string kMinPMIDParamKey = "Histo_MinPMID_int";
  const std::string kMaxPMIDParamKey = "Histo_MaxPMID_int";

  const int kNumOfThresholds = 2;
  double fSigChLeadTrailMaxTime = 23000.0;
  double fSigChEdgeMaxTime = 5000.0;
  bool fSaveControlHistos = true;
  double fScalingFactor = 1.0;
  int fMinPMID = 0;
  int fMaxPMID = 1;

  void initialiseHistograms();
};

#endif /* !SIGNALFINDER_H */
