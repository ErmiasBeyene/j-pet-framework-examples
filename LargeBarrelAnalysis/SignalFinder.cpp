/**
 *  @copyright Copyright 2021 The J-PET Framework Authors. All rights reserved.
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
 *  @file SignalFinder.cpp
 */

#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetTimeWindow/JPetTimeWindow.h>
#include <JPetWriter/JPetWriter.h>

#include "SignalFinder.h"
#include <string>

#include <utility>
#include <vector>

using namespace jpet_options_tools;
using namespace std;

SignalFinder::SignalFinder(const char* name) : JPetUserTask(name) {}

SignalFinder::~SignalFinder() {}

bool SignalFinder::init()
{
  INFO("Signal finding started.");
  fOutputEvents = new JPetTimeWindow("JPetRawSignal");

  // Reading values from the user options if available
  // Time window parameter for leading edge
  if (isOptionSet(fParams.getOptions(), kEdgeMaxTimeParamKey))
  {
    fSigChEdgeMaxTime = getOptionAsDouble(fParams.getOptions(), kEdgeMaxTimeParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kEdgeMaxTimeParamKey.c_str(), fSigChEdgeMaxTime));
  }

  // Time window parameter for leading-trailing comparison
  if (isOptionSet(fParams.getOptions(), kLeadTrailMaxTimeParamKey))
  {
    fSigChLeadTrailMaxTime = getOptionAsDouble(fParams.getOptions(), kLeadTrailMaxTimeParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kLeadTrailMaxTimeParamKey.c_str(),
                 fSigChLeadTrailMaxTime));
  }

  // Get bool for using corrupted Signal Channels
  if (isOptionSet(fParams.getOptions(), kUseCorruptedSigChParamKey))
  {
    fUseCorruptedSigCh = getOptionAsBool(fParams.getOptions(), kUseCorruptedSigChParamKey);
    if (fUseCorruptedSigCh)
    {
      WARNING("Signal Finder is using Corrupted Signal Channels, as set by the user");
    }
    else
    {
      WARNING("Signal Finder is NOT using Corrupted Signal Channels, as set by the user");
    }
  }
  else
  {
    WARNING("Signal Finder is not using Corrupted Signal Channels (default option)");
  }
  // Reference detector photomultiplier identifier
  if (isOptionSet(fParams.getOptions(), kRefPMIDParamKey))
  {
    fRefPMID = getOptionAsInt(fParams.getOptions(), kRefPMIDParamKey);
    WARNING("Signal Finder is ignoring Corrupted Signal Channels from the Refference detector as specified by the user (photomultiplier number + "
            "std::to_string(fRefPMID)");
  }
  else
  {
    WARNING(
        "Signal Finder is ignoring Corrupted Signal Channels from the Refference detector(default photomultiplier number + std::to_string(fRefPMID)");
  }
  // Getting flags for saving histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey))
  {
    fSaveControlHistos = getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
  }
  // Check if the user requested ordering of thresholds by value
  if (isOptionSet(fParams.getOptions(), kOrderThresholdsByValueKey))
  {
    fOrderThresholdsByValue = getOptionAsBool(fParams.getOptions(), kOrderThresholdsByValueKey);
  }
  if (fOrderThresholdsByValue)
  {
    INFO("Threshold reordering was requested. Thresholds will be ordered by their values according to provided detector setup file.");
    fThresholdOrderings = SignalFinderTools::findThresholdOrders(getParamBank());
  }

  // Creating control histograms
  if (fSaveControlHistos)
  {
    initialiseHistograms();
  }
  return true;
}

bool SignalFinder::exec()
{
  // Getting the data from event in an apropriate format
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent))
  {
    // Distribute signal channels by PM IDs and filter out Corrupted SigChs if requested
    auto& sigChByPM = SignalFinderTools::getSigChByPM(timeWindow, fUseCorruptedSigCh, fRefPMID);
    // Building signals
    auto allSignals = SignalFinderTools::buildAllSignals(sigChByPM, fSigChEdgeMaxTime, fSigChLeadTrailMaxTime, getStatistics(), fSaveControlHistos,
                                                         fThresholdOrderings);
    // Saving method invocation
    saveRawSignals(allSignals);
  }
  else
  {
    return false;
  }
  return true;
}

bool SignalFinder::terminate()
{
  INFO("Signal finding ended.");
  return true;
}

void SignalFinder::saveRawSignals(const vector<JPetRawSignal>& rawSigVec)
{
  for (auto& rawSig : rawSigVec)
  {
    auto leads = rawSig.getPoints(JPetSigCh::Leading);
    auto trails = rawSig.getPoints(JPetSigCh::Trailing);
    if (leads.size() == trails.size())
    {
      fOutputEvents->add<JPetRawSignal>(rawSig);
      if (fSaveControlHistos)
      {
        auto pmID = rawSig.getPM().getID();
        for (int thr_i = 0; thr_i < leads.size() && thr_i < trails.size(); ++thr_i)
        {
          double tot_thr = trails.at(thr_i).getValue() - leads.at(thr_i).getValue();
          getStatistics().fillHistogram(Form("thr_tot_%d_pm", thr_i + 1), pmID, tot_thr);

          for (int thr_j = 0; thr_j < leads.size() && thr_j < trails.size(); ++thr_j)
          {
            if (thr_j == thr_i)
            {
              continue;
            }
            double tdiff = leads.at(thr_j).getValue() - leads.at(thr_i).getValue();
            getStatistics().fillHistogram(Form("thr_tdiff_%d_%d_pm", max(thr_j + 1, thr_i + 1), min(thr_j + 1, thr_i + 1)), pmID, tdiff);
          }
        }
      }
    }
  }
}

void SignalFinder::initialiseHistograms()
{
  getStatistics().createHistogramWithAxes(new TH1D("unused_sigch_all", "Unused Signal Channels", 8, 0.5, 8.5), "Signal label", "Number of SigChs");
  vector<pair<unsigned, string>> binLabels1 = {make_pair(1, "THR 1 Lead"),  make_pair(2, "THR 1 Trail"), make_pair(3, "THR 2 Lead"),
                                               make_pair(4, "THR 2 Trail"), make_pair(5, "THR 3 Lead"),  make_pair(6, "THR 3 Trail"),
                                               make_pair(7, "THR 4 Lead"),  make_pair(8, "THR 4 Trail")};
  getStatistics().setHistogramBinLabel("unused_sigch_all", getStatistics().AxisLabel::kXaxis, binLabels1);

  getStatistics().createHistogramWithAxes(new TH1D("unused_sigch_good", "Unused Signal Channels with GOOD flag", 8, 0.5, 8.5), "Signal label",
                                          "Number of SigChs");
  getStatistics().setHistogramBinLabel("unused_sigch_good", getStatistics().AxisLabel::kXaxis, binLabels1);

  getStatistics().createHistogramWithAxes(new TH1D("unused_sigch_corr", "Unused Signal Channels with CORRUPTED flag", 8, 0.5, 8.5), "Signal label",
                                          "Number of SigChs");
  getStatistics().setHistogramBinLabel("unused_sigch_corr", getStatistics().AxisLabel::kXaxis, binLabels1);

  vector<pair<unsigned, string>> binLabels2 = {make_pair(1, "GOOD"), make_pair(2, "CORRUPTED"), make_pair(3, "UNKNOWN")};
  getStatistics().createHistogramWithAxes(new TH1D("good_v_bad_raw_sigs", "Number of good and corrupted signals created", 3, 0.5, 3.5), "Flag",
                                          "Number of Raw Signals");
  getStatistics().setHistogramBinLabel("good_v_bad_raw_sigs", getStatistics().AxisLabel::kXaxis, binLabels2);

  auto minPMID = getParamBank().getPMs().begin()->first;
  auto maxPMID = getParamBank().getPMs().rbegin()->first;

  for (int i = 1; i <= kNumOfThresholds; ++i)
  {
    getStatistics().createHistogramWithAxes(new TH2D(Form("thr_tot_%d_pm", i), Form("TOT on THR%d per PM", i), maxPMID - minPMID + 1, minPMID - 0.5,
                                                     maxPMID + 0.5, 200, 0.0, 1.1 * fSigChLeadTrailMaxTime),
                                            "PM ID", "TOT [ps]");

    for (int j = i + 1; j <= kNumOfThresholds; ++j)
    {
      if (i == j)
      {
        continue;
      }
      getStatistics().createHistogramWithAxes(new TH2D(Form("thr_tdiff_%d_%d_pm", j, i), Form("Offsets THR%d THR%d per PM", j, i),
                                                       maxPMID - minPMID + 1, minPMID - 0.5, maxPMID + 0.5, 200, -1.1 * fSigChEdgeMaxTime,
                                                       1.1 * fSigChEdgeMaxTime),
                                              "PM ID", "THR signals time diff [ps]");
    }
  }
}
