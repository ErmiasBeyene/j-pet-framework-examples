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
 *  @file SignalTransformer.cpp
 */

#include "JPetWriter/JPetWriter.h"

#include "HitFinderTools.h"
#include "SignalTransformer.h"

using namespace jpet_options_tools;
using namespace std;

SignalTransformer::SignalTransformer(const char* name) : JPetUserTask(name) {}

SignalTransformer::~SignalTransformer() {}

bool SignalTransformer::init()
{
  INFO("Signal transforming started: Raw to Reco and Phys");
  fOutputEvents = new JPetTimeWindow("JPetPhysSignal");
  // Getting bool for using bad signals
  if (isOptionSet(fParams.getOptions(), kUseCorruptedSignalsParamKey))
  {
    fUseCorruptedSignals = getOptionAsBool(fParams.getOptions(), kUseCorruptedSignalsParamKey);
    if (fUseCorruptedSignals)
    {
      WARNING("Signal Transformer is using Corrupted Signals, as set by the user");
    }
    else
    {
      WARNING("Signal Transformer is NOT using Corrupted Signals, as set by the user");
    }
  }
  else
  {
    WARNING("Signal Transformer is not using Corrupted Signals (default option)");
  }
  // Walk correction constants (for each threshold separately)
  if (isOptionSet(fParams.getOptions(), kWalkCorrConst1ParamKey))
  {
    fWalkCorrConst[0] = getOptionAsFloat(fParams.getOptions(), kWalkCorrConst1ParamKey);
  }
  if (isOptionSet(fParams.getOptions(), kWalkCorrConst2ParamKey))
  {
    fWalkCorrConst[1] = getOptionAsFloat(fParams.getOptions(), kWalkCorrConst2ParamKey);
  }
  if (isOptionSet(fParams.getOptions(), kWalkCorrConst3ParamKey))
  {
    fWalkCorrConst[2] = getOptionAsFloat(fParams.getOptions(), kWalkCorrConst3ParamKey);
  }
  if (isOptionSet(fParams.getOptions(), kWalkCorrConst4ParamKey))
  {
    fWalkCorrConst[3] = getOptionAsFloat(fParams.getOptions(), kWalkCorrConst4ParamKey);
  }
  // Getting bool for saving histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey))
  {
    fSaveControlHistos = getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
  }

  // Control histograms
  if (fSaveControlHistos)
  {
    initialiseHistograms();
  }
  return true;
}

bool SignalTransformer::exec()
{
  if (auto& timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent))
  {
    uint n = timeWindow->getNumberOfEvents();
    for (uint i = 0; i < n; ++i)
    {
      auto& rawSignal = dynamic_cast<const JPetRawSignal&>(timeWindow->operator[](i));
      if (!fUseCorruptedSignals && rawSignal.getRecoFlag() == JPetBaseSignal::Corrupted)
      {
        continue;
      }
      if (fSaveControlHistos)
      {
        auto leads = rawSignal.getPoints(JPetSigCh::Leading, JPetRawSignal::ByThrNum);
        auto trails = rawSignal.getPoints(JPetSigCh::Trailing, JPetRawSignal::ByThrNum);
        for (unsigned int i = 0; i < leads.size(); i++)
        {
          getStatistics().fillHistogram("raw_sigs_multi", 2 * i + 1);
        }
        for (unsigned int i = 0; i < trails.size(); i++)
        {
          getStatistics().fillHistogram("raw_sigs_multi", 2 * (i + 1));
        }
        if (rawSignal.getRecoFlag() == JPetBaseSignal::Good)
        {
          getStatistics().fillHistogram("good_vs_bad_signals", 1);
          for (unsigned int i = 0; i < leads.size(); i++)
          {
            getStatistics().fillHistogram("raw_sigs_multi_good", 2 * i + 1);
          }
          for (unsigned int i = 0; i < trails.size(); i++)
          {
            getStatistics().fillHistogram("raw_sigs_multi_good", 2 * (i + 1));
          }
        }
        else if (rawSignal.getRecoFlag() == JPetBaseSignal::Corrupted)
        {
          getStatistics().fillHistogram("good_vs_bad_signals", 2);
          for (unsigned int i = 0; i < leads.size(); i++)
          {
            getStatistics().fillHistogram("raw_sigs_multi_corr", 2 * i + 1);
            if (leads.at(i).getRecoFlag() == JPetSigCh::Good)
            {
              getStatistics().fillHistogram("raw_sigs_multi_corr_sigch_good", 2 * i + 1);
            }
            else if (leads.at(i).getRecoFlag() == JPetSigCh::Corrupted)
            {
              getStatistics().fillHistogram("raw_sigs_multi_corr_sigch_corr", 2 * i + 1);
            }
          }
          for (unsigned int i = 0; i < trails.size(); i++)
          {
            getStatistics().fillHistogram("raw_sigs_multi_corr", 2 * (i + 1));
            if (trails.at(i).getRecoFlag() == JPetSigCh::Good)
            {
              getStatistics().fillHistogram("raw_sigs_multi_corr_sigch_good", 2 * (i + 1));
            }
            else if (trails.at(i).getRecoFlag() == JPetSigCh::Corrupted)
            {
              getStatistics().fillHistogram("raw_sigs_multi_corr_sigch_corr", 2 * (i + 1));
            }
          }
        }
        else if (rawSignal.getRecoFlag() == JPetBaseSignal::Unknown)
        {
          getStatistics().fillHistogram("good_vs_bad_signals", 3);
        }
      }
      // Make Reco Signal from Raw Signal
      auto recoSignal = createRecoSignal(rawSignal);
      // Make Phys Signal from Reco Signal and save
      auto physSignal = createPhysSignal(recoSignal);
      fOutputEvents->add<JPetPhysSignal>(physSignal);
    }
  }
  else
  {
    return false;
  }
  return true;
}

bool SignalTransformer::terminate()
{
  INFO("Signal transforming finished");
  return true;
}

/**
 * Method rewrites Raw Signal to Reco Signal. All fields set to -1.
 */
JPetRecoSignal SignalTransformer::createRecoSignal(const JPetRawSignal& rawSignal)
{
  JPetRecoSignal recoSignal;
  recoSignal.setRawSignal(rawSignal);
  recoSignal.setAmplitude(-1.0);
  recoSignal.setOffset(-1.0);
  recoSignal.setCharge(-1.0);
  recoSignal.setDelay(-1.0);
  recoSignal.setRecoFlag(rawSignal.getRecoFlag());
  recoSignal.setPM(rawSignal.getPM());
  return recoSignal;
}

/**
 * Method rewrites Reco Signal to Phys Signal.
 * Time of Signal set to time of the Leading Signal Channel at the lowest threshold.
 * Other fields are set to -1, quality fields set to 0.
 */
JPetPhysSignal SignalTransformer::createPhysSignal(const JPetRecoSignal& recoSignal)
{
  JPetPhysSignal physSignal;
  correctForWalk(recoSignal);
  vector<JPetSigCh> leadingSigChVec = recoSignal.getRawSignal().getPoints(JPetSigCh::Leading, JPetRawSignal::ByThrValue);
  physSignal.setPM(recoSignal.getPM());
  physSignal.setRecoSignal(recoSignal);
  physSignal.setPhe(-1.0);
  physSignal.setQualityOfPhe(0.0);
  physSignal.setQualityOfTime(0.0);
  physSignal.setRecoFlag(recoSignal.getRecoFlag());
  physSignal.setTime(leadingSigChVec.at(0).getValue());

  if (fSaveControlHistos)
  {
    auto tot1 = HitFinderTools::calculateSignalTOT(physSignal, HitFinderTools::kSimplified);
    auto tot2 = HitFinderTools::calculateSignalTOT(physSignal, HitFinderTools::kThresholdRectangular);
    auto tot3 = HitFinderTools::calculateSignalTOT(physSignal, HitFinderTools::kThresholdTrapeze);

    if (tot1 > 0.0)
    {
      getStatistics().fillHistogram("phys_sig_tot_pm_simple", physSignal.getPM().getID(), tot1);
    }
    if (tot2 > 0.0)
    {
      getStatistics().fillHistogram("phys_sig_tot_pm_rectangular", physSignal.getPM().getID(), tot2);
    }
    if (tot3 > 0.0)
    {
      getStatistics().fillHistogram("phys_sig_tot_pm_trapeze", physSignal.getPM().getID(), tot3);
    }
  }

  return physSignal;
}

/**
 * Walk correction applied to the SigCh times on both edges
 */
void SignalTransformer::correctForWalk(const JPetRecoSignal& recoSignal)
{
  vector<JPetSigCh> leadingSigChVec = recoSignal.getRawSignal().getPoints(JPetSigCh::Leading, JPetRawSignal::ByThrValue);
  vector<JPetSigCh> trailingSigChVec = recoSignal.getRawSignal().getPoints(JPetSigCh::Trailing, JPetRawSignal::ByThrValue);
  double TOT = 0.;
  for (unsigned i = 0; i < leadingSigChVec.size() && i < trailingSigChVec.size(); i++)
  {
    TOT += trailingSigChVec.at(i).getValue() - leadingSigChVec.at(i).getValue();
  }
  for (unsigned i = 0; i < leadingSigChVec.size(); i++)
  {
    if (TOT > 0. && fWalkCorrConst[i] > 0.)
    {
      double WalkCorr = fWalkCorrConst[i] / sqrt(TOT);
      leadingSigChVec.at(i).setValue(leadingSigChVec.at(i).getValue() - WalkCorr);
      getStatistics().fillHistogram("WalkCorrLead", WalkCorr);
    }
    for (unsigned i = 0; i < trailingSigChVec.size(); i++)
    {
      if (TOT > 0. && fWalkCorrConst[i] > 0.)
      {
        double WalkCorr = fWalkCorrConst[i] / sqrt(TOT);
        trailingSigChVec.at(i).setValue(trailingSigChVec.at(i).getValue() - WalkCorr);
        getStatistics().fillHistogram("WalkCorrTrail", WalkCorr);
      }
    }
  }
}
void SignalTransformer::initialiseHistograms()
{
  getStatistics().createHistogramWithAxes(new TH1D("good_vs_bad_signals", "Number of good and corrupted signals created", 3, 0.5, 3.5), "Flag",
                                          "Number of Signals");
  vector<pair<unsigned, string>> binLabels1 = {make_pair(1, "GOOD"), make_pair(2, "CORRUPTED"), make_pair(3, "UNKNOWN")};
  getStatistics().setHistogramBinLabel("good_vs_bad_signals", getStatistics().AxisLabel::kXaxis, binLabels1);

  getStatistics().createHistogramWithAxes(new TH1D("raw_sigs_multi", "Multiplicity of created Raw Signals", 8, 0.5, 8.5), "Signal label",
                                          "Number of SigChs");

  vector<pair<unsigned, string>> binLabels2 = {make_pair(1, "THR 1 Lead"),  make_pair(2, "THR 1 Trail"), make_pair(3, "THR 2 Lead"),
                                               make_pair(4, "THR 2 Trail"), make_pair(5, "THR 3 Lead"),  make_pair(6, "THR 3 Trail"),
                                               make_pair(7, "THR 4 Lead"),  make_pair(8, "THR 4 Trail")};
  getStatistics().setHistogramBinLabel("raw_sigs_multi", getStatistics().AxisLabel::kXaxis, binLabels2);

  getStatistics().createHistogramWithAxes(new TH1D("raw_sigs_multi_good", "Multiplicity of created Raw Signals with GOOD flag", 8, 0.5, 8.5),
                                          "Signal label", "Number of SigChs");
  getStatistics().setHistogramBinLabel("raw_sigs_multi_good", getStatistics().AxisLabel::kXaxis, binLabels2);

  getStatistics().createHistogramWithAxes(new TH1D("raw_sigs_multi_corr", "Multiplicity of created Raw Signals with CORRUPTED flag", 8, 0.5, 8.5),
                                          "Signal label", "Number of SigChs");
  getStatistics().setHistogramBinLabel("raw_sigs_multi_corr", getStatistics().AxisLabel::kXaxis, binLabels2);

  getStatistics().createHistogramWithAxes(
      new TH1D("raw_sigs_multi_corr_sigch_good", "Multiplicity of created Raw Signals with CORRUPTED flag - GOOD SigCh only", 8, 0.5, 8.5),
      "Signal label", "Number of GOOD SigChs");
  getStatistics().setHistogramBinLabel("raw_sigs_multi_corr_sigch_good", getStatistics().AxisLabel::kXaxis, binLabels2);

  getStatistics().createHistogramWithAxes(
      new TH1D("raw_sigs_multi_corr_sigch_corr", "Multiplicity of created Raw Signals with CORRUPTED flag - CORRUPTED SigCh only", 8, 0.5, 8.5),
      "Signal label", "Number of CORRUPTED SigChs");
  getStatistics().setHistogramBinLabel("raw_sigs_multi_corr_sigch_corr", getStatistics().AxisLabel::kXaxis, binLabels2);

  getStatistics().createHistogramWithAxes(new TH1D("WalkCorrLead", "Walk Correction applied on the leading edge", 1000, -0.025, 49.975),
                                          "Walk Correction [ps]", "Walk Correction applied on the leading edge");

  getStatistics().createHistogramWithAxes(new TH1D("WalkCorrTrail", "Walk Correction applied on the trailing edge", 1000, -0.025, 49.975),
                                          "Walk Correction [ps]", "Walk Correction applied on the trailing edge");

  auto minPMID = getParamBank().getPMs().begin()->first;
  auto maxPMID = getParamBank().getPMs().rbegin()->first;

  getStatistics().createHistogramWithAxes(new TH2D("phys_sig_tot_pm_simple", "TOT of created signals, calculated with simplified method",
                                                   maxPMID - minPMID + 1, minPMID - 0.5, maxPMID + 0.5, 200, 0.0, 50000.0),
                                          "PM ID", "Signal TOT [ps]");

  getStatistics().createHistogramWithAxes(new TH2D("phys_sig_tot_pm_rectangular", "TOT of created signals, calculated with rectangular method",
                                                   maxPMID - minPMID + 1, minPMID - 0.5, maxPMID + 0.5, 200, 0.0, 50000.0),
                                          "PM ID", "Signal TOT [ps]");

  getStatistics().createHistogramWithAxes(new TH2D("phys_sig_tot_pm_trapeze", "TOT of created signals, calculated with trapeze method",
                                                   maxPMID - minPMID + 1, minPMID - 0.5, maxPMID + 0.5, 200, 0.0, 50000.0),
                                          "PM ID", "Signal TOT [ps]");
}
