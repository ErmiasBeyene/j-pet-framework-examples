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
 *  @file SignalFinder.cpp
 */

using namespace std;

#include "SignalFinder.h"
#include "SignalFinderTools.h"
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetTimeWindow/JPetTimeWindow.h>
#include <JPetWriter/JPetWriter.h>
#include <TRandom.h>
#include <string>
#include <utility>
#include <vector>

using namespace jpet_options_tools;

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

  if (isOptionSet(fParams.getOptions(), kMinPMIDParamKey))
  {
    fMinPMID = getOptionAsInt(fParams.getOptions(), kMinPMIDParamKey);
  }
  else
  {
    fMinPMID = getParamBank().getPMs().begin()->first;
  }

  if (isOptionSet(fParams.getOptions(), kMaxPMIDParamKey))
  {
    fMaxPMID = getOptionAsInt(fParams.getOptions(), kMaxPMIDParamKey);
  }
  else
  {
    fMaxPMID = getParamBank().getPMs().rbegin()->first;
  }

  // Getting bool for saving histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey))
  {
    fSaveControlHistos = getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
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
    // Distribute signal channels by PM IDs
    auto& sigChByPM = SignalFinderTools::getSigChByPM(timeWindow);
    // Building signals
    auto allSignals = SignalFinderTools::buildAllSignals(sigChByPM, fSigChEdgeMaxTime, fSigChLeadTrailMaxTime, kNumOfThresholds, getStatistics(),
                                                         fSaveControlHistos);
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

/**
 * Saving Raw Signals that have leading-trailing pairs,
 * otherwise filling histogram with incomplete signals
 */
void SignalFinder::saveRawSignals(const vector<JPetRawSignal>& rawSigVec)
{
  if (fSaveControlHistos && rawSigVec.size() > 0)
  {
    getStatistics().getHisto1D("rawsig_tslot")->Fill(rawSigVec.size());
  }

  for (auto& rawSig : rawSigVec)
  {

    auto leads = rawSig.getPoints(JPetSigCh::Leading, JPetRawSignal::ByThrValue);
    auto trails = rawSig.getPoints(JPetSigCh::Trailing, JPetRawSignal::ByThrValue);

    // Saving only signals with lead-trail pair on threshold
    if (!(leads.size() == kNumOfThresholds && trails.size() == kNumOfThresholds))
    {
      continue;
    }
    fOutputEvents->add<JPetRawSignal>(rawSig);

    if (fSaveControlHistos)
    {

      auto pmID = rawSig.getPM().getID();

      if (pmID < fMinPMID || pmID > fMaxPMID)
      {
        continue;
      }

      // Filling histo with average TOT (sum divided by number of thresholds (2))
      getStatistics().getHisto2D("tot_sipm")->Fill(pmID, rawSig.getTOT() / kNumOfThresholds);

      if (rawSig.getPM().getPosition() > -99.0)
      {
        getStatistics().getHisto1D("wls_sig_pos")->Fill(rawSig.getPM().getPosition());
      }

      // Multiplicities and occupancies downscaled
      if (gRandom->Uniform() < fScalingFactor)
      {
        getStatistics().getHisto1D("rawsig_per_pm")->Fill(pmID);
        getStatistics().getHisto1D("rawsig_multi")->Fill(leads.size() + trails.size());

        for (auto& sigCh : leads)
        {
          getStatistics().getHisto1D("rawsig_thr_occ")->Fill(sigCh.getChannel().getThresholdNumber());
        }

        for (auto& sigCh : trails)
        {
          getStatistics().getHisto1D("rawsig_thr_occ")->Fill(sigCh.getChannel().getThresholdNumber());
        }
      }
    }
  }
}

void SignalFinder::initialiseHistograms()
{

  // Unused objects stats
  getStatistics().createHistogram(new TH1F("unused_sigch_thr", "Unused Signal Channels per THR (downscaled)", 5, 0.5, 5.5));
  getStatistics().getHisto1D("unused_sigch_thr")->GetXaxis()->SetBinLabel(1, "THR 1 Lead");
  getStatistics().getHisto1D("unused_sigch_thr")->GetXaxis()->SetBinLabel(2, "THR 1 Trail");
  getStatistics().getHisto1D("unused_sigch_thr")->GetXaxis()->SetBinLabel(3, "THR 2 Lead");
  getStatistics().getHisto1D("unused_sigch_thr")->GetXaxis()->SetBinLabel(4, "THR 2 Trail");
  getStatistics().getHisto1D("unused_sigch_thr")->GetXaxis()->SetBinLabel(5, "  ");
  getStatistics().getHisto1D("unused_sigch_thr")->GetYaxis()->SetTitle("Number of SigChs");

  getStatistics().createHistogram(
      new TH1F("unused_sigch_pm", "Unused Signal Channels per SiPM", fMaxPMID - fMinPMID + 1, fMinPMID - 0.5, fMaxPMID + 0.5));
  getStatistics().getHisto1D("unused_sigch_pm")->GetXaxis()->SetTitle("SiPM ID");
  getStatistics().getHisto1D("unused_sigch_pm")->GetYaxis()->SetTitle("Number of Signal Channels");

  // Occupancies and multiplicities
  getStatistics().createHistogram(new TH1F("rawsig_per_pm", "Raw Signals per SiPM", fMaxPMID - fMinPMID + 1, fMinPMID - 0.5, fMaxPMID + 0.5));
  getStatistics().getHisto1D("rawsig_per_pm")->GetXaxis()->SetTitle("SiPM ID");
  getStatistics().getHisto1D("rawsig_per_pm")->GetYaxis()->SetTitle("Number of Raw Signals");

  getStatistics().createHistogram(new TH1F("rawsig_thr_occ", "Thresholds occupation in createds Raw Signals", 3, 0.5, 3.5));
  getStatistics().getHisto1D("rawsig_thr_occ")->GetXaxis()->SetTitle("Threshold number");
  getStatistics().getHisto1D("rawsig_thr_occ")->GetYaxis()->SetTitle("Number of Signal Channels");

  getStatistics().createHistogram(new TH1F("rawsig_multi", "Raw Signal Multiplicity", 6, 0.5, 6.5));
  getStatistics().getHisto1D("rawsig_multi")->GetXaxis()->SetTitle("Total number of SigChs in RawSig");
  getStatistics().getHisto1D("rawsig_multi")->GetYaxis()->SetTitle("Number of Signal Channels");

  getStatistics().createHistogram(new TH1F("rawsig_tslot", "Number of Raw Signals in Time Window", 70, 0.5, 71.5));
  getStatistics().getHisto1D("rawsig_tslot")->GetXaxis()->SetTitle("Number of Raw Signal in Time Window");
  getStatistics().getHisto1D("rawsig_tslot")->GetYaxis()->SetTitle("Number of Time Windows");

  getStatistics().createHistogram(new TH1F("wls_sig_pos", "Signal occupancy of WLS SiPMs in position along Z axis", 64, -20.38, 20.38));
  getStatistics().getHisto1D("wls_sig_pos")->GetXaxis()->SetTitle("SiPM position [cm]");
  getStatistics().getHisto1D("wls_sig_pos")->GetYaxis()->SetTitle("Number of Raw Signals");

  getStatistics().createHistogram(new TH2F("tot_sipm", "Signal Time over Threshold per SiPM", fMaxPMID - fMinPMID + 1, fMinPMID - 0.5, fMaxPMID + 0.5,
                                           200, 0.0, 1.1 * fSigChLeadTrailMaxTime));
  getStatistics().getHisto2D("tot_sipm")->GetXaxis()->SetTitle("SiPM ID");
  getStatistics().getHisto2D("tot_sipm")->GetYaxis()->SetTitle("TOT [ps]");
}
