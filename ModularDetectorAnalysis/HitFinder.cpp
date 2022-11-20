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
 *  @file HitFinder.cpp
 */

using namespace std;

#include <JPetAnalysisTools/JPetAnalysisTools.h>
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "HitFinder.h"
#include "HitFinderTools.h"
#include <map>
#include <string>
#include <vector>

using namespace jpet_options_tools;

HitFinder::HitFinder(const char* name) : JPetUserTask(name) {}

HitFinder::~HitFinder() {}

bool HitFinder::init()
{
  INFO("Hit finding Started");
  fOutputEvents = new JPetTimeWindow("JPetHit");

  // Reading values from the user options if available
  // Allowed time difference between signals on A and B sides - min and max
  if (isOptionSet(fParams.getOptions(), kMinABTimeDiffParamKey))
  {
    fMinABTimeDiff = getOptionAsDouble(fParams.getOptions(), kMinABTimeDiffParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kMinABTimeDiffParamKey.c_str(), fMinABTimeDiff));
  }

  if (isOptionSet(fParams.getOptions(), kMaxABTimeDiffParamKey))
  {
    fMaxABTimeDiff = getOptionAsDouble(fParams.getOptions(), kMaxABTimeDiffParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kMaxABTimeDiffParamKey.c_str(), fMaxABTimeDiff));
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

bool HitFinder::exec()
{
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent))
  {
    auto mappedSignals = HitFinderTools::getMappedSignals(timeWindow);
    auto allHits = HitFinderTools::matchAllSignals(mappedSignals, fMinABTimeDiff, fMaxABTimeDiff, getStatistics(), fSaveControlHistos);
    saveHits(allHits);
  }
  else
    return false;
  return true;
}

bool HitFinder::terminate()
{
  INFO("Hit finding ended");
  return true;
}

void HitFinder::saveHits(const std::vector<JPetHit>& hits)
{
  if (fSaveControlHistos)
  {
    getStatistics().getHisto1D("hits_tslot")->Fill(hits.size());
  }

  auto sortedHits = JPetAnalysisTools::getHitsOrderedByTime(hits);
  for (auto& hit : sortedHits)
  {
    fOutputEvents->add<JPetHit>(hit);
    if (fSaveControlHistos)
    {
      if (!hit.getScin().isNullObject())
      {

        auto minScin = min(hit.getSignalA().getMatrix().getScin().getID(), hit.getSignalB().getMatrix().getScin().getID());
        auto maxScin = max(hit.getSignalA().getMatrix().getScin().getID(), hit.getSignalB().getMatrix().getScin().getID());

        auto minMtx = min(hit.getSignalA().getMatrix().getID(), hit.getSignalB().getMatrix().getID());
        auto maxMtx = max(hit.getSignalA().getMatrix().getID(), hit.getSignalB().getMatrix().getID());

        getStatistics().getHisto2D("hit_scin_scin_id_1")->Fill(minScin, maxScin);
        getStatistics().getHisto2D("hit_scin_scin_id_2")->Fill(minScin, maxScin);
        getStatistics().getHisto2D("hit_scin_scin_id_3")->Fill(minScin, maxScin);

        getStatistics().getHisto2D("hit_mtx_mtx_id_1")->Fill(minMtx, maxMtx);
        getStatistics().getHisto2D("hit_mtx_mtx_id_2")->Fill(minMtx, maxMtx);
        getStatistics().getHisto2D("hit_mtx_mtx_id_3")->Fill(minMtx, maxMtx);
      }
    }
    //   auto multi = ((int)hit.getQualityOfEnergy());
    //
    //   getStatistics().getHisto2D("time_diff_per_scin")->Fill(hit.getTimeDiff(), hit.getScin().getID());
    //   getStatistics().getHisto1D("hit_sig_multi")->Fill(multi);
    //   if (hit.getPosZ() != -99.0)
    //   {
    //     getStatistics().getHisto1D("hit_pos_z_wls")->Fill(hit.getPosZ());
    //   }

    // auto scinID = hit.getScin().getID();
    // getStatistics().getHisto1D(Form("hit_sig_multi_scin_%d", scinID))->Fill(hit.getQualityOfEnergy());
    // getStatistics().getHisto1D(Form("hit_tdiff_scin_%d_m_%d", hit.getScin().getID(), multi))->Fill(hit.getTimeDiff());
    // if (hit.getEnergy() != 0.0)
    // {
    //   getStatistics().getHisto1D(Form("hit_tot_scin_%d_m_%d", hit.getScin().getID(), multi))->Fill(hit.getEnergy() / ((double)multi));
    // }
    // stats.getHisto2D("hit_pos_per_scin")->Fill(hit.getPosZ(), hit.getScin().getID());
    // }
  }
}

void HitFinder::initialiseHistograms()
{
  getStatistics().createHistogram(new TH1F("hits_tslot", "Number of Hits in Time Window", 30, 0.5, 31.5));
  getStatistics().getHisto1D("hits_tslot")->GetXaxis()->SetTitle("Hits in Time Slot");
  getStatistics().getHisto1D("hits_tslot")->GetYaxis()->SetTitle("Number of Time Slots");

  auto minScinID = getParamBank().getScins().begin()->first;
  auto maxScinID = getParamBank().getScins().rbegin()->first;

  auto minMtxID = getParamBank().getMatrices().begin()->first;
  auto maxMtxID = getParamBank().getMatrices().rbegin()->first;

  // Checking signal building
  getStatistics().createHistogram(new TH2F("hit_scin_scin_id_1", "Scins ID of AB Signals in Hits", 13, 0.5, 13.5, 13, 0.5, 13.5));
  getStatistics().getHisto2D("hit_scin_scin_id_1")->GetXaxis()->SetTitle("Scintillator ID");
  getStatistics().getHisto2D("hit_scin_scin_id_1")->GetYaxis()->SetTitle("Scintillator ID");

  getStatistics().createHistogram(new TH2F("hit_scin_scin_id_2", "Scins ID of AB Signals in Hits", 13, 13.5, 26.5, 13, 13.5, 26.5));
  getStatistics().getHisto2D("hit_scin_scin_id_2")->GetXaxis()->SetTitle("Scintillator ID");
  getStatistics().getHisto2D("hit_scin_scin_id_2")->GetYaxis()->SetTitle("Scintillator ID");

  getStatistics().createHistogram(new TH2F("hit_scin_scin_id_3", "Scins ID of AB Signals in Hits", 13, 26.5, 39.5, 13, 26.5, 39.5));
  getStatistics().getHisto2D("hit_scin_scin_id_3")->GetXaxis()->SetTitle("Scintillator ID");
  getStatistics().getHisto2D("hit_scin_scin_id_3")->GetYaxis()->SetTitle("Scintillator ID");

  getStatistics().createHistogram(new TH2F("hit_mtx_mtx_id_1", "Matrix ID of AB Signals in Hits", 26, 0.5, 26.5, 26, 0.5, 26.5));
  getStatistics().getHisto2D("hit_mtx_mtx_id_1")->GetXaxis()->SetTitle("Matrix ID");
  getStatistics().getHisto2D("hit_mtx_mtx_id_1")->GetYaxis()->SetTitle("Matrix ID");

  getStatistics().createHistogram(new TH2F("hit_mtx_mtx_id_2", "Matrix ID of AB Signals in Hits", 26, 26.5, 52.5, 26, 26.5, 52.5));
  getStatistics().getHisto2D("hit_mtx_mtx_id_2")->GetXaxis()->SetTitle("Matrix ID");
  getStatistics().getHisto2D("hit_mtx_mtx_id_2")->GetYaxis()->SetTitle("Matrix ID");

  getStatistics().createHistogram(new TH2F("hit_mtx_mtx_id_3", "Matrix ID of AB Signals in Hits", 26, 52.5, 78.5, 26, 52.5, 78.5));
  getStatistics().getHisto2D("hit_mtx_mtx_id_3")->GetXaxis()->SetTitle("Matrix ID");
  getStatistics().getHisto2D("hit_mtx_mtx_id_3")->GetYaxis()->SetTitle("Matrix ID");

  getStatistics().createHistogram(new TH1F("hit_tdiff", "A-B Signals Time Difference", 200, -1.1 * fMaxABTimeDiff, 1.1 * fMaxABTimeDiff));
  getStatistics().getHisto1D("hit_tdiff")->GetXaxis()->SetTitle("A-B time difference [ps]");
  getStatistics().getHisto1D("hit_tdiff")->GetYaxis()->SetTitle("Number of Hits");

  getStatistics().createHistogram(
      new TH1F("hit_per_scin", "Number of Hits in Scintillator", maxScinID - minScinID + 1, minScinID - 0.5, maxScinID + 0.5));
  getStatistics().getHisto1D("hit_per_scin")->GetXaxis()->SetTitle("Scintillator ID");
  getStatistics().getHisto1D("hit_per_scin")->GetYaxis()->SetTitle("Number of Hits");

  getStatistics().createHistogram(new TH2F("time_diff_per_scin", "Signals Time Difference per Scintillator ID", 200, -1.1 * fMaxABTimeDiff,
                                           1.1 * fMaxABTimeDiff, maxScinID - minScinID + 1, minScinID - 0.5, maxScinID + 0.5));
  getStatistics().getHisto2D("time_diff_per_scin")->GetXaxis()->SetTitle("A-B time difference [ps]");
  getStatistics().getHisto2D("time_diff_per_scin")->GetYaxis()->SetTitle("ID of Scintillator");

  getStatistics().createHistogram(new TH2F("hit_pos_XY", "Hit Position XY projection", 31, -15.5, 15.5, 21, -10.5, 10.5));
  getStatistics().getHisto2D("hit_pos_XY")->GetXaxis()->SetTitle("Y [cm]");
  getStatistics().getHisto2D("hit_pos_XY")->GetYaxis()->SetTitle("X [cm]");

  getStatistics().createHistogram(new TH1F("hit_pos_z", "Hit Z axis position", 100, -25.0, 25.0));
  getStatistics().getHisto1D("hit_pos_z")->GetXaxis()->SetTitle("z [cm]");
  getStatistics().getHisto1D("hit_pos_z")->GetYaxis()->SetTitle("Number of Hits");

  getStatistics().createHistogram(new TH1F("hit_pos_z_wls", "Hit Z axis position based on WLS position", 100, -25.0, 25.0));
  getStatistics().getHisto1D("hit_pos_z_wls")->GetXaxis()->SetTitle("z [cm]");
  getStatistics().getHisto1D("hit_pos_z_wls")->GetYaxis()->SetTitle("Number of Hits");

  // getStatistics().createHistogram(new TH2F(
  //   "hit_pos_per_scin", "Hit Position per Scintillator ID",
  //   200, -50.0, 50.0, maxScinID-minScinID+1, minScinID-0.5, maxScinID+0.5
  // ));
  // getStatistics().getHisto2D("hit_pos_per_scin")
  // ->GetXaxis()->SetTitle("Hit z position [cm]");
  // getStatistics().getHisto2D("hit_pos_per_scin")
  // ->GetYaxis()->SetTitle("ID of Scintillator");

  // Multiplicity of signals in Hits
  getStatistics().createHistogram(new TH1F("hit_sig_multi", "Number of signals from SiPMs in created hit", 11, -0.5, 10.5));
  getStatistics().getHisto1D("hit_sig_multi")->GetXaxis()->SetTitle("Number of signals");
  getStatistics().getHisto1D("hit_sig_multi")->GetYaxis()->SetTitle("Number of Hits");

  getStatistics().createHistogram(
      new TH1F("remain_signals_scin", "Number of Unused Signals in Scintillator", maxScinID - minScinID + 1, minScinID - 0.5, maxScinID + 0.5));
  getStatistics().getHisto1D("remain_signals_scin")->GetXaxis()->SetTitle("Scintillator ID");
  getStatistics().getHisto1D("remain_signals_scin")->GetYaxis()->SetTitle("Number of Unused Signals in Scintillator");

  getStatistics().createHistogram(
      new TH1F("remain_signals_tdiff", "Time Diff of an unused signal and the consecutive one", 200, fMaxABTimeDiff, 5.0 * fMaxABTimeDiff));
  getStatistics().getHisto1D("remain_signals_tdiff")->GetXaxis()->SetTitle("Time difference [ps]");
  getStatistics().getHisto1D("remain_signals_tdiff")->GetYaxis()->SetTitle("Number of Signals");

  for (int mtxID_i = minMtxID; mtxID_i <= maxMtxID; ++mtxID_i)
  {
    for (int mtxID_j = mtxID_i; mtxID_j <= maxMtxID; ++mtxID_j)
    {
      getStatistics().createHistogram(new TH1F(Form("mtx_%d_%d_time_diff", mtxID_i, mtxID_j),
                                               Form("Time difference of signals in matrices IDs %d and %d", mtxID_i, mtxID_j), 100, -20000.0,
                                               20000.0));
      getStatistics().getHisto1D(Form("mtx_%d_%d_time_diff", mtxID_i, mtxID_j))->GetXaxis()->SetTitle("time difference [ps]");
      getStatistics().getHisto1D(Form("mtx_%d_%d_time_diff", mtxID_i, mtxID_j))->GetYaxis()->SetTitle("Number of Hits");
    }
  }

  // Time diff and TOT per scin per multi
  // for (int scinID = minScinID; scinID <= maxScinID; scinID++)
  // {
  //
  //   getStatistics().createHistogram(
  //       new TH1F(Form("hit_sig_multi_scin_%d", scinID), Form("Number of signals from SiPMs in created hit for scin %d", scinID), 10, -0.5, 9.5));
  //   getStatistics().getHisto1D(Form("hit_sig_multi_scin_%d", scinID))->GetXaxis()->SetTitle("Number of signals");
  //   getStatistics().getHisto1D(Form("hit_sig_multi_scin_%d", scinID))->GetYaxis()->SetTitle("Number of Hits");
  //
  //   for (int multi = 2; multi <= 8; multi++)
  //   {
  //     getStatistics().createHistogram(new TH1F(Form("hit_tdiff_scin_%d_m_%d", scinID, multi),
  //                                              Form("Hit time difference, scin %d,  multiplicity %d", scinID, multi), 200, -1.1 * fMaxABTimeDiff,
  //                                              1.1 * fMaxABTimeDiff));
  //     getStatistics().getHisto1D(Form("hit_tdiff_scin_%d_m_%d", scinID, multi))->GetXaxis()->SetTitle("Time difference [ps]");
  //     getStatistics().getHisto1D(Form("hit_tdiff_scin_%d_m_%d", scinID, multi))->GetYaxis()->SetTitle("Number of Hits");
  //
  //     getStatistics().createHistogram(new TH1F(Form("hit_tot_scin_%d_m_%d", scinID, multi),
  //                                              Form("Hit TOT divided by multiplicity, scin %d multi %d", scinID, multi), 200, 0.0, 300000.0));
  //     getStatistics().getHisto1D(Form("hit_tot_scin_%d_m_%d", scinID, multi))->GetXaxis()->SetTitle("Time over Threshold [ps]");
  //     getStatistics().getHisto1D(Form("hit_tot_scin_%d_m_%d", scinID, multi))->GetYaxis()->SetTitle("Number of Hits");
  //   }
  // }

  // Time difference spectra - histograms upper limit is time window length
  // getStatistics().createHistogram(new TH1F("ab_tdiff_all", "All time differences of two A-B signals", 200, 0.0, 50000000.0));
  // getStatistics().getHisto1D("ab_tdiff_all")->GetXaxis()->SetTitle("Time difference [ps]");
  // getStatistics().getHisto1D("ab_tdiff_all")->GetYaxis()->SetTitle("Number of pairs");
  //
  // getStatistics().createHistogram(new TH1F("ab_tdiff_rej", "Time Diff of two signals not in coincidence", 200, 0.0, 50000000.0));
  // getStatistics().getHisto1D("ab_tdiff_rej")->GetXaxis()->SetTitle("Time difference [ps]");
  // getStatistics().getHisto1D("ab_tdiff_rej")->GetYaxis()->SetTitle("Number of pairs");
  //
  // getStatistics().createHistogram(new TH1F("ab_tdiff_acc", "Time Diff of two signals in coincidence", 200, fMinABTimeDiff, fMaxABTimeDiff));
  // getStatistics().getHisto1D("ab_tdiff_acc")->GetXaxis()->SetTitle("Time difference [ps]");
  // getStatistics().getHisto1D("ab_tdiff_acc")->GetYaxis()->SetTitle("Number of pairs");
  //
  // getStatistics().createHistogram(new TH1F("hit_wls_tdiff_all", "All time differences of hit and WLS signal", 200, 0.0, 50000000.0));
  // getStatistics().getHisto1D("hit_wls_tdiff_all")->GetXaxis()->SetTitle("Time difference [ps]");
  // getStatistics().getHisto1D("hit_wls_tdiff_all")->GetYaxis()->SetTitle("Number of pairs");
  //
  // getStatistics().createHistogram(new TH1F("hit_wls_tdiff_rej", "Time differences of hit and WLS signal not in conincidece", 200, 0.0,
  // 50000000.0)); getStatistics().getHisto1D("hit_wls_tdiff_rej")->GetXaxis()->SetTitle("Time difference [ps]");
  // getStatistics().getHisto1D("hit_wls_tdiff_rej")->GetYaxis()->SetTitle("Number of pairs");
  //
  // getStatistics().createHistogram(
  //     new TH1F("hit_wls_tdiff_acc", "Time differences of hit and WLS signal not in conincidece", 200, fMinABTimeDiff, fMaxABTimeDiff));
  // getStatistics().getHisto1D("hit_wls_tdiff_acc")->GetXaxis()->SetTitle("Time difference [ps]");
  // getStatistics().getHisto1D("hit_wls_tdiff_acc")->GetYaxis()->SetTitle("Number of pairs");

  // Unused sigals stats
  // getStatistics().createHistogram(new TH1F(
  //   "remain_signals_tdiff", "Time Diff of an unused signal and the consecutive one",
  //   200, 0.0, 5*fMaxABTimeDiff
  // ));
  // getStatistics().getHisto1D("remain_signals_tdiff")->GetXaxis()->SetTitle("Time difference [ps]");
  // getStatistics().getHisto1D("remain_signals_tdiff")->GetYaxis()->SetTitle("Number of Signals");

  // getStatistics().createHistogram(new TH1F(
  //   "remain_signals_per_scin", "Number of Unused Signals in Scintillator",
  //   maxScinID-minScinID+1, minScinID-0.5, maxScinID+0.5
  // ));
  // getStatistics().getHisto1D("remain_signals_per_scin")->GetXaxis()->SetTitle("ID of Scintillator");
  // getStatistics().getHisto1D("remain_signals_per_scin")->GetYaxis()->SetTitle("Number of Unused Signals in Scintillator");
}
