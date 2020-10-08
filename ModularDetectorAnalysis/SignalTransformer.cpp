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

#include "SignalTransformerTools.h"
#include "JPetWriter/JPetWriter.h"
#include "SignalTransformer.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

using namespace jpet_options_tools;
using namespace std;
namespace pt = boost::property_tree;

SignalTransformer::SignalTransformer(const char* name): JPetUserTask(name) {}

SignalTransformer::~SignalTransformer() {}

bool SignalTransformer::init()
{
  INFO("Signal transforming started: Raw to Matrix Signal");
  fOutputEvents = new JPetTimeWindow("JPetMatrixSignal");

  // Getting bool for saving histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey)) {
    fSaveControlHistos = getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
  }

  // Signal merging time parameter
  if (isOptionSet(fParams.getOptions(), kMergeSignalsTimeParamKey)) {
    fMergingTime = getOptionAsDouble(fParams.getOptions(), kMergeSignalsTimeParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kMergeSignalsTimeParamKey.c_str(), fMergingTime
    ));
  }

  // Reading file with offsets to property tree - scintillator synchronization
  // if (isOptionSet(fParams.getOptions(), kScinSyncFileParamKey)) {
  //   auto scinSyncFileName = getOptionAsString(fParams.getOptions(), kScinSyncFileParamKey);
  //   pt::read_json(scinSyncFileName, fScinSyncTree);
  // }

  // if (isOptionSet(fParams.getOptions(), kMinScinIDParamKey)) {
  //   fMinScinID = getOptionAsInt(fParams.getOptions(), kMinScinIDParamKey);
  // } else {
  //   fMinScinID = getParamBank().getScins().begin()->first;
  // }
  //
  // if (isOptionSet(fParams.getOptions(), kMaxScinIDParamKey)) {
  //   fMaxScinID = getOptionAsInt(fParams.getOptions(), kMaxScinIDParamKey);
  // } else {
  //   fMaxScinID = getParamBank().getScins().rbegin()->first;
  // }

  // if (isOptionSet(fParams.getOptions(), kMinPMIDParamKey)) {
  //   fMinPMID = getOptionAsInt(fParams.getOptions(), kMinPMIDParamKey);
  // } else {
  //   fMinPMID = getParamBank().getPMs().begin()->first;
  // }
  //
  // if (isOptionSet(fParams.getOptions(), kMaxPMIDParamKey)) {
  //   fMaxPMID = getOptionAsInt(fParams.getOptions(), kMaxPMIDParamKey);
  // } else {
  //   fMaxPMID = getParamBank().getPMs().rbegin()->first;
  // }

  // Control histograms
  if(fSaveControlHistos) { initialiseHistograms(); }
  return true;
}

bool SignalTransformer::exec()
{
  // Getting the data from event in an apropriate format
  if(auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

    // Distribute Raw Signals per Matrices
    auto rawSigPMMap = SignalTransformerTools::getRawSigPMMap(timeWindow);

    auto mergedMtxSignals = SignalTransformerTools::mergeSignalsAllMtx(
      getParamBank(), rawSigPMMap, fMergingTime, getStatistics(), fSaveControlHistos
    );

    saveMatrixSignals(mergedMtxSignals);

  } else { return false; }
  return true;
}

bool SignalTransformer::terminate()
{
  // if(fSaveControlHistos){
  //   namespace pt = boost::property_tree;
  //   using namespace std;
  //
  //   if (isOptionSet(fParams.getOptions(), kOutputOffsetsFileParamKey)) {
  //     INFO("Signal transforming - printing out offsets for SiPMs in matrices");
  //     auto outputOffsetsFile = getOptionAsString(fParams.getOptions(), kOutputOffsetsFileParamKey);
  //
  //     pt::ptree root;
  //     pt::ptree sipm_node;
  //
  //     for(int pmID=fMinPMID; pmID<=fMaxPMID; pmID++){
  //       auto mean = getStatistics().getHisto1D(Form("offset_sipm_%d", pmID))->GetMean();
  //       sipm_node.put(to_string(pmID), mean);
  //     }
  //     root.add_child("sipm_offsets", sipm_node);
  //
  //     // Merging used calibration with new one - iteration alike
  //     if (isOptionSet(fParams.getOptions(), kInputOffsetsFileParamKey)) {
  //
  //       auto inputOffsetsFile = getOptionAsString(fParams.getOptions(), kInputOffsetsFileParamKey);
  //
  //       pt::ptree rootOld;
  //       pt::read_json(inputOffsetsFile, rootOld);
  //
  //       pt::ptree new_root;
  //       pt::ptree new_sipm_node;
  //
  //       for(int pmID=fMinPMID; pmID<=fMaxPMID; pmID++){
  //         double oldOffset = rootOld.get("sipm_offsets."+to_string(pmID), 0.0);
  //         double newOffset = root.get("sipm_offsets."+to_string(pmID), 0.0);
  //         new_sipm_node.put(to_string(pmID), oldOffset+newOffset);
  //       }
  //
  //       new_root.add_child("sipm_offsets", new_sipm_node);
  //       pt::write_json(outputOffsetsFile, new_root);
  //
  //     } else {
  //       pt::write_json(outputOffsetsFile, root);
  //     }
  //   }
  // }

  INFO("Signal transforming finished");
  return true;
}

/**
* Save objects and make some histograms
*/
void SignalTransformer::saveMatrixSignals(const std::vector<JPetMatrixSignal>& mtxSigVec)
{
  if(fSaveControlHistos && mtxSigVec.size() > 0){
    getStatistics().getHisto1D("mtxsig_tslot")->Fill(mtxSigVec.size());
  }

  for (auto& mtxSig : mtxSigVec) {
    fOutputEvents->add<JPetMatrixSignal>(mtxSig);
    if(fSaveControlHistos){
      getStatistics().getHisto1D("mtxsig_per_matrix")->Fill(mtxSig.getMatrix().getID());
    }
  }
}

void SignalTransformer::initialiseHistograms()
{
  getStatistics().createHistogram(new TH1F(
    "mtxsig_tslot", "Number of Matrix Signals in Time Window", 50, 0.5, 51.5
  ));
  getStatistics().getHisto1D("mtxsig_tslot")->GetXaxis()->SetTitle("Number of Matrix Signal in Time Window");
  getStatistics().getHisto1D("mtxsig_tslot")->GetYaxis()->SetTitle("Number of Time Windows");

  auto minMatrixID = getParamBank().getMatrices().begin()->first;
  auto maxMatrixID = getParamBank().getMatrices().rbegin()->first;

  getStatistics().createHistogram(new TH1F(
    "mtxsig_per_matrix", "Number of Matrix Signals per matrix",
    maxMatrixID-minMatrixID+1, minMatrixID-0.5, maxMatrixID+0.5
  ));
  getStatistics().getHisto1D("mtxsig_per_matrix")->GetXaxis()->SetTitle("Matrix ID");
  getStatistics().getHisto1D("mtxsig_per_matrix")->GetYaxis()->SetTitle("Number of Matrix Signals");

  // Same offsets but with SiPM ID in the name
  // for(int pmID=fMinPMID; pmID<=fMaxPMID; pmID++){
  //   getStatistics().createHistogram(new TH1F(
  //     Form("offset_sipm_%d", pmID), Form("Offset for SiPM %d", pmID),
  //     200, -1.1*fMergingTime, 1.1*fMergingTime
  //   ));
  //   getStatistics().getHisto1D(Form("offset_sipm_%d", pmID))->GetXaxis()->SetTitle("Time difference [ps]");
  //   getStatistics().getHisto1D(Form("offset_sipm_%d", pmID))->GetYaxis()->SetTitle("Number of Raw Signal pairs");
  // }

  // WLS
  auto minWLSID = getParamBank().getWLSs().begin()->first;
  auto maxWLSID = getParamBank().getWLSs().rbegin()->first;

  getStatistics().createHistogram(new TH1F(
    "wls_sig_occ", "WLS occupancy", maxWLSID-minWLSID+1, minWLSID-0.5, maxWLSID+0.5
  ));
  getStatistics().getHisto1D("wls_sig_occ")->GetXaxis()->SetTitle("WLS ID");
  getStatistics().getHisto1D("wls_sig_occ")->GetYaxis()->SetTitle("Number of Matrix signals");

  getStatistics().createHistogram(new TH1F(
    "wls_sig_z_pos", "WLS matrix signal position based on TOT", 200, -25.0, 25.0
  ));
  getStatistics().getHisto1D("wls_sig_z_pos")->GetXaxis()->SetTitle("z [cm]");
  getStatistics().getHisto1D("wls_sig_z_pos")->GetYaxis()->SetTitle("Number of Matrix signals");

  for(auto mtx_p : getParamBank().getMatrices()){
    auto mtxID = mtx_p.second->getID();

    if(mtx_p.second->getType()=="WLS") {
      auto wlsID =  mtx_p.second->getWLS().getID();
      getStatistics().createHistogram(new TH1F(
        Form("wls_%d_tot", wlsID), Form("Average ToT on WLS ID %d", wlsID),
        200, 0.0, 400000.0
      ));
      getStatistics().getHisto1D(Form("wls_%d_tot", wlsID))->GetXaxis()->SetTitle("TOT [ps]");
      getStatistics().getHisto1D(Form("wls_%d_tot", wlsID))->GetYaxis()->SetTitle("Number of Mtx Signals");

      for(auto pmID : mtx_p.second->getPMIDs()) {
        if(pmID!=-1){
          getStatistics().createHistogram(new TH1F(
            Form("wls_%d_sipm_%d_tot", wlsID, pmID),
            Form("ToT of signals SiPM %d on WLS ID %d", pmID, wlsID),
            200, 0.0, 400000.0
          ));
          getStatistics().getHisto1D(Form("wls_%d_sipm_%d_tot", wlsID, pmID))->GetXaxis()->SetTitle("TOT [ps]");
          getStatistics().getHisto1D(Form("wls_%d_sipm_%d_tot", wlsID, pmID))->GetYaxis()->SetTitle("Number of Signals");
        }
      }
    } else {
      auto scinID = mtx_p.second->getScin().getID();
      auto hName = Form("scin_%d_%s_tot", scinID, mtx_p.second->getType().c_str());
      getStatistics().createHistogram(new TH1F(
        hName, Form("Average ToT on Scin ID %d %s", scinID, mtx_p.second->getType().c_str()),
        200, 0.0, 400000.0
      ));
      getStatistics().getHisto1D(hName)->GetXaxis()->SetTitle("TOT [ps]");
      getStatistics().getHisto1D(hName)->GetYaxis()->SetTitle("Number of Mtx Signals");
    }
  }
}
