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
  if (isOptionSet(fParams.getOptions(), kScinSyncFileParamKey)) {
    auto scinSyncFileName = getOptionAsString(fParams.getOptions(), kScinSyncFileParamKey);
    pt::read_json(scinSyncFileName, fScinSyncTree);
  }

  if (isOptionSet(fParams.getOptions(), kMinScinIDParamKey)) {
    fMinScinID = getOptionAsInt(fParams.getOptions(), kMinScinIDParamKey);
  } else {
    fMinScinID = getParamBank().getScins().begin()->first;
  }

  if (isOptionSet(fParams.getOptions(), kMaxScinIDParamKey)) {
    fMaxScinID = getOptionAsInt(fParams.getOptions(), kMaxScinIDParamKey);
  } else {
    fMaxScinID = getParamBank().getScins().rbegin()->first;
  }

  if (isOptionSet(fParams.getOptions(), kMinPMIDParamKey)) {
    fMinPMID = getOptionAsInt(fParams.getOptions(), kMinPMIDParamKey);
  } else {
    fMinPMID = getParamBank().getPMs().begin()->first;
  }

  if (isOptionSet(fParams.getOptions(), kMaxPMIDParamKey)) {
    fMaxPMID = getOptionAsInt(fParams.getOptions(), kMaxPMIDParamKey);
  } else {
    fMaxPMID = getParamBank().getPMs().rbegin()->first;
  }

  // Control histograms
  if(fSaveControlHistos) { initialiseHistograms(); }
  return true;
}

bool SignalTransformer::exec()
{
  // Getting the data from event in an apropriate format
  if(auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

    // Distribute Raw Signals per Matrices
    auto rawSigMtxWLSMap = SignalTransformerTools::getRawSigMtxWLSMap(timeWindow);

    // Merging max. 4 Raw Signals into a MatrixSignal
    auto mergedScinSignals = SignalTransformerTools::mergeScinSiPMSignals(
      rawSigMtxWLSMap.first, fMergingTime, fScinSyncTree, getStatistics(), fSaveControlHistos
    );

    auto mergedWLSSignals = SignalTransformerTools::mergeWLSSiPMSignals(
      getParamBank().getWLSs(), rawSigMtxWLSMap.second, fMergingTime,
      getStatistics(), fSaveControlHistos
    );

    mergedScinSignals.insert(mergedScinSignals.end(), mergedWLSSignals.begin(), mergedWLSSignals.end());

    if (mergedScinSignals.size()>0) { saveMatrixSignals(mergedScinSignals); }

  } else { return false; }
  return true;
}

bool SignalTransformer::terminate()
{
  if(fSaveControlHistos){
    namespace pt = boost::property_tree;
    using namespace std;

    if (isOptionSet(fParams.getOptions(), kOutputOffsetsFileParamKey)) {
      INFO("Signal transforming - printing out offsets for SiPMs in matrices");
      auto outputOffsetsFile = getOptionAsString(fParams.getOptions(), kOutputOffsetsFileParamKey);

      pt::ptree root;
      pt::ptree sipm_node;

      for(int pmID=fMinPMID; pmID<=fMaxPMID; pmID++){
        auto mean = getStatistics().getHisto1D(Form("offset_sipm_%d", pmID))->GetMean();
        sipm_node.put(to_string(pmID), mean);
      }
      root.add_child("sipm_offsets", sipm_node);

      // Merging used calibration with new one - iteration alike
      if (isOptionSet(fParams.getOptions(), kInputOffsetsFileParamKey)) {

        auto inputOffsetsFile = getOptionAsString(fParams.getOptions(), kInputOffsetsFileParamKey);

        pt::ptree rootOld;
        pt::read_json(inputOffsetsFile, rootOld);

        pt::ptree new_root;
        pt::ptree new_sipm_node;

        for(int pmID=fMinPMID; pmID<=fMaxPMID; pmID++){
          double oldOffset = rootOld.get("sipm_offsets."+to_string(pmID), 0.0);
          double newOffset = root.get("sipm_offsets."+to_string(pmID), 0.0);
          new_sipm_node.put(to_string(pmID), oldOffset+newOffset);
        }

        new_root.add_child("sipm_offsets", new_sipm_node);
        pt::write_json(outputOffsetsFile, new_root);

      } else {
        pt::write_json(outputOffsetsFile, root);
      }
    }
  }

  INFO("Signal transforming finished");
  return true;
}

/**
* Save objects and make some histograms
*/
void SignalTransformer::saveMatrixSignals(const std::vector<JPetMatrixSignal>& mtxSigVec)
{
  if(fSaveControlHistos){
    getStatistics().getHisto1D("mtxsig_tslot")->Fill(mtxSigVec.size());
  }

  for (auto& mtxSig : mtxSigVec) { fOutputEvents->add<JPetMatrixSignal>(mtxSig); }
}

void SignalTransformer::initialiseHistograms()
{
  // MatrixSignal multiplicity
  getStatistics().createHistogram(new TH1F(
    "mtxsig_multi", "Multiplicity of MatrixSignals matched in scintillators", 5, 0.5, 5.5
  ));
  getStatistics().getHisto1D("mtxsig_multi")->GetXaxis()->SetTitle("Number of Raw Signals in Matrix Signal");
  getStatistics().getHisto1D("mtxsig_multi")->GetYaxis()->SetTitle("Number of Matrix Signals");

  getStatistics().createHistogram(new TH1F(
    "mtxsig_tslot", "Number of Matrix Signals in Time Window", 50, -0.5, 49.5
  ));
  getStatistics().getHisto1D("mtxsig_tslot")->GetXaxis()->SetTitle("Number of Matrix Signal in Time Window");
  getStatistics().getHisto1D("mtxsig_tslot")->GetYaxis()->SetTitle("Number of Time Windows");

  getStatistics().createHistogram(new TH1F(
    "mtxsig_per_scin_sideA", "Number of MatrixSignals per scintillator side A",
    fMaxScinID-fMinScinID+1, fMinScinID-0.5, fMaxScinID+0.5
  ));
  getStatistics().getHisto1D("mtxsig_per_scin_sideA")->GetXaxis()->SetTitle("Scin ID");
  getStatistics().getHisto1D("mtxsig_per_scin_sideA")->GetYaxis()->SetTitle("Number of Matrix Signals");

  getStatistics().createHistogram(new TH1F(
    "mtxsig_per_scin_sideB", "Number of MatrixSignals per scintillator side B",
    fMaxScinID-fMinScinID+1, fMinScinID-0.5, fMaxScinID+0.5
  ));
  getStatistics().getHisto1D("mtxsig_per_scin_sideB")->GetXaxis()->SetTitle("Scin ID");
  getStatistics().getHisto1D("mtxsig_per_scin_sideB")->GetYaxis()->SetTitle("Number of Matrix Signals");

  // Same offsets but with SiPM ID in the name
  for(int pmID=fMinPMID; pmID<=fMaxPMID; pmID++){
    getStatistics().createHistogram(new TH1F(
      Form("offset_sipm_%d", pmID), Form("Offset for SiPM %d", pmID),
      200, -1.1*fMergingTime, 1.1*fMergingTime
    ));
    getStatistics().getHisto1D(Form("offset_sipm_%d", pmID))->GetXaxis()->SetTitle("Time difference [ps]");
    getStatistics().getHisto1D(Form("offset_sipm_%d", pmID))->GetYaxis()->SetTitle("Number of Raw Signal pairs");
  }

  // For each WLS
  for(auto wlsEle : getParamBank().getWLSs()){
    auto wlsID = wlsEle.second->getID();

    getStatistics().createHistogram(new TH1F(
      Form("wls_%d_tot", wlsID), Form("ToT on WLS ID %d", wlsID),
      200, 0.0, 200000.0
    ));
    getStatistics().getHisto1D(Form("wls_%d_tot", wlsID))->GetXaxis()->SetTitle("TOT [ps]");
    getStatistics().getHisto1D(Form("wls_%d_tot", wlsID))->GetYaxis()->SetTitle("Number of Mtx Signals");

    for(auto pmID : wlsEle.second->getPMIDs()) {
      if(pmID!=-1){
        getStatistics().createHistogram(new TH1F(
          Form("wls_%d_sipm_%d_tot", wlsID, pmID),
          Form("ToT of signals SiPM %d on WLS ID %d", pmID, wlsID),
          200, 0.0, 200000.0
        ));
        getStatistics().getHisto1D(Form("wls_%d_sipm_%d_tot", wlsID, pmID))->GetXaxis()->SetTitle("TOT [ps]");
        getStatistics().getHisto1D(Form("wls_%d_sipm_%d_tot", wlsID, pmID))->GetYaxis()->SetTitle("Number of Signals");
      }
    }
  }
}
