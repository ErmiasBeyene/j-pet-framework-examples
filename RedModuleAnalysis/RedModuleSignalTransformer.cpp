/**
 *  @copyright Copyright 2022 The J-PET Framework Authors. All rights reserved.
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
 *  @file RedModuleSignalTransformer.cpp
 */

#include "RedModuleSignalTransformer.h"
#include "RedModuleSignalTransformerTools.h"
#include <boost/property_tree/json_parser.hpp>

using namespace jpet_options_tools;

RedModuleSignalTransformer::RedModuleSignalTransformer(const char* name) : JPetUserTask(name) {}

RedModuleSignalTransformer::~RedModuleSignalTransformer() {}

bool RedModuleSignalTransformer::init()
{
  INFO("Signal Transformer started: PM to Matrix Signal");
  fOutputEvents = new JPetTimeWindow("JPetMatrixSignal");

  // Getting bools for saving control and calibration histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey))
  {
    fSaveControlHistos = getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
  }
  if (isOptionSet(fParams.getOptions(), kSaveCalibHistosParamKey))
  {
    fSaveCalibHistos = getOptionAsBool(fParams.getOptions(), kSaveCalibHistosParamKey);
  }

  // Reading file with Side B signals correction to property tree
  if (isOptionSet(fParams.getOptions(), kConstantsFileParamKey))
  {
    boost::property_tree::read_json(getOptionAsString(fParams.getOptions(), kConstantsFileParamKey), fConstansTree);
  }

  // Reading WLS config file
  if (isOptionSet(fParams.getOptions(), kWLSConfigFileParamKey))
  {
    boost::property_tree::read_json(getOptionAsString(fParams.getOptions(), kWLSConfigFileParamKey), fWLSConfigTree);
  }

  // Signal merging time parameter
  if (isOptionSet(fParams.getOptions(), kMergeSignalsTimeParamKey))
  {
    fMergingTime = getOptionAsDouble(fParams.getOptions(), kMergeSignalsTimeParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kMergeSignalsTimeParamKey.c_str(), fMergingTime));
  }

  // For plotting ToT histograms
  if (isOptionSet(fParams.getOptions(), kToTHistoUpperLimitParamKey))
  {
    fToTHistoUpperLimit = getOptionAsDouble(fParams.getOptions(), kToTHistoUpperLimitParamKey);
  }

  if (isOptionSet(fParams.getOptions(), kWLSSlotIDParamKey))
  {
    fWLSSlotID = getOptionAsInt(fParams.getOptions(), kWLSSlotIDParamKey);
  }
  INFO(Form("Using slot with ID %d as WLS set.", fWLSSlotID));

  // Control histograms
  if (fSaveControlHistos)
  {
    initialiseHistograms();
  }
  return true;
}

bool RedModuleSignalTransformer::exec()
{
  // Getting the data from event in an apropriate format
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent))
  {
    // Distribute PM Signals per Matrices
    auto pmSigMtxMap = RedModuleSignalTransformerTools::getPMSigMtxMap(timeWindow);

    if (fSaveCalibHistos)
    {
      RedModuleSignalTransformerTools::plotWLSSignalsTimeDiffs(pmSigMtxMap[JPetMatrix::WLS], getStatistics(), 401, 464);
    }

    // Merging max. 4 PM Signals into a MatrixSignal and separately signals on WLS SiPMs
    auto mergedSignals =
        RedModuleSignalTransformerTools::mergeSignalsAllSiPMs(pmSigMtxMap, fMergingTime, fConstansTree, fWLSConfigTree, getParamBank());

    // Saving method invocation
    if (mergedSignals.size() > 0)
    {
      saveMatrixSignals(mergedSignals);
    }
  }
  else
  {
    return false;
  }
  return true;
}

bool RedModuleSignalTransformer::terminate()
{
  INFO("Signal Transformer finished");
  return true;
}

/**
 * Save objects and make some histograms
 */
void RedModuleSignalTransformer::saveMatrixSignals(const std::vector<JPetMatrixSignal>& mtxSigVec)
{
  if (mtxSigVec.size() > 0 && fSaveControlHistos)
  {
    getStatistics().fillHistogram("mtxsig_tslot", mtxSigVec.size());
  }
  for (auto& mtxSig : mtxSigVec)
  {
    fOutputEvents->add<JPetMatrixSignal>(mtxSig);

    if (fSaveControlHistos)
    {
      auto scinID = mtxSig.getMatrix().getScin().getID();
      getStatistics().fillHistogram("mtxsig_multi", mtxSig.getPMSignals().size());
      if (mtxSig.getMatrix().getSide() == JPetMatrix::SideA)
      {
        getStatistics().fillHistogram("mtxsig_scin_sideA", scinID);
        getStatistics().fillHistogram("mtxsig_sideA_tot", scinID, mtxSig.getToT());
      }
      else if (mtxSig.getMatrix().getSide() == JPetMatrix::SideB)
      {
        getStatistics().fillHistogram("mtxsig_scin_sideB", scinID);
        getStatistics().fillHistogram("mtxsig_sideB_tot", scinID, mtxSig.getToT());
      }
      else if (mtxSig.getMatrix().getSide() == JPetMatrix::WLS)
      {
        getStatistics().fillHistogram("mtxsig_wls", scinID);
        getStatistics().fillHistogram("mtxsig_wls_multi", scinID, mtxSig.getPMSignals().size());

        getStatistics().fillHistogram("mtxsig_wls_tot", scinID, mtxSig.getToT());
        if (mtxSig.getPMSignals().size() == 1)
        {
          getStatistics().fillHistogram("mtxsig_wls_tot_multi1", scinID, mtxSig.getToT());
        }
        else if (mtxSig.getPMSignals().size() > 1)
        {
          getStatistics().fillHistogram("mtxsig_wls_tot_multi2p", scinID, mtxSig.getToT());
        }
      }
    }

    if (fSaveCalibHistos)
    {
      // Filling histograms gor each channel in Matrix SiPMs to produce
      // channel offsets with the respect to channel on 1st THR of SiPM mtx pos 1
      auto sigMap = mtxSig.getPMSignals();

      // Check only for matrices attached to scintillators
      if (mtxSig.getMatrix().getScin().getSlot().getType() == JPetSlot::Module)
      {
        if (sigMap.find(1) != sigMap.end())
        {
          auto t_1_1 = sigMap.at(1).getLeadTrailPairs().at(0).first.getTime();

          for (auto pmSig : sigMap)
          {
            auto pairs = pmSig.second.getLeadTrailPairs();
            for (auto pair : pairs)
            {
              auto t_ch_i = pair.first.getTime();
              auto channelID = pair.first.getChannel().getID();
              if (t_1_1 == t_ch_i)
              {
                continue;
              }
              getStatistics().fillHistogram("mtx_channel_offsets", channelID, t_ch_i - t_1_1);
            }
          }
        }
      }
      else if (mtxSig.getMatrix().getScin().getSlot().getType() == JPetSlot::WLS)
      {
        int prevSiPMID = 999;
        double t_prev_1 = -9999.0;

        for (auto pmSig : sigMap)
        {
          if (prevSiPMID == 999)
          {
            prevSiPMID = pmSig.second.getPM().getID();
            t_prev_1 = pmSig.second.getLeadTrailPairs().at(0).first.getTime();
            if (pmSig.second.getLeadTrailPairs().size() > 1)
            {
              auto t_prev_2 = pmSig.second.getLeadTrailPairs().at(1).first.getTime();
              auto channelID = pmSig.second.getLeadTrailPairs().at(1).first.getChannel().getID();
              getStatistics().fillHistogram("mtx_channel_offsets", channelID, t_prev_2 - t_prev_1);
            }
          }
          else
          {
            auto pairs = pmSig.second.getLeadTrailPairs();
            for (auto pair : pairs)
            {
              auto t_ch_i = pair.first.getTime();
              auto channelID = pair.first.getChannel().getID();
              getStatistics().fillHistogram("mtx_channel_offsets", channelID, t_ch_i - t_prev_1);
            }
          }
        }
      }
    }
  }
}

void RedModuleSignalTransformer::initialiseHistograms()
{
  auto minScinID = getParamBank().getScins().begin()->first;
  auto maxScinID = getParamBank().getScins().rbegin()->first;

  // MatrixSignal multiplicity
  getStatistics().createHistogramWithAxes(new TH1D("mtxsig_multi", "Multiplicity of matched MatrixSignals", 5, 0.5, 5.5),
                                          "Number of PM Signals in Matrix Signal", "Number of Matrix Signals");

  getStatistics().createHistogramWithAxes(new TH1D("mtxsig_tslot", "Number of Matrix Signals in Time Window", 100, 0.5, 100.5),
                                          "Number of Matrix Signals in Time Window", "Number of Time Windows");

  getStatistics().createHistogramWithAxes(
      new TH1D("mtxsig_scin_sideA", "Number of Matrix Signals per scintillator side A", maxScinID - minScinID + 1, minScinID - 0.5, maxScinID + 0.5),
      "Scin ID", "Number of Matrix Signals");

  getStatistics().createHistogramWithAxes(
      new TH1D("mtxsig_scin_sideB", "Number of Matrix Signals per scintillator side B", maxScinID - minScinID + 1, minScinID - 0.5, maxScinID + 0.5),
      "Scin ID", "Number of Matrix Signals");

  getStatistics().createHistogramWithAxes(new TH2D("mtxsig_sideA_tot", "Matrix Signal ToT - Side A per scintillator", maxScinID - minScinID + 1,
                                                   minScinID - 0.5, maxScinID + 0.5, 200, 0.0, fToTHistoUpperLimit),
                                          "Scin ID", "ToT [ps]");

  getStatistics().createHistogramWithAxes(new TH2D("mtxsig_sideB_tot", "Matrix Signal ToT - Side B per scintillator", maxScinID - minScinID + 1,
                                                   minScinID - 0.5, maxScinID + 0.5, 200, 0.0, fToTHistoUpperLimit),
                                          "Scin ID", "ToT [ps]");

  getStatistics().createHistogramWithAxes(
      new TH1D("mtxsig_wls", "Number of Matrix Signals per WLS", maxScinID - minScinID + 1, minScinID - 0.5, maxScinID + 0.5), "Scin ID",
      "Number of Matrix Signals");

  getStatistics().createHistogramWithAxes(new TH2D("mtxsig_wls_tot", "Matrix Signal ToT - WLS layer", maxScinID - minScinID + 1, minScinID - 0.5,
                                                   maxScinID + 0.5, 200, 0.0, fToTHistoUpperLimit),
                                          "Scin ID", "ToT [ps]");

  getStatistics().createHistogramWithAxes(new TH2D("mtxsig_wls_tot_multi1", "Matrix Signal ToT - WLS layer", maxScinID - minScinID + 1,
                                                   minScinID - 0.5, maxScinID + 0.5, 200, 0.0, fToTHistoUpperLimit),
                                          "Scin ID", "ToT [ps]");

  getStatistics().createHistogramWithAxes(new TH2D("mtxsig_wls_tot_multi2p", "Matrix Signal ToT - WLS layer", maxScinID - minScinID + 1,
                                                   minScinID - 0.5, maxScinID + 0.5, 200, 0.0, fToTHistoUpperLimit),
                                          "Scin ID", "ToT [ps]");

  getStatistics().createHistogramWithAxes(
      new TH2D("mtxsig_wls_multi", "WLS Matrix Signal Multiplicity", maxScinID - minScinID + 1, minScinID - 0.5, maxScinID + 0.5, 3, 0.5, 3.5),
      "Scin ID", "Number of PM signals merged into WLS signal");

  // SiPM calibrations
  if (fSaveCalibHistos)
  {
    auto minChannelID = getParamBank().getChannels().begin()->first;
    auto maxChannelID = getParamBank().getChannels().rbegin()->first;

    getStatistics().createHistogramWithAxes(new TH2D("mtx_channel_offsets", "Offset of Channel in Matrix vs. Channel ID",
                                                     maxChannelID - minChannelID + 1, minChannelID - 0.5, maxChannelID + 0.5, 200, -fMergingTime,
                                                     fMergingTime),
                                            "Channel ID", "Offset");

    auto minSiPMID = 401;
    auto maxSiPMID = 464;

    getStatistics().createHistogramWithAxes(new TH2D("wls_sipm_calib", "Time differences between consecutive SiPM signals in WLS layer",
                                                     maxSiPMID - minSiPMID + 1, minSiPMID - 0.5, maxSiPMID + 0.5, 200, -fMergingTime, fMergingTime),
                                            "SiPM ID", "Signal time difference");
  }
}
