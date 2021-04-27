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
 *  @file threshold_offests.C
 *
 *  @brief Script for reading histograms with threshold offsets and producing calibraiton json file
 *
 *  This script uses histograms, that are produced by task SignalFinder:
 *  "thr_tdiff_2_1_pm", "thr_tdiff_3_1_pm" and "thr_tdiff_4_1_pm".
 *  Channels that belong to the same PM are synchronized to THR1 channel based on
 *  time differences between leading channel signals THR_i-THR1.
 *
 *  Basic usage:
 *  root> .L threshold_offests.C
 *  root> threshold_offests("file_with_calib_histos.root")
 *  -- this will produce file "calibration_constants.json" with the results. If the
 *  file exists, the result of this calibration will be appended to the existing tree.
 */

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TMath.h>

#include <fstream>
#include <iostream>
#include <vector>

namespace bpt = boost::property_tree;

const int kNumberOfThresholds = 4;

void threshold_offests(std::string fileName, std::string calibJSONFileName = "calibration_constants.json", bool saveResult = false,
                       std::string resultDir = "./", int minPMID = 1, int maxPMID = 384)
{
  TFile* inputFile = new TFile(fileName.c_str(), "READ");

  bpt::ptree tree;
  ifstream calibFile(calibJSONFileName.c_str());
  if (calibFile.good())
  {
    bpt::read_json(calibJSONFileName, tree);
  }

  if (inputFile->IsOpen())
  {
    for (int thr = 2; thr <= kNumberOfThresholds; ++thr)
    {
      TH2D* thrTimeDiffs = dynamic_cast<TH2D*>(inputFile->Get(Form("thr_tdiff_%d_1_pm", thr)));

      for (int pmID = minPMID; pmID <= maxPMID; ++pmID)
      {
        TH1D* offsetHist = thrTimeDiffs->ProjectionY(Form("offset_pm_%d", pmID), pmID - minPMID + 1, pmID - minPMID + 1);
        offsetHist->SetLineWidth(2);
        offsetHist->SetLineColor(kBlue);

        if (offsetHist->GetEntries() < 100)
        {
          continue;
        }

        // Offset is a time indicated by bin with highest number of counts
        double offset = offsetHist->GetBinCenter(offsetHist->GetMaximumBin());
        tree.put("pm." + to_string(pmID) + ".offset_thr_" + to_string(thr), offset);

        if (saveResult)
        {
          auto name = Form("offset_pm_%d_thr_%d", pmID, thr);
          TCanvas* can = new TCanvas(name, name, 900, 720);
          offsetHist->Draw();
          TLine* line = new TLine(offset, offsetHist->GetMinimum(), offset, offsetHist->GetMaximum());
          line->SetLineWidth(2);
          line->SetLineColor(kRed);
          line->Draw("same");
          can->SaveAs(Form("%s/offset_pm_%d_thr_%d.png", resultDir.c_str(), pmID, thr));
        }
      }
    }
  }

  // Saving tree into json file
  bpt::write_json(calibJSONFileName, tree);
}
