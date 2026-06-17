#ifndef CREATEDCODE_HISTOGRAMMAKERS_HELPERS_TWOELECTRONSELECTEDPARTICLEPLOTS_HH
#define CREATEDCODE_HISTOGRAMMAKERS_HELPERS_TWOELECTRONSELECTEDPARTICLEPLOTS_HH

//----------------------------------------------------------------------------------
//
// twoElectronSelectedParticlePlots.hh
// Written by Andrew Boldy University of South Carolina, 2026
// Assisted by Codex
//
// Purpose:
//   Read the histogram ROOT file written by the two-electron truth/reco
//   comparer and save diagnostic PDF plots into a dedicated output directory.
//
//   The helper deliberately keeps plotting separate from the histogram-filling
//   code.  The comparer writes a ROOT file first; this helper opens that file
//   afterwards and produces human-readable PDFs under categorized
//   subdirectories in:
//
//       Plots/TruthVsRecoPlots/
//
//   That keeps the event loop simple and makes it easy to regenerate plots
//   without rerunning the full ntuple scan once the histogram file exists.
//
//----------------------------------------------------------------------------------

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TAxis.h>
#include <TBox.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TSystem.h>

namespace twoelectronplots
{
  // These bounds match the stopping-target geometry overlay used by the
  // original vertexer macro.
  constexpr double kStoppingTargetOuterRadius = 75.0;
  constexpr double kStoppingTargetZMin = -4700.053;
  constexpr double kStoppingTargetZMax = -3899.947;
  constexpr double kFoilPointTruthDistancePlotMin = 0.0;
  constexpr double kFoilPointTruthDistancePlotMax = 1000.0;
  constexpr double kAverageAbsLineParameterPlotMin = 0.0;
  constexpr double kAverageAbsLineParameterPlotMax = 1200.0;
  constexpr double kRecoTruthDeltaZPlotMin = -1000.0;
  constexpr double kRecoTruthDeltaZPlotMax = 1000.0;
  constexpr double kSharedFoilIndexPlotMin = -0.5;
  constexpr double kSharedFoilIndexPlotMax = 36.5;
  constexpr double kFoilCountPlotMin = -0.5;
  constexpr double kFoilCountPlotMax = 37.5;
  constexpr double kOpeningAngleDegreesPlotMin = 0.0;
  constexpr double kOpeningAngleDegreesPlotMax = 180.0;
  constexpr int kSelectedFoilCount = 37;
  constexpr int kDeltaZByFoilPlotsPerCanvas = 9;

  inline std::string selectedFoilDeltaZHistogramName(int foilIndex)
  {
    return "hRecoSpaceSelectedDeltaZSelectedFoil" +
           std::to_string(foilIndex);
  }

  inline void ensureDirectoryPath(const std::string& directoryPath)
  {
    if (!directoryPath.empty() && gSystem != nullptr)
    {
      gSystem->mkdir(directoryPath.c_str(), true);
    }
  }

  inline void ensureOutputDirectory(const std::string& outputPath)
  {
    const std::string::size_type slashPosition = outputPath.find_last_of("/\\");
    if (slashPosition == std::string::npos)
    {
      return;
    }

    const std::string outputDirectory = outputPath.substr(0, slashPosition);
    ensureDirectoryPath(outputDirectory);
  }

  inline TH1F* get1DHistogram(TFile& file, const std::string& name)
  {
    return dynamic_cast<TH1F*>(file.Get(name.c_str()));
  }

  inline TH2F* get2DHistogram(TFile& file, const std::string& name)
  {
    return dynamic_cast<TH2F*>(file.Get(name.c_str()));
  }

  inline TGraph* getGraph(TFile& file, const std::string& name)
  {
    return dynamic_cast<TGraph*>(file.Get(name.c_str()));
  }

  inline void style1DHistogram(TH1* histogram)
  {
    if (histogram == nullptr)
    {
      return;
    }

    histogram->SetLineWidth(2);
  }

  inline void draw2DHistogram(TH2* histogram)
  {
    if (histogram == nullptr)
    {
      return;
    }

    if (gPad != nullptr)
    {
      gPad->SetRightMargin(0.15);
      gPad->SetGridx();
      gPad->SetGridy();
    }
    histogram->SetStats(false);
    histogram->Draw("COLZ");
  }

  inline void styleScatterGraph(TGraph* graph, Color_t color)
  {
    if (graph == nullptr)
    {
      return;
    }

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.65);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
  }

  inline void drawScatterGraph(TGraph* graph,
                               double xMin,
                               double xMax,
                               double yMin,
                               double yMax)
  {
    if (graph == nullptr)
    {
      return;
    }

    if (gPad != nullptr)
    {
      gPad->SetGridx();
      gPad->SetGridy();
    }

    graph->SetMinimum(yMin);
    graph->SetMaximum(yMax);
    graph->Draw("AP");
    if (graph->GetXaxis() != nullptr)
    {
      graph->GetXaxis()->SetLimits(xMin, xMax);
    }
    if (gPad != nullptr)
    {
      gPad->Modified();
      gPad->Update();
    }
  }

  inline void drawDeltaZBySelectedFoilHistogram(TH1* histogram, int foilIndex)
  {
    if (histogram == nullptr)
    {
      return;
    }

    if (gPad != nullptr)
    {
      gPad->SetGridx();
      gPad->SetGridy();
    }

    const Long64_t entries =
      static_cast<Long64_t>(histogram->GetEntries());
    const std::string title =
      "selected foil sindex " + std::to_string(foilIndex) +
      ", N = " + std::to_string(entries);
    histogram->SetTitle(title.c_str());
    histogram->SetLineColor(kBlue + 1);
    histogram->SetLineWidth(2);
    histogram->SetStats(false);
    histogram->Draw("HIST E");
  }

  inline void saveDeltaZBySelectedFoilCanvases(
    const std::vector<TH1F*>& deltaZBySelectedFoilHistograms,
    const std::string& outputDirectory)
  {
    std::vector<int> selectedFoilIndices;
    selectedFoilIndices.reserve(deltaZBySelectedFoilHistograms.size());
    for (size_t foilIndex = 0;
         foilIndex < deltaZBySelectedFoilHistograms.size();
         ++foilIndex)
    {
      TH1F* histogram = deltaZBySelectedFoilHistograms.at(foilIndex);
      if (histogram != nullptr && histogram->GetEntries() > 0.0)
      {
        selectedFoilIndices.push_back(static_cast<int>(foilIndex));
      }
    }

    for (size_t firstFoilInCanvas = 0;
         firstFoilInCanvas < selectedFoilIndices.size();
         firstFoilInCanvas += kDeltaZByFoilPlotsPerCanvas)
    {
      const int pageNumber =
        static_cast<int>(firstFoilInCanvas / kDeltaZByFoilPlotsPerCanvas) + 1;
      const std::string canvasName =
        "cRecoSpaceSelectedDeltaZBySelectedFoilPage" +
        std::to_string(pageNumber);
      const std::string canvasTitle =
        "Space-selected #Delta z by selected foil, page " +
        std::to_string(pageNumber);
      TCanvas canvas(canvasName.c_str(), canvasTitle.c_str(), 1800, 1400);
      canvas.Divide(3, 3);

      for (int padIndex = 0;
           padIndex < kDeltaZByFoilPlotsPerCanvas;
           ++padIndex)
      {
        const size_t selectedIndex =
          firstFoilInCanvas + static_cast<size_t>(padIndex);
        if (selectedIndex >= selectedFoilIndices.size())
        {
          break;
        }

        const int foilIndex = selectedFoilIndices.at(selectedIndex);
        canvas.cd(padIndex + 1);
        drawDeltaZBySelectedFoilHistogram(
          deltaZBySelectedFoilHistograms.at(foilIndex),
          foilIndex);
      }

      canvas.SaveAs(
        (outputDirectory + "/RecoSpaceSelectedDeltaZBySelectedFoil_Page" +
         std::to_string(pageNumber) + ".pdf").c_str());
    }
  }

  inline void drawStoppingTargetBoxXY()
  {
    TBox stoppingTargetBox(-kStoppingTargetOuterRadius,
                           -kStoppingTargetOuterRadius,
                           kStoppingTargetOuterRadius,
                           kStoppingTargetOuterRadius);
    stoppingTargetBox.SetFillStyle(0);
    stoppingTargetBox.SetLineColor(kRed + 1);
    stoppingTargetBox.SetLineWidth(3);
    stoppingTargetBox.DrawClone("same");
  }

  inline void drawStoppingTargetBoxXZ()
  {
    TBox stoppingTargetBox(-kStoppingTargetOuterRadius,
                           kStoppingTargetZMin,
                           kStoppingTargetOuterRadius,
                           kStoppingTargetZMax);
    stoppingTargetBox.SetFillStyle(0);
    stoppingTargetBox.SetLineColor(kRed + 1);
    stoppingTargetBox.SetLineWidth(3);
    stoppingTargetBox.DrawClone("same");
  }

  inline void drawStoppingTargetBoxYZ()
  {
    TBox stoppingTargetBox(-kStoppingTargetOuterRadius,
                           kStoppingTargetZMin,
                           kStoppingTargetOuterRadius,
                           kStoppingTargetZMax);
    stoppingTargetBox.SetFillStyle(0);
    stoppingTargetBox.SetLineColor(kRed + 1);
    stoppingTargetBox.SetLineWidth(3);
    stoppingTargetBox.DrawClone("same");
  }

  inline bool saveTruthVsRecoPlotsFromFile(const std::string& histogramFileName,
                                           const std::string& outputDirectory =
                                             "Plots/TruthVsRecoPlots")
  {
    TFile* histogramFile = TFile::Open(histogramFileName.c_str(), "READ");
    if (histogramFile == nullptr || histogramFile->IsZombie())
    {
      std::cerr << "ERROR: could not open histogram file: "
                << histogramFileName << std::endl;
      return false;
    }

    ensureDirectoryPath(outputDirectory);
    const std::string monteCarloTruthDirectory =
      outputDirectory + "/MonteCarloTruth";
    const std::string rawRecoDirectory =
      outputDirectory + "/RawReco";
    const std::string anglesDirectory =
      outputDirectory + "/Angles";
    const std::string deltaZVsFoilDirectory =
      outputDirectory + "/DeltaZVsFoil";
    ensureDirectoryPath(monteCarloTruthDirectory);
    ensureDirectoryPath(rawRecoDirectory);
    ensureDirectoryPath(anglesDirectory);
    ensureDirectoryPath(deltaZVsFoilDirectory);

    TH1F* hMCTruthOriginT = get1DHistogram(*histogramFile,
                                           "hMCTruthRank0DownstreamElectronOriginT");
    TH1F* hMCTruthOriginX = get1DHistogram(*histogramFile,
                                           "hMCTruthRank0DownstreamElectronOriginX");
    TH1F* hMCTruthOriginY = get1DHistogram(*histogramFile,
                                           "hMCTruthRank0DownstreamElectronOriginY");
    TH1F* hMCTruthOriginZ = get1DHistogram(*histogramFile,
                                           "hMCTruthRank0DownstreamElectronOriginZ");
    TH2F* hMCTruthOriginXY = get2DHistogram(*histogramFile,
                                            "hMCTruthRank0DownstreamElectronOriginXY");
    TH2F* hMCTruthOriginXZ = get2DHistogram(*histogramFile,
                                            "hMCTruthRank0DownstreamElectronOriginXZ");
    TH2F* hMCTruthOriginYZ = get2DHistogram(*histogramFile,
                                            "hMCTruthRank0DownstreamElectronOriginYZ");
    TH1F* hMCTruthMomentum = get1DHistogram(*histogramFile,
                                            "hMCTruthRank0DownstreamElectronMomentum");

    TH1F* hRecoMomentum = get1DHistogram(*histogramFile,
                                         "hRecoSelectedDownstreamElectronMomentum");

    TH1F* hRecoVertexX = get1DHistogram(*histogramFile, "hRecoTwoElectronVertexX");
    TH1F* hRecoVertexY = get1DHistogram(*histogramFile, "hRecoTwoElectronVertexY");
    TH1F* hRecoVertexZ = get1DHistogram(*histogramFile, "hRecoTwoElectronVertexZ");
    TH2F* hRecoVertexXY = get2DHistogram(*histogramFile, "hRecoTwoElectronVertexXY");
    TH2F* hRecoVertexXZ = get2DHistogram(*histogramFile, "hRecoTwoElectronVertexXZ");
    TH2F* hRecoVertexYZ = get2DHistogram(*histogramFile, "hRecoTwoElectronVertexYZ");
    TH2F* hRecoVertexSpaceSelectedMomentumTheta1Theta2 = get2DHistogram(
      *histogramFile,
      "hRecoTwoElectronVertexSpaceSelectedMomentumTheta1Theta2");
    TH1F* hRecoVertexSpaceSelectedMomentumOpeningAngle = get1DHistogram(
      *histogramFile,
      "hRecoTwoElectronVertexSpaceSelectedMomentumOpeningAngle");
    TH2F* hRecoVertexFoilIndexMatchedXY = get2DHistogram(
      *histogramFile,
      "hRecoTwoElectronVertexFoilIndexMatchedXY");
    TH2F* hRecoVertexFoilIndexMatchedXZ = get2DHistogram(
      *histogramFile,
      "hRecoTwoElectronVertexFoilIndexMatchedXZ");
    TH2F* hRecoVertexFoilIndexMatchedYZ = get2DHistogram(
      *histogramFile,
      "hRecoTwoElectronVertexFoilIndexMatchedYZ");
    TH1F* hRecoVertexDistance = get1DHistogram(*histogramFile,
                                                "hRecoTwoElectronVertexLineDistance");
    TH1F* hRecoVertexLineParameterS = get1DHistogram(
      *histogramFile,
      "hRecoTwoElectronVertexLineParameterS");
    TH1F* hRecoVertexLineParameterT = get1DHistogram(
      *histogramFile,
      "hRecoTwoElectronVertexLineParameterT");
    TH1F* hRecoVertexAbsLineParameterS = get1DHistogram(
      *histogramFile,
      "hRecoTwoElectronVertexAbsLineParameterS");
    TH1F* hRecoVertexAbsLineParameterT = get1DHistogram(
      *histogramFile,
      "hRecoTwoElectronVertexAbsLineParameterT");
    TH2F* hRecoVertexLineParameterST = get2DHistogram(
      *histogramFile,
      "hRecoTwoElectronVertexLineParameterST");
    TH2F* hRecoVertexAbsLineParameterST = get2DHistogram(
      *histogramFile,
      "hRecoTwoElectronVertexAbsLineParameterST");
    TH2F* hRecoTrackMultiplicityVsDeltaZ = get2DHistogram(
      *histogramFile,
      "hRecoTrackMultiplicityVsDeltaZ");
    TH1F* hRecoVertexTruthDeltaX = get1DHistogram(*histogramFile,
                                                   "hRecoTruthVertexDeltaX");
    TH1F* hRecoVertexTruthDeltaY = get1DHistogram(*histogramFile,
                                                   "hRecoTruthVertexDeltaY");
    TH1F* hRecoVertexTruthDeltaZ = get1DHistogram(*histogramFile,
                                                   "hRecoTruthVertexDeltaZ");
    TH1F* hRecoVertexTruthDistance = get1DHistogram(*histogramFile,
                                                     "hRecoTruthVertexDistance");
    TH1F* hRecoVertexFoilIndexMatchedTruthDeltaX = get1DHistogram(
      *histogramFile,
      "hRecoTruthVertexFoilIndexMatchedDeltaX");
    TH1F* hRecoVertexFoilIndexMatchedTruthDeltaY = get1DHistogram(
      *histogramFile,
      "hRecoTruthVertexFoilIndexMatchedDeltaY");
    TH1F* hRecoVertexFoilIndexMatchedTruthDeltaZ = get1DHistogram(
      *histogramFile,
      "hRecoTruthVertexFoilIndexMatchedDeltaZ");
    TH1F* hRecoVertexFoilIndexMatchedTruthDistance = get1DHistogram(
      *histogramFile,
      "hRecoTruthVertexFoilIndexMatchedDistance");
    TH1F* hRecoVertexSelectedSegmentDeltaTTest = get1DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexSelectedSegmentDeltaT");
    TH1F* hRecoVertexMinTimeDifferenceTest = get1DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeDifference");
    TGraph* gRecoAllSharedFoilCandidateMaxLvsDeltaZ = getGraph(
      *histogramFile,
      "gRecoAllSharedFoilCandidateMaxLvsDeltaZ");
    TGraph* gRecoSpaceSelectedSharedFoilMaxLvsDeltaZ = getGraph(
      *histogramFile,
      "gRecoSpaceSelectedSharedFoilMaxLvsDeltaZ");
    TGraph* gRecoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZ =
      getGraph(*histogramFile,
               "gRecoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZ");
    TGraph* gRecoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZ =
      getGraph(*histogramFile,
               "gRecoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZ");
    TGraph* gRecoSpaceSelectedSharedFoilNumberVsDeltaZ =
      getGraph(*histogramFile,
               "gRecoSpaceSelectedSharedFoilNumberVsDeltaZ");
    std::vector<TH1F*> hRecoSpaceSelectedDeltaZBySelectedFoil;
    hRecoSpaceSelectedDeltaZBySelectedFoil.reserve(kSelectedFoilCount);
    bool allSelectedFoilDeltaZHistogramsPresent = true;
    for (int foilIndex = 0; foilIndex < kSelectedFoilCount; ++foilIndex)
    {
      TH1F* histogram =
        get1DHistogram(*histogramFile,
                       selectedFoilDeltaZHistogramName(foilIndex));
      hRecoSpaceSelectedDeltaZBySelectedFoil.push_back(histogram);
      if (histogram == nullptr)
      {
        allSelectedFoilDeltaZHistogramsPresent = false;
      }
    }
    TGraph* gRecoSpaceSelectedSharedFoilCountVsDeltaZ =
      getGraph(*histogramFile,
               "gRecoSpaceSelectedSharedFoilCountVsDeltaZ");
    TGraph* gRecoSpaceSelectedMaxFoilsHitVsDeltaZ =
      getGraph(*histogramFile,
               "gRecoSpaceSelectedMaxFoilsHitVsDeltaZ");
    TGraph* gRecoSpaceSelectedAbsDeltaFoilsHitVsDeltaZ =
      getGraph(*histogramFile,
               "gRecoSpaceSelectedAbsDeltaFoilsHitVsDeltaZ");
    TGraph* gRecoSpaceSelectedOpeningAngleVsDeltaZ =
      getGraph(*histogramFile,
               "gRecoSpaceSelectedOpeningAngleVsDeltaZ");
    TH1F* hTestRecoVertexMinTimeX = get1DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeX");
    TH1F* hTestRecoVertexMinTimeY = get1DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeY");
    TH1F* hTestRecoVertexMinTimeZ = get1DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeZ");
    TH2F* hTestRecoVertexMinTimeXY = get2DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeXY");
    TH2F* hTestRecoVertexMinTimeXZ = get2DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeXZ");
    TH2F* hTestRecoVertexMinTimeYZ = get2DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeYZ");
    TH2F* hTestRecoVertexMinTimeMomentumTheta1Theta2 = get2DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeMomentumTheta1Theta2");
    TH1F* hTestRecoVertexMinTimeMomentumOpeningAngle = get1DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeMomentumOpeningAngle");
    TH2F* hTestRecoVertexMinTimeFoilIndexMatchedXY = get2DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeFoilIndexMatchedXY");
    TH2F* hTestRecoVertexMinTimeFoilIndexMatchedXZ = get2DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeFoilIndexMatchedXZ");
    TH2F* hTestRecoVertexMinTimeFoilIndexMatchedYZ = get2DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeFoilIndexMatchedYZ");
    TGraph* gTestRecoMinTimeSharedFoilMaxLvsDeltaZ = getGraph(
      *histogramFile,
      "gTESTRecoMinTimeSharedFoilMaxLvsDeltaZ");
    TGraph* gTestRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZ =
      getGraph(*histogramFile,
               "gTESTRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZ");
    TH1F* hTestRecoVertexMinTimeDistance = get1DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeLineDistance");
    TH1F* hTestRecoVertexMinTimeTruthDeltaX = get1DHistogram(
      *histogramFile,
      "hTESTRecoTruthVertexMinTimeDeltaX");
    TH1F* hTestRecoVertexMinTimeTruthDeltaY = get1DHistogram(
      *histogramFile,
      "hTESTRecoTruthVertexMinTimeDeltaY");
    TH1F* hTestRecoVertexMinTimeTruthDeltaZ = get1DHistogram(
      *histogramFile,
      "hTESTRecoTruthVertexMinTimeDeltaZ");
    TH1F* hTestRecoVertexMinTimeTruthDistance = get1DHistogram(
      *histogramFile,
      "hTESTRecoTruthVertexMinTimeDistance");
    TH1F* hTestRecoVertexMinTimeFoilIndexMatchedTruthDeltaX = get1DHistogram(
      *histogramFile,
      "hTESTRecoTruthVertexMinTimeFoilIndexMatchedDeltaX");
    TH1F* hTestRecoVertexMinTimeFoilIndexMatchedTruthDeltaY = get1DHistogram(
      *histogramFile,
      "hTESTRecoTruthVertexMinTimeFoilIndexMatchedDeltaY");
    TH1F* hTestRecoVertexMinTimeFoilIndexMatchedTruthDeltaZ = get1DHistogram(
      *histogramFile,
      "hTESTRecoTruthVertexMinTimeFoilIndexMatchedDeltaZ");
    TH1F* hTestRecoVertexMinTimeFoilIndexMatchedTruthDistance = get1DHistogram(
      *histogramFile,
      "hTESTRecoTruthVertexMinTimeFoilIndexMatchedDistance");

    if (hMCTruthOriginT == nullptr || hMCTruthOriginX == nullptr ||
        hMCTruthOriginY == nullptr || hMCTruthOriginZ == nullptr ||
        hMCTruthOriginXY == nullptr || hMCTruthOriginXZ == nullptr ||
        hMCTruthOriginYZ == nullptr || hMCTruthMomentum == nullptr ||
        hRecoMomentum == nullptr || hRecoVertexX == nullptr ||
        hRecoVertexY == nullptr || hRecoVertexZ == nullptr ||
        hRecoVertexXY == nullptr || hRecoVertexXZ == nullptr ||
        hRecoVertexYZ == nullptr ||
        hRecoVertexSpaceSelectedMomentumTheta1Theta2 == nullptr ||
        hRecoVertexSpaceSelectedMomentumOpeningAngle == nullptr ||
        hRecoVertexDistance == nullptr ||
        hRecoVertexLineParameterS == nullptr ||
        hRecoVertexLineParameterT == nullptr ||
        hRecoVertexAbsLineParameterS == nullptr ||
        hRecoVertexAbsLineParameterT == nullptr ||
        hRecoVertexLineParameterST == nullptr ||
        hRecoVertexAbsLineParameterST == nullptr ||
        hRecoTrackMultiplicityVsDeltaZ == nullptr ||
        hRecoVertexFoilIndexMatchedXY == nullptr ||
        hRecoVertexFoilIndexMatchedXZ == nullptr ||
        hRecoVertexFoilIndexMatchedYZ == nullptr ||
        hRecoVertexTruthDeltaX == nullptr || hRecoVertexTruthDeltaY == nullptr ||
        hRecoVertexTruthDeltaZ == nullptr || hRecoVertexTruthDistance == nullptr ||
        hRecoVertexFoilIndexMatchedTruthDeltaX == nullptr ||
        hRecoVertexFoilIndexMatchedTruthDeltaY == nullptr ||
        hRecoVertexFoilIndexMatchedTruthDeltaZ == nullptr ||
        hRecoVertexFoilIndexMatchedTruthDistance == nullptr ||
        hRecoVertexSelectedSegmentDeltaTTest == nullptr ||
        hRecoVertexMinTimeDifferenceTest == nullptr ||
        gRecoAllSharedFoilCandidateMaxLvsDeltaZ == nullptr ||
        gRecoSpaceSelectedSharedFoilMaxLvsDeltaZ == nullptr ||
        gRecoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZ == nullptr ||
        gRecoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZ == nullptr ||
        gRecoSpaceSelectedSharedFoilNumberVsDeltaZ == nullptr ||
        !allSelectedFoilDeltaZHistogramsPresent ||
        gRecoSpaceSelectedSharedFoilCountVsDeltaZ == nullptr ||
        gRecoSpaceSelectedMaxFoilsHitVsDeltaZ == nullptr ||
        gRecoSpaceSelectedAbsDeltaFoilsHitVsDeltaZ == nullptr ||
        gRecoSpaceSelectedOpeningAngleVsDeltaZ == nullptr ||
        hTestRecoVertexMinTimeX == nullptr ||
        hTestRecoVertexMinTimeY == nullptr ||
        hTestRecoVertexMinTimeZ == nullptr ||
        hTestRecoVertexMinTimeXY == nullptr ||
        hTestRecoVertexMinTimeXZ == nullptr ||
        hTestRecoVertexMinTimeYZ == nullptr ||
        hTestRecoVertexMinTimeMomentumTheta1Theta2 == nullptr ||
        hTestRecoVertexMinTimeMomentumOpeningAngle == nullptr ||
        hTestRecoVertexMinTimeFoilIndexMatchedXY == nullptr ||
        hTestRecoVertexMinTimeFoilIndexMatchedXZ == nullptr ||
        hTestRecoVertexMinTimeFoilIndexMatchedYZ == nullptr ||
        gTestRecoMinTimeSharedFoilMaxLvsDeltaZ == nullptr ||
        gTestRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZ == nullptr ||
        hTestRecoVertexMinTimeDistance == nullptr ||
        hTestRecoVertexMinTimeTruthDeltaX == nullptr ||
        hTestRecoVertexMinTimeTruthDeltaY == nullptr ||
        hTestRecoVertexMinTimeTruthDeltaZ == nullptr ||
        hTestRecoVertexMinTimeTruthDistance == nullptr ||
        hTestRecoVertexMinTimeFoilIndexMatchedTruthDeltaX == nullptr ||
        hTestRecoVertexMinTimeFoilIndexMatchedTruthDeltaY == nullptr ||
        hTestRecoVertexMinTimeFoilIndexMatchedTruthDeltaZ == nullptr ||
        hTestRecoVertexMinTimeFoilIndexMatchedTruthDistance == nullptr)
    {
      std::cerr << "ERROR: one or more expected plot objects are missing from "
                << histogramFileName << std::endl;
      histogramFile->Close();
      delete histogramFile;
      return false;
    }

    style1DHistogram(hMCTruthOriginT);
    style1DHistogram(hMCTruthOriginX);
    style1DHistogram(hMCTruthOriginY);
    style1DHistogram(hMCTruthOriginZ);
    style1DHistogram(hMCTruthMomentum);
    style1DHistogram(hRecoMomentum);
    style1DHistogram(hRecoVertexX);
    style1DHistogram(hRecoVertexY);
    style1DHistogram(hRecoVertexZ);
    style1DHistogram(hRecoVertexSpaceSelectedMomentumOpeningAngle);
    style1DHistogram(hRecoVertexDistance);
    style1DHistogram(hRecoVertexLineParameterS);
    style1DHistogram(hRecoVertexLineParameterT);
    style1DHistogram(hRecoVertexAbsLineParameterS);
    style1DHistogram(hRecoVertexAbsLineParameterT);
    style1DHistogram(hRecoVertexTruthDeltaX);
    style1DHistogram(hRecoVertexTruthDeltaY);
    style1DHistogram(hRecoVertexTruthDeltaZ);
    style1DHistogram(hRecoVertexTruthDistance);
    style1DHistogram(hRecoVertexFoilIndexMatchedTruthDeltaX);
    style1DHistogram(hRecoVertexFoilIndexMatchedTruthDeltaY);
    style1DHistogram(hRecoVertexFoilIndexMatchedTruthDeltaZ);
    style1DHistogram(hRecoVertexFoilIndexMatchedTruthDistance);
    style1DHistogram(hRecoVertexSelectedSegmentDeltaTTest);
    style1DHistogram(hRecoVertexMinTimeDifferenceTest);
    styleScatterGraph(gRecoAllSharedFoilCandidateMaxLvsDeltaZ, kBlue + 1);
    styleScatterGraph(gRecoSpaceSelectedSharedFoilMaxLvsDeltaZ, kGreen + 2);
    styleScatterGraph(gRecoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZ,
                      kBlue + 1);
    styleScatterGraph(gRecoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZ,
                      kGreen + 2);
    styleScatterGraph(gRecoSpaceSelectedSharedFoilNumberVsDeltaZ, kAzure + 2);
    styleScatterGraph(gRecoSpaceSelectedSharedFoilCountVsDeltaZ, kTeal + 2);
    styleScatterGraph(gRecoSpaceSelectedMaxFoilsHitVsDeltaZ, kOrange + 7);
    styleScatterGraph(gRecoSpaceSelectedAbsDeltaFoilsHitVsDeltaZ, kViolet + 1);
    styleScatterGraph(gRecoSpaceSelectedOpeningAngleVsDeltaZ, kRed + 1);
    style1DHistogram(hTestRecoVertexMinTimeX);
    style1DHistogram(hTestRecoVertexMinTimeY);
    style1DHistogram(hTestRecoVertexMinTimeZ);
    style1DHistogram(hTestRecoVertexMinTimeMomentumOpeningAngle);
    style1DHistogram(hTestRecoVertexMinTimeDistance);
    style1DHistogram(hTestRecoVertexMinTimeTruthDeltaX);
    style1DHistogram(hTestRecoVertexMinTimeTruthDeltaY);
    style1DHistogram(hTestRecoVertexMinTimeTruthDeltaZ);
    style1DHistogram(hTestRecoVertexMinTimeTruthDistance);
    style1DHistogram(hTestRecoVertexMinTimeFoilIndexMatchedTruthDeltaX);
    style1DHistogram(hTestRecoVertexMinTimeFoilIndexMatchedTruthDeltaY);
    style1DHistogram(hTestRecoVertexMinTimeFoilIndexMatchedTruthDeltaZ);
    style1DHistogram(hTestRecoVertexMinTimeFoilIndexMatchedTruthDistance);
    styleScatterGraph(gTestRecoMinTimeSharedFoilMaxLvsDeltaZ, kMagenta + 2);
    styleScatterGraph(gTestRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZ,
                      kMagenta + 2);

    TCanvas cTruthOrigin("cTruthOrigin",
                         "MC truth rank-0 downstream electron origin",
                         1200,
                         800);
    cTruthOrigin.Divide(2, 2);
    cTruthOrigin.cd(1); hMCTruthOriginT->Draw("HIST E");
    cTruthOrigin.cd(2); hMCTruthOriginX->Draw("HIST E");
    cTruthOrigin.cd(3); hMCTruthOriginY->Draw("HIST E");
    cTruthOrigin.cd(4); hMCTruthOriginZ->Draw("HIST E");
    cTruthOrigin.SaveAs(
      (monteCarloTruthDirectory + "/TruthOrigin_1D.pdf").c_str());

    TCanvas cTruthOriginMaps("cTruthOriginMaps",
                             "MC truth rank-0 downstream electron origin maps",
                             1500,
                             500);
    cTruthOriginMaps.Divide(3, 1);
    cTruthOriginMaps.cd(1); draw2DHistogram(hMCTruthOriginXY); drawStoppingTargetBoxXY();
    cTruthOriginMaps.cd(2); draw2DHistogram(hMCTruthOriginXZ); drawStoppingTargetBoxXZ();
    cTruthOriginMaps.cd(3); draw2DHistogram(hMCTruthOriginYZ); drawStoppingTargetBoxYZ();
    cTruthOriginMaps.SaveAs(
      (monteCarloTruthDirectory + "/TruthOrigin_2DMaps.pdf").c_str());

    TCanvas cTruthMomentum("cTruthMomentum",
                           "MC truth rank-0 downstream electron momentum",
                           900,
                           700);
    hMCTruthMomentum->Draw("HIST E");
    cTruthMomentum.SaveAs(
      (monteCarloTruthDirectory + "/TruthMomentum.pdf").c_str());

    TCanvas cRecoMomentum("cRecoMomentum",
                          "Selected reconstructed downstream electron momentum",
                          900,
                          700);
    hRecoMomentum->Draw("HIST E");
    cRecoMomentum.SaveAs((rawRecoDirectory + "/RecoMomentum.pdf").c_str());

    TCanvas cRecoVertex("cRecoVertex",
                        "Selected reconstructed two-electron vertex position",
                        1200,
                        800);
    cRecoVertex.Divide(2, 2);
    cRecoVertex.cd(1); hRecoVertexX->Draw("HIST E");
    cRecoVertex.cd(2); hRecoVertexY->Draw("HIST E");
    cRecoVertex.cd(3); hRecoVertexZ->Draw("HIST E");
    cRecoVertex.cd(4); hRecoVertexDistance->Draw("HIST E");
    cRecoVertex.SaveAs((rawRecoDirectory + "/RecoVertex_1D.pdf").c_str());

    TCanvas cRecoVertexLineParameters(
      "cRecoVertexLineParameters",
      "Selected reconstructed two-electron vertex line parameters",
      1200,
      800);
    cRecoVertexLineParameters.Divide(2, 2);
    cRecoVertexLineParameters.cd(1);
    hRecoVertexLineParameterS->Draw("HIST E");
    cRecoVertexLineParameters.cd(2);
    hRecoVertexLineParameterT->Draw("HIST E");
    cRecoVertexLineParameters.cd(3);
    hRecoVertexAbsLineParameterS->Draw("HIST E");
    cRecoVertexLineParameters.cd(4);
    hRecoVertexAbsLineParameterT->Draw("HIST E");
    cRecoVertexLineParameters.SaveAs(
      (rawRecoDirectory + "/RecoVertexLineParameters_1D.pdf").c_str());

    TCanvas cRecoVertexLineParameterST(
      "cRecoVertexLineParameterST",
      "Selected reconstructed two-electron vertex s vs t",
      900,
      700);
    draw2DHistogram(hRecoVertexLineParameterST);
    cRecoVertexLineParameterST.SaveAs(
      (rawRecoDirectory + "/RecoVertexLineParameterST_2D.pdf").c_str());

    TCanvas cRecoVertexAbsLineParameterST(
      "cRecoVertexAbsLineParameterST",
      "Selected reconstructed two-electron vertex |s| vs |t|",
      900,
      700);
    draw2DHistogram(hRecoVertexAbsLineParameterST);
    cRecoVertexAbsLineParameterST.SaveAs(
      (rawRecoDirectory + "/RecoVertexAbsLineParameterST_2D.pdf").c_str());

    TCanvas cRecoTrackMultiplicityVsDeltaZ(
      "cRecoTrackMultiplicityVsDeltaZ",
      "Selected events: reconstructed track multiplicity vs #Delta z",
      900,
      700);
    draw2DHistogram(hRecoTrackMultiplicityVsDeltaZ);
    cRecoTrackMultiplicityVsDeltaZ.SaveAs(
      (rawRecoDirectory + "/RecoTrackMultiplicityVsDeltaZ_2D.pdf").c_str());

    TCanvas cRecoVertexResidual("cRecoVertexResidual",
                                "Selected reconstructed minus truth vertex residual",
                                1600,
                                800);
    cRecoVertexResidual.Divide(2, 2);
    cRecoVertexResidual.cd(1); hRecoVertexTruthDeltaX->Draw("HIST E");
    cRecoVertexResidual.cd(2); hRecoVertexTruthDeltaY->Draw("HIST E");
    cRecoVertexResidual.cd(3); hRecoVertexTruthDeltaZ->Draw("HIST E");
    cRecoVertexResidual.cd(4); hRecoVertexTruthDistance->Draw("HIST E");
    cRecoVertexResidual.SaveAs(
      (rawRecoDirectory + "/RecoVertexTruthResidualXYZ.pdf").c_str());

    TCanvas cRecoVertexMaps("cRecoVertexMaps",
                            "Selected reconstructed two-electron vertex maps",
                            1500,
                            500);
    cRecoVertexMaps.Divide(3, 1);
    cRecoVertexMaps.cd(1); draw2DHistogram(hRecoVertexXY); drawStoppingTargetBoxXY();
    cRecoVertexMaps.cd(2); draw2DHistogram(hRecoVertexXZ); drawStoppingTargetBoxXZ();
    cRecoVertexMaps.cd(3); draw2DHistogram(hRecoVertexYZ); drawStoppingTargetBoxYZ();
    cRecoVertexMaps.SaveAs(
      (rawRecoDirectory + "/RecoVertex_2DMaps.pdf").c_str());

    TCanvas cRecoVertexAngularDiagnostics(
      "cRecoVertexAngularDiagnostics",
      "Selected shared ST_Foils angular diagnostics",
      1600,
      1200);
    cRecoVertexAngularDiagnostics.Divide(2, 2);
    cRecoVertexAngularDiagnostics.cd(1);
    draw2DHistogram(hRecoVertexSpaceSelectedMomentumTheta1Theta2);
    cRecoVertexAngularDiagnostics.cd(2);
    draw2DHistogram(hTestRecoVertexMinTimeMomentumTheta1Theta2);
    cRecoVertexAngularDiagnostics.cd(3);
    hRecoVertexSpaceSelectedMomentumOpeningAngle->SetLineColor(kGreen + 2);
    hTestRecoVertexMinTimeMomentumOpeningAngle->SetLineColor(kMagenta + 2);
    const double openingAngleMaximum =
      std::max(hRecoVertexSpaceSelectedMomentumOpeningAngle->GetMaximum(),
               hTestRecoVertexMinTimeMomentumOpeningAngle->GetMaximum());
    if (openingAngleMaximum > 0.0)
    {
      hRecoVertexSpaceSelectedMomentumOpeningAngle->SetMaximum(
        1.15 * openingAngleMaximum);
    }
    hRecoVertexSpaceSelectedMomentumOpeningAngle->Draw("HIST E");
    hTestRecoVertexMinTimeMomentumOpeningAngle->Draw("HIST E SAME");
    TLegend openingAngleLegend(0.48, 0.72, 0.88, 0.88);
    openingAngleLegend.SetBorderSize(0);
    openingAngleLegend.SetFillStyle(0);
    openingAngleLegend.AddEntry(hRecoVertexSpaceSelectedMomentumOpeningAngle,
                                "space-selected",
                                "l");
    openingAngleLegend.AddEntry(hTestRecoVertexMinTimeMomentumOpeningAngle,
                                "TEST min-|#Delta t|",
                                "l");
    openingAngleLegend.Draw();
    cRecoVertexAngularDiagnostics.cd(4);
    drawScatterGraph(gRecoSpaceSelectedOpeningAngleVsDeltaZ,
                     kOpeningAngleDegreesPlotMin,
                     kOpeningAngleDegreesPlotMax,
                     kRecoTruthDeltaZPlotMin,
                     kRecoTruthDeltaZPlotMax);
    cRecoVertexAngularDiagnostics.SaveAs(
      (anglesDirectory + "/RecoVertexAngularDiagnostics.pdf").c_str());

    TCanvas cRecoVertexFoilIndexMatchedMaps(
      "cRecoVertexFoilIndexMatchedMaps",
      "Foil-index matched reconstructed two-electron vertex maps",
      1500,
      500);
    cRecoVertexFoilIndexMatchedMaps.Divide(3, 1);
    cRecoVertexFoilIndexMatchedMaps.cd(1);
    draw2DHistogram(hRecoVertexFoilIndexMatchedXY); drawStoppingTargetBoxXY();
    cRecoVertexFoilIndexMatchedMaps.cd(2);
    draw2DHistogram(hRecoVertexFoilIndexMatchedXZ); drawStoppingTargetBoxXZ();
    cRecoVertexFoilIndexMatchedMaps.cd(3);
    draw2DHistogram(hRecoVertexFoilIndexMatchedYZ); drawStoppingTargetBoxYZ();
    cRecoVertexFoilIndexMatchedMaps.SaveAs(
      (deltaZVsFoilDirectory + "/RecoVertexFoilIndexMatched_2DMaps.pdf").c_str());

    TCanvas cRecoVertexFoilIndexMatchedResidual(
      "cRecoVertexFoilIndexMatchedResidual",
      "Foil-index matched reconstructed minus truth vertex residual",
      1600,
      800);
    cRecoVertexFoilIndexMatchedResidual.Divide(2, 2);
    cRecoVertexFoilIndexMatchedResidual.cd(1);
    hRecoVertexFoilIndexMatchedTruthDeltaX->Draw("HIST E");
    cRecoVertexFoilIndexMatchedResidual.cd(2);
    hRecoVertexFoilIndexMatchedTruthDeltaY->Draw("HIST E");
    cRecoVertexFoilIndexMatchedResidual.cd(3);
    hRecoVertexFoilIndexMatchedTruthDeltaZ->Draw("HIST E");
    cRecoVertexFoilIndexMatchedResidual.cd(4);
    hRecoVertexFoilIndexMatchedTruthDistance->Draw("HIST E");
    cRecoVertexFoilIndexMatchedResidual.SaveAs(
      (deltaZVsFoilDirectory + "/RecoVertexFoilIndexMatchedTruthResidualXYZ.pdf").c_str());

    TCanvas cRecoVertexDeltaTTest("cRecoVertexDeltaTTest",
                                  "Selected reconstructed two-electron segment time differences",
                                  900,
                                  700);
    hRecoVertexSelectedSegmentDeltaTTest->Draw("HIST E");
    cRecoVertexDeltaTTest.SaveAs(
      (rawRecoDirectory + "/RecoVertexSelectedSegmentDeltaT_TEST.pdf").c_str());

    TCanvas cRecoVertexMinTimeDifferenceTest("cRecoVertexMinTimeDifferenceTest",
                                             "Minimum-|#Delta t| shared ST_Foils pair time difference",
                                             900,
                                             700);
    hRecoVertexMinTimeDifferenceTest->Draw("HIST E");
    cRecoVertexMinTimeDifferenceTest.SaveAs(
      (rawRecoDirectory + "/RecoVertexMinTimeDifference_TEST.pdf").c_str());

    TCanvas cRecoSharedFoilMaxLvsDeltaZ(
      "cRecoSharedFoilMaxLvsDeltaZ",
      "Shared same-foil pair choices: #Delta z vs max(L_{1}, L_{2})",
      1800,
      600);
    cRecoSharedFoilMaxLvsDeltaZ.Divide(3, 1);
    cRecoSharedFoilMaxLvsDeltaZ.cd(1);
    drawScatterGraph(gRecoAllSharedFoilCandidateMaxLvsDeltaZ,
                     kFoilPointTruthDistancePlotMin,
                     kFoilPointTruthDistancePlotMax,
                     kRecoTruthDeltaZPlotMin,
                     kRecoTruthDeltaZPlotMax);
    cRecoSharedFoilMaxLvsDeltaZ.cd(2);
    drawScatterGraph(gRecoSpaceSelectedSharedFoilMaxLvsDeltaZ,
                     kFoilPointTruthDistancePlotMin,
                     kFoilPointTruthDistancePlotMax,
                     kRecoTruthDeltaZPlotMin,
                     kRecoTruthDeltaZPlotMax);
    cRecoSharedFoilMaxLvsDeltaZ.cd(3);
    drawScatterGraph(gTestRecoMinTimeSharedFoilMaxLvsDeltaZ,
                     kFoilPointTruthDistancePlotMin,
                     kFoilPointTruthDistancePlotMax,
                     kRecoTruthDeltaZPlotMin,
                     kRecoTruthDeltaZPlotMax);
    cRecoSharedFoilMaxLvsDeltaZ.SaveAs(
      (deltaZVsFoilDirectory + "/RecoSharedFoilMaxLvsDeltaZ_Scatter.pdf").c_str());

    TCanvas cRecoSharedFoilAvgAbsLineParameterVsDeltaZ(
      "cRecoSharedFoilAvgAbsLineParameterVsDeltaZ",
      "Shared same-foil pair choices: #Delta z vs (|s|+|t|)/2",
      1800,
      600);
    cRecoSharedFoilAvgAbsLineParameterVsDeltaZ.Divide(3, 1);
    cRecoSharedFoilAvgAbsLineParameterVsDeltaZ.cd(1);
    drawScatterGraph(gRecoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZ,
                     kAverageAbsLineParameterPlotMin,
                     kAverageAbsLineParameterPlotMax,
                     kRecoTruthDeltaZPlotMin,
                     kRecoTruthDeltaZPlotMax);
    cRecoSharedFoilAvgAbsLineParameterVsDeltaZ.cd(2);
    drawScatterGraph(gRecoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZ,
                     kAverageAbsLineParameterPlotMin,
                     kAverageAbsLineParameterPlotMax,
                     kRecoTruthDeltaZPlotMin,
                     kRecoTruthDeltaZPlotMax);
    cRecoSharedFoilAvgAbsLineParameterVsDeltaZ.cd(3);
    drawScatterGraph(gTestRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZ,
                     kAverageAbsLineParameterPlotMin,
                     kAverageAbsLineParameterPlotMax,
                     kRecoTruthDeltaZPlotMin,
                     kRecoTruthDeltaZPlotMax);
    cRecoSharedFoilAvgAbsLineParameterVsDeltaZ.SaveAs(
      (deltaZVsFoilDirectory + "/RecoSharedFoilAvgAbsLineParameterVsDeltaZ_Scatter.pdf").c_str());

    TCanvas cRecoSpaceSelectedDeltaZFoilRelations(
      "cRecoSpaceSelectedDeltaZFoilRelations",
      "Space-selected shared same-foil pair: #Delta z foil diagnostics",
      1600,
      1200);
    cRecoSpaceSelectedDeltaZFoilRelations.Divide(2, 2);
    cRecoSpaceSelectedDeltaZFoilRelations.cd(1);
    drawScatterGraph(gRecoSpaceSelectedSharedFoilNumberVsDeltaZ,
                     kSharedFoilIndexPlotMin,
                     kSharedFoilIndexPlotMax,
                     kRecoTruthDeltaZPlotMin,
                     kRecoTruthDeltaZPlotMax);
    cRecoSpaceSelectedDeltaZFoilRelations.cd(2);
    drawScatterGraph(gRecoSpaceSelectedSharedFoilCountVsDeltaZ,
                     kFoilCountPlotMin,
                     kFoilCountPlotMax,
                     kRecoTruthDeltaZPlotMin,
                     kRecoTruthDeltaZPlotMax);
    cRecoSpaceSelectedDeltaZFoilRelations.cd(3);
    drawScatterGraph(gRecoSpaceSelectedMaxFoilsHitVsDeltaZ,
                     kFoilCountPlotMin,
                     kFoilCountPlotMax,
                     kRecoTruthDeltaZPlotMin,
                     kRecoTruthDeltaZPlotMax);
    cRecoSpaceSelectedDeltaZFoilRelations.cd(4);
    drawScatterGraph(gRecoSpaceSelectedAbsDeltaFoilsHitVsDeltaZ,
                     kFoilCountPlotMin,
                     kFoilCountPlotMax,
                     kRecoTruthDeltaZPlotMin,
                     kRecoTruthDeltaZPlotMax);
    cRecoSpaceSelectedDeltaZFoilRelations.SaveAs(
      (deltaZVsFoilDirectory + "/RecoSpaceSelectedDeltaZFoilRelations_Scatter.pdf").c_str());

    saveDeltaZBySelectedFoilCanvases(
      hRecoSpaceSelectedDeltaZBySelectedFoil,
      deltaZVsFoilDirectory);

    TCanvas cTestRecoVertexMinTime("cTestRecoVertexMinTime",
                                   "TEST: minimum-|#Delta t| reconstructed two-electron vertex position",
                                   1200,
                                   800);
    cTestRecoVertexMinTime.Divide(2, 2);
    cTestRecoVertexMinTime.cd(1); hTestRecoVertexMinTimeX->Draw("HIST E");
    cTestRecoVertexMinTime.cd(2); hTestRecoVertexMinTimeY->Draw("HIST E");
    cTestRecoVertexMinTime.cd(3); hTestRecoVertexMinTimeZ->Draw("HIST E");
    cTestRecoVertexMinTime.cd(4); hTestRecoVertexMinTimeDistance->Draw("HIST E");
    cTestRecoVertexMinTime.SaveAs(
      (rawRecoDirectory + "/RecoVertexMinTimeChoice_1D_TEST.pdf").c_str());

    TCanvas cTestRecoVertexMinTimeResidual("cTestRecoVertexMinTimeResidual",
                                           "TEST: minimum-|#Delta t| reconstructed minus truth vertex residual",
                                           1600,
                                           800);
    cTestRecoVertexMinTimeResidual.Divide(2, 2);
    cTestRecoVertexMinTimeResidual.cd(1);
    hTestRecoVertexMinTimeTruthDeltaX->Draw("HIST E");
    cTestRecoVertexMinTimeResidual.cd(2);
    hTestRecoVertexMinTimeTruthDeltaY->Draw("HIST E");
    cTestRecoVertexMinTimeResidual.cd(3);
    hTestRecoVertexMinTimeTruthDeltaZ->Draw("HIST E");
    cTestRecoVertexMinTimeResidual.cd(4);
    hTestRecoVertexMinTimeTruthDistance->Draw("HIST E");
    cTestRecoVertexMinTimeResidual.SaveAs(
      (rawRecoDirectory + "/RecoVertexMinTimeChoiceTruthResidualXYZ_TEST.pdf").c_str());

    TCanvas cTestRecoVertexMinTimeMaps("cTestRecoVertexMinTimeMaps",
                                       "TEST: minimum-|#Delta t| reconstructed two-electron vertex maps",
                                       1500,
                                       500);
    cTestRecoVertexMinTimeMaps.Divide(3, 1);
    cTestRecoVertexMinTimeMaps.cd(1);
    draw2DHistogram(hTestRecoVertexMinTimeXY); drawStoppingTargetBoxXY();
    cTestRecoVertexMinTimeMaps.cd(2);
    draw2DHistogram(hTestRecoVertexMinTimeXZ); drawStoppingTargetBoxXZ();
    cTestRecoVertexMinTimeMaps.cd(3);
    draw2DHistogram(hTestRecoVertexMinTimeYZ); drawStoppingTargetBoxYZ();
    cTestRecoVertexMinTimeMaps.SaveAs(
      (rawRecoDirectory + "/RecoVertexMinTimeChoice_2DMaps_TEST.pdf").c_str());

    TCanvas cTestRecoVertexMinTimeFoilIndexMatchedMaps(
      "cTestRecoVertexMinTimeFoilIndexMatchedMaps",
      "TEST: foil-index matched minimum-|#Delta t| reconstructed vertex maps",
      1500,
      500);
    cTestRecoVertexMinTimeFoilIndexMatchedMaps.Divide(3, 1);
    cTestRecoVertexMinTimeFoilIndexMatchedMaps.cd(1);
    draw2DHistogram(hTestRecoVertexMinTimeFoilIndexMatchedXY);
    drawStoppingTargetBoxXY();
    cTestRecoVertexMinTimeFoilIndexMatchedMaps.cd(2);
    draw2DHistogram(hTestRecoVertexMinTimeFoilIndexMatchedXZ);
    drawStoppingTargetBoxXZ();
    cTestRecoVertexMinTimeFoilIndexMatchedMaps.cd(3);
    draw2DHistogram(hTestRecoVertexMinTimeFoilIndexMatchedYZ);
    drawStoppingTargetBoxYZ();
    cTestRecoVertexMinTimeFoilIndexMatchedMaps.SaveAs(
      (deltaZVsFoilDirectory + "/RecoVertexMinTimeFoilIndexMatched_2DMaps_TEST.pdf").c_str());

    TCanvas cTestRecoVertexMinTimeFoilIndexMatchedResidual(
      "cTestRecoVertexMinTimeFoilIndexMatchedResidual",
      "TEST: foil-index matched minimum-|#Delta t| reconstructed minus truth vertex residual",
      1600,
      800);
    cTestRecoVertexMinTimeFoilIndexMatchedResidual.Divide(2, 2);
    cTestRecoVertexMinTimeFoilIndexMatchedResidual.cd(1);
    hTestRecoVertexMinTimeFoilIndexMatchedTruthDeltaX->Draw("HIST E");
    cTestRecoVertexMinTimeFoilIndexMatchedResidual.cd(2);
    hTestRecoVertexMinTimeFoilIndexMatchedTruthDeltaY->Draw("HIST E");
    cTestRecoVertexMinTimeFoilIndexMatchedResidual.cd(3);
    hTestRecoVertexMinTimeFoilIndexMatchedTruthDeltaZ->Draw("HIST E");
    cTestRecoVertexMinTimeFoilIndexMatchedResidual.cd(4);
    hTestRecoVertexMinTimeFoilIndexMatchedTruthDistance->Draw("HIST E");
    cTestRecoVertexMinTimeFoilIndexMatchedResidual.SaveAs(
      (deltaZVsFoilDirectory + "/RecoVertexMinTimeFoilIndexMatchedTruthResidualXYZ_TEST.pdf").c_str());

    histogramFile->Close();
    delete histogramFile;
    return true;
  }
}

#endif

