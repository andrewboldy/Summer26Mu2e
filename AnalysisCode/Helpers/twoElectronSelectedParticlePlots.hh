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
//   afterwards and produces human-readable PDFs under:
//
//       Plots/TruthVsRecoPlots/
//
//   That keeps the event loop simple and makes it easy to regenerate plots
//   without rerunning the full ntuple scan once the histogram file exists.
//
//----------------------------------------------------------------------------------

#include <string>

#include <iostream>

#include <TCanvas.h>
#include <TBox.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TSystem.h>

namespace twoelectronplots
{
  // These bounds match the stopping-target geometry overlay used by the
  // original vertexer macro.
  constexpr double kStoppingTargetOuterRadius = 75.0;
  constexpr double kStoppingTargetZMin = -4700.053;
  constexpr double kStoppingTargetZMax = -3899.947;

  inline void ensureOutputDirectory(const std::string& outputPath)
  {
    const std::string::size_type slashPosition = outputPath.find_last_of("/\\");
    if (slashPosition == std::string::npos)
    {
      return;
    }

    const std::string outputDirectory = outputPath.substr(0, slashPosition);
    if (!outputDirectory.empty() && gSystem != nullptr)
    {
      gSystem->mkdir(outputDirectory.c_str(), true);
    }
  }

  inline TH1F* get1DHistogram(TFile& file, const std::string& name)
  {
    return dynamic_cast<TH1F*>(file.Get(name.c_str()));
  }

  inline TH2F* get2DHistogram(TFile& file, const std::string& name)
  {
    return dynamic_cast<TH2F*>(file.Get(name.c_str()));
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

    ensureOutputDirectory(outputDirectory + "/placeholder.pdf");

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
    TH1F* hRecoVertexDistance = get1DHistogram(*histogramFile,
                                                "hRecoTwoElectronVertexLineDistance");
    TH1F* hRecoVertexTruthDeltaX = get1DHistogram(*histogramFile,
                                                   "hRecoTruthVertexDeltaX");
    TH1F* hRecoVertexTruthDeltaY = get1DHistogram(*histogramFile,
                                                   "hRecoTruthVertexDeltaY");
    TH1F* hRecoVertexTruthDeltaZ = get1DHistogram(*histogramFile,
                                                   "hRecoTruthVertexDeltaZ");
    TH1F* hRecoVertexTruthDistance = get1DHistogram(*histogramFile,
                                                     "hRecoTruthVertexDistance");
    TH1F* hRecoVertexSelectedSegmentDeltaTTest = get1DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexSelectedSegmentDeltaT");
    TH1F* hRecoVertexMinTimeDifferenceTest = get1DHistogram(
      *histogramFile,
      "hTESTRecoTwoElectronVertexMinTimeDifference");

    if (hMCTruthOriginT == nullptr || hMCTruthOriginX == nullptr ||
        hMCTruthOriginY == nullptr || hMCTruthOriginZ == nullptr ||
        hMCTruthOriginXY == nullptr || hMCTruthOriginXZ == nullptr ||
        hMCTruthOriginYZ == nullptr || hMCTruthMomentum == nullptr ||
        hRecoMomentum == nullptr || hRecoVertexX == nullptr ||
        hRecoVertexY == nullptr || hRecoVertexZ == nullptr ||
        hRecoVertexXY == nullptr || hRecoVertexXZ == nullptr ||
        hRecoVertexYZ == nullptr || hRecoVertexDistance == nullptr ||
        hRecoVertexTruthDeltaX == nullptr || hRecoVertexTruthDeltaY == nullptr ||
        hRecoVertexTruthDeltaZ == nullptr || hRecoVertexTruthDistance == nullptr ||
        hRecoVertexSelectedSegmentDeltaTTest == nullptr ||
        hRecoVertexMinTimeDifferenceTest == nullptr)
    {
      std::cerr << "ERROR: one or more expected histograms are missing from "
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
    style1DHistogram(hRecoVertexDistance);
    style1DHistogram(hRecoVertexTruthDeltaX);
    style1DHistogram(hRecoVertexTruthDeltaY);
    style1DHistogram(hRecoVertexTruthDeltaZ);
    style1DHistogram(hRecoVertexTruthDistance);
    style1DHistogram(hRecoVertexSelectedSegmentDeltaTTest);
    style1DHistogram(hRecoVertexMinTimeDifferenceTest);

    TCanvas cTruthOrigin("cTruthOrigin",
                         "MC truth rank-0 downstream electron origin",
                         1200,
                         800);
    cTruthOrigin.Divide(2, 2);
    cTruthOrigin.cd(1); hMCTruthOriginT->Draw("HIST E");
    cTruthOrigin.cd(2); hMCTruthOriginX->Draw("HIST E");
    cTruthOrigin.cd(3); hMCTruthOriginY->Draw("HIST E");
    cTruthOrigin.cd(4); hMCTruthOriginZ->Draw("HIST E");
    cTruthOrigin.SaveAs((outputDirectory + "/TruthOrigin_1D.pdf").c_str());

    TCanvas cTruthOriginMaps("cTruthOriginMaps",
                             "MC truth rank-0 downstream electron origin maps",
                             1500,
                             500);
    cTruthOriginMaps.Divide(3, 1);
    cTruthOriginMaps.cd(1); draw2DHistogram(hMCTruthOriginXY); drawStoppingTargetBoxXY();
    cTruthOriginMaps.cd(2); draw2DHistogram(hMCTruthOriginXZ); drawStoppingTargetBoxXZ();
    cTruthOriginMaps.cd(3); draw2DHistogram(hMCTruthOriginYZ); drawStoppingTargetBoxYZ();
    cTruthOriginMaps.SaveAs((outputDirectory + "/TruthOrigin_2DMaps.pdf").c_str());

    TCanvas cTruthMomentum("cTruthMomentum",
                           "MC truth rank-0 downstream electron momentum",
                           900,
                           700);
    hMCTruthMomentum->Draw("HIST E");
    cTruthMomentum.SaveAs((outputDirectory + "/TruthMomentum.pdf").c_str());

    TCanvas cRecoMomentum("cRecoMomentum",
                          "Selected reconstructed downstream electron momentum",
                          900,
                          700);
    hRecoMomentum->Draw("HIST E");
    cRecoMomentum.SaveAs((outputDirectory + "/RecoMomentum.pdf").c_str());

    TCanvas cRecoVertex("cRecoVertex",
                        "Selected reconstructed two-electron vertex position",
                        1200,
                        800);
    cRecoVertex.Divide(2, 2);
    cRecoVertex.cd(1); hRecoVertexX->Draw("HIST E");
    cRecoVertex.cd(2); hRecoVertexY->Draw("HIST E");
    cRecoVertex.cd(3); hRecoVertexZ->Draw("HIST E");
    cRecoVertex.cd(4); hRecoVertexDistance->Draw("HIST E");
    cRecoVertex.SaveAs((outputDirectory + "/RecoVertex_1D.pdf").c_str());

    TCanvas cRecoVertexResidual("cRecoVertexResidual",
                                "Selected reconstructed minus truth vertex residual",
                                1600,
                                800);
    cRecoVertexResidual.Divide(2, 2);
    cRecoVertexResidual.cd(1); hRecoVertexTruthDeltaX->Draw("HIST E");
    cRecoVertexResidual.cd(2); hRecoVertexTruthDeltaY->Draw("HIST E");
    cRecoVertexResidual.cd(3); hRecoVertexTruthDeltaZ->Draw("HIST E");
    cRecoVertexResidual.cd(4); hRecoVertexTruthDistance->Draw("HIST E");
    cRecoVertexResidual.SaveAs((outputDirectory + "/RecoVertexTruthResidualXYZ.pdf").c_str());

    TCanvas cRecoVertexMaps("cRecoVertexMaps",
                            "Selected reconstructed two-electron vertex maps",
                            1500,
                            500);
    cRecoVertexMaps.Divide(3, 1);
    cRecoVertexMaps.cd(1); draw2DHistogram(hRecoVertexXY); drawStoppingTargetBoxXY();
    cRecoVertexMaps.cd(2); draw2DHistogram(hRecoVertexXZ); drawStoppingTargetBoxXZ();
    cRecoVertexMaps.cd(3); draw2DHistogram(hRecoVertexYZ); drawStoppingTargetBoxYZ();
    cRecoVertexMaps.SaveAs((outputDirectory + "/RecoVertex_2DMaps.pdf").c_str());

    TCanvas cRecoVertexDeltaTTest("cRecoVertexDeltaTTest",
                                  "Selected reconstructed two-electron segment time differences",
                                  900,
                                  700);
    hRecoVertexSelectedSegmentDeltaTTest->Draw("HIST E");
    cRecoVertexDeltaTTest.SaveAs(
      (outputDirectory + "/RecoVertexSelectedSegmentDeltaT_TEST.pdf").c_str());

    TCanvas cRecoVertexMinTimeDifferenceTest("cRecoVertexMinTimeDifferenceTest",
                                             "Minimum-|#Delta t| shared ST_Foils pair time difference",
                                             900,
                                             700);
    hRecoVertexMinTimeDifferenceTest->Draw("HIST E");
    cRecoVertexMinTimeDifferenceTest.SaveAs(
      (outputDirectory + "/RecoVertexMinTimeDifference_TEST.pdf").c_str());

    histogramFile->Close();
    delete histogramFile;
    return true;
  }
}

#endif

