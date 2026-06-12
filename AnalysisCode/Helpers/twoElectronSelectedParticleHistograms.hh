#ifndef CREATEDCODE_HISTOGRAMMAKERS_HELPERS_TWOELECTRONSELECTEDPARTICLEHISTOGRAMS_HH
#define CREATEDCODE_HISTOGRAMMAKERS_HELPERS_TWOELECTRONSELECTEDPARTICLEHISTOGRAMS_HH

//----------------------------------------------------------------------------------
//
// twoElectronSelectedParticleHistograms.hh
// Written by Andrew Boldy University of South Carolina, 2026
// Assisted by Codex
//
// Purpose:
//   Header-only histogram helper for the ntuple-based two-electron truth/reco
//   comparison workflow.
//
//   The comparer macro owns the event loop and decides when an event passes the
//   reconstructed-event selection.  This helper owns the repetitive histogram
//   mechanics:
//
//     - book MC truth origin position histograms
//     - book MC truth momentum histograms
//     - book reconstructed momentum histograms
//     - book reconstructed two-track vertex position histograms
//     - book TEST reconstructed vertex histograms for the minimum-|Delta t|
//       shared ST_Foils crossing choice
//     - fill only rank-0 downstream electron truth particles
//     - write the full histogram set to a ROOT file
//
// Selection contract:
//   This helper is intentionally strict about MC truth particle selection:
//
//     sim.valid
//     sim.rank == 0
//     sim.pdg == 11
//     sim.mom.z() > 0
//
//   The caller is responsible for enforcing the event-level reconstruction
//   selection before calling the fill functions.  In the current comparer that
//   means exactly two reconstructed downstream electron tracks, each with an
//   associated calorimeter hit and reconstructed momentum in 50-53 MeV/c.
//   The momentum histogram range is intentionally wider: 30-55 MeV/c.
//
// Coordinate and unit conventions:
//   - positions are in mm
//   - momenta are in MeV/c
//   - times are in ns
//   - EventNtuple SimInfo and TrkSegInfo positions are in tracker/detector
//     coordinates, consistent with the existing vertexer macros
//
//----------------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>

#include "EventNtuple/inc/SimInfo.hh"
#include "twoParticleVertexer.hh"

namespace twoelectronhist
{
  // Histogram ranges are kept in one config object so later tuning does not
  // touch the event-selection code.  The default position bounds match the
  // original vertexer macro so the truth and reco position plots are directly
  // comparable.
  struct HistogramConfig
  {
    int momentumBins = 140;
    double momentumMin = 30.0;
    double momentumMax = 55.0;

    int timeBins = 200;
    double timeMin = 0.0;
    double timeMax = 2000.0;

    int transverseBins = 200;
    double transverseMin = -200.0;
    double transverseMax = 200.0;

    int zBins = 405;
    double zMin = -5000.0;
    double zMax = -3700.0;

    int vertexDistanceBins = 200;
    double vertexDistanceMin = 0.0;
    double vertexDistanceMax = 250.0;

    int lineParameterBins = 240;
    double lineParameterMin = -1200.0;
    double lineParameterMax = 1200.0;
    int absLineParameterBins = 240;
    double absLineParameterMin = 0.0;
    double absLineParameterMax = 1200.0;

    int timeDifferenceBins = 200;
    double timeDifferenceMin = -250.0;
    double timeDifferenceMax = 250.0;

    double foilPointTruthDistanceMin = 0.0;
    double foilPointTruthDistanceMax = 1000.0;
    double recoTruthDeltaZMin = -1000.0;
    double recoTruthDeltaZMax = 1000.0;
  };

  struct HistogramBook
  {
    HistogramConfig config;

    TH1F* mcTruthOriginT = nullptr;
    TH1F* mcTruthOriginX = nullptr;
    TH1F* mcTruthOriginY = nullptr;
    TH1F* mcTruthOriginZ = nullptr;
    TH2F* mcTruthOriginXY = nullptr;
    TH2F* mcTruthOriginXZ = nullptr;
    TH2F* mcTruthOriginYZ = nullptr;
    TH1F* mcTruthMomentum = nullptr;

    TH1F* recoMomentum = nullptr;

    TH1F* recoVertexX = nullptr;
    TH1F* recoVertexY = nullptr;
    TH1F* recoVertexZ = nullptr;
    TH2F* recoVertexXY = nullptr;
    TH2F* recoVertexXZ = nullptr;
    TH2F* recoVertexYZ = nullptr;
    TH2F* recoVertexFoilIndexMatchedXY = nullptr;
    TH2F* recoVertexFoilIndexMatchedXZ = nullptr;
    TH2F* recoVertexFoilIndexMatchedYZ = nullptr;
    TH1F* recoVertexLineDistance = nullptr;
    TH1F* recoVertexLineParameterS = nullptr;
    TH1F* recoVertexLineParameterT = nullptr;
    TH1F* recoVertexAbsLineParameterS = nullptr;
    TH1F* recoVertexAbsLineParameterT = nullptr;
    TH2F* recoVertexAbsLineParameterST = nullptr;
    TH1F* recoVertexSelectedSegmentDeltaTTest = nullptr;
    TH1F* recoVertexMinTimeDifferenceTest = nullptr;
    TGraph* recoAllSharedFoilCandidateMaxLvsDeltaZ = nullptr;
    TGraph* recoSpaceSelectedSharedFoilMaxLvsDeltaZ = nullptr;
    TH1F* testRecoVertexMinTimeX = nullptr;
    TH1F* testRecoVertexMinTimeY = nullptr;
    TH1F* testRecoVertexMinTimeZ = nullptr;
    TH2F* testRecoVertexMinTimeXY = nullptr;
    TH2F* testRecoVertexMinTimeXZ = nullptr;
    TH2F* testRecoVertexMinTimeYZ = nullptr;
    TH2F* testRecoVertexMinTimeFoilIndexMatchedXY = nullptr;
    TH2F* testRecoVertexMinTimeFoilIndexMatchedXZ = nullptr;
    TH2F* testRecoVertexMinTimeFoilIndexMatchedYZ = nullptr;
    TGraph* testRecoMinTimeSharedFoilMaxLvsDeltaZ = nullptr;
    TH1F* testRecoVertexMinTimeLineDistance = nullptr;
    TH1F* recoVertexTruthDeltaX = nullptr;
    TH1F* recoVertexTruthDeltaY = nullptr;
    TH1F* recoVertexTruthDeltaZ = nullptr;
    TH1F* recoVertexTruthDistance = nullptr;
    TH1F* recoVertexFoilIndexMatchedTruthDeltaX = nullptr;
    TH1F* recoVertexFoilIndexMatchedTruthDeltaY = nullptr;
    TH1F* recoVertexFoilIndexMatchedTruthDeltaZ = nullptr;
    TH1F* recoVertexFoilIndexMatchedTruthDistance = nullptr;
    TH1F* testRecoVertexMinTimeTruthDeltaX = nullptr;
    TH1F* testRecoVertexMinTimeTruthDeltaY = nullptr;
    TH1F* testRecoVertexMinTimeTruthDeltaZ = nullptr;
    TH1F* testRecoVertexMinTimeTruthDistance = nullptr;
    TH1F* testRecoVertexMinTimeFoilIndexMatchedTruthDeltaX = nullptr;
    TH1F* testRecoVertexMinTimeFoilIndexMatchedTruthDeltaY = nullptr;
    TH1F* testRecoVertexMinTimeFoilIndexMatchedTruthDeltaZ = nullptr;
    TH1F* testRecoVertexMinTimeFoilIndexMatchedTruthDistance = nullptr;

    Long64_t selectedEvents = 0;
    Long64_t mcTruthEntries = 0;
    Long64_t recoMomentumEntries = 0;
    Long64_t recoVertexEntries = 0;
    Long64_t recoVertexFoilIndexMatchedEntries = 0;
    Long64_t recoVertexLineParameterEntries = 0;
    Long64_t recoVertexSelectedSegmentDeltaTEntries = 0;
    Long64_t recoVertexMinTimeDifferenceEntries = 0;
    Long64_t recoAllSharedFoilCandidateMaxLvsDeltaZEntries = 0;
    Long64_t recoSpaceSelectedSharedFoilMaxLvsDeltaZEntries = 0;
    Long64_t testRecoVertexMinTimeEntries = 0;
    Long64_t testRecoVertexMinTimeFoilIndexMatchedEntries = 0;
    Long64_t testRecoMinTimeSharedFoilMaxLvsDeltaZEntries = 0;
    Long64_t recoVertexTruthResidualEntries = 0;
    Long64_t recoVertexFoilIndexMatchedTruthResidualEntries = 0;
    Long64_t testRecoVertexMinTimeTruthResidualEntries = 0;
    Long64_t testRecoVertexMinTimeFoilIndexMatchedTruthResidualEntries = 0;
  };

  inline TH1F* make1D(const std::string& name,
                      const std::string& title,
                      int bins,
                      double min,
                      double max)
  {
    TH1F* histogram = new TH1F(name.c_str(), title.c_str(), bins, min, max);
    histogram->SetDirectory(nullptr);
    return histogram;
  }

  inline TH2F* make2D(const std::string& name,
                      const std::string& title,
                      int xBins,
                      double xMin,
                      double xMax,
                      int yBins,
                      double yMin,
                      double yMax)
  {
    TH2F* histogram = new TH2F(name.c_str(),
                              title.c_str(),
                              xBins,
                              xMin,
                              xMax,
                              yBins,
                              yMin,
                              yMax);
    histogram->SetDirectory(nullptr);
    return histogram;
  }

  inline TGraph* makeGraph(const std::string& name,
                           const std::string& title)
  {
    TGraph* graph = new TGraph();
    graph->SetName(name.c_str());
    graph->SetTitle(title.c_str());
    return graph;
  }

  inline HistogramBook bookHistograms(const HistogramConfig& config = HistogramConfig())
  {
    HistogramBook book;
    book.config = config;

    book.mcTruthOriginT =
      make1D("hMCTruthRank0DownstreamElectronOriginT",
             "MC truth rank-0 downstream e^{-} origin time;t_{0} [ns];entries",
             config.timeBins,
             config.timeMin,
             config.timeMax);
    book.mcTruthOriginX =
      make1D("hMCTruthRank0DownstreamElectronOriginX",
             "MC truth rank-0 downstream e^{-} origin x;x_{0} [mm];entries",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax);
    book.mcTruthOriginY =
      make1D("hMCTruthRank0DownstreamElectronOriginY",
             "MC truth rank-0 downstream e^{-} origin y;y_{0} [mm];entries",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax);
    book.mcTruthOriginZ =
      make1D("hMCTruthRank0DownstreamElectronOriginZ",
             "MC truth rank-0 downstream e^{-} origin z;z_{0} [mm];entries",
             config.zBins,
             config.zMin,
             config.zMax);
    book.mcTruthOriginXY =
      make2D("hMCTruthRank0DownstreamElectronOriginXY",
             "MC truth rank-0 downstream e^{-} origin;x_{0} [mm];y_{0} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.transverseBins,
             config.transverseMin,
             config.transverseMax);
    book.mcTruthOriginXZ =
      make2D("hMCTruthRank0DownstreamElectronOriginXZ",
             "MC truth rank-0 downstream e^{-} origin;x_{0} [mm];z_{0} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.zBins,
             config.zMin,
             config.zMax);
    book.mcTruthOriginYZ =
      make2D("hMCTruthRank0DownstreamElectronOriginYZ",
             "MC truth rank-0 downstream e^{-} origin;y_{0} [mm];z_{0} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.zBins,
             config.zMin,
             config.zMax);
    book.mcTruthMomentum =
      make1D("hMCTruthRank0DownstreamElectronMomentum",
             "MC truth rank-0 downstream e^{-} origin momentum;|p_{true}| [MeV/c];entries",
             config.momentumBins,
             config.momentumMin,
             config.momentumMax);

    book.recoMomentum =
      make1D("hRecoSelectedDownstreamElectronMomentum",
             "Selected reconstructed downstream e^{-} momentum;|p_{reco}| [MeV/c];entries",
             config.momentumBins,
             config.momentumMin,
             config.momentumMax);

    book.recoVertexX =
      make1D("hRecoTwoElectronVertexX",
             "Selected two-track reconstructed vertex x;x_{vtx} [mm];entries",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax);
    book.recoVertexY =
      make1D("hRecoTwoElectronVertexY",
             "Selected two-track reconstructed vertex y;y_{vtx} [mm];entries",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax);
    book.recoVertexZ =
      make1D("hRecoTwoElectronVertexZ",
             "Selected two-track reconstructed vertex z;z_{vtx} [mm];entries",
             config.zBins,
             config.zMin,
             config.zMax);
    book.recoVertexXY =
      make2D("hRecoTwoElectronVertexXY",
             "Selected two-track reconstructed vertex;x_{vtx} [mm];y_{vtx} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.transverseBins,
             config.transverseMin,
             config.transverseMax);
    book.recoVertexXZ =
      make2D("hRecoTwoElectronVertexXZ",
             "Selected two-track reconstructed vertex;x_{vtx} [mm];z_{vtx} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.zBins,
             config.zMin,
             config.zMax);
    book.recoVertexYZ =
      make2D("hRecoTwoElectronVertexYZ",
             "Selected two-track reconstructed vertex;y_{vtx} [mm];z_{vtx} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.zBins,
             config.zMin,
             config.zMax);
    book.recoVertexFoilIndexMatchedXY =
      make2D("hRecoTwoElectronVertexFoilIndexMatchedXY",
             "Foil-index matched selected two-track reconstructed vertex;x_{vtx} [mm];y_{vtx} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.transverseBins,
             config.transverseMin,
             config.transverseMax);
    book.recoVertexFoilIndexMatchedXZ =
      make2D("hRecoTwoElectronVertexFoilIndexMatchedXZ",
             "Foil-index matched selected two-track reconstructed vertex;x_{vtx} [mm];z_{vtx} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.zBins,
             config.zMin,
             config.zMax);
    book.recoVertexFoilIndexMatchedYZ =
      make2D("hRecoTwoElectronVertexFoilIndexMatchedYZ",
             "Foil-index matched selected two-track reconstructed vertex;y_{vtx} [mm];z_{vtx} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.zBins,
             config.zMin,
             config.zMax);
    book.recoVertexLineDistance =
      make1D("hRecoTwoElectronVertexLineDistance",
             "Selected two-track closest-line distance;line-line distance [mm];entries",
             config.vertexDistanceBins,
             config.vertexDistanceMin,
             config.vertexDistanceMax);
    book.recoVertexLineParameterS =
      make1D("hRecoTwoElectronVertexLineParameterS",
             "Selected two-track closest-approach line parameter s;s [mm];entries",
             config.lineParameterBins,
             config.lineParameterMin,
             config.lineParameterMax);
    book.recoVertexLineParameterT =
      make1D("hRecoTwoElectronVertexLineParameterT",
             "Selected two-track closest-approach line parameter t;t [mm];entries",
             config.lineParameterBins,
             config.lineParameterMin,
             config.lineParameterMax);
    book.recoVertexAbsLineParameterS =
      make1D("hRecoTwoElectronVertexAbsLineParameterS",
             "Selected two-track closest-approach absolute line parameter |s|;|s| [mm];entries",
             config.absLineParameterBins,
             config.absLineParameterMin,
             config.absLineParameterMax);
    book.recoVertexAbsLineParameterT =
      make1D("hRecoTwoElectronVertexAbsLineParameterT",
             "Selected two-track closest-approach absolute line parameter |t|;|t| [mm];entries",
             config.absLineParameterBins,
             config.absLineParameterMin,
             config.absLineParameterMax);
    book.recoVertexAbsLineParameterST =
      make2D("hRecoTwoElectronVertexAbsLineParameterST",
             "Selected two-track closest-approach absolute line parameters;|s| [mm];|t| [mm]",
             config.absLineParameterBins,
             config.absLineParameterMin,
             config.absLineParameterMax,
             config.absLineParameterBins,
             config.absLineParameterMin,
             config.absLineParameterMax);
    book.recoVertexSelectedSegmentDeltaTTest =
      make1D("hTESTRecoTwoElectronVertexSelectedSegmentDeltaT",
             "TEST: selected two-track segment time difference;#Delta t [ns];entries",
             config.timeDifferenceBins,
             config.timeDifferenceMin,
             config.timeDifferenceMax);
    book.recoVertexMinTimeDifferenceTest =
      make1D("hTESTRecoTwoElectronVertexMinTimeDifference",
             "TEST: minimum-|#Delta t| shared ST_Foils pair time difference;#Delta t [ns];entries",
             config.timeDifferenceBins,
             config.timeDifferenceMin,
             config.timeDifferenceMax);
    book.recoAllSharedFoilCandidateMaxLvsDeltaZ =
      makeGraph("gRecoAllSharedFoilCandidateMaxLvsDeltaZ",
                "All shared same-foil candidate pairs;max(L_{1}, L_{2}) [mm];#Delta z_{reco-truth} [mm]");
    book.recoSpaceSelectedSharedFoilMaxLvsDeltaZ =
      makeGraph("gRecoSpaceSelectedSharedFoilMaxLvsDeltaZ",
                "Space-selected shared same-foil pair;max(L_{1}, L_{2}) [mm];#Delta z_{reco-truth} [mm]");
    book.testRecoVertexMinTimeX =
      make1D("hTESTRecoTwoElectronVertexMinTimeX",
             "TEST: minimum-|#Delta t| shared ST_Foils reconstructed vertex x;x_{vtx} [mm];entries",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax);
    book.testRecoVertexMinTimeY =
      make1D("hTESTRecoTwoElectronVertexMinTimeY",
             "TEST: minimum-|#Delta t| shared ST_Foils reconstructed vertex y;y_{vtx} [mm];entries",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax);
    book.testRecoVertexMinTimeZ =
      make1D("hTESTRecoTwoElectronVertexMinTimeZ",
             "TEST: minimum-|#Delta t| shared ST_Foils reconstructed vertex z;z_{vtx} [mm];entries",
             config.zBins,
             config.zMin,
             config.zMax);
    book.testRecoVertexMinTimeXY =
      make2D("hTESTRecoTwoElectronVertexMinTimeXY",
             "TEST: minimum-|#Delta t| shared ST_Foils reconstructed vertex;x_{vtx} [mm];y_{vtx} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.transverseBins,
             config.transverseMin,
             config.transverseMax);
    book.testRecoVertexMinTimeXZ =
      make2D("hTESTRecoTwoElectronVertexMinTimeXZ",
             "TEST: minimum-|#Delta t| shared ST_Foils reconstructed vertex;x_{vtx} [mm];z_{vtx} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.zBins,
             config.zMin,
             config.zMax);
    book.testRecoVertexMinTimeYZ =
      make2D("hTESTRecoTwoElectronVertexMinTimeYZ",
             "TEST: minimum-|#Delta t| shared ST_Foils reconstructed vertex;y_{vtx} [mm];z_{vtx} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.zBins,
             config.zMin,
             config.zMax);
    book.testRecoVertexMinTimeFoilIndexMatchedXY =
      make2D("hTESTRecoTwoElectronVertexMinTimeFoilIndexMatchedXY",
             "TEST: foil-index matched minimum-|#Delta t| reconstructed vertex;x_{vtx} [mm];y_{vtx} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.transverseBins,
             config.transverseMin,
             config.transverseMax);
    book.testRecoVertexMinTimeFoilIndexMatchedXZ =
      make2D("hTESTRecoTwoElectronVertexMinTimeFoilIndexMatchedXZ",
             "TEST: foil-index matched minimum-|#Delta t| reconstructed vertex;x_{vtx} [mm];z_{vtx} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.zBins,
             config.zMin,
             config.zMax);
    book.testRecoVertexMinTimeFoilIndexMatchedYZ =
      make2D("hTESTRecoTwoElectronVertexMinTimeFoilIndexMatchedYZ",
             "TEST: foil-index matched minimum-|#Delta t| reconstructed vertex;y_{vtx} [mm];z_{vtx} [mm]",
             config.transverseBins,
             config.transverseMin,
             config.transverseMax,
             config.zBins,
             config.zMin,
             config.zMax);
    book.testRecoMinTimeSharedFoilMaxLvsDeltaZ =
      makeGraph("gTESTRecoMinTimeSharedFoilMaxLvsDeltaZ",
                "TEST: minimum-|#Delta t| shared same-foil pair;max(L_{1}, L_{2}) [mm];#Delta z_{reco-truth} [mm]");
    book.testRecoVertexMinTimeLineDistance =
      make1D("hTESTRecoTwoElectronVertexMinTimeLineDistance",
             "TEST: minimum-|#Delta t| shared ST_Foils closest-line distance;line-line distance [mm];entries",
             config.vertexDistanceBins,
             config.vertexDistanceMin,
             config.vertexDistanceMax);
    book.recoVertexTruthDeltaX =
      make1D("hRecoTruthVertexDeltaX",
             "Selected two-track vertex minus truth origin x difference;#Delta x [mm];entries",
             config.transverseBins,
             -1000.0,
             1000.0);
    book.recoVertexTruthDeltaY =
      make1D("hRecoTruthVertexDeltaY",
             "Selected two-track vertex minus truth origin y difference;#Delta y [mm];entries",
             config.transverseBins,
             -1000.0,
             1000.0);
    book.recoVertexTruthDeltaZ =
      make1D("hRecoTruthVertexDeltaZ",
             "Selected two-track vertex minus truth origin z difference;#Delta z [mm];entries",
             config.transverseBins,
             -1000.0,
             1000.0);
    book.recoVertexTruthDistance =
      make1D("hRecoTruthVertexDistance",
             "Selected two-track vertex minus truth origin residual magnitude;|#Delta r| [mm];entries",
             config.vertexDistanceBins,
             config.vertexDistanceMin,
             1000.0);
    book.recoVertexFoilIndexMatchedTruthDeltaX =
      make1D("hRecoTruthVertexFoilIndexMatchedDeltaX",
             "Foil-index matched vertex minus truth origin x difference;#Delta x [mm];entries",
             config.transverseBins,
             -1000.0,
             1000.0);
    book.recoVertexFoilIndexMatchedTruthDeltaY =
      make1D("hRecoTruthVertexFoilIndexMatchedDeltaY",
             "Foil-index matched vertex minus truth origin y difference;#Delta y [mm];entries",
             config.transverseBins,
             -1000.0,
             1000.0);
    book.recoVertexFoilIndexMatchedTruthDeltaZ =
      make1D("hRecoTruthVertexFoilIndexMatchedDeltaZ",
             "Foil-index matched vertex minus truth origin z difference;#Delta z [mm];entries",
             config.transverseBins,
             -1000.0,
             1000.0);
    book.recoVertexFoilIndexMatchedTruthDistance =
      make1D("hRecoTruthVertexFoilIndexMatchedDistance",
             "Foil-index matched vertex minus truth origin residual magnitude;|#Delta r| [mm];entries",
             config.vertexDistanceBins,
             config.vertexDistanceMin,
             1000.0);
    book.testRecoVertexMinTimeTruthDeltaX =
      make1D("hTESTRecoTruthVertexMinTimeDeltaX",
             "TEST: minimum-|#Delta t| vertex minus truth origin x difference;#Delta x [mm];entries",
             config.transverseBins,
             -1000.0,
             1000.0);
    book.testRecoVertexMinTimeTruthDeltaY =
      make1D("hTESTRecoTruthVertexMinTimeDeltaY",
             "TEST: minimum-|#Delta t| vertex minus truth origin y difference;#Delta y [mm];entries",
             config.transverseBins,
             -1000.0,
             1000.0);
    book.testRecoVertexMinTimeTruthDeltaZ =
      make1D("hTESTRecoTruthVertexMinTimeDeltaZ",
             "TEST: minimum-|#Delta t| vertex minus truth origin z difference;#Delta z [mm];entries",
             config.transverseBins,
             -1000.0,
             1000.0);
    book.testRecoVertexMinTimeTruthDistance =
      make1D("hTESTRecoTruthVertexMinTimeDistance",
             "TEST: minimum-|#Delta t| vertex minus truth origin residual magnitude;|#Delta r| [mm];entries",
             config.vertexDistanceBins,
             config.vertexDistanceMin,
             1000.0);
    book.testRecoVertexMinTimeFoilIndexMatchedTruthDeltaX =
      make1D("hTESTRecoTruthVertexMinTimeFoilIndexMatchedDeltaX",
             "TEST: foil-index matched minimum-|#Delta t| vertex minus truth origin x difference;#Delta x [mm];entries",
             config.transverseBins,
             -1000.0,
             1000.0);
    book.testRecoVertexMinTimeFoilIndexMatchedTruthDeltaY =
      make1D("hTESTRecoTruthVertexMinTimeFoilIndexMatchedDeltaY",
             "TEST: foil-index matched minimum-|#Delta t| vertex minus truth origin y difference;#Delta y [mm];entries",
             config.transverseBins,
             -1000.0,
             1000.0);
    book.testRecoVertexMinTimeFoilIndexMatchedTruthDeltaZ =
      make1D("hTESTRecoTruthVertexMinTimeFoilIndexMatchedDeltaZ",
             "TEST: foil-index matched minimum-|#Delta t| vertex minus truth origin z difference;#Delta z [mm];entries",
             config.transverseBins,
             -1000.0,
             1000.0);
    book.testRecoVertexMinTimeFoilIndexMatchedTruthDistance =
      make1D("hTESTRecoTruthVertexMinTimeFoilIndexMatchedDistance",
             "TEST: foil-index matched minimum-|#Delta t| vertex minus truth origin residual magnitude;|#Delta r| [mm];entries",
             config.vertexDistanceBins,
             config.vertexDistanceMin,
             1000.0);

    return book;
  }

  inline bool isRankZeroDownstreamElectron(const mu2e::SimInfo& sim,
                                           int electronPdg = 11)
  {
    return sim.valid &&
           sim.rank == 0 &&
           sim.pdg == electronPdg &&
           sim.mom.z() > 0.0;
  }

  inline void fillMCTruthParticle(HistogramBook& book, const mu2e::SimInfo& sim)
  {
    book.mcTruthOriginT->Fill(sim.time);
    book.mcTruthOriginX->Fill(sim.pos.x());
    book.mcTruthOriginY->Fill(sim.pos.y());
    book.mcTruthOriginZ->Fill(sim.pos.z());
    book.mcTruthOriginXY->Fill(sim.pos.x(), sim.pos.y());
    book.mcTruthOriginXZ->Fill(sim.pos.x(), sim.pos.z());
    book.mcTruthOriginYZ->Fill(sim.pos.y(), sim.pos.z());
    book.mcTruthMomentum->Fill(sim.mom.R());
    ++book.mcTruthEntries;
  }

  inline int fillRankZeroDownstreamElectronTruthFromSelectedTracks(
    HistogramBook& book,
    const std::vector<std::vector<mu2e::SimInfo>>* truthSimByTrack,
    const std::vector<size_t>& selectedTrackIndices,
    int electronPdg = 11)
  {
    if (truthSimByTrack == nullptr)
    {
      return 0;
    }

    int filled = 0;
    std::set<int> filledSimIds;

    for (const size_t trackIndex : selectedTrackIndices)
    {
      if (trackIndex >= truthSimByTrack->size())
      {
        continue;
      }

      const auto& truthSims = truthSimByTrack->at(trackIndex);
      for (const auto& sim : truthSims)
      {
        if (!isRankZeroDownstreamElectron(sim, electronPdg))
        {
          continue;
        }

        // trkmcsim is organized per reconstructed track.  If two selected
        // tracks point back to the same SimParticle, avoid double-counting that
        // truth particle in event-level truth histograms.
        if (sim.id >= 0 && filledSimIds.find(sim.id) != filledSimIds.end())
        {
          continue;
        }
        if (sim.id >= 0)
        {
          filledSimIds.insert(sim.id);
        }

        fillMCTruthParticle(book, sim);
        ++filled;
      }
    }

    return filled;
  }

  inline void fillRecoMomentum(HistogramBook& book, double momentum)
  {
    if (momentum <= 0.0)
    {
      return;
    }

    book.recoMomentum->Fill(momentum);
    ++book.recoMomentumEntries;
  }

  inline void fillRecoVertex(HistogramBook& book,
                             const twoparticlevertexer::VertexResult& vertex)
  {
    if (!vertex.valid)
    {
      return;
    }

    book.recoVertexX->Fill(vertex.vertex.x());
    book.recoVertexY->Fill(vertex.vertex.y());
    book.recoVertexZ->Fill(vertex.vertex.z());
    book.recoVertexXY->Fill(vertex.vertex.x(), vertex.vertex.y());
    book.recoVertexXZ->Fill(vertex.vertex.x(), vertex.vertex.z());
    book.recoVertexYZ->Fill(vertex.vertex.y(), vertex.vertex.z());
    book.recoVertexLineDistance->Fill(vertex.distance);

    const double lineParameterS = vertex.firstLineParameter;
    const double lineParameterT = vertex.secondLineParameter;
    if (std::isfinite(lineParameterS) && std::isfinite(lineParameterT))
    {
      const double absLineParameterS = std::fabs(lineParameterS);
      const double absLineParameterT = std::fabs(lineParameterT);
      book.recoVertexLineParameterS->Fill(lineParameterS);
      book.recoVertexLineParameterT->Fill(lineParameterT);
      book.recoVertexAbsLineParameterS->Fill(absLineParameterS);
      book.recoVertexAbsLineParameterT->Fill(absLineParameterT);
      book.recoVertexAbsLineParameterST->Fill(absLineParameterS,
                                              absLineParameterT);
      ++book.recoVertexLineParameterEntries;
    }

    ++book.recoVertexEntries;
  }

  inline void fillRecoVertexFoilIndexMatchedMaps(
    HistogramBook& book,
    const twoparticlevertexer::VertexResult& vertex)
  {
    if (!vertex.valid)
    {
      return;
    }

    book.recoVertexFoilIndexMatchedXY->Fill(vertex.vertex.x(), vertex.vertex.y());
    book.recoVertexFoilIndexMatchedXZ->Fill(vertex.vertex.x(), vertex.vertex.z());
    book.recoVertexFoilIndexMatchedYZ->Fill(vertex.vertex.y(), vertex.vertex.z());
    ++book.recoVertexFoilIndexMatchedEntries;
  }

  inline void fillRecoVertexSelectedSegmentTimeDifference(
    HistogramBook& book,
    const twoparticlevertexer::VertexResult& vertex)
  {
    if (!vertex.valid)
    {
      return;
    }

    book.recoVertexSelectedSegmentDeltaTTest->Fill(vertex.deltaInputTime);
    ++book.recoVertexSelectedSegmentDeltaTEntries;
  }

  inline void fillRecoVertexMinTimeDifferenceTest(
    HistogramBook& book,
    const twoparticlevertexer::VertexResult& vertex)
  {
    if (!vertex.valid)
    {
      return;
    }

    book.recoVertexMinTimeDifferenceTest->Fill(vertex.deltaInputTime);
    ++book.recoVertexMinTimeDifferenceEntries;
  }

  inline bool validMaxLvsDeltaZInputs(double maxPointTruthDistance,
                                      double recoMinusTruthZ)
  {
    return std::isfinite(maxPointTruthDistance) &&
           std::isfinite(recoMinusTruthZ) &&
           maxPointTruthDistance >= 0.0;
  }

  inline void fillRecoAllSharedFoilCandidateMaxLvsDeltaZ(
    HistogramBook& book,
    double maxPointTruthDistance,
    double recoMinusTruthZ)
  {
    if (!validMaxLvsDeltaZInputs(maxPointTruthDistance, recoMinusTruthZ))
    {
      return;
    }

    book.recoAllSharedFoilCandidateMaxLvsDeltaZ->SetPoint(
      book.recoAllSharedFoilCandidateMaxLvsDeltaZ->GetN(),
      maxPointTruthDistance,
      recoMinusTruthZ);
    ++book.recoAllSharedFoilCandidateMaxLvsDeltaZEntries;
  }

  inline void fillRecoSpaceSelectedSharedFoilMaxLvsDeltaZ(
    HistogramBook& book,
    double maxPointTruthDistance,
    double recoMinusTruthZ)
  {
    if (!validMaxLvsDeltaZInputs(maxPointTruthDistance, recoMinusTruthZ))
    {
      return;
    }

    book.recoSpaceSelectedSharedFoilMaxLvsDeltaZ->SetPoint(
      book.recoSpaceSelectedSharedFoilMaxLvsDeltaZ->GetN(),
      maxPointTruthDistance,
      recoMinusTruthZ);
    ++book.recoSpaceSelectedSharedFoilMaxLvsDeltaZEntries;
  }

  inline void fillRecoMinTimeSharedFoilMaxLvsDeltaZTest(
    HistogramBook& book,
    double maxPointTruthDistance,
    double recoMinusTruthZ)
  {
    if (!validMaxLvsDeltaZInputs(maxPointTruthDistance, recoMinusTruthZ))
    {
      return;
    }

    book.testRecoMinTimeSharedFoilMaxLvsDeltaZ->SetPoint(
      book.testRecoMinTimeSharedFoilMaxLvsDeltaZ->GetN(),
      maxPointTruthDistance,
      recoMinusTruthZ);
    ++book.testRecoMinTimeSharedFoilMaxLvsDeltaZEntries;
  }

  inline void fillRecoVertexMinTimeChoiceTest(
    HistogramBook& book,
    const twoparticlevertexer::VertexResult& vertex)
  {
    if (!vertex.valid)
    {
      return;
    }

    book.testRecoVertexMinTimeX->Fill(vertex.vertex.x());
    book.testRecoVertexMinTimeY->Fill(vertex.vertex.y());
    book.testRecoVertexMinTimeZ->Fill(vertex.vertex.z());
    book.testRecoVertexMinTimeXY->Fill(vertex.vertex.x(), vertex.vertex.y());
    book.testRecoVertexMinTimeXZ->Fill(vertex.vertex.x(), vertex.vertex.z());
    book.testRecoVertexMinTimeYZ->Fill(vertex.vertex.y(), vertex.vertex.z());
    book.testRecoVertexMinTimeLineDistance->Fill(vertex.distance);
    ++book.testRecoVertexMinTimeEntries;
  }

  inline void fillRecoVertexMinTimeFoilIndexMatchedMapsTest(
    HistogramBook& book,
    const twoparticlevertexer::VertexResult& vertex)
  {
    if (!vertex.valid)
    {
      return;
    }

    book.testRecoVertexMinTimeFoilIndexMatchedXY->Fill(vertex.vertex.x(),
                                                       vertex.vertex.y());
    book.testRecoVertexMinTimeFoilIndexMatchedXZ->Fill(vertex.vertex.x(),
                                                       vertex.vertex.z());
    book.testRecoVertexMinTimeFoilIndexMatchedYZ->Fill(vertex.vertex.y(),
                                                       vertex.vertex.z());
    ++book.testRecoVertexMinTimeFoilIndexMatchedEntries;
  }

  inline void fillRecoVertexTruthResidual(HistogramBook& book, const XYZVectorF& delta)
  {
    if (!std::isfinite(delta.x()) || !std::isfinite(delta.y()) || !std::isfinite(delta.z()))
    {
      return;
    }

    book.recoVertexTruthDeltaX->Fill(delta.x());
    book.recoVertexTruthDeltaY->Fill(delta.y());
    book.recoVertexTruthDeltaZ->Fill(delta.z());
    book.recoVertexTruthDistance->Fill(delta.R());
    ++book.recoVertexTruthResidualEntries;
  }

  inline void fillRecoVertexFoilIndexMatchedTruthResidual(HistogramBook& book,
                                                          const XYZVectorF& delta)
  {
    if (!std::isfinite(delta.x()) || !std::isfinite(delta.y()) || !std::isfinite(delta.z()))
    {
      return;
    }

    book.recoVertexFoilIndexMatchedTruthDeltaX->Fill(delta.x());
    book.recoVertexFoilIndexMatchedTruthDeltaY->Fill(delta.y());
    book.recoVertexFoilIndexMatchedTruthDeltaZ->Fill(delta.z());
    book.recoVertexFoilIndexMatchedTruthDistance->Fill(delta.R());
    ++book.recoVertexFoilIndexMatchedTruthResidualEntries;
  }

  inline void fillRecoVertexMinTimeTruthResidualTest(HistogramBook& book,
                                                     const XYZVectorF& delta)
  {
    if (!std::isfinite(delta.x()) || !std::isfinite(delta.y()) || !std::isfinite(delta.z()))
    {
      return;
    }

    book.testRecoVertexMinTimeTruthDeltaX->Fill(delta.x());
    book.testRecoVertexMinTimeTruthDeltaY->Fill(delta.y());
    book.testRecoVertexMinTimeTruthDeltaZ->Fill(delta.z());
    book.testRecoVertexMinTimeTruthDistance->Fill(delta.R());
    ++book.testRecoVertexMinTimeTruthResidualEntries;
  }

  inline void fillRecoVertexMinTimeFoilIndexMatchedTruthResidualTest(
    HistogramBook& book,
    const XYZVectorF& delta)
  {
    if (!std::isfinite(delta.x()) || !std::isfinite(delta.y()) || !std::isfinite(delta.z()))
    {
      return;
    }

    book.testRecoVertexMinTimeFoilIndexMatchedTruthDeltaX->Fill(delta.x());
    book.testRecoVertexMinTimeFoilIndexMatchedTruthDeltaY->Fill(delta.y());
    book.testRecoVertexMinTimeFoilIndexMatchedTruthDeltaZ->Fill(delta.z());
    book.testRecoVertexMinTimeFoilIndexMatchedTruthDistance->Fill(delta.R());
    ++book.testRecoVertexMinTimeFoilIndexMatchedTruthResidualEntries;
  }

  inline void writeOne(TH1* histogram)
  {
    if (histogram != nullptr)
    {
      histogram->Write();
    }
  }

  inline void writeOne(TGraph* graph)
  {
    if (graph != nullptr)
    {
      graph->Write();
    }
  }

  inline bool writeHistograms(HistogramBook& book, const std::string& outputName)
  {
    TFile* outputFile = TFile::Open(outputName.c_str(), "RECREATE");
    if (outputFile == nullptr || outputFile->IsZombie())
    {
      std::cerr << "ERROR: could not create histogram output file: "
                << outputName << std::endl;
      return false;
    }

    outputFile->cd();
    writeOne(book.mcTruthOriginT);
    writeOne(book.mcTruthOriginX);
    writeOne(book.mcTruthOriginY);
    writeOne(book.mcTruthOriginZ);
    writeOne(book.mcTruthOriginXY);
    writeOne(book.mcTruthOriginXZ);
    writeOne(book.mcTruthOriginYZ);
    writeOne(book.mcTruthMomentum);
    writeOne(book.recoMomentum);
    writeOne(book.recoVertexX);
    writeOne(book.recoVertexY);
    writeOne(book.recoVertexZ);
    writeOne(book.recoVertexXY);
    writeOne(book.recoVertexXZ);
    writeOne(book.recoVertexYZ);
    writeOne(book.recoVertexFoilIndexMatchedXY);
    writeOne(book.recoVertexFoilIndexMatchedXZ);
    writeOne(book.recoVertexFoilIndexMatchedYZ);
    writeOne(book.recoVertexLineDistance);
    writeOne(book.recoVertexLineParameterS);
    writeOne(book.recoVertexLineParameterT);
    writeOne(book.recoVertexAbsLineParameterS);
    writeOne(book.recoVertexAbsLineParameterT);
    writeOne(book.recoVertexAbsLineParameterST);
    writeOne(book.recoVertexSelectedSegmentDeltaTTest);
    writeOne(book.recoVertexMinTimeDifferenceTest);
    writeOne(book.recoAllSharedFoilCandidateMaxLvsDeltaZ);
    writeOne(book.recoSpaceSelectedSharedFoilMaxLvsDeltaZ);
    writeOne(book.testRecoVertexMinTimeX);
    writeOne(book.testRecoVertexMinTimeY);
    writeOne(book.testRecoVertexMinTimeZ);
    writeOne(book.testRecoVertexMinTimeXY);
    writeOne(book.testRecoVertexMinTimeXZ);
    writeOne(book.testRecoVertexMinTimeYZ);
    writeOne(book.testRecoVertexMinTimeFoilIndexMatchedXY);
    writeOne(book.testRecoVertexMinTimeFoilIndexMatchedXZ);
    writeOne(book.testRecoVertexMinTimeFoilIndexMatchedYZ);
    writeOne(book.testRecoMinTimeSharedFoilMaxLvsDeltaZ);
    writeOne(book.testRecoVertexMinTimeLineDistance);
    writeOne(book.recoVertexTruthDeltaX);
    writeOne(book.recoVertexTruthDeltaY);
    writeOne(book.recoVertexTruthDeltaZ);
    writeOne(book.recoVertexTruthDistance);
    writeOne(book.recoVertexFoilIndexMatchedTruthDeltaX);
    writeOne(book.recoVertexFoilIndexMatchedTruthDeltaY);
    writeOne(book.recoVertexFoilIndexMatchedTruthDeltaZ);
    writeOne(book.recoVertexFoilIndexMatchedTruthDistance);
    writeOne(book.testRecoVertexMinTimeTruthDeltaX);
    writeOne(book.testRecoVertexMinTimeTruthDeltaY);
    writeOne(book.testRecoVertexMinTimeTruthDeltaZ);
    writeOne(book.testRecoVertexMinTimeTruthDistance);
    writeOne(book.testRecoVertexMinTimeFoilIndexMatchedTruthDeltaX);
    writeOne(book.testRecoVertexMinTimeFoilIndexMatchedTruthDeltaY);
    writeOne(book.testRecoVertexMinTimeFoilIndexMatchedTruthDeltaZ);
    writeOne(book.testRecoVertexMinTimeFoilIndexMatchedTruthDistance);
    outputFile->Close();
    delete outputFile;
    return true;
  }
}

#endif

