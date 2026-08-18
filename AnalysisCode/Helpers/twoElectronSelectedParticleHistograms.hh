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
//   means the event contains exactly two reconstructed downstream electron
//   tracks total, and both have associated calorimeter hits and reconstructed
//   momentum in 50-53 MeV/c.
//   The momentum histogram range is intentionally wider: 30-55 MeV/c.
//
// Coordinate and unit conventions:
//   - positions are in mm
//   - momenta are in MeV/c
//   - times are in ns
//   - angles are in degrees
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
#include "EventNtuple/inc/TrkSegInfo.hh"
#include "twoParticleVertexer.hh"

namespace twoelectronhist
{
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kRadiansToDegrees = 180.0 / kPi;

  // Histogram ranges are kept in one config object so later tuning does not
  // touch the event-selection code.  The default position bounds match the
  // original vertexer macro so the truth and reco position plots are directly
  // comparable.
  struct HistogramConfig
  {
    int momentumBins = 140;
    double momentumMin = 30.0;
    double momentumMax = 55.0;
    int momentumComponentBins = 240;
    double momentumComponentMin = -60.0;
    double momentumComponentMax = 60.0;

    int angleBins = 180;
    double thetaMin = 0.0;
    double thetaMax = 180.0;

    int timeBins = 200;
    double timeMin = 0.0;
    double timeMax = 2000.0;
    int firstSTFoilTimeDifferenceBins = 400;
    double firstSTFoilTimeDifferenceMin = 0.0;
    double firstSTFoilTimeDifferenceMax = 2000.0;
    int timingResolutionBins = 400;
    double timingResolutionMin = -2000.0;
    double timingResolutionMax = 2000.0;

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

    int recoTrackMultiplicityBins = 20;
    double recoTrackMultiplicityMin = -0.5;
    double recoTrackMultiplicityMax = 19.5;

    double foilPointTruthDistanceMin = 0.0;
    double foilPointTruthDistanceMax = 1000.0;
    double recoTruthDeltaZMin = -1000.0;
    double recoTruthDeltaZMax = 1000.0;

    int stoppingTargetFoils = 37;
    int deltaZByFoilBins = 200;
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
    TH1F* timingSelectedEventTruthOriginTime = nullptr;
    TH1F* timingSelectedTrackFirstSTFoilTime = nullptr;
    TH1F* timingSelectedTrackFirstSTFoilDeltaT = nullptr;
    TH1F* timingTruthMinusEarlierFirstSTFoilTime = nullptr;
    TH1F* timingTruthMinusLaterFirstSTFoilTime = nullptr;
    TH1F* timingTruthMinusAverageFirstSTFoilTime = nullptr;
    TH2F* timingTruthMinusEarlierVsLaterFirstSTFoilTime = nullptr;

    TH1F* recoVertexX = nullptr;
    TH1F* recoVertexY = nullptr;
    TH1F* recoVertexZ = nullptr;
    TH2F* recoVertexXY = nullptr;
    TH2F* recoVertexXZ = nullptr;
    TH2F* recoVertexYZ = nullptr;
    TH2F* recoVertexSpaceSelectedMomentumTheta1Theta2 = nullptr;
    TH1F* recoVertexSpaceSelectedMomentumOpeningAngle = nullptr;
    TH2F* recoVertexFoilIndexMatchedXY = nullptr;
    TH2F* recoVertexFoilIndexMatchedXZ = nullptr;
    TH2F* recoVertexFoilIndexMatchedYZ = nullptr;
    TH1F* recoVertexLineDistance = nullptr;
    TH1F* recoVertexLineParameterS = nullptr;
    TH1F* recoVertexLineParameterT = nullptr;
    TH1F* recoVertexAbsLineParameterS = nullptr;
    TH1F* recoVertexAbsLineParameterT = nullptr;
    TH2F* recoVertexLineParameterST = nullptr;
    TH2F* recoVertexAbsLineParameterST = nullptr;
    TH2F* recoTrackMultiplicityVsDeltaZ = nullptr;
    TH1F* recoVertexSelectedSegmentDeltaTTest = nullptr;
    TH1F* recoVertexMinTimeDifferenceTest = nullptr;
    TH2F* testRecoVertexMinTimeDeltaTVsTruthDeltaZ = nullptr;
    TGraph* recoAllSharedFoilCandidateMaxLvsDeltaZ = nullptr;
    TGraph* recoSpaceSelectedSharedFoilMaxLvsDeltaZ = nullptr;
    TGraph* recoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZ = nullptr;
    TGraph* recoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZ = nullptr;
    TGraph* recoSpaceSelectedSharedFoilNumberVsDeltaZ = nullptr;
    std::vector<TH1F*> recoSpaceSelectedDeltaZBySelectedFoil;
    TH1F* recoSelectedDownstreamElectronTrackCountByFoil = nullptr;
    std::vector<TH1F*> recoSelectedTrackFoilIntersectionZByFoil;
    std::vector<TH1F*> recoSelectedTrackFoilIntersectionPxByFoil;
    std::vector<TH1F*> recoSelectedTrackFoilIntersectionPyByFoil;
    std::vector<TH1F*> recoSelectedTrackFoilIntersectionPzByFoil;
    std::vector<TH1F*> recoSelectedTrackFoilIntersectionPByFoil;
    TGraph* recoSpaceSelectedSharedFoilCountVsDeltaZ = nullptr;
    TGraph* recoSpaceSelectedMaxFoilsHitVsDeltaZ = nullptr;
    TGraph* recoSpaceSelectedAbsDeltaFoilsHitVsDeltaZ = nullptr;
    TGraph* recoSpaceSelectedOpeningAngleVsDeltaZ = nullptr;
    TH1F* testRecoVertexMinTimeX = nullptr;
    TH1F* testRecoVertexMinTimeY = nullptr;
    TH1F* testRecoVertexMinTimeZ = nullptr;
    TH2F* testRecoVertexMinTimeXY = nullptr;
    TH2F* testRecoVertexMinTimeXZ = nullptr;
    TH2F* testRecoVertexMinTimeYZ = nullptr;
    TH2F* testRecoVertexMinTimeMomentumTheta1Theta2 = nullptr;
    TH1F* testRecoVertexMinTimeMomentumOpeningAngle = nullptr;
    TH2F* testRecoVertexMinTimeFoilIndexMatchedXY = nullptr;
    TH2F* testRecoVertexMinTimeFoilIndexMatchedXZ = nullptr;
    TH2F* testRecoVertexMinTimeFoilIndexMatchedYZ = nullptr;
    TGraph* testRecoMinTimeSharedFoilMaxLvsDeltaZ = nullptr;
    TGraph* testRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZ = nullptr;
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
    Long64_t timingSelectedEventTruthOriginTimeEntries = 0;
    Long64_t timingSelectedTrackFirstSTFoilTimeEntries = 0;
    Long64_t timingSelectedTrackFirstSTFoilDeltaTEntries = 0;
    Long64_t timingTruthMinusFirstSTFoilEntries = 0;
    Long64_t recoVertexEntries = 0;
    Long64_t recoVertexMomentumThetaEntries = 0;
    Long64_t recoVertexMomentumOpeningAngleEntries = 0;
    Long64_t recoVertexFoilIndexMatchedEntries = 0;
    Long64_t recoVertexLineParameterEntries = 0;
    Long64_t recoTrackMultiplicityVsDeltaZEntries = 0;
    Long64_t recoVertexSelectedSegmentDeltaTEntries = 0;
    Long64_t recoVertexMinTimeDifferenceEntries = 0;
    Long64_t testRecoVertexMinTimeDeltaTVsTruthDeltaZEntries = 0;
    Long64_t recoAllSharedFoilCandidateMaxLvsDeltaZEntries = 0;
    Long64_t recoSpaceSelectedSharedFoilMaxLvsDeltaZEntries = 0;
    Long64_t recoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZEntries = 0;
    Long64_t recoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZEntries = 0;
    Long64_t recoSpaceSelectedSharedFoilNumberVsDeltaZEntries = 0;
    std::vector<Long64_t> recoSpaceSelectedDeltaZBySelectedFoilEntries;
    Long64_t recoSelectedDownstreamElectronTrackCountByFoilEntries = 0;
    std::vector<Long64_t> recoSelectedTrackFoilIntersectionZByFoilEntries;
    Long64_t recoSelectedTrackFoilIntersectionZEntries = 0;
    std::vector<Long64_t> recoSelectedTrackFoilIntersectionMomentumByFoilEntries;
    Long64_t recoSelectedTrackFoilIntersectionMomentumEntries = 0;
    Long64_t recoSpaceSelectedSharedFoilCountVsDeltaZEntries = 0;
    Long64_t recoSpaceSelectedMaxFoilsHitVsDeltaZEntries = 0;
    Long64_t recoSpaceSelectedAbsDeltaFoilsHitVsDeltaZEntries = 0;
    Long64_t recoSpaceSelectedOpeningAngleVsDeltaZEntries = 0;
    Long64_t testRecoVertexMinTimeEntries = 0;
    Long64_t testRecoVertexMinTimeMomentumThetaEntries = 0;
    Long64_t testRecoVertexMinTimeMomentumOpeningAngleEntries = 0;
    Long64_t testRecoVertexMinTimeFoilIndexMatchedEntries = 0;
    Long64_t testRecoMinTimeSharedFoilMaxLvsDeltaZEntries = 0;
    Long64_t testRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZEntries = 0;
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
    book.timingSelectedEventTruthOriginTime =
      make1D("hTimingSelectedEventTruthOriginTime",
             "Selected event MC truth origin time;t_{truth} [ns];events",
             config.timeBins,
             config.timeMin,
             config.timeMax);
    book.timingSelectedTrackFirstSTFoilTime =
      make1D("hTimingSelectedTrackFirstSTFoilTime",
             "Selected track earliest ST_Foils intersection time;t_{first ST} [ns];tracks",
             config.timeBins,
             config.timeMin,
             config.timeMax);
    book.timingSelectedTrackFirstSTFoilDeltaT =
      make1D("hTimingSelectedTrackFirstSTFoilDeltaT",
             "Selected event earliest ST_Foils time separation;t_{later first ST} - t_{earlier first ST} [ns];events",
             config.firstSTFoilTimeDifferenceBins,
             config.firstSTFoilTimeDifferenceMin,
             config.firstSTFoilTimeDifferenceMax);
    book.timingTruthMinusEarlierFirstSTFoilTime =
      make1D("hTimingTruthMinusEarlierFirstSTFoilTime",
             "Selected event timing residual for earlier first ST_Foils track;t_{truth} - t_{earlier first ST} [ns];events",
             config.timingResolutionBins,
             config.timingResolutionMin,
             config.timingResolutionMax);
    book.timingTruthMinusLaterFirstSTFoilTime =
      make1D("hTimingTruthMinusLaterFirstSTFoilTime",
             "Selected event timing residual for later first ST_Foils track;t_{truth} - t_{later first ST} [ns];events",
             config.timingResolutionBins,
             config.timingResolutionMin,
             config.timingResolutionMax);
    book.timingTruthMinusAverageFirstSTFoilTime =
      make1D("hTimingTruthMinusAverageFirstSTFoilTime",
             "Selected event timing residual using average first ST_Foils time;t_{truth} - #LTt_{first ST}#GT [ns];events",
             config.timingResolutionBins,
             config.timingResolutionMin,
             config.timingResolutionMax);
    book.timingTruthMinusEarlierVsLaterFirstSTFoilTime =
      make2D("hTimingTruthMinusEarlierVsLaterFirstSTFoilTime",
             "Selected event timing residuals;t_{truth} - t_{earlier first ST} [ns];t_{truth} - t_{later first ST} [ns]",
             config.timingResolutionBins,
             config.timingResolutionMin,
             config.timingResolutionMax,
             config.timingResolutionBins,
             config.timingResolutionMin,
             config.timingResolutionMax);

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
    book.recoVertexSpaceSelectedMomentumTheta1Theta2 =
      make2D("hRecoTwoElectronVertexSpaceSelectedMomentumTheta1Theta2",
             "Space-selected shared ST_Foils momentum polar angles;#theta_{1} [deg];#theta_{2} [deg]",
             config.angleBins,
             config.thetaMin,
             config.thetaMax,
             config.angleBins,
             config.thetaMin,
             config.thetaMax);
    book.recoVertexSpaceSelectedMomentumOpeningAngle =
      make1D("hRecoTwoElectronVertexSpaceSelectedMomentumOpeningAngle",
             "Space-selected shared ST_Foils momentum opening angle;opening angle [deg];entries",
             config.angleBins,
             config.thetaMin,
             config.thetaMax);
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
    book.recoVertexLineParameterST =
      make2D("hRecoTwoElectronVertexLineParameterST",
             "Selected two-track closest-approach line parameters;s [mm];t [mm]",
             config.lineParameterBins,
             config.lineParameterMin,
             config.lineParameterMax,
             config.lineParameterBins,
             config.lineParameterMin,
             config.lineParameterMax);
    book.recoVertexAbsLineParameterST =
      make2D("hRecoTwoElectronVertexAbsLineParameterST",
             "Selected two-track closest-approach absolute line parameters;|s| [mm];|t| [mm]",
             config.absLineParameterBins,
             config.absLineParameterMin,
             config.absLineParameterMax,
             config.absLineParameterBins,
             config.absLineParameterMin,
             config.absLineParameterMax);
    book.recoTrackMultiplicityVsDeltaZ =
      make2D("hRecoTrackMultiplicityVsDeltaZ",
             "Selected events;reconstructed tracks in event;#Delta z_{reco-truth} [mm]",
             config.recoTrackMultiplicityBins,
             config.recoTrackMultiplicityMin,
             config.recoTrackMultiplicityMax,
             config.deltaZByFoilBins,
             config.recoTruthDeltaZMin,
             config.recoTruthDeltaZMax);
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
    book.testRecoVertexMinTimeDeltaTVsTruthDeltaZ =
      make2D("hTESTRecoTwoElectronVertexMinTimeDeltaTVsTruthDeltaZ",
             "TEST: time-selected shared ST_Foils pair;#Delta t_{selected} [ns];z_{reco}-z_{truth} [mm]",
             config.timeDifferenceBins,
             config.timeDifferenceMin,
             config.timeDifferenceMax,
             config.deltaZByFoilBins,
             config.recoTruthDeltaZMin,
             config.recoTruthDeltaZMax);
    book.recoAllSharedFoilCandidateMaxLvsDeltaZ =
      makeGraph("gRecoAllSharedFoilCandidateMaxLvsDeltaZ",
                "All shared same-foil candidate pairs;max(L_{1}, L_{2}) [mm];#Delta z_{reco-truth} [mm]");
    book.recoSpaceSelectedSharedFoilMaxLvsDeltaZ =
      makeGraph("gRecoSpaceSelectedSharedFoilMaxLvsDeltaZ",
                "Space-selected shared same-foil pair;max(L_{1}, L_{2}) [mm];#Delta z_{reco-truth} [mm]");
    book.recoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZ =
      makeGraph("gRecoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZ",
                "All shared same-foil candidate pairs;(|s|+|t|)/2 [mm];#Delta z_{reco-truth} [mm]");
    book.recoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZ =
      makeGraph("gRecoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZ",
                "Space-selected shared same-foil pair;(|s|+|t|)/2 [mm];#Delta z_{reco-truth} [mm]");
    book.recoSpaceSelectedSharedFoilNumberVsDeltaZ =
      makeGraph("gRecoSpaceSelectedSharedFoilNumberVsDeltaZ",
                "Space-selected shared same-foil pair;shared foil sindex;#Delta z_{reco-truth} [mm]");
    book.recoSpaceSelectedDeltaZBySelectedFoil.reserve(config.stoppingTargetFoils);
    book.recoSpaceSelectedDeltaZBySelectedFoilEntries.assign(
      config.stoppingTargetFoils,
      0);
    for (int foilIndex = 0; foilIndex < config.stoppingTargetFoils; ++foilIndex)
    {
      book.recoSpaceSelectedDeltaZBySelectedFoil.push_back(
        make1D("hRecoSpaceSelectedDeltaZSelectedFoil" +
                 std::to_string(foilIndex),
               "Space-selected shared foil sindex " +
                 std::to_string(foilIndex) +
                 ";#Delta z_{reco-truth} [mm];entries",
               config.deltaZByFoilBins,
               config.recoTruthDeltaZMin,
               config.recoTruthDeltaZMax));
    }
    book.recoSelectedDownstreamElectronTrackCountByFoil =
      make1D("hRecoSelectedDownstreamElectronTrackCountByFoil",
             "Selected downstream e^{-} tracks by ST_Foils sindex;foil sindex;selected track count",
             config.stoppingTargetFoils,
             -0.5,
             static_cast<double>(config.stoppingTargetFoils) - 0.5);
    book.recoSelectedTrackFoilIntersectionZByFoil.reserve(
      config.stoppingTargetFoils);
    book.recoSelectedTrackFoilIntersectionZByFoilEntries.assign(
      config.stoppingTargetFoils,
      0);
    for (int foilIndex = 0; foilIndex < config.stoppingTargetFoils; ++foilIndex)
    {
      book.recoSelectedTrackFoilIntersectionZByFoil.push_back(
        make1D("hRecoSelectedTrackFoilIntersectionZFoil" +
                 std::to_string(foilIndex),
               "Selected downstream e^{-} track ST_Foils intersection z, foil sindex " +
                 std::to_string(foilIndex) +
                 ";z_{intersection} [mm];entries",
               config.zBins,
               config.zMin,
               config.zMax));
    }
    book.recoSelectedTrackFoilIntersectionPxByFoil.reserve(
      config.stoppingTargetFoils);
    book.recoSelectedTrackFoilIntersectionPyByFoil.reserve(
      config.stoppingTargetFoils);
    book.recoSelectedTrackFoilIntersectionPzByFoil.reserve(
      config.stoppingTargetFoils);
    book.recoSelectedTrackFoilIntersectionPByFoil.reserve(
      config.stoppingTargetFoils);
    book.recoSelectedTrackFoilIntersectionMomentumByFoilEntries.assign(
      config.stoppingTargetFoils,
      0);
    for (int foilIndex = 0; foilIndex < config.stoppingTargetFoils; ++foilIndex)
    {
      book.recoSelectedTrackFoilIntersectionPxByFoil.push_back(
        make1D("hRecoSelectedTrackFoilIntersectionPxFoil" +
                 std::to_string(foilIndex),
               "Selected downstream e^{-} track ST_Foils intersection p_{x}, foil sindex " +
                 std::to_string(foilIndex) +
                 ";p_{x} [MeV/c];entries",
               config.momentumComponentBins,
               config.momentumComponentMin,
               config.momentumComponentMax));
      book.recoSelectedTrackFoilIntersectionPyByFoil.push_back(
        make1D("hRecoSelectedTrackFoilIntersectionPyFoil" +
                 std::to_string(foilIndex),
               "Selected downstream e^{-} track ST_Foils intersection p_{y}, foil sindex " +
                 std::to_string(foilIndex) +
                 ";p_{y} [MeV/c];entries",
               config.momentumComponentBins,
               config.momentumComponentMin,
               config.momentumComponentMax));
      book.recoSelectedTrackFoilIntersectionPzByFoil.push_back(
        make1D("hRecoSelectedTrackFoilIntersectionPzFoil" +
                 std::to_string(foilIndex),
               "Selected downstream e^{-} track ST_Foils intersection p_{z}, foil sindex " +
                 std::to_string(foilIndex) +
                 ";p_{z} [MeV/c];entries",
               config.momentumComponentBins,
               config.momentumComponentMin,
               config.momentumComponentMax));
      book.recoSelectedTrackFoilIntersectionPByFoil.push_back(
        make1D("hRecoSelectedTrackFoilIntersectionPFoil" +
                 std::to_string(foilIndex),
               "Selected downstream e^{-} track ST_Foils intersection |p|, foil sindex " +
                 std::to_string(foilIndex) +
                 ";|p| [MeV/c];entries",
               config.momentumBins,
               config.momentumMin,
               config.momentumMax));
    }
    book.recoSpaceSelectedSharedFoilCountVsDeltaZ =
      makeGraph("gRecoSpaceSelectedSharedFoilCountVsDeltaZ",
                "Space-selected track pair;number of shared ST_Foils indices;#Delta z_{reco-truth} [mm]");
    book.recoSpaceSelectedMaxFoilsHitVsDeltaZ =
      makeGraph("gRecoSpaceSelectedMaxFoilsHitVsDeltaZ",
                "Space-selected shared same-foil pair;max unique ST_Foils hit by either track;#Delta z_{reco-truth} [mm]");
    book.recoSpaceSelectedAbsDeltaFoilsHitVsDeltaZ =
      makeGraph("gRecoSpaceSelectedAbsDeltaFoilsHitVsDeltaZ",
                "Space-selected shared same-foil pair;|N_{foils,1}-N_{foils,2}|;#Delta z_{reco-truth} [mm]");
    book.recoSpaceSelectedOpeningAngleVsDeltaZ =
      makeGraph("gRecoSpaceSelectedOpeningAngleVsDeltaZ",
                "Space-selected shared same-foil pair;opening angle [deg];#Delta z_{reco-truth} [mm]");
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
    book.testRecoVertexMinTimeMomentumTheta1Theta2 =
      make2D("hTESTRecoTwoElectronVertexMinTimeMomentumTheta1Theta2",
             "TEST: minimum-|#Delta t| shared ST_Foils momentum polar angles;#theta_{1} [deg];#theta_{2} [deg]",
             config.angleBins,
             config.thetaMin,
             config.thetaMax,
             config.angleBins,
             config.thetaMin,
             config.thetaMax);
    book.testRecoVertexMinTimeMomentumOpeningAngle =
      make1D("hTESTRecoTwoElectronVertexMinTimeMomentumOpeningAngle",
             "TEST: minimum-|#Delta t| shared ST_Foils momentum opening angle;opening angle [deg];entries",
             config.angleBins,
             config.thetaMin,
             config.thetaMax);
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
    book.testRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZ =
      makeGraph("gTESTRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZ",
                "TEST: minimum-|#Delta t| shared same-foil pair;(|s|+|t|)/2 [mm];#Delta z_{reco-truth} [mm]");
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

  inline void fillSelectedEventTimingDiagnostics(
    HistogramBook& book,
    double truthOriginTime,
    double firstTrackFirstSTFoilTime,
    double secondTrackFirstSTFoilTime)
  {
    if (!std::isfinite(truthOriginTime) ||
        !std::isfinite(firstTrackFirstSTFoilTime) ||
        !std::isfinite(secondTrackFirstSTFoilTime))
    {
      return;
    }

    const double earlierFirstSTFoilTime =
      std::min(firstTrackFirstSTFoilTime, secondTrackFirstSTFoilTime);
    const double laterFirstSTFoilTime =
      std::max(firstTrackFirstSTFoilTime, secondTrackFirstSTFoilTime);
    const double averageFirstSTFoilTime =
      0.5 * (firstTrackFirstSTFoilTime + secondTrackFirstSTFoilTime);
    const double firstSTFoilDeltaT =
      laterFirstSTFoilTime - earlierFirstSTFoilTime;
    const double truthMinusEarlierFirstSTFoilTime =
      truthOriginTime - earlierFirstSTFoilTime;
    const double truthMinusLaterFirstSTFoilTime =
      truthOriginTime - laterFirstSTFoilTime;
    const double truthMinusAverageFirstSTFoilTime =
      truthOriginTime - averageFirstSTFoilTime;

    book.timingSelectedEventTruthOriginTime->Fill(truthOriginTime);
    ++book.timingSelectedEventTruthOriginTimeEntries;

    book.timingSelectedTrackFirstSTFoilTime->Fill(firstTrackFirstSTFoilTime);
    book.timingSelectedTrackFirstSTFoilTime->Fill(secondTrackFirstSTFoilTime);
    book.timingSelectedTrackFirstSTFoilTimeEntries += 2;

    book.timingSelectedTrackFirstSTFoilDeltaT->Fill(firstSTFoilDeltaT);
    ++book.timingSelectedTrackFirstSTFoilDeltaTEntries;

    book.timingTruthMinusEarlierFirstSTFoilTime->Fill(
      truthMinusEarlierFirstSTFoilTime);
    book.timingTruthMinusLaterFirstSTFoilTime->Fill(
      truthMinusLaterFirstSTFoilTime);
    book.timingTruthMinusAverageFirstSTFoilTime->Fill(
      truthMinusAverageFirstSTFoilTime);
    book.timingTruthMinusEarlierVsLaterFirstSTFoilTime->Fill(
      truthMinusEarlierFirstSTFoilTime,
      truthMinusLaterFirstSTFoilTime);
    ++book.timingTruthMinusFirstSTFoilEntries;
  }

  inline double momentumPolarTheta(const twoparticlevertexer::Line3D& line)
  {
    if (!line.valid)
    {
      return -1.0;
    }

    const double transverseMomentumDirection =
      std::sqrt(line.unitDirection.x() * line.unitDirection.x() +
                line.unitDirection.y() * line.unitDirection.y());
    return kRadiansToDegrees *
           std::atan2(transverseMomentumDirection, line.unitDirection.z());
  }

  inline bool fillMomentumThetaPair(
    TH2F* histogram,
    const twoparticlevertexer::VertexResult& vertex)
  {
    if (histogram == nullptr || !vertex.valid)
    {
      return false;
    }

    const double theta1 = momentumPolarTheta(vertex.firstLine);
    const double theta2 = momentumPolarTheta(vertex.secondLine);
    if (!std::isfinite(theta1) || !std::isfinite(theta2) ||
        theta1 < 0.0 || theta2 < 0.0)
    {
      return false;
    }

    histogram->Fill(theta1, theta2);
    return true;
  }

  inline double momentumOpeningAngle(
    const twoparticlevertexer::VertexResult& vertex)
  {
    if (!vertex.valid || !vertex.firstLine.valid || !vertex.secondLine.valid)
    {
      return -1.0;
    }

    const double directionDotProduct =
      vertex.firstLine.unitDirection.Dot(vertex.secondLine.unitDirection);
    const double clampedDotProduct =
      std::max(-1.0, std::min(1.0, directionDotProduct));
    return kRadiansToDegrees * std::acos(clampedDotProduct);
  }

  inline bool fillMomentumOpeningAngle(
    TH1F* histogram,
    const twoparticlevertexer::VertexResult& vertex)
  {
    if (histogram == nullptr || !vertex.valid)
    {
      return false;
    }

    const double openingAngle = momentumOpeningAngle(vertex);
    if (!std::isfinite(openingAngle) || openingAngle < 0.0)
    {
      return false;
    }

    histogram->Fill(openingAngle);
    return true;
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
    if (fillMomentumThetaPair(book.recoVertexSpaceSelectedMomentumTheta1Theta2,
                              vertex))
    {
      ++book.recoVertexMomentumThetaEntries;
    }
    if (fillMomentumOpeningAngle(
          book.recoVertexSpaceSelectedMomentumOpeningAngle,
          vertex))
    {
      ++book.recoVertexMomentumOpeningAngleEntries;
    }

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
      book.recoVertexLineParameterST->Fill(lineParameterS, lineParameterT);
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

  inline void fillRecoTrackMultiplicityVsDeltaZ(HistogramBook& book,
                                                size_t recoTrackCount,
                                                double recoMinusTruthZ)
  {
    if (book.recoTrackMultiplicityVsDeltaZ == nullptr ||
        !std::isfinite(recoMinusTruthZ))
    {
      return;
    }

    book.recoTrackMultiplicityVsDeltaZ->Fill(
      static_cast<double>(recoTrackCount),
      recoMinusTruthZ);
    ++book.recoTrackMultiplicityVsDeltaZEntries;
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

  // Both inputs are from the same time-selected vertex: the shared ST_Foils
  // pair with minimum |deltaInputTime|.  This deliberately does not use the
  // space-selected pair.
  inline void fillRecoVertexMinTimeDeltaTVsTruthDeltaZTest(
    HistogramBook& book,
    const twoparticlevertexer::VertexResult& vertex,
    double recoMinusTruthZ)
  {
    if (!vertex.valid || !std::isfinite(vertex.deltaInputTime) ||
        !std::isfinite(recoMinusTruthZ))
    {
      return;
    }

    book.testRecoVertexMinTimeDeltaTVsTruthDeltaZ->Fill(
      vertex.deltaInputTime,
      recoMinusTruthZ);
    ++book.testRecoVertexMinTimeDeltaTVsTruthDeltaZEntries;
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

  inline bool validAverageAbsLineParameterVsDeltaZInputs(
    double averageAbsLineParameter,
    double recoMinusTruthZ)
  {
    return std::isfinite(averageAbsLineParameter) &&
           std::isfinite(recoMinusTruthZ) &&
           averageAbsLineParameter >= 0.0;
  }

  inline void fillRecoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZ(
    HistogramBook& book,
    double averageAbsLineParameter,
    double recoMinusTruthZ)
  {
    if (!validAverageAbsLineParameterVsDeltaZInputs(averageAbsLineParameter,
                                                    recoMinusTruthZ))
    {
      return;
    }

    book.recoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZ->SetPoint(
      book.recoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZ->GetN(),
      averageAbsLineParameter,
      recoMinusTruthZ);
    ++book.recoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZEntries;
  }

  inline void fillRecoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZ(
    HistogramBook& book,
    double averageAbsLineParameter,
    double recoMinusTruthZ)
  {
    if (!validAverageAbsLineParameterVsDeltaZInputs(averageAbsLineParameter,
                                                    recoMinusTruthZ))
    {
      return;
    }

    book.recoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZ->SetPoint(
      book.recoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZ->GetN(),
      averageAbsLineParameter,
      recoMinusTruthZ);
    ++book.recoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZEntries;
  }

  inline void fillRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZTest(
    HistogramBook& book,
    double averageAbsLineParameter,
    double recoMinusTruthZ)
  {
    if (!validAverageAbsLineParameterVsDeltaZInputs(averageAbsLineParameter,
                                                    recoMinusTruthZ))
    {
      return;
    }

    book.testRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZ->SetPoint(
      book.testRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZ->GetN(),
      averageAbsLineParameter,
      recoMinusTruthZ);
    ++book.testRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZEntries;
  }

  inline bool validDeltaZRelationInputs(double xValue, double recoMinusTruthZ)
  {
    return std::isfinite(xValue) && std::isfinite(recoMinusTruthZ);
  }

  inline void fillDeltaZRelationGraph(TGraph* graph,
                                      Long64_t& entries,
                                      double xValue,
                                      double recoMinusTruthZ)
  {
    if (graph == nullptr ||
        !validDeltaZRelationInputs(xValue, recoMinusTruthZ))
    {
      return;
    }

    graph->SetPoint(graph->GetN(), xValue, recoMinusTruthZ);
    ++entries;
  }

  inline void fillRecoSpaceSelectedDeltaZForSelectedFoil(
    HistogramBook& book,
    int sharedFoilIndex,
    double recoMinusTruthZ)
  {
    if (sharedFoilIndex < 0 ||
        sharedFoilIndex >=
          static_cast<int>(book.recoSpaceSelectedDeltaZBySelectedFoil.size()) ||
        !std::isfinite(recoMinusTruthZ))
    {
      return;
    }

    TH1F* histogram =
      book.recoSpaceSelectedDeltaZBySelectedFoil.at(sharedFoilIndex);
    if (histogram == nullptr)
    {
      return;
    }

    histogram->Fill(recoMinusTruthZ);
    ++book.recoSpaceSelectedDeltaZBySelectedFoilEntries.at(sharedFoilIndex);
  }

  inline void fillRecoSpaceSelectedSharedFoilNumberVsDeltaZ(
    HistogramBook& book,
    int sharedFoilIndex,
    double recoMinusTruthZ)
  {
    if (sharedFoilIndex < 0)
    {
      return;
    }

    fillDeltaZRelationGraph(book.recoSpaceSelectedSharedFoilNumberVsDeltaZ,
                            book.recoSpaceSelectedSharedFoilNumberVsDeltaZEntries,
                            static_cast<double>(sharedFoilIndex),
                            recoMinusTruthZ);
    fillRecoSpaceSelectedDeltaZForSelectedFoil(book,
                                               sharedFoilIndex,
                                               recoMinusTruthZ);
  }

  inline void fillRecoSelectedTrackFoilIntersectionZForFoil(
    HistogramBook& book,
    int foilIndex,
    double intersectionZ)
  {
    if (foilIndex < 0 ||
        foilIndex >=
          static_cast<int>(
            book.recoSelectedTrackFoilIntersectionZByFoil.size()) ||
        !std::isfinite(intersectionZ))
    {
      return;
    }

    TH1F* histogram =
      book.recoSelectedTrackFoilIntersectionZByFoil.at(foilIndex);
    if (histogram == nullptr)
    {
      return;
    }

    histogram->Fill(intersectionZ);
    ++book.recoSelectedTrackFoilIntersectionZByFoilEntries.at(foilIndex);
    ++book.recoSelectedTrackFoilIntersectionZEntries;
  }

  inline void fillRecoSelectedDownstreamElectronTrackCountForFoil(
    HistogramBook& book,
    int foilIndex)
  {
    if (book.recoSelectedDownstreamElectronTrackCountByFoil == nullptr ||
        foilIndex < 0 ||
        foilIndex >= book.config.stoppingTargetFoils)
    {
      return;
    }

    book.recoSelectedDownstreamElectronTrackCountByFoil->Fill(
      static_cast<double>(foilIndex));
    ++book.recoSelectedDownstreamElectronTrackCountByFoilEntries;
  }

  inline void fillRecoSelectedTrackFoilIntersectionMomentumForFoil(
    HistogramBook& book,
    int foilIndex,
    const XYZVectorF& momentum)
  {
    if (foilIndex < 0 ||
        foilIndex >=
          static_cast<int>(
            book.recoSelectedTrackFoilIntersectionPxByFoil.size()) ||
        foilIndex >=
          static_cast<int>(
            book.recoSelectedTrackFoilIntersectionPyByFoil.size()) ||
        foilIndex >=
          static_cast<int>(
            book.recoSelectedTrackFoilIntersectionPzByFoil.size()) ||
        foilIndex >=
          static_cast<int>(
            book.recoSelectedTrackFoilIntersectionPByFoil.size()))
    {
      return;
    }

    const double px = momentum.x();
    const double py = momentum.y();
    const double pz = momentum.z();
    const double p = momentum.R();
    if (!std::isfinite(px) || !std::isfinite(py) ||
        !std::isfinite(pz) || !std::isfinite(p))
    {
      return;
    }

    TH1F* pxHistogram =
      book.recoSelectedTrackFoilIntersectionPxByFoil.at(foilIndex);
    TH1F* pyHistogram =
      book.recoSelectedTrackFoilIntersectionPyByFoil.at(foilIndex);
    TH1F* pzHistogram =
      book.recoSelectedTrackFoilIntersectionPzByFoil.at(foilIndex);
    TH1F* pHistogram =
      book.recoSelectedTrackFoilIntersectionPByFoil.at(foilIndex);
    if (pxHistogram == nullptr || pyHistogram == nullptr ||
        pzHistogram == nullptr || pHistogram == nullptr)
    {
      return;
    }

    pxHistogram->Fill(px);
    pyHistogram->Fill(py);
    pzHistogram->Fill(pz);
    pHistogram->Fill(p);
    ++book.recoSelectedTrackFoilIntersectionMomentumByFoilEntries.at(foilIndex);
    ++book.recoSelectedTrackFoilIntersectionMomentumEntries;
  }

  inline void fillRecoSelectedTrackFoilIntersectionZByFoil(
    HistogramBook& book,
    const std::vector<std::vector<mu2e::TrkSegInfo>>* trackSegments,
    const std::vector<size_t>& selectedTrackIndices)
  {
    if (trackSegments == nullptr)
    {
      return;
    }

    for (const size_t trackIndex : selectedTrackIndices)
    {
      if (trackIndex >= trackSegments->size())
      {
        continue;
      }

      std::set<int> foilIndicesForTrack;
      for (const auto& segment : trackSegments->at(trackIndex))
      {
        if (segment.sid != mu2e::SurfaceIdDetail::ST_Foils)
        {
          continue;
        }

        fillRecoSelectedTrackFoilIntersectionZForFoil(book,
                                                      segment.sindex,
                                                      segment.pos.z());
        fillRecoSelectedTrackFoilIntersectionMomentumForFoil(book,
                                                             segment.sindex,
                                                             segment.mom);
        if (segment.sindex >= 0 &&
            segment.sindex < book.config.stoppingTargetFoils)
        {
          foilIndicesForTrack.insert(segment.sindex);
        }
      }

      for (const int foilIndex : foilIndicesForTrack)
      {
        fillRecoSelectedDownstreamElectronTrackCountForFoil(book,
                                                            foilIndex);
      }
    }
  }

  inline void fillRecoSpaceSelectedSharedFoilCountVsDeltaZ(
    HistogramBook& book,
    int sharedFoilCount,
    double recoMinusTruthZ)
  {
    if (sharedFoilCount < 0)
    {
      return;
    }

    fillDeltaZRelationGraph(book.recoSpaceSelectedSharedFoilCountVsDeltaZ,
                            book.recoSpaceSelectedSharedFoilCountVsDeltaZEntries,
                            static_cast<double>(sharedFoilCount),
                            recoMinusTruthZ);
  }

  inline void fillRecoSpaceSelectedMaxFoilsHitVsDeltaZ(
    HistogramBook& book,
    int maxFoilsHit,
    double recoMinusTruthZ)
  {
    if (maxFoilsHit < 0)
    {
      return;
    }

    fillDeltaZRelationGraph(book.recoSpaceSelectedMaxFoilsHitVsDeltaZ,
                            book.recoSpaceSelectedMaxFoilsHitVsDeltaZEntries,
                            static_cast<double>(maxFoilsHit),
                            recoMinusTruthZ);
  }

  inline void fillRecoSpaceSelectedAbsDeltaFoilsHitVsDeltaZ(
    HistogramBook& book,
    int absDeltaFoilsHit,
    double recoMinusTruthZ)
  {
    if (absDeltaFoilsHit < 0)
    {
      return;
    }

    fillDeltaZRelationGraph(book.recoSpaceSelectedAbsDeltaFoilsHitVsDeltaZ,
                            book.recoSpaceSelectedAbsDeltaFoilsHitVsDeltaZEntries,
                            static_cast<double>(absDeltaFoilsHit),
                            recoMinusTruthZ);
  }

  inline void fillRecoSpaceSelectedOpeningAngleVsDeltaZ(
    HistogramBook& book,
    double openingAngleDegrees,
    double recoMinusTruthZ)
  {
    if (openingAngleDegrees < 0.0)
    {
      return;
    }

    fillDeltaZRelationGraph(book.recoSpaceSelectedOpeningAngleVsDeltaZ,
                            book.recoSpaceSelectedOpeningAngleVsDeltaZEntries,
                            openingAngleDegrees,
                            recoMinusTruthZ);
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
    if (fillMomentumThetaPair(book.testRecoVertexMinTimeMomentumTheta1Theta2,
                              vertex))
    {
      ++book.testRecoVertexMinTimeMomentumThetaEntries;
    }
    if (fillMomentumOpeningAngle(
          book.testRecoVertexMinTimeMomentumOpeningAngle,
          vertex))
    {
      ++book.testRecoVertexMinTimeMomentumOpeningAngleEntries;
    }
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
    writeOne(book.timingSelectedEventTruthOriginTime);
    writeOne(book.timingSelectedTrackFirstSTFoilTime);
    writeOne(book.timingSelectedTrackFirstSTFoilDeltaT);
    writeOne(book.timingTruthMinusEarlierFirstSTFoilTime);
    writeOne(book.timingTruthMinusLaterFirstSTFoilTime);
    writeOne(book.timingTruthMinusAverageFirstSTFoilTime);
    writeOne(book.timingTruthMinusEarlierVsLaterFirstSTFoilTime);
    writeOne(book.recoVertexX);
    writeOne(book.recoVertexY);
    writeOne(book.recoVertexZ);
    writeOne(book.recoVertexXY);
    writeOne(book.recoVertexXZ);
    writeOne(book.recoVertexYZ);
    writeOne(book.recoVertexSpaceSelectedMomentumTheta1Theta2);
    writeOne(book.recoVertexSpaceSelectedMomentumOpeningAngle);
    writeOne(book.recoVertexFoilIndexMatchedXY);
    writeOne(book.recoVertexFoilIndexMatchedXZ);
    writeOne(book.recoVertexFoilIndexMatchedYZ);
    writeOne(book.recoVertexLineDistance);
    writeOne(book.recoVertexLineParameterS);
    writeOne(book.recoVertexLineParameterT);
    writeOne(book.recoVertexAbsLineParameterS);
    writeOne(book.recoVertexAbsLineParameterT);
    writeOne(book.recoVertexLineParameterST);
    writeOne(book.recoVertexAbsLineParameterST);
    writeOne(book.recoTrackMultiplicityVsDeltaZ);
    writeOne(book.recoVertexSelectedSegmentDeltaTTest);
    writeOne(book.recoVertexMinTimeDifferenceTest);
    writeOne(book.testRecoVertexMinTimeDeltaTVsTruthDeltaZ);
    writeOne(book.recoAllSharedFoilCandidateMaxLvsDeltaZ);
    writeOne(book.recoSpaceSelectedSharedFoilMaxLvsDeltaZ);
    writeOne(book.recoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZ);
    writeOne(book.recoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZ);
    writeOne(book.recoSpaceSelectedSharedFoilNumberVsDeltaZ);
    for (TH1F* histogram : book.recoSpaceSelectedDeltaZBySelectedFoil)
    {
      writeOne(histogram);
    }
    writeOne(book.recoSelectedDownstreamElectronTrackCountByFoil);
    for (TH1F* histogram : book.recoSelectedTrackFoilIntersectionZByFoil)
    {
      writeOne(histogram);
    }
    for (TH1F* histogram : book.recoSelectedTrackFoilIntersectionPxByFoil)
    {
      writeOne(histogram);
    }
    for (TH1F* histogram : book.recoSelectedTrackFoilIntersectionPyByFoil)
    {
      writeOne(histogram);
    }
    for (TH1F* histogram : book.recoSelectedTrackFoilIntersectionPzByFoil)
    {
      writeOne(histogram);
    }
    for (TH1F* histogram : book.recoSelectedTrackFoilIntersectionPByFoil)
    {
      writeOne(histogram);
    }
    writeOne(book.recoSpaceSelectedSharedFoilCountVsDeltaZ);
    writeOne(book.recoSpaceSelectedMaxFoilsHitVsDeltaZ);
    writeOne(book.recoSpaceSelectedAbsDeltaFoilsHitVsDeltaZ);
    writeOne(book.recoSpaceSelectedOpeningAngleVsDeltaZ);
    writeOne(book.testRecoVertexMinTimeX);
    writeOne(book.testRecoVertexMinTimeY);
    writeOne(book.testRecoVertexMinTimeZ);
    writeOne(book.testRecoVertexMinTimeXY);
    writeOne(book.testRecoVertexMinTimeXZ);
    writeOne(book.testRecoVertexMinTimeYZ);
    writeOne(book.testRecoVertexMinTimeMomentumTheta1Theta2);
    writeOne(book.testRecoVertexMinTimeMomentumOpeningAngle);
    writeOne(book.testRecoVertexMinTimeFoilIndexMatchedXY);
    writeOne(book.testRecoVertexMinTimeFoilIndexMatchedXZ);
    writeOne(book.testRecoVertexMinTimeFoilIndexMatchedYZ);
    writeOne(book.testRecoMinTimeSharedFoilMaxLvsDeltaZ);
    writeOne(book.testRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZ);
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

