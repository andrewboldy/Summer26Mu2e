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
    double momentumMin = 0.0;
    double momentumMax = 140.0;

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
    double vertexDistanceMax = 500.0;
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
    TH1F* recoVertexLineDistance = nullptr;
    TH1F* recoVertexTruthDeltaZ = nullptr;
    TH1F* recoVertexTruthAbsDeltaZ = nullptr;

    Long64_t selectedEvents = 0;
    Long64_t mcTruthEntries = 0;
    Long64_t recoMomentumEntries = 0;
    Long64_t recoVertexEntries = 0;
    Long64_t recoVertexTruthResidualEntries = 0;
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
    book.recoVertexLineDistance =
      make1D("hRecoTwoElectronVertexLineDistance",
             "Selected two-track closest-line distance;line-line distance [mm];entries",
             config.vertexDistanceBins,
             config.vertexDistanceMin,
             config.vertexDistanceMax);
    book.recoVertexTruthDeltaZ =
      make1D("hRecoTruthVertexDeltaZ",
             "Selected two-track vertex minus truth origin z difference;#Delta z [mm];entries",
             config.transverseBins,
             -1000.0,
             1000.0);
    book.recoVertexTruthAbsDeltaZ =
      make1D("hRecoTruthVertexAbsDeltaZ",
             "Absolute selected two-track vertex minus truth origin z difference;|#Delta z| [mm];entries",
             config.transverseBins,
             0.0,
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
    ++book.recoVertexEntries;
  }

  inline void fillRecoVertexTruthResidual(HistogramBook& book, double deltaZ)
  {
    if (!std::isfinite(deltaZ))
    {
      return;
    }

    book.recoVertexTruthDeltaZ->Fill(deltaZ);
    book.recoVertexTruthAbsDeltaZ->Fill(std::fabs(deltaZ));
    ++book.recoVertexTruthResidualEntries;
  }

  inline void writeOne(TH1* histogram)
  {
    if (histogram != nullptr)
    {
      histogram->Write();
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
    writeOne(book.recoVertexLineDistance);
    writeOne(book.recoVertexTruthDeltaZ);
    writeOne(book.recoVertexTruthAbsDeltaZ);
    outputFile->Close();
    delete outputFile;
    return true;
  }
}

#endif

