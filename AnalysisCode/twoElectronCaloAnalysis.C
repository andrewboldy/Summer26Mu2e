//----------------------------------------------------------------------------------
//
// twoElectronCaloAnalysis.C
// Written by Andrew Boldy University of South Carolina, 2026
// Assisted by Codex
//
// Purpose:
//   1. Loop over an EventNtuple ROOT file or filelist.
//   2. Find events with exactly two reconstructed e-minus tracks whose
//      trkmcsim MC match is a valid rank-0 electron.
//   3. For those selected events, print calorimeter information:
//        - track-associated calorimeter cluster information from trkcalohit
//        - event-level calorimeter clusters
//        - crystal-hit energies belonging to each cluster
//   4. Make calorimeter-energy histograms specifically for the two selected
//      reconstructed electron tracks.
//   5. Count whether the two selected track-calo matches are on the same disk
//      or on different disks.
//   6. Save CaloHitter overlays for the first selected Front/Front, Back/Back,
//      and Front/Back track-calo events.
//   7. Also draw and save a blank two-disk calorimeter map through CaloHitter.
//
// Coordinate note:
//   EventNtuple stores caloclusters.cog_ in the calorimeter disk front-face
//   frame.  Individual calohits store crystalId_ and eDep_, but do not store
//   individual crystal xyz coordinates.  The xyz printed for crystal-hit lines
//   below is therefore the parent cluster COG, not a per-crystal center.
//
//----------------------------------------------------------------------------------

// Standard C++ includes used for file I/O, formatted printing, string assembly,
// and temporary containers.
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// ROOT plotting and output utilities.  TH1F/TH2F hold the histograms, TCanvas
// saves them as PDFs, TFile writes a ROOT histogram file, and TSystem creates
// the plot directory if needed.
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TSystem.h>

// RooUtil opens EventNtuple files and exposes the event content through simple
// C++ objects.  common_cuts.hh gives helper functions like is_e_minus(track).
#include "EventNtuple/rooutil/inc/RooUtil.hh"
#include "EventNtuple/rooutil/inc/common_cuts.hh"

// Side helper for drawing the blank calorimeter disk geometry.  This macro does
// not fill that geometry yet, but the include keeps the calorimeter drawing data
// structures available for the next analysis layer.
#include "CaloHitter.hh"

// Helper for matching the two selected trkcalohit objects back to event-level
// caloclusters and calculating the distance between their cluster COGs.  The
// implementation is included here so this ROOT macro can run without a separate
// build-system target for the helper.
#include "TwoElectronCaloClusterDistance.C"

using namespace std;
using namespace rooutil;

void twoElectronCaloAnalysis(const string& generatorName,
                             const string& fileName,
                             const int maxSelectedEventsToPrint = -1,
                             // Denominator for generator-level efficiency
                             // percentages.  The current production sample was
                             // thrown with 100,000 events, but this can be
                             // overridden from the ROOT call if needed.
                             const long long totalThrownEvents = 100000,
                             // Per-crystal text lines are very large.  Keep
                             // them off by default so the text output stays
                             // below common upload limits while the histograms
                             // still receive every crystal hit.
                             const bool printCrystalHitDetails = false)
{
  // First do a simple existence/open check.  RooUtil will handle the detailed
  // ROOT parsing later, but this gives a clear error for a bad path or filelist.
  ifstream file(fileName);
  if (!file.is_open())
  {
    cerr << "ERROR: could not open input file or filelist: " << fileName << endl;
    return;
  }
  file.close();

  TStopwatch timer;
  timer.Start();

  // Keep ROOT plotting non-interactive.  This is important when running on a VM
  // or in batch because canvases otherwise try to open GUI windows.
  const bool wasBatchMode = gROOT->IsBatch();
  const bool oldAddDirectoryStatus = TH1::AddDirectoryStatus();
  gROOT->SetBatch(true);
  TH1::AddDirectory(false);

  // RooUtil accepts either a ROOT file path or a filelist.  From here on, the
  // code treats every entry as one EventNtuple event.
  RooUtil util(fileName);
  const int numEvents = util.GetNEvents();
  cout << "There are " << numEvents << " entries in the filelist." << endl;

  // The text output is intentionally flat and searchable.  Later plotting or
  // filtering scripts can grep for EVENT, RANK0_ELECTRON, TRACK_CALO, etc.
  const string outputFileName = "twoElectronCaloAnalysis_" + generatorName + ".txt";
  ofstream outputFile(outputFileName);
  if (!outputFile.is_open())
  {
    cerr << "ERROR: could not create output text file: " << outputFileName << endl;
    TH1::AddDirectory(oldAddDirectoryStatus);
    gROOT->SetBatch(wasBatchMode);
    return;
  }

  outputFile << "# twoElectronCaloAnalysis output\n"
             << "# Selected event definition: exactly two valid rank-0 trkmcsim electrons among reconstructed e-minus tracks.\n"
             << "# Cluster coordinates are caloclusters.cog_ in the calorimeter disk front-face frame.\n"
             << "# Cluster COG distance matches each selected trkcalohit back to caloclusters by disk, energy, time, and size.\n"
             << "# Cross-disk COG distances are computed in the stored disk front-face coordinates, not global detector coordinates.\n"
             << "# Crystal-hit lines use parent cluster COG xyz because EventNtuple calohits do not store per-crystal xyz.\n"
             << "# Per-crystal CRYSTAL_HIT detail lines are disabled unless printCrystalHitDetails is true.\n"
             << "# Crystal-hit histograms are still filled even when text detail lines are disabled.\n"
             << "# Momentum plots use reconstructed track momentum from trkcalohit.mom at the track-calo association.\n"
             << "# This is not a calorimeter-only momentum measurement.\n"
             << "# ALL histograms use reconstructed e-minus tracks with no trkmcsim/MC-truth requirement.\n"
             << "# Disk labels: raw disk 0 = Front, raw disk 1 = Back.\n"
             << "# Units: energy in MeV, momentum in MeV/c, position in mm.\n";

  // Print every analysis line both to the terminal and to the output text file.
  // Keeping this in one lambda prevents the cout and file output from drifting.
  auto printLine = [&outputFile](const string& line) {
    cout << line << endl;
    outputFile << line << '\n';
  };

  // All calorimeter plots for this analysis live in the requested subdirectory.
  // The CaloHitter module also uses this directory for geometry-only disk plots.
  const string caloPlotsDirectory = "Plots/CaloHitPlots";
  if (gSystem != nullptr)
  {
    gSystem->mkdir(caloPlotsDirectory.c_str(), true);
  }

  // Common histogram bounds for the reconstructed calorimeter-energy plots.
  // The requested analysis range is 0-70 MeV.  With 140 bins, each bin is
  // 0.5 MeV wide, which keeps the spectra readable while preserving useful
  // shape information.
  const double caloEnergyMin = 0.0;
  const double caloEnergyMax = 70.0;
  const int caloEnergyBins = 140;

  // The summed-energy plot contains E0 + E1 for the two selected reconstructed
  // electrons, so the upper edge is twice the single-electron upper edge.
  const double twoElectronCaloEnergySumMax = 2.0 * caloEnergyMax;

  // Crystal-hit sums can represent all reconstructed hit energy in one cluster
  // or in an entire selected event.  Give them a wider range than individual
  // crystal hits because many crystals can contribute to one shower.
  const double clusterCrystalHitEnergySumMax = twoElectronCaloEnergySumMax;
  const double eventCrystalHitEnergySumMax = 4.0 * caloEnergyMax;

  // The matched-track-momentum plots intentionally use the same 0-70 visual
  // range as the calorimeter-energy plots.  The value is the reconstructed
  // track momentum stored in trkcalohit.mom at the track-calo association, not
  // a calorimeter-only momentum measurement.
  const double caloMomentumMin = 0.0;
  const double caloMomentumMax = 70.0;
  const int caloMomentumBins = 140;

  // The summed-momentum plot also needs space for both reconstructed electrons.
  const double twoElectronCaloMomentumSumMax = 2.0 * caloMomentumMax;

  // Track-associated calorimeter histograms.
  //
  // These use trkcalohit, which is the reconstructed track-to-calo match.  The
  // energy is the matched calorimeter energy.  The momentum and POCA position
  // are the reconstructed track state stored at the track-calo association.
  TH1F* hTrackCaloEnergy = new TH1F(
    "hTrackCaloEnergy",
    "Track-associated calorimeter energy;E_{calo} [MeV];Matched rank-0 electron tracks",
    caloEnergyBins, caloEnergyMin, caloEnergyMax);
  TH1F* hTrackCaloMomentum = new TH1F(
    "hTrackCaloMomentum",
    "Reconstructed track momentum at track-calo association;p_{reco track at calo assoc} [MeV/c];Matched rank-0 electron tracks",
    caloMomentumBins, caloMomentumMin, caloMomentumMax);
  TH2F* hTrackCaloPOCAXY = new TH2F(
    "hTrackCaloPOCAXY",
    "Track-calo POCA position;x [mm];y [mm]",
    200, -1000.0, 1000.0, 200, -1000.0, 1000.0);
  TH1F* hTrackCaloPOCAZ = new TH1F(
    "hTrackCaloPOCAZ",
    "Track-calo POCA z position;z [mm];Matched rank-0 electron tracks",
    240, -12000.0, 12000.0);

  // Reconstruction-only "ALL" track-calo histograms.
  //
  // These are intentionally filled before any trkmcsim/rank-0 requirement is
  // applied.  The only physics-object requirement is a reconstructed e-minus
  // track with a usable track-to-calorimeter association.  They provide the
  // no-MC-truth reference sample for comparison with the exact two-rank0-electron
  // selected sample.
  TH1F* hAllRecoTrackCaloEnergy = new TH1F(
    "hAllRecoTrackCaloEnergy",
    "ALL reco e^{-} tracks: track-associated calorimeter energy;E_{calo} [MeV];Reconstructed e^{-} tracks",
    caloEnergyBins, caloEnergyMin, caloEnergyMax);
  TH1F* hAllRecoTrackCaloMomentum = new TH1F(
    "hAllRecoTrackCaloMomentum",
    "ALL reco e^{-} tracks: reconstructed track momentum at track-calo association;p_{reco track at calo assoc} [MeV/c];Reconstructed e^{-} tracks",
    caloMomentumBins, caloMomentumMin, caloMomentumMax);
  TH2F* hAllRecoTrackCaloPOCAXY = new TH2F(
    "hAllRecoTrackCaloPOCAXY",
    "ALL reco e^{-} tracks: track-calo POCA position;x [mm];y [mm]",
    200, -1000.0, 1000.0, 200, -1000.0, 1000.0);
  TH1F* hAllRecoTrackCaloPOCAZ = new TH1F(
    "hAllRecoTrackCaloPOCAZ",
    "ALL reco e^{-} tracks: track-calo POCA z position;z [mm];Reconstructed e^{-} tracks",
    240, -12000.0, 12000.0);

  // Dedicated two-electron calorimeter-energy histograms.
  //
  // These are the electron-count-like plots: they are filled only from the two
  // selected reconstructed rank-0 electron tracks in the event, using the
  // calorimeter energy in the track-matched trkcalohit object.
  //
  // "Electron 0" and "electron 1" are the two selected tracks in event.tracks
  // order.  They are not sorted by energy, momentum, disk number, time, or
  // leading/subleading rank.
  TH1F* hTwoElectronTrackCaloEnergyAll = new TH1F(
    "hTwoElectronTrackCaloEnergyAll",
    "Two-track events: electron calorimeter energy;E_{calo} [MeV];Selected reconstructed electrons",
    caloEnergyBins, caloEnergyMin, caloEnergyMax);
  TH1F* hTwoElectronTrackCaloEnergyElectron0 = new TH1F(
    "hTwoElectronTrackCaloEnergyElectron0",
    "Two-track events: electron 0 calorimeter energy;E_{calo} [MeV];Events",
    caloEnergyBins, caloEnergyMin, caloEnergyMax);
  TH1F* hTwoElectronTrackCaloEnergyElectron1 = new TH1F(
    "hTwoElectronTrackCaloEnergyElectron1",
    "Two-track events: electron 1 calorimeter energy;E_{calo} [MeV];Events",
    caloEnergyBins, caloEnergyMin, caloEnergyMax);
  TH1F* hTwoElectronTrackCaloEnergySum = new TH1F(
    "hTwoElectronTrackCaloEnergySum",
    "Two-track events: summed electron calorimeter energy;E_{calo,0}+E_{calo,1} [MeV];Events",
    2 * caloEnergyBins, caloEnergyMin, twoElectronCaloEnergySumMax);
  TH2F* hTwoElectronTrackCaloEnergyPair = new TH2F(
    "hTwoElectronTrackCaloEnergyPair",
    "Two-track events: electron calorimeter energy pair;electron 0 E_{calo} [MeV];electron 1 E_{calo} [MeV]",
    caloEnergyBins, caloEnergyMin, caloEnergyMax, caloEnergyBins, caloEnergyMin, caloEnergyMax);

  // Dedicated two-electron reconstructed track-momentum-at-calo histograms.
  //
  // These mirror the energy plots above, but use the reconstructed track
  // momentum vector stored in trkcalohit at the track-calo association.  The
  // same event.tracks ordering defines electron 0 and electron 1 here.
  TH1F* hTwoElectronTrackCaloMomentumAll = new TH1F(
    "hTwoElectronTrackCaloMomentumAll",
    "Two-track events: reconstructed track momentum at track-calo association;p_{reco track at calo assoc} [MeV/c];Selected reconstructed electrons",
    caloMomentumBins, caloMomentumMin, caloMomentumMax);
  TH1F* hTwoElectronTrackCaloMomentumElectron0 = new TH1F(
    "hTwoElectronTrackCaloMomentumElectron0",
    "Two-track events: electron 0 reconstructed track momentum at track-calo association;p_{reco track at calo assoc} [MeV/c];Events",
    caloMomentumBins, caloMomentumMin, caloMomentumMax);
  TH1F* hTwoElectronTrackCaloMomentumElectron1 = new TH1F(
    "hTwoElectronTrackCaloMomentumElectron1",
    "Two-track events: electron 1 reconstructed track momentum at track-calo association;p_{reco track at calo assoc} [MeV/c];Events",
    caloMomentumBins, caloMomentumMin, caloMomentumMax);
  TH1F* hTwoElectronTrackCaloMomentumSum = new TH1F(
    "hTwoElectronTrackCaloMomentumSum",
    "Two-track events: summed reconstructed track momentum at track-calo association;p_{reco track,0}+p_{reco track,1} at calo assoc [MeV/c];Events",
    2 * caloMomentumBins, caloMomentumMin, twoElectronCaloMomentumSumMax);
  TH2F* hTwoElectronTrackCaloMomentumPair = new TH2F(
    "hTwoElectronTrackCaloMomentumPair",
    "Two-track events: reconstructed track momentum pair at track-calo association;electron 0 p_{reco track at calo assoc} [MeV/c];electron 1 p_{reco track at calo assoc} [MeV/c]",
    caloMomentumBins, caloMomentumMin, caloMomentumMax, caloMomentumBins, caloMomentumMin, caloMomentumMax);

  // Event-level reconstructed calorimeter cluster histograms.
  //
  // caloclusters stores reconstructed energy and timing.  We still print COG in
  // the text dump, but we no longer make COG plots in this macro.
  //
  // These cluster histograms are broader diagnostics: they include all
  // reconstructed calorimeter clusters found in selected two-electron events,
  // not only clusters that are matched to the two selected tracks.
  TH1F* hClusterEnergy = new TH1F(
    "hClusterEnergy",
    "Reconstructed calorimeter cluster energy;E_{cluster} [MeV];Clusters",
    caloEnergyBins, caloEnergyMin, caloEnergyMax);
  TH1F* hClusterTime = new TH1F(
    "hClusterTime",
    "Reconstructed calorimeter cluster time;t_{cluster} [ns];Clusters",
    240, 0.0, 2400.0);

  // Crystal-hit histograms.  These use event.calohits, which provides crystal ID
  // and deposited energy.  The energy-by-ID histogram is the bridge to the
  // future CaloHitter overlay: each x bin corresponds to one calorimeter crystal.
  //
  // A single electron shower can produce many crystal-hit entries.  Therefore
  // this section will usually have many more entries than the two-electron
  // track-calo histograms.
  TH1F* hCrystalHitEnergy = new TH1F(
    "hCrystalHitEnergy",
    "Calorimeter crystal-hit energy;E_{hit} [MeV];Crystal hits",
    caloEnergyBins, caloEnergyMin, caloEnergyMax);
  TH1F* hCrystalEnergyById = new TH1F(
    "hCrystalEnergyById",
    "Summed crystal-hit energy by crystal ID;crystal ID;#Sigma E_{hit} [MeV]",
    1348, -0.5, 1347.5);
  TH1F* hClusterCrystalHitEnergySum = new TH1F(
    "hClusterCrystalHitEnergySum",
    "Summed crystal-hit energy per reconstructed cluster;#Sigma E_{hit in cluster} [MeV];Clusters",
    2 * caloEnergyBins, caloEnergyMin, clusterCrystalHitEnergySumMax);
  TH1F* hSelectedEventCrystalHitEnergySum = new TH1F(
    "hSelectedEventCrystalHitEnergySum",
    "Summed crystal-hit energy per selected event;#Sigma E_{hit in selected event} [MeV];Selected events",
    4 * caloEnergyBins, caloEnergyMin, eventCrystalHitEnergySumMax);
  TH2F* hClusterEnergyVsCrystalHitEnergySum = new TH2F(
    "hClusterEnergyVsCrystalHitEnergySum",
    "Cluster energy vs summed crystal-hit energy;E_{cluster} [MeV];#Sigma E_{hit in cluster} [MeV]",
    2 * caloEnergyBins, caloEnergyMin, clusterCrystalHitEnergySumMax,
    2 * caloEnergyBins, caloEnergyMin, clusterCrystalHitEnergySumMax);

  // Two-electron matched-cluster COG distance histograms.  These are filled
  // only for selected two-rank0-electron events where both selected trkcalohit
  // objects can be matched back to event-level caloclusters.
  TH1F* hTwoElectronClusterCogDistance3D =
    twoelectroncalocog::makeClusterCogDistance3DHistogram();
  TH1F* hTwoElectronClusterCogDistanceXY =
    twoelectroncalocog::makeClusterCogDistanceXYHistogram();

  // Small local record for the two electrons we are trying to identify.  It
  // keeps the reconstructed track index, its rank-0 MC truth match, and the
  // optional track-to-calorimeter association in one object.
  struct Rank0ElectronTrack {
    size_t trackIndex = 0;
    const mu2e::SimInfo* sim = nullptr;
    const mu2e::TrkCaloHitInfo* trkcalohit = nullptr;
  };

  // Per-electron display coordinates for the selected track-calo association.
  // The CaloHitter overlay expects disk-local x/y positions in mm.  Here we use
  // the x/y values carried by the track-calo POCA object because that is the
  // available reconstructed track-to-calorimeter impact information.
  struct SelectedElectronCaloPoint {
    bool valid = false;
    int disk = -1;
    size_t trackIndex = 0;
    double x = 0.0;
    double y = 0.0;
    double energy = -1.0;
  };

  // One saved example event display.  The analysis code fills only generic
  // calohitter::SelectedHit objects so CaloHitter remains a reusable drawing
  // helper instead of hard-coding this analysis' front/back categories.
  struct FirstDiskTopologyDisplay {
    bool found = false;
    int entry = -1;
    int run = -1;
    int subrun = -1;
    int event = -1;
    string eventLabel;
    vector<calohitter::SelectedHit> selectedHits;
  };

  FirstDiskTopologyDisplay firstFrontFrontDisplay;
  FirstDiskTopologyDisplay firstBackBackDisplay;
  FirstDiskTopologyDisplay firstFrontBackDisplay;

  // Counters are separated into "selected" and "printed" because the optional
  // maxSelectedEventsToPrint limit can stop printing before the scan is done.
  // Histogram counters are kept separately from print counters so the final
  // summary can answer two different questions:
  //   - how many physics objects were selected from the whole ntuple?
  //   - how many text lines were actually written under the print limit?
  long long selectedEventCount = 0;
  long long printedSelectedEventCount = 0;
  long long printedRank0ElectronCount = 0;
  long long printedTrackCaloCount = 0;
  long long printedClusterCount = 0;
  long long printedCrystalHitCount = 0;
  long long printedCrystalHitSummaryCount = 0;
  long long suppressedCrystalHitDetailCount = 0;
  long long allRecoEMinusTrackCount = 0;
  long long allRecoEMinusTracksWithoutCalo = 0;
  long long allRecoEMinusTrackCaloFills = 0;
  long long selectedEventsWithoutCaloClusters = 0;
  long long selectedRank0TracksWithoutCalo = 0;
  long long twoElectronTrackCaloEnergyFills = 0;
  long long selectedEventsWithBothElectronCaloEnergies = 0;
  long long twoElectronTrackCaloMomentumFills = 0;
  long long selectedEventsWithBothElectronCaloMomenta = 0;
  long long selectedEventsWithBothElectronCaloDisks = 0;
  long long selectedEventsWithSameElectronCaloDisk = 0;
  long long selectedEventsWithDifferentElectronCaloDisks = 0;
  long long selectedEventsWithBothElectronsOnFrontDisk = 0;
  long long selectedEventsWithBothElectronsOnBackDisk = 0;
  long long selectedEventsWithFrontBackElectronCaloDisks = 0;
  long long selectedEventsWithClusterCogDistance = 0;
  long long selectedEventsWithoutClusterCogDistance = 0;

  // Main EventNtuple loop.  Each event is examined independently, then rejected
  // unless it has exactly two reconstructed e- tracks with valid rank-0 truth.
  for (int i_event = 0; i_event < numEvents; ++i_event)
  {
    auto& event = util.GetEvent(i_event);

    // This vector will contain only reconstructed e- tracks whose trkmcsim list
    // contains a valid rank-0 electron.  The event passes only if this ends at 2.
    vector<Rank0ElectronTrack> rank0ElectronTracks;

    // Loop over reconstructed tracks.  We do the reconstructed-particle cut
    // first so non-electron tracks never enter the MC-truth matching logic.
    for (size_t i_track = 0; i_track < event.tracks.size(); ++i_track)
    {
      auto& track = event.tracks.at(i_track);
      if (!is_e_minus(track))
      {
        continue;
      }

      ++allRecoEMinusTrackCount;

      // Fill the reconstruction-only ALL reference sample before any MC-truth
      // requirement.  This keeps the comparison sample independent of trkmcsim.
      const auto* allRecoTrkCaloHit = track.trkcalohit;
      if (allRecoTrkCaloHit == nullptr ||
          allRecoTrkCaloHit->did < 0 ||
          allRecoTrkCaloHit->edep < 0.0)
      {
        ++allRecoEMinusTracksWithoutCalo;
      }
      else
      {
        hAllRecoTrackCaloEnergy->Fill(allRecoTrkCaloHit->edep);
        hAllRecoTrackCaloMomentum->Fill(allRecoTrkCaloHit->mom.R());
        hAllRecoTrackCaloPOCAXY->Fill(
          allRecoTrkCaloHit->poca.x(), allRecoTrkCaloHit->poca.y());
        hAllRecoTrackCaloPOCAZ->Fill(allRecoTrkCaloHit->poca.z());
        ++allRecoEMinusTrackCaloFills;
      }

      // Some tracks have no truth-match vector.  Those cannot be rank-0
      // electrons by this analysis definition.
      if (track.trkmcsim == nullptr)
      {
        continue;
      }

      // trkmcsim can contain multiple possible truth matches.  rank == 0 is
      // EventNtuple's best match to this reconstructed track.
      for (const auto& mctrack : *(track.trkmcsim))
      {
        if (!(mctrack.valid && mctrack.pdg == 11 && mctrack.rank == 0))
        {
          continue;
        }

        Rank0ElectronTrack electronTrack;
        electronTrack.trackIndex = i_track;
        electronTrack.sim = &mctrack;
        electronTrack.trkcalohit = track.trkcalohit;
        rank0ElectronTracks.push_back(electronTrack);
      }
    }

    // The requested sample is exactly two rank-0 electrons.  Events with one,
    // three, or more such tracks are deliberately skipped.
    if (rank0ElectronTracks.size() != 2)
    {
      continue;
    }

    ++selectedEventCount;

    // A negative maxSelectedEventsToPrint means "print all selected events".
    // Nonnegative values limit only the text dump; histograms still use every
    // selected event so plotted distributions are not biased by the print limit.
    const bool printThisEvent =
      (maxSelectedEventsToPrint < 0 || printedSelectedEventCount < maxSelectedEventsToPrint);
    auto printSelectedLine = [&printLine, printThisEvent](const string& line) {
      if (printThisEvent)
      {
        printLine(line);
      }
    };

    if (printThisEvent)
    {
      ++printedSelectedEventCount;
      printedRank0ElectronCount += rank0ElectronTracks.size();
    }

    int run = -1;
    int subrun = -1;
    int eventNumber = -1;
    // evtinfo should normally be present, but keep sentinel values if the ntuple
    // was made without that branch.
    if (event.evtinfo != nullptr)
    {
      run = event.evtinfo->run;
      subrun = event.evtinfo->subrun;
      eventNumber = event.evtinfo->event;
    }

    auto eventDisplayLabel = [&](const string& topologyName) {
      ostringstream label;
      label << topologyName
            << " selected two-electron track-calo event"
            << " | entry=" << i_event
            << " run=" << run
            << " subrun=" << subrun
            << " event=" << eventNumber;
      return label.str();
    };

    {
      ostringstream line;
      line << "\nEVENT entry=" << i_event
           << " run=" << run
           << " subrun=" << subrun
           << " event=" << eventNumber
           << " rank0_electrons=" << rank0ElectronTracks.size();
      printSelectedLine(line.str());
    }

    // Print one block per selected rank-0 electron.  This is the track-level
    // view: MC truth momentum plus the optional reconstructed track-calo match.
    //
    // The arrays below hold the two per-electron values until the end of the
    // event.  The boolean flags carry the real validity information; the -1
    // placeholders are just nonphysical sentinels for unfilled values.
    double twoElectronCaloEnergies[2] = {-1.0, -1.0};
    bool hasTwoElectronCaloEnergy[2] = {false, false};
    double twoElectronCaloMomenta[2] = {-1.0, -1.0};
    bool hasTwoElectronCaloMomentum[2] = {false, false};
    SelectedElectronCaloPoint selectedElectronCaloPoints[2];
    bool currentEventIsFrontFront = false;
    bool currentEventIsBackBack = false;
    bool currentEventIsFrontBack = false;
    vector<calohitter::SelectedHit> currentEventCrystalHitHighlights;
    vector<int> currentEventHighlightedCrystalIds;
    double selectedEventCrystalHitEnergySum = 0.0;

    for (size_t i_electron = 0; i_electron < rank0ElectronTracks.size(); ++i_electron)
    {
      const auto& electronTrack = rank0ElectronTracks.at(i_electron);
      const auto* sim = electronTrack.sim;
      const auto* trkcalohit = electronTrack.trkcalohit;

      {
        ostringstream line;
        line << "  RANK0_ELECTRON electron_index=" << i_electron
             << " trk_index=" << electronTrack.trackIndex
             << " sim_id=" << (sim != nullptr ? sim->id : -1)
             << " sim_mom=" << fixed << setprecision(6)
             << (sim != nullptr ? sim->mom.R() : -1.0)
             << " sim_p=("
             << (sim != nullptr ? sim->mom.x() : 0.0) << ", "
             << (sim != nullptr ? sim->mom.y() : 0.0) << ", "
             << (sim != nullptr ? sim->mom.z() : 0.0) << ")";
        printSelectedLine(line.str());
      }

      // trkcalohit is the track-associated calorimeter result.  It is not the
      // full event-level cluster collection; it is the calorimeter information
      // matched back to this reconstructed track.
      //
      // Negative disk IDs or negative energies are treated as unusable matches.
      // The track still helped the event pass the two-electron selection, but
      // it cannot contribute to the electron calo-energy/momentum histograms.
      if (trkcalohit == nullptr || trkcalohit->did < 0 || trkcalohit->edep < 0.0)
      {
        ++selectedRank0TracksWithoutCalo;
        ostringstream line;
        line << "    TRACK_CALO no associated calorimeter cluster for trk_index="
             << electronTrack.trackIndex;
        printSelectedLine(line.str());
        continue;
      }

      // did is the disk ID.  edep is the matched calorimeter energy.  poca is
      // the point of closest approach information carried by the track-calo
      // association, not a per-crystal center position.  Filling happens here,
      // after the validity checks, so every entry in the following track-calo
      // plots corresponds to one usable selected electron.
      hTrackCaloEnergy->Fill(trkcalohit->edep);
      const double trackCaloMomentum = trkcalohit->mom.R();
      hTrackCaloMomentum->Fill(trackCaloMomentum);
      hTrackCaloPOCAXY->Fill(trkcalohit->poca.x(), trkcalohit->poca.y());
      hTrackCaloPOCAZ->Fill(trkcalohit->poca.z());

      // Fill the electron-specific calorimeter-energy plots.  These are not
      // crystal-hit counts; each fill corresponds to one of the two selected
      // reconstructed electron tracks.  The "all" histogram receives both
      // electrons; the split histograms keep the first and second selected
      // electron separate for comparison.
      hTwoElectronTrackCaloEnergyAll->Fill(trkcalohit->edep);
      if (i_electron == 0)
      {
        hTwoElectronTrackCaloEnergyElectron0->Fill(trkcalohit->edep);
      }
      else if (i_electron == 1)
      {
        hTwoElectronTrackCaloEnergyElectron1->Fill(trkcalohit->edep);
      }
      if (i_electron < 2)
      {
        twoElectronCaloEnergies[i_electron] = trkcalohit->edep;
        hasTwoElectronCaloEnergy[i_electron] = true;
      }
      ++twoElectronTrackCaloEnergyFills;

      // Momentum plots are filled with the reconstructed track momentum stored
      // in the track-calo association, again one fill per selected electron.
      // This mirrors the energy logic exactly so energy and momentum entry
      // counts can be compared directly in the final summary.
      hTwoElectronTrackCaloMomentumAll->Fill(trackCaloMomentum);
      if (i_electron == 0)
      {
        hTwoElectronTrackCaloMomentumElectron0->Fill(trackCaloMomentum);
      }
      else if (i_electron == 1)
      {
        hTwoElectronTrackCaloMomentumElectron1->Fill(trackCaloMomentum);
      }
      if (i_electron < 2)
      {
        twoElectronCaloMomenta[i_electron] = trackCaloMomentum;
        hasTwoElectronCaloMomentum[i_electron] = true;
      }
      ++twoElectronTrackCaloMomentumFills;

      // Keep the two selected track-calo hit positions for the front/back
      // topology counters and for the first-event CaloHitter overlays.  The
      // drawing helper currently knows the two Mu2e disk IDs, 0 and 1, so other
      // disk values are still printed above but are not used for these overlays.
      if (i_electron < 2 && (trkcalohit->did == 0 || trkcalohit->did == 1))
      {
        selectedElectronCaloPoints[i_electron].valid = true;
        selectedElectronCaloPoints[i_electron].disk = trkcalohit->did;
        selectedElectronCaloPoints[i_electron].trackIndex = electronTrack.trackIndex;
        selectedElectronCaloPoints[i_electron].x = trkcalohit->poca.x();
        selectedElectronCaloPoints[i_electron].y = trkcalohit->poca.y();
        selectedElectronCaloPoints[i_electron].energy = trkcalohit->edep;
      }

      if (printThisEvent)
      {
        ++printedTrackCaloCount;
      }
      ostringstream line;
      line << "    TRACK_CALO"
           << " disk=" << trkcalohit->did
           << " disk_label=" << calohitter::diskShortLabel(trkcalohit->did)
           << " energy=" << fixed << setprecision(6) << trkcalohit->edep
           << " energyErr=" << trkcalohit->edeperr
           << " active=" << trkcalohit->active
           << " poca_xyz=(" << trkcalohit->poca.x()
           << ", " << trkcalohit->poca.y()
           << ", " << trkcalohit->poca.z() << ")"
           << " reco_track_mom_at_calo_assoc=" << trkcalohit->mom.R()
           << " reco_track_mom_xyz_at_calo_assoc=(" << trkcalohit->mom.x()
           << ", " << trkcalohit->mom.y()
           << ", " << trkcalohit->mom.z() << ")"
           << " doca=" << trkcalohit->doca
           << " dt=" << trkcalohit->dt;
      printSelectedLine(line.str());
    }

    // Pair-level two-electron calorimeter observables are filled only when both
    // selected reconstructed electrons have valid track-associated calo energy.
    // If one track lacks a usable calorimeter match, the event still contributes
    // to the single-electron distributions above, but not to the pair sum or the
    // two-dimensional pair plot.
    if (hasTwoElectronCaloEnergy[0] && hasTwoElectronCaloEnergy[1])
    {
      hTwoElectronTrackCaloEnergySum->Fill(twoElectronCaloEnergies[0] + twoElectronCaloEnergies[1]);
      hTwoElectronTrackCaloEnergyPair->Fill(twoElectronCaloEnergies[0], twoElectronCaloEnergies[1]);
      ++selectedEventsWithBothElectronCaloEnergies;
    }
    // Apply the same "both electrons must be valid" rule to the momentum pair
    // plots, because the pair sum is undefined if either reconstructed track
    // momentum at the track-calo association is missing.
    if (hasTwoElectronCaloMomentum[0] && hasTwoElectronCaloMomentum[1])
    {
      hTwoElectronTrackCaloMomentumSum->Fill(twoElectronCaloMomenta[0] + twoElectronCaloMomenta[1]);
      hTwoElectronTrackCaloMomentumPair->Fill(twoElectronCaloMomenta[0], twoElectronCaloMomenta[1]);
      ++selectedEventsWithBothElectronCaloMomenta;
    }

    // Count the disk topology of the two selected electron track-calo matches.
    // This is deliberately based on trkcalohit.did, not on MC truth, because the
    // question is where the reconstructed track-calo associations landed.
    if (selectedElectronCaloPoints[0].valid && selectedElectronCaloPoints[1].valid)
    {
      ++selectedEventsWithBothElectronCaloDisks;

      const int disk0 = selectedElectronCaloPoints[0].disk;
      const int disk1 = selectedElectronCaloPoints[1].disk;

      string diskTopologyCategory;
      if (disk0 == disk1)
      {
        ++selectedEventsWithSameElectronCaloDisk;
        diskTopologyCategory = "same";

        if (disk0 == 0)
        {
          ++selectedEventsWithBothElectronsOnFrontDisk;
          currentEventIsFrontFront = true;
        }
        else if (disk0 == 1)
        {
          ++selectedEventsWithBothElectronsOnBackDisk;
          currentEventIsBackBack = true;
        }
      }
      else
      {
        ++selectedEventsWithDifferentElectronCaloDisks;
        diskTopologyCategory = "different";

        // With the current two-disk calorimeter, any valid different-disk pair
        // is one Front hit and one Back hit, regardless of electron order.
        if ((disk0 == 0 && disk1 == 1) || (disk0 == 1 && disk1 == 0))
        {
          ++selectedEventsWithFrontBackElectronCaloDisks;
          currentEventIsFrontBack = true;
        }
      }

      {
        ostringstream line;
        line << "  TRACK_CALO_DISK_TOPOLOGY"
             << " electron0_disk=" << disk0
             << " electron0_disk_label=" << calohitter::diskShortLabel(disk0)
             << " electron1_disk=" << disk1
             << " electron1_disk_label=" << calohitter::diskShortLabel(disk1)
             << " relation=" << diskTopologyCategory;
        printSelectedLine(line.str());
      }
    }

    // Match the two selected track-associated calorimeter objects back to the
    // event-level caloclusters collection and measure the distance between the
    // two cluster COGs.  The COGs are caloclusters.cog_ values in the disk
    // front-face coordinate frame.
    const auto cogDistanceResult =
      twoelectroncalocog::calculateTwoElectronClusterCogDistance(
        rank0ElectronTracks.at(0).trkcalohit,
        rank0ElectronTracks.at(1).trkcalohit,
        event.caloclusters);
    twoelectroncalocog::fillClusterCogDistanceHistograms(
      cogDistanceResult,
      hTwoElectronClusterCogDistance3D,
      hTwoElectronClusterCogDistanceXY);
    if (cogDistanceResult.valid)
    {
      ++selectedEventsWithClusterCogDistance;
    }
    else
    {
      ++selectedEventsWithoutClusterCogDistance;
    }
    printSelectedLine(
      twoelectroncalocog::formatClusterCogDistanceLine(cogDistanceResult));

    // The event-level caloclusters collection contains calorimeter clusters
    // independent of whether a particular track was matched to one.  A selected
    // two-electron event can still have no reconstructed calorimeter cluster.
    //
    // This block is intentionally broader than the electron-only plots above.
    // It is the diagnostic view of all reconstructed calorimeter activity in
    // the selected event, which is why its entries can outnumber the two tracks.
    if (event.caloclusters == nullptr)
    {
      ++selectedEventsWithoutCaloClusters;
      printSelectedLine("  CALO_CLUSTER branch missing or disabled for this event.");
    }
    else if (event.caloclusters->empty())
    {
      ++selectedEventsWithoutCaloClusters;
      printSelectedLine("  CALO_CLUSTER no reconstructed calorimeter clusters in this event.");
    }
    else
    {
      // Loop over all reconstructed calorimeter clusters in this selected event.
      for (size_t i_cluster = 0; i_cluster < event.caloclusters->size(); ++i_cluster)
      {
        const auto& cluster = event.caloclusters->at(i_cluster);
        hClusterEnergy->Fill(cluster.energyDep_);
        hClusterTime->Fill(cluster.time_);

        if (printThisEvent)
        {
          ++printedClusterCount;
        }

        {
          // Cluster COG is stored in the calorimeter disk front-face coordinate
          // frame in EventNtuple.  It is not a global Mu2e coordinate.
          ostringstream line;
          line << "  CALO_CLUSTER"
               << " cluster_index=" << i_cluster
               << " disk=" << cluster.diskID_
               << " disk_label=" << calohitter::diskShortLabel(cluster.diskID_)
               << " energy=" << fixed << setprecision(6) << cluster.energyDep_
               << " energyErr=" << cluster.energyDepErr_
               << " time=" << cluster.time_
               << " size=" << cluster.size_
               << " isSplit=" << cluster.isSplit_
               << " cog_xyz=(" << cluster.cog_.x()
               << ", " << cluster.cog_.y()
               << ", " << cluster.cog_.z() << ")";
          printSelectedLine(line.str());
        }

        // The cluster owns a list of indices into event.calohits.  Without the
        // calohits branch, the macro can print cluster information but cannot
        // descend to individual crystal-hit energies.
        if (event.calohits == nullptr)
        {
          printSelectedLine("    CRYSTAL_HIT calohits branch missing or disabled; cannot print per-crystal energies.");
          continue;
        }

        // Each hit index is resolved back into the event-level calohits vector.
        // These entries carry crystalId_ and eDep_, which are the key pieces for
        // connecting deposited energy to individual calorimeter crystals.
        //
        // The crystal-hit histograms are always filled.  The per-hit text dump
        // is optional because those lines dominate the output file size.  When
        // detailed lines are disabled, print one compact summary line per
        // cluster instead.
        int validClusterCrystalHitCount = 0;
        int invalidClusterCrystalHitCount = 0;
        double validClusterCrystalHitEnergySum = 0.0;

        // Guard every index before dereferencing.  This keeps the macro useful
        // even if a file has an unusual cluster-hit reference: the bad index is
        // printed and the rest of the event can still be processed.
        for (const int hitIndex : cluster.hits_)
        {
          if (hitIndex < 0 || static_cast<size_t>(hitIndex) >= event.calohits->size())
          {
            ++invalidClusterCrystalHitCount;
            ostringstream line;
            line << "    CRYSTAL_HIT invalid hit index " << hitIndex
                 << " for cluster_index=" << i_cluster;
            printSelectedLine(line.str());
            continue;
          }

          const auto& hit = event.calohits->at(hitIndex);
          ++validClusterCrystalHitCount;
          validClusterCrystalHitEnergySum += hit.eDep_;
          selectedEventCrystalHitEnergySum += hit.eDep_;
          hCrystalHitEnergy->Fill(hit.eDep_);
          if (hit.crystalId_ >= 0 && hit.crystalId_ < 1348)
          {
            // Fill by weight, not by count.  The bin content is summed deposited
            // energy for that crystal ID across all selected events.
            hCrystalEnergyById->Fill(hit.crystalId_, hit.eDep_);

            // Save one highlight per hit crystal for this selected event.  The
            // first-event displays later reuse this vector so the CaloHitter PDF
            // shows every reconstructed crystal assigned to the event's clusters.
            bool alreadyHighlighted = false;
            for (const int highlightedCrystalId : currentEventHighlightedCrystalIds)
            {
              if (highlightedCrystalId == hit.crystalId_)
              {
                alreadyHighlighted = true;
                break;
              }
            }
            if (!alreadyHighlighted)
            {
              currentEventHighlightedCrystalIds.push_back(hit.crystalId_);
              currentEventCrystalHitHighlights.push_back(
                calohitter::SelectedHit::highlightCrystal(hit.crystalId_));
            }
          }

          if (printThisEvent && printCrystalHitDetails)
          {
            ++printedCrystalHitCount;
          }
          else if (printThisEvent)
          {
            ++suppressedCrystalHitDetailCount;
          }

          if (printCrystalHitDetails)
          {
            // EventNtuple calohits do not include the crystal center xyz directly.
            // For now, the printed xyz is the parent cluster COG.  CaloHitter will
            // let us convert crystal_id to a drawn crystal location in the next pass.
            ostringstream line;
            line << "    CRYSTAL_HIT"
                 << " hit_index=" << hitIndex
                 << " crystal_id=" << hit.crystalId_
                 << " parent_cluster=" << hit.clusterIdx_
                 << " disk=" << cluster.diskID_
                 << " disk_label=" << calohitter::diskShortLabel(cluster.diskID_)
                 << " energy=" << fixed << setprecision(6) << hit.eDep_
                 << " energyErr=" << hit.eDepErr_
                 << " time=" << hit.time_
                 << " nSiPMs=" << hit.nSiPMs_
                 << " parent_cluster_cog_xyz=(" << cluster.cog_.x()
                 << ", " << cluster.cog_.y()
                 << ", " << cluster.cog_.z() << ")";
            printSelectedLine(line.str());
          }
        }

        if (!printCrystalHitDetails)
        {
          if (printThisEvent)
          {
            ++printedCrystalHitSummaryCount;
          }

          ostringstream line;
          line << "    CRYSTAL_HIT_SUMMARY"
               << " cluster_index=" << i_cluster
               << " valid_hits=" << validClusterCrystalHitCount
               << " invalid_hits=" << invalidClusterCrystalHitCount
               << " summed_energy=" << fixed << setprecision(6)
               << validClusterCrystalHitEnergySum
               << " detail_lines_suppressed=" << validClusterCrystalHitCount;
          printSelectedLine(line.str());
        }

        hClusterCrystalHitEnergySum->Fill(validClusterCrystalHitEnergySum);
        hClusterEnergyVsCrystalHitEnergySum->Fill(
          cluster.energyDep_, validClusterCrystalHitEnergySum);
      }
    }

    hSelectedEventCrystalHitEnergySum->Fill(selectedEventCrystalHitEnergySum);

    // After the cluster loop, the event-level crystal highlight vector is
    // complete.  Save the first example of each requested disk topology with
    // both the selected track-calo markers and all reconstructed hit crystals.
    auto makeTrackSelectedHit = [](const SelectedElectronCaloPoint& point,
                                   const int electronIndex) {
      ostringstream label;
      label << "e" << electronIndex
            << " trk " << point.trackIndex
            << " E=" << fixed << setprecision(1) << point.energy;

      const Color_t markerColor = electronIndex == 0 ? kRed + 1 : kBlue + 1;
      const Style_t markerStyle = electronIndex == 0 ? 20 : 21;
      return calohitter::SelectedHit::fromXY(
        point.disk, point.x, point.y, label.str(), markerColor, markerStyle, 1.45);
    };

    auto saveFirstTopologyDisplay = [&](FirstDiskTopologyDisplay& display,
                                        const string& topologyName) {
      if (display.found)
      {
        return;
      }

      display.found = true;
      display.entry = i_event;
      display.run = run;
      display.subrun = subrun;
      display.event = eventNumber;
      display.eventLabel = eventDisplayLabel(topologyName);
      display.selectedHits.clear();
      display.selectedHits.insert(
        display.selectedHits.end(),
        currentEventCrystalHitHighlights.begin(),
        currentEventCrystalHitHighlights.end());
      display.selectedHits.push_back(makeTrackSelectedHit(selectedElectronCaloPoints[0], 0));
      display.selectedHits.push_back(makeTrackSelectedHit(selectedElectronCaloPoints[1], 1));
    };

    if (currentEventIsFrontFront)
    {
      saveFirstTopologyDisplay(firstFrontFrontDisplay, "Front/Front");
    }
    if (currentEventIsBackBack)
    {
      saveFirstTopologyDisplay(firstBackBackDisplay, "Back/Back");
    }
    if (currentEventIsFrontBack)
    {
      saveFirstTopologyDisplay(firstFrontBackDisplay, "Front/Back");
    }
  }

  // Write a machine-readable summary at the end of the text file.  The same
  // information is also printed below for quick terminal checks.
  //
  // The percentage table uses several denominators because they answer
  // different questions:
  //   - selected: fraction of the exact two-rank0-electron selected sample
  //   - valid_disk: fraction of selected events where both electrons have a
  //     usable Front/Back track-calo disk ID
  //   - ntuple: fraction of all events stored in this EventNtuple
  //   - thrown: fraction of the original generated/thrown sample
  auto percentOf = [](const long long count, const long long denominator) {
    if (denominator <= 0)
    {
      return -1.0;
    }
    return 100.0 * static_cast<double>(count) / static_cast<double>(denominator);
  };

  auto percentText = [&percentOf](const long long count, const long long denominator) {
    ostringstream text;
    if (denominator <= 0)
    {
      text << "n/a";
    }
    else
    {
      text << fixed << setprecision(6) << percentOf(count, denominator) << "%";
    }
    return text.str();
  };

  auto writeTopologyPercentLine = [&](const string& label, const long long count) {
    outputFile << "# disk_topology_percent " << label
               << " count " << count
               << " pct_of_selected_two_rank0 " << percentOf(count, selectedEventCount)
               << " pct_of_valid_disk_topology " << percentOf(count, selectedEventsWithBothElectronCaloDisks)
               << " pct_of_ntuple_entries " << percentOf(count, numEvents)
               << " pct_of_thrown_events " << percentOf(count, totalThrownEvents)
               << '\n';
  };

  outputFile << "\n# Summary\n"
             << "# selected_events " << selectedEventCount << '\n'
             << "# ntuple_entries " << numEvents << '\n'
             << "# total_thrown_events " << totalThrownEvents << '\n'
             << "# ntuple_entries_pct_of_thrown_events " << percentOf(numEvents, totalThrownEvents) << '\n'
             << "# printed_selected_events " << printedSelectedEventCount << '\n'
             << "# printed_rank0_electrons " << printedRank0ElectronCount << '\n'
             << "# printed_track_calo_entries " << printedTrackCaloCount << '\n'
             << "# selected_rank0_tracks_without_track_calo " << selectedRank0TracksWithoutCalo << '\n'
             << "# printed_clusters " << printedClusterCount << '\n'
             << "# printed_crystal_hits " << printedCrystalHitCount << '\n'
             << "# printed_crystal_hit_summaries " << printedCrystalHitSummaryCount << '\n'
             << "# suppressed_crystal_hit_detail_lines " << suppressedCrystalHitDetailCount << '\n'
             << "# print_crystal_hit_details " << (printCrystalHitDetails ? 1 : 0) << '\n'
             << "# all_reco_e_minus_tracks_no_mc_truth_requirement " << allRecoEMinusTrackCount << '\n'
             << "# all_reco_e_minus_tracks_without_track_calo " << allRecoEMinusTracksWithoutCalo << '\n'
             << "# all_reco_e_minus_track_calo_fills " << allRecoEMinusTrackCaloFills << '\n'
             << "# selected_events_without_calo_clusters " << selectedEventsWithoutCaloClusters << '\n'
             << "# two_electron_track_calo_energy_fills " << twoElectronTrackCaloEnergyFills << '\n'
             << "# selected_events_with_both_electron_calo_energies " << selectedEventsWithBothElectronCaloEnergies << '\n'
             << "# two_electron_track_calo_momentum_fills " << twoElectronTrackCaloMomentumFills << '\n'
             << "# two_electron_reco_track_momentum_at_calo_assoc_fills " << twoElectronTrackCaloMomentumFills << '\n'
             << "# selected_events_with_both_electron_calo_momenta " << selectedEventsWithBothElectronCaloMomenta << '\n'
             << "# selected_events_with_both_reco_track_momenta_at_calo_assoc " << selectedEventsWithBothElectronCaloMomenta << '\n'
             << "# selected_events_with_both_electron_calo_disks " << selectedEventsWithBothElectronCaloDisks << '\n'
             << "# selected_events_with_same_electron_calo_disk " << selectedEventsWithSameElectronCaloDisk << '\n'
             << "# selected_events_with_different_electron_calo_disks " << selectedEventsWithDifferentElectronCaloDisks << '\n'
             << "# selected_events_with_both_electrons_on_front_disk " << selectedEventsWithBothElectronsOnFrontDisk << '\n'
             << "# selected_events_with_both_electrons_on_back_disk " << selectedEventsWithBothElectronsOnBackDisk << '\n'
             << "# selected_events_with_front_back_electron_calo_disks " << selectedEventsWithFrontBackElectronCaloDisks << '\n'
             << "# selected_events_with_cluster_cog_distance " << selectedEventsWithClusterCogDistance << '\n'
             << "# selected_events_without_cluster_cog_distance " << selectedEventsWithoutClusterCogDistance << '\n'
             << "# first_front_front_display_found " << (firstFrontFrontDisplay.found ? 1 : 0) << '\n'
             << "# first_front_front_display_entry " << firstFrontFrontDisplay.entry
             << " run " << firstFrontFrontDisplay.run
             << " subrun " << firstFrontFrontDisplay.subrun
             << " event " << firstFrontFrontDisplay.event << '\n'
             << "# first_back_back_display_found " << (firstBackBackDisplay.found ? 1 : 0) << '\n'
             << "# first_back_back_display_entry " << firstBackBackDisplay.entry
             << " run " << firstBackBackDisplay.run
             << " subrun " << firstBackBackDisplay.subrun
             << " event " << firstBackBackDisplay.event << '\n'
             << "# first_front_back_display_found " << (firstFrontBackDisplay.found ? 1 : 0) << '\n'
             << "# first_front_back_display_entry " << firstFrontBackDisplay.entry
             << " run " << firstFrontBackDisplay.run
             << " subrun " << firstFrontBackDisplay.subrun
             << " event " << firstFrontBackDisplay.event << '\n'
             << "# histogram_all_reco_track_calo_entries " << hAllRecoTrackCaloEnergy->GetEntries() << '\n'
             << "# histogram_all_reco_track_calo_momentum_entries " << hAllRecoTrackCaloMomentum->GetEntries() << '\n'
             << "# histogram_track_calo_entries " << hTrackCaloEnergy->GetEntries() << '\n'
             << "# histogram_two_electron_track_calo_entries " << hTwoElectronTrackCaloEnergyAll->GetEntries() << '\n'
             << "# histogram_two_electron_track_calo_momentum_entries " << hTwoElectronTrackCaloMomentumAll->GetEntries() << '\n'
             << "# histogram_two_electron_reco_track_momentum_at_calo_assoc_entries " << hTwoElectronTrackCaloMomentumAll->GetEntries() << '\n'
             << "# histogram_cluster_entries " << hClusterEnergy->GetEntries() << '\n'
             << "# histogram_crystal_hit_entries " << hCrystalHitEnergy->GetEntries() << '\n'
             << "# histogram_cluster_crystal_hit_energy_sum_entries " << hClusterCrystalHitEnergySum->GetEntries() << '\n'
             << "# histogram_selected_event_crystal_hit_energy_sum_entries " << hSelectedEventCrystalHitEnergySum->GetEntries() << '\n'
             << "# histogram_two_electron_cluster_cog_distance_entries " << hTwoElectronClusterCogDistance3D->GetEntries() << '\n';

  outputFile << "# selected_events_pct_of_ntuple_entries "
             << percentOf(selectedEventCount, numEvents) << '\n'
             << "# selected_events_pct_of_thrown_events "
             << percentOf(selectedEventCount, totalThrownEvents) << '\n';
  writeTopologyPercentLine(
    "both_valid_disk_ids", selectedEventsWithBothElectronCaloDisks);
  writeTopologyPercentLine(
    "same_disk", selectedEventsWithSameElectronCaloDisk);
  writeTopologyPercentLine(
    "different_disks", selectedEventsWithDifferentElectronCaloDisks);
  writeTopologyPercentLine(
    "both_front", selectedEventsWithBothElectronsOnFrontDisk);
  writeTopologyPercentLine(
    "both_back", selectedEventsWithBothElectronsOnBackDisk);
  writeTopologyPercentLine(
    "front_back", selectedEventsWithFrontBackElectronCaloDisks);

  outputFile.close();

  // Terminal summary for a quick sanity check after a long run.
  cout << "\nSummary:" << endl;
  cout << "  selected two-rank0-electron events found: " << selectedEventCount << endl;
  cout << "  selected two-rank0-electron events printed: " << printedSelectedEventCount << endl;
  cout << "  printed rank-0 electrons: " << printedRank0ElectronCount << endl;
  cout << "  track-associated calo entries printed: " << printedTrackCaloCount << endl;
  cout << "  rank-0 tracks without track calo association: " << selectedRank0TracksWithoutCalo << endl;
  cout << "  event-level calo clusters printed: " << printedClusterCount << endl;
  cout << "  detailed crystal-hit energy lines printed: " << printedCrystalHitCount << endl;
  cout << "  crystal-hit summary lines printed: " << printedCrystalHitSummaryCount << endl;
  cout << "  detailed crystal-hit lines suppressed: " << suppressedCrystalHitDetailCount << endl;
  cout << "  print detailed crystal-hit lines: " << (printCrystalHitDetails ? "yes" : "no") << endl;
  cout << "  ALL reco e-minus tracks without MC-truth requirement: "
       << allRecoEMinusTrackCount << endl;
  cout << "  ALL reco e-minus tracks without track calo association: "
       << allRecoEMinusTracksWithoutCalo << endl;
  cout << "  ALL reco e-minus track-calo histogram fills: "
       << allRecoEMinusTrackCaloFills << endl;
  cout << "  selected events without calo clusters: " << selectedEventsWithoutCaloClusters << endl;
  cout << "  two-electron track-calo energy histogram fills: " << twoElectronTrackCaloEnergyFills << endl;
  cout << "  selected events with both electron calo energies: " << selectedEventsWithBothElectronCaloEnergies << endl;
  cout << "  two-electron reconstructed track momentum-at-calo-association histogram fills: "
       << twoElectronTrackCaloMomentumFills << endl;
  cout << "  selected events with both reconstructed track momenta at calo association: "
       << selectedEventsWithBothElectronCaloMomenta << endl;
  cout << "  selected events with both electron calo disk IDs: " << selectedEventsWithBothElectronCaloDisks << endl;
  cout << "  selected events with both electrons on the same disk: " << selectedEventsWithSameElectronCaloDisk << endl;
  cout << "  selected events with electrons on different disks: " << selectedEventsWithDifferentElectronCaloDisks << endl;
  cout << "  selected events with both electrons on the Front disk: " << selectedEventsWithBothElectronsOnFrontDisk << endl;
  cout << "  selected events with both electrons on the Back disk: " << selectedEventsWithBothElectronsOnBackDisk << endl;
  cout << "  selected events with one Front electron and one Back electron: " << selectedEventsWithFrontBackElectronCaloDisks << endl;
  cout << "  selected events with matched cluster COG distance: "
       << selectedEventsWithClusterCogDistance << endl;
  cout << "  selected events without matched cluster COG distance: "
       << selectedEventsWithoutClusterCogDistance << endl;
  cout << "\nDisk topology percentages:" << endl;
  cout << "  denominator selected two-rank0-electron events: " << selectedEventCount << endl;
  cout << "  denominator valid two-electron calo disk IDs: " << selectedEventsWithBothElectronCaloDisks << endl;
  cout << "  denominator ntuple entries: " << numEvents << endl;
  cout << "  denominator thrown events: " << totalThrownEvents << endl;
  cout << "  ntuple entries as percent of thrown events: "
       << percentText(numEvents, totalThrownEvents) << endl;
  cout << "  selected two-rank0-electron events as percent of ntuple entries: "
       << percentText(selectedEventCount, numEvents) << endl;
  cout << "  selected two-rank0-electron events as percent of thrown events: "
       << percentText(selectedEventCount, totalThrownEvents) << endl;

  auto printTopologyPercentLine = [&](const string& label, const long long count) {
    cout << "  " << label << ": " << count
         << " | selected=" << percentText(count, selectedEventCount)
         << " | valid_disk=" << percentText(count, selectedEventsWithBothElectronCaloDisks)
         << " | ntuple=" << percentText(count, numEvents)
         << " | thrown=" << percentText(count, totalThrownEvents)
         << endl;
  };

  printTopologyPercentLine(
    "both electrons have valid Front/Back disk IDs",
    selectedEventsWithBothElectronCaloDisks);
  printTopologyPercentLine(
    "both electrons on the same disk",
    selectedEventsWithSameElectronCaloDisk);
  printTopologyPercentLine(
    "electrons on different disks",
    selectedEventsWithDifferentElectronCaloDisks);
  printTopologyPercentLine(
    "both electrons on the Front disk",
    selectedEventsWithBothElectronsOnFrontDisk);
  printTopologyPercentLine(
    "both electrons on the Back disk",
    selectedEventsWithBothElectronsOnBackDisk);
  printTopologyPercentLine(
    "one Front electron and one Back electron",
    selectedEventsWithFrontBackElectronCaloDisks);

  cout << "  first Front/Front display entry: " << firstFrontFrontDisplay.entry
       << " run=" << firstFrontFrontDisplay.run
       << " subrun=" << firstFrontFrontDisplay.subrun
       << " event=" << firstFrontFrontDisplay.event << endl;
  cout << "  first Back/Back display entry: " << firstBackBackDisplay.entry
       << " run=" << firstBackBackDisplay.run
       << " subrun=" << firstBackBackDisplay.subrun
       << " event=" << firstBackBackDisplay.event << endl;
  cout << "  first Front/Back display entry: " << firstFrontBackDisplay.entry
       << " run=" << firstFrontBackDisplay.run
       << " subrun=" << firstFrontBackDisplay.subrun
       << " event=" << firstFrontBackDisplay.event << endl;
  cout << "Wrote text output to " << outputFileName << endl;

  // Write the ROOT histograms to a file so the plotted distributions can be
  // reopened, rebinned, or overlaid later without rerunning the ntuple loop.
  const string histogramRootFileName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_RecoCaloHistograms.root";
  TFile histogramFile(histogramRootFileName.c_str(), "RECREATE");
  if (!histogramFile.IsZombie())
  {
    // Reconstruction-only ALL reference sample.  These histograms use
    // reconstructed e-minus tracks and do not require trkmcsim/MC truth.
    hAllRecoTrackCaloEnergy->Write();
    hAllRecoTrackCaloMomentum->Write();
    hAllRecoTrackCaloPOCAXY->Write();
    hAllRecoTrackCaloPOCAZ->Write();

    // Track-matched objects for the selected rank-0 electrons.
    hTrackCaloEnergy->Write();
    hTrackCaloMomentum->Write();
    hTrackCaloPOCAXY->Write();
    hTrackCaloPOCAZ->Write();

    // Electron-only calorimeter energy spectra for events that pass the exact
    // two-rank0-electron selection.
    hTwoElectronTrackCaloEnergyAll->Write();
    hTwoElectronTrackCaloEnergyElectron0->Write();
    hTwoElectronTrackCaloEnergyElectron1->Write();
    hTwoElectronTrackCaloEnergySum->Write();
    hTwoElectronTrackCaloEnergyPair->Write();

    // Electron-only reconstructed track momentum-at-calo-association spectra
    // with the same selection and pair logic as the energy histograms.
    hTwoElectronTrackCaloMomentumAll->Write();
    hTwoElectronTrackCaloMomentumElectron0->Write();
    hTwoElectronTrackCaloMomentumElectron1->Write();
    hTwoElectronTrackCaloMomentumSum->Write();
    hTwoElectronTrackCaloMomentumPair->Write();

    // Broader event-level calorimeter diagnostics for the selected events.
    hClusterEnergy->Write();
    hClusterTime->Write();
    hCrystalHitEnergy->Write();
    hCrystalEnergyById->Write();
    hClusterCrystalHitEnergySum->Write();
    hSelectedEventCrystalHitEnergySum->Write();
    hClusterEnergyVsCrystalHitEnergySum->Write();
    hTwoElectronClusterCogDistance3D->Write();
    hTwoElectronClusterCogDistanceXY->Write();
    histogramFile.Close();
    cout << "Wrote calorimeter histogram ROOT file to " << histogramRootFileName << endl;
  }
  else
  {
    cerr << "ERROR: could not create calorimeter histogram ROOT file: "
         << histogramRootFileName << endl;
  }

  auto styleLineHistogram = [](TH1* histogram,
                               const Color_t lineColor,
                               const Style_t lineStyle = 1) {
    if (histogram == nullptr)
    {
      return;
    }
    histogram->SetLineColor(lineColor);
    histogram->SetMarkerColor(lineColor);
    histogram->SetLineStyle(lineStyle);
    histogram->SetLineWidth(2);
    histogram->SetStats(false);
  };

  auto drawAllVsTwoElectronOverlay = [&](TH1* allHistogram,
                                         TH1* selectedHistogram,
                                         const string& title,
                                         const bool logY) {
    if (allHistogram == nullptr || selectedHistogram == nullptr)
    {
      return;
    }

    styleLineHistogram(allHistogram, kBlack, 1);
    styleLineHistogram(selectedHistogram, kRed + 1, 2);

    const double allMaximum = allHistogram->GetMaximum();
    const double selectedMaximum = selectedHistogram->GetMaximum();
    const double largerMaximum = allMaximum > selectedMaximum ? allMaximum : selectedMaximum;
    const double yMaximum = largerMaximum > 0.0 ? largerMaximum * (logY ? 20.0 : 1.25) : 1.0;

    allHistogram->SetTitle(title.c_str());
    allHistogram->SetMinimum(logY ? 0.5 : 0.0);
    allHistogram->SetMaximum(yMaximum);
    allHistogram->Draw("HIST");
    selectedHistogram->Draw("HIST SAME");

    TLegend* legend = new TLegend(0.48, 0.74, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(allHistogram, "ALL reco e^{-}, no MC-truth requirement", "l");
    legend->AddEntry(selectedHistogram, "Two-electron MC-truth-selected reco tracks", "l");
    legend->Draw();
  };

  // Save a compact set of PDF plots for the reconstructed calorimeter
  // observables now available in the ntuple.
  //
  // PDF group 1: basic track-calo quantities.  This is the first place to look
  // when checking whether the two selected reconstructed electrons actually
  // have usable calorimeter associations.
  const string trackCaloPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_TrackMatchedCalo.pdf";
  TCanvas* cTrackCalo = new TCanvas(
    "cTrackCalo",
    "Track-associated calorimeter quantities",
    1400, 1000);
  cTrackCalo->Divide(2, 2);
  cTrackCalo->cd(1);
  hTrackCaloEnergy->Draw("HIST");
  cTrackCalo->cd(2);
  hTrackCaloMomentum->Draw("HIST");
  cTrackCalo->cd(3);
  hTrackCaloPOCAXY->Draw("COLZ");
  cTrackCalo->cd(4);
  hTrackCaloPOCAZ->Draw("HIST");
  cTrackCalo->SaveAs(trackCaloPdfName.c_str());

  // PDF group 2: reconstruction-only ALL track-calo quantities.  These plots
  // use every reconstructed e-minus track with a usable track-calo association,
  // with no trkmcsim/MC-truth requirement.
  const string allRecoTrackCaloPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_ALLRecoTrackMatchedCalo.pdf";
  TCanvas* cAllRecoTrackCalo = new TCanvas(
    "cAllRecoTrackCalo",
    "ALL reconstructed e-minus track-calo quantities",
    1400, 1000);
  cAllRecoTrackCalo->Divide(2, 2);
  cAllRecoTrackCalo->cd(1);
  styleLineHistogram(hAllRecoTrackCaloEnergy, kBlack, 1);
  hAllRecoTrackCaloEnergy->Draw("HIST");
  cAllRecoTrackCalo->cd(2);
  styleLineHistogram(hAllRecoTrackCaloMomentum, kBlack, 1);
  hAllRecoTrackCaloMomentum->Draw("HIST");
  cAllRecoTrackCalo->cd(3);
  hAllRecoTrackCaloPOCAXY->SetStats(false);
  hAllRecoTrackCaloPOCAXY->Draw("COLZ");
  cAllRecoTrackCalo->cd(4);
  styleLineHistogram(hAllRecoTrackCaloPOCAZ, kBlack, 1);
  hAllRecoTrackCaloPOCAZ->Draw("HIST");
  cAllRecoTrackCalo->SaveAs(allRecoTrackCaloPdfName.c_str());

  // PDF group 3: raw-entry overlays of the reconstruction-only ALL sample and
  // the two-electron MC-truth-selected reconstructed tracks.  These quantities
  // are still reconstructed track/calo quantities; the MC truth is only used to
  // define the selected comparison sample.
  const string allVsTwoElectronTrackCaloOverlayPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_ALLVsTwoElectronMCTruthTrackCaloOverlay.pdf";
  TCanvas* cAllVsTwoElectronTrackCaloOverlay = new TCanvas(
    "cAllVsTwoElectronTrackCaloOverlay",
    "ALL reco tracks vs two-electron MC-truth-selected reco tracks",
    1500, 500);
  cAllVsTwoElectronTrackCaloOverlay->Divide(3, 1);
  cAllVsTwoElectronTrackCaloOverlay->cd(1);
  drawAllVsTwoElectronOverlay(
    hAllRecoTrackCaloEnergy, hTwoElectronTrackCaloEnergyAll,
    "Track-associated calorimeter energy: ALL vs two-electron MC truth;E_{calo} [MeV];Reconstructed e^{-} tracks",
    false);
  cAllVsTwoElectronTrackCaloOverlay->cd(2);
  drawAllVsTwoElectronOverlay(
    hAllRecoTrackCaloMomentum, hTwoElectronTrackCaloMomentumAll,
    "Reconstructed track momentum at track-calo association: ALL vs two-electron MC truth;p_{reco track at calo assoc} [MeV/c];Reconstructed e^{-} tracks",
    false);
  cAllVsTwoElectronTrackCaloOverlay->cd(3);
  drawAllVsTwoElectronOverlay(
    hAllRecoTrackCaloPOCAZ, hTrackCaloPOCAZ,
    "Track-calo POCA z: ALL vs two-electron MC truth;z [mm];Reconstructed e^{-} tracks",
    false);
  cAllVsTwoElectronTrackCaloOverlay->SaveAs(allVsTwoElectronTrackCaloOverlayPdfName.c_str());

  // PDF group 4: log-y version of the same overlays.  This keeps the
  // two-electron selected distributions visible when the ALL sample is much
  // larger in raw entry count.
  const string allVsTwoElectronTrackCaloOverlayLogPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_ALLVsTwoElectronMCTruthTrackCaloOverlayLogY.pdf";
  TCanvas* cAllVsTwoElectronTrackCaloOverlayLog = new TCanvas(
    "cAllVsTwoElectronTrackCaloOverlayLog",
    "ALL reco tracks vs two-electron MC-truth-selected reco tracks log y",
    1500, 500);
  cAllVsTwoElectronTrackCaloOverlayLog->Divide(3, 1);
  cAllVsTwoElectronTrackCaloOverlayLog->cd(1);
  cAllVsTwoElectronTrackCaloOverlayLog->GetPad(1)->SetLogy();
  drawAllVsTwoElectronOverlay(
    hAllRecoTrackCaloEnergy, hTwoElectronTrackCaloEnergyAll,
    "Track-associated calorimeter energy: ALL vs two-electron MC truth log y;E_{calo} [MeV];Reconstructed e^{-} tracks",
    true);
  cAllVsTwoElectronTrackCaloOverlayLog->cd(2);
  cAllVsTwoElectronTrackCaloOverlayLog->GetPad(2)->SetLogy();
  drawAllVsTwoElectronOverlay(
    hAllRecoTrackCaloMomentum, hTwoElectronTrackCaloMomentumAll,
    "Reconstructed track momentum at track-calo association: ALL vs two-electron MC truth log y;p_{reco track at calo assoc} [MeV/c];Reconstructed e^{-} tracks",
    true);
  cAllVsTwoElectronTrackCaloOverlayLog->cd(3);
  cAllVsTwoElectronTrackCaloOverlayLog->GetPad(3)->SetLogy();
  drawAllVsTwoElectronOverlay(
    hAllRecoTrackCaloPOCAZ, hTrackCaloPOCAZ,
    "Track-calo POCA z: ALL vs two-electron MC truth log y;z [mm];Reconstructed e^{-} tracks",
    true);
  cAllVsTwoElectronTrackCaloOverlayLog->SaveAs(allVsTwoElectronTrackCaloOverlayLogPdfName.c_str());

  // PDF group 5: linear-scale calorimeter-energy spectra for the two selected
  // electrons.  The first pad combines both electrons, the next two pads split
  // them by event.tracks order, and the final pad keeps the pair information.
  const string twoElectronTrackCaloPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_TwoElectronTrackCaloEnergy.pdf";
  TCanvas* cTwoElectronTrackCalo = new TCanvas(
    "cTwoElectronTrackCalo",
    "Two selected reconstructed electrons: calorimeter energy",
    1400, 1000);
  cTwoElectronTrackCalo->Divide(2, 2);
  cTwoElectronTrackCalo->cd(1);
  styleLineHistogram(hTwoElectronTrackCaloEnergyAll, kBlack, 1);
  hTwoElectronTrackCaloEnergyAll->Draw("HIST");
  cTwoElectronTrackCalo->cd(2);
  hTwoElectronTrackCaloEnergyElectron0->Draw("HIST");
  cTwoElectronTrackCalo->cd(3);
  hTwoElectronTrackCaloEnergyElectron1->Draw("HIST");
  cTwoElectronTrackCalo->cd(4);
  hTwoElectronTrackCaloEnergyPair->Draw("COLZ");
  cTwoElectronTrackCalo->SaveAs(twoElectronTrackCaloPdfName.c_str());

  // PDF group 6: the same combined electron calorimeter-energy spectrum on a
  // log y-axis.  The histogram content is unchanged; only the drawing scale is
  // different so tails and low-count bins are easier to see.
  const string twoElectronTrackCaloLogPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_TwoElectronTrackCaloEnergyLogY.pdf";
  TCanvas* cTwoElectronTrackCaloLog = new TCanvas(
    "cTwoElectronTrackCaloLog",
    "Two selected reconstructed electrons: calorimeter energy log y",
    900, 700);
  cTwoElectronTrackCaloLog->SetLogy();
  styleLineHistogram(hTwoElectronTrackCaloEnergyAll, kBlack, 1);
  hTwoElectronTrackCaloEnergyAll->SetMinimum(0.5);
  hTwoElectronTrackCaloEnergyAll->Draw("HIST");
  cTwoElectronTrackCaloLog->SaveAs(twoElectronTrackCaloLogPdfName.c_str());

  // PDF group 7: linear-scale reconstructed track momentum-at-calo-association
  // spectra.  These plots are the momentum counterparts of the energy plots
  // above, using trkcalohit.mom.R().
  const string twoElectronTrackCaloMomentumPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_TwoElectronRecoTrackMomentumAtCaloAssociation.pdf";
  TCanvas* cTwoElectronTrackCaloMomentum = new TCanvas(
    "cTwoElectronTrackCaloMomentum",
    "Two selected reconstructed electrons: reconstructed track momentum at track-calo association",
    1400, 1000);
  cTwoElectronTrackCaloMomentum->Divide(2, 2);
  cTwoElectronTrackCaloMomentum->cd(1);
  styleLineHistogram(hTwoElectronTrackCaloMomentumAll, kBlack, 1);
  hTwoElectronTrackCaloMomentumAll->Draw("HIST");
  cTwoElectronTrackCaloMomentum->cd(2);
  hTwoElectronTrackCaloMomentumElectron0->Draw("HIST");
  cTwoElectronTrackCaloMomentum->cd(3);
  hTwoElectronTrackCaloMomentumElectron1->Draw("HIST");
  cTwoElectronTrackCaloMomentum->cd(4);
  hTwoElectronTrackCaloMomentumPair->Draw("COLZ");
  cTwoElectronTrackCaloMomentum->SaveAs(twoElectronTrackCaloMomentumPdfName.c_str());

  // PDF group 8: log-scale version of the combined reconstructed track
  // momentum-at-calo-association spectrum.  This is useful for checking small
  // tails without changing the selected sample.
  const string twoElectronTrackCaloMomentumLogPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_TwoElectronRecoTrackMomentumAtCaloAssociationLogY.pdf";
  TCanvas* cTwoElectronTrackCaloMomentumLog = new TCanvas(
    "cTwoElectronTrackCaloMomentumLog",
    "Two selected reconstructed electrons: reconstructed track momentum at track-calo association log y",
    900, 700);
  cTwoElectronTrackCaloMomentumLog->SetLogy();
  styleLineHistogram(hTwoElectronTrackCaloMomentumAll, kBlack, 1);
  hTwoElectronTrackCaloMomentumAll->SetMinimum(0.5);
  hTwoElectronTrackCaloMomentumAll->Draw("HIST");
  cTwoElectronTrackCaloMomentumLog->SaveAs(twoElectronTrackCaloMomentumLogPdfName.c_str());

  // PDF group 9: event-level reconstructed clusters.  These plots are not
  // restricted to the exact two matched clusters; they show all clusters present
  // in events that passed the two-electron selection.
  const string clusterPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_RecoCaloClusters.pdf";
  TCanvas* cClusters = new TCanvas(
    "cClusters",
    "Reconstructed calorimeter clusters",
    1000, 500);
  cClusters->Divide(2, 1);
  cClusters->cd(1);
  hClusterEnergy->Draw("HIST");
  cClusters->cd(2);
  hClusterTime->Draw("HIST");
  cClusters->SaveAs(clusterPdfName.c_str());

  // PDF group 10: matched cluster COG distance for the two selected electrons.
  // The matching is from each selected trkcalohit back to caloclusters, then
  // the distance is calculated using caloclusters.cog_ in disk front-face
  // coordinates.
  const string clusterCogDistancePdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_TwoElectronClusterCogDistance.pdf";
  TCanvas* cClusterCogDistance = new TCanvas(
    "cClusterCogDistance",
    "Two selected reconstructed electrons: matched calo cluster COG distance",
    1200, 500);
  cClusterCogDistance->Divide(2, 1);
  cClusterCogDistance->cd(1);
  styleLineHistogram(hTwoElectronClusterCogDistance3D, kBlack, 1);
  hTwoElectronClusterCogDistance3D->Draw("HIST");
  cClusterCogDistance->cd(2);
  styleLineHistogram(hTwoElectronClusterCogDistanceXY, kBlack, 1);
  hTwoElectronClusterCogDistanceXY->Draw("HIST");
  cClusterCogDistance->SaveAs(clusterCogDistancePdfName.c_str());
  if (gSystem != nullptr && gSystem->AccessPathName(clusterCogDistancePdfName.c_str()))
  {
    cerr << "ERROR: cluster COG distance PDF was not found after SaveAs: "
         << clusterCogDistancePdfName << endl;
  }
  else
  {
    cout << "Verified cluster COG distance PDF: "
         << clusterCogDistancePdfName << endl;
  }

  // PDF group 11: crystal-hit diagnostics.  The left plot counts individual
  // crystal hits by deposited energy.  The right plot accumulates deposited
  // energy by global crystal ID, which is the information that can later be
  // painted onto the CaloHitter geometry.
  const string crystalPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_CrystalHitEnergy.pdf";
  TCanvas* cCrystalHits = new TCanvas(
    "cCrystalHits",
    "Calorimeter crystal-hit energy",
    1400, 600);
  cCrystalHits->Divide(2, 1);
  cCrystalHits->cd(1);
  hCrystalHitEnergy->Draw("HIST");
  cCrystalHits->cd(2);
  hCrystalEnergyById->Draw("HIST");
  cCrystalHits->SaveAs(crystalPdfName.c_str());

  // PDF group 12: summed crystal-hit energies.  These plots answer the question
  // "how much energy is carried by all reconstructed hit crystals together?"
  // rather than treating each crystal as one independent entry.
  const string crystalSumPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_CrystalHitEnergySums.pdf";
  TCanvas* cCrystalHitSums = new TCanvas(
    "cCrystalHitSums",
    "Summed calorimeter crystal-hit energy",
    1400, 1000);
  cCrystalHitSums->Divide(2, 2);
  cCrystalHitSums->cd(1);
  hClusterCrystalHitEnergySum->Draw("HIST");
  cCrystalHitSums->cd(2);
  hSelectedEventCrystalHitEnergySum->Draw("HIST");
  cCrystalHitSums->cd(3);
  hClusterEnergyVsCrystalHitEnergySum->Draw("COLZ");
  cCrystalHitSums->SaveAs(crystalSumPdfName.c_str());

  cout << "Wrote calorimeter PDF plots to:" << endl;
  cout << "  " << trackCaloPdfName << endl;
  cout << "  " << allRecoTrackCaloPdfName << endl;
  cout << "  " << allVsTwoElectronTrackCaloOverlayPdfName << endl;
  cout << "  " << allVsTwoElectronTrackCaloOverlayLogPdfName << endl;
  cout << "  " << twoElectronTrackCaloPdfName << endl;
  cout << "  " << twoElectronTrackCaloLogPdfName << endl;
  cout << "  " << twoElectronTrackCaloMomentumPdfName << endl;
  cout << "  " << twoElectronTrackCaloMomentumLogPdfName << endl;
  cout << "  " << clusterPdfName << endl;
  cout << "  " << clusterCogDistancePdfName << endl;
  cout << "  " << crystalPdfName << endl;
  cout << "  " << crystalSumPdfName << endl;

  // Draw the calorimeter after the event loop.  Right now CaloHitter only knows
  // the blank crystal geometry, but this location is intentional: later we can
  // use the accumulated event/crystal information from the analysis above to
  // alter crystal colors, labels, or hit markers before saving the diagram.
  const string blankCaloDiskPdfName =
    "Plots/CaloHitPlots/twoElectronCaloAnalysis_BlankCaloDisks_" + generatorName + ".pdf";
  cout << "Writing blank calorimeter disk PDF: " << blankCaloDiskPdfName << endl;
  const vector<calohitter::SelectedHit> noSelectedHits;
  calohitter::saveCalorimeterPdf(
    blankCaloDiskPdfName, noSelectedHits, "cTwoElectronCaloAnalysisBlankDisks");

  // Save the first selected examples for each requested disk topology.  The
  // overlay vectors were filled during the event loop, but the actual drawing is
  // done here after all histograms and standard PDFs have been written.
  auto saveFirstTopologyPdf = [](const FirstDiskTopologyDisplay& display,
                                 const string& outputPdf,
                                 const string& topologyName,
                                 const string& canvasName) {
    if (!display.found)
    {
      cout << "No " << topologyName
           << " selected two-electron track-calo event found; not writing "
           << outputPdf << endl;
      return;
    }

    cout << "Writing " << topologyName
         << " selected two-electron track-calo display: " << outputPdf
         << " from entry=" << display.entry
         << " run=" << display.run
         << " subrun=" << display.subrun
         << " event=" << display.event
         << " overlays=" << display.selectedHits.size() << endl;
    calohitter::saveCalorimeterPdf(
      outputPdf, display.selectedHits, canvasName, display.eventLabel);
  };

  const string firstFrontFrontCaloDiskPdfName =
    "Plots/CaloHitPlots/twoElectronCaloAnalysis_" + generatorName +
    "_FirstFrontFrontTrackCaloEvent.pdf";
  const string firstBackBackCaloDiskPdfName =
    "Plots/CaloHitPlots/twoElectronCaloAnalysis_" + generatorName +
    "_FirstBackBackTrackCaloEvent.pdf";
  const string firstFrontBackCaloDiskPdfName =
    "Plots/CaloHitPlots/twoElectronCaloAnalysis_" + generatorName +
    "_FirstFrontBackTrackCaloEvent.pdf";

  saveFirstTopologyPdf(
    firstFrontFrontDisplay, firstFrontFrontCaloDiskPdfName,
    "Front/Front", "cTwoElectronCaloAnalysisFirstFrontFront");
  saveFirstTopologyPdf(
    firstBackBackDisplay, firstBackBackCaloDiskPdfName,
    "Back/Back", "cTwoElectronCaloAnalysisFirstBackBack");
  saveFirstTopologyPdf(
    firstFrontBackDisplay, firstFrontBackCaloDiskPdfName,
    "Front/Back", "cTwoElectronCaloAnalysisFirstFrontBack");

  timer.Stop();
  // Timing is useful when deciding whether to print all selected events or only
  // the first few during development.
  cout << "CPU time: " << timer.CpuTime() << " s, real time: " << timer.RealTime() << " s" << endl;

  TH1::AddDirectory(oldAddDirectoryStatus);
  gROOT->SetBatch(wasBatchMode);
}

