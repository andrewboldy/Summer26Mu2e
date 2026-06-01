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

using namespace std;
using namespace rooutil;

void twoElectronCaloAnalysis(const string& generatorName,
                             const string& fileName,
                             const int maxSelectedEventsToPrint = -1,
                             // Denominator for generator-level efficiency
                             // percentages.  The current production sample was
                             // thrown with 100,000 events, but this can be
                             // overridden from the ROOT call if needed.
                             const long long totalThrownEvents = 100000)
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
             << "# Crystal-hit lines use parent cluster COG xyz because EventNtuple calohits do not store per-crystal xyz.\n"
             << "# Calorimeter momentum is represented by the matched reconstructed track momentum from trkcalohit.mom.\n"
             << "# Disk labels: raw disk 0 = Front, raw disk 1 = Back.\n"
             << "# Units: energy in MeV, position in mm.\n";

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

  // The matched-momentum plots intentionally use the same 0-70 visual range as
  // the calorimeter-energy plots.  The value is the track momentum stored in
  // trkcalohit.mom, not a calorimeter-only momentum measurement.
  const double caloMomentumMin = 0.0;
  const double caloMomentumMax = 70.0;
  const int caloMomentumBins = 140;

  // The summed-momentum plot also needs space for both reconstructed electrons.
  const double twoElectronCaloMomentumSumMax = 2.0 * caloMomentumMax;

  // Track-associated calorimeter histograms.
  //
  // These use trkcalohit, which is the reconstructed track-to-calo match.  The
  // energy is the matched calorimeter energy.  The momentum and POCA position
  // are the reconstructed track state at the calorimeter association.
  TH1F* hTrackCaloEnergy = new TH1F(
    "hTrackCaloEnergy",
    "Track-associated calorimeter energy;E_{calo} [MeV];Matched rank-0 electron tracks",
    caloEnergyBins, caloEnergyMin, caloEnergyMax);
  TH1F* hTrackCaloMomentum = new TH1F(
    "hTrackCaloMomentum",
    "Track momentum at calorimeter association;p_{track} [MeV/c];Matched rank-0 electron tracks",
    caloMomentumBins, caloMomentumMin, caloMomentumMax);
  TH2F* hTrackCaloPOCAXY = new TH2F(
    "hTrackCaloPOCAXY",
    "Track-calo POCA position;x [mm];y [mm]",
    200, -1000.0, 1000.0, 200, -1000.0, 1000.0);
  TH1F* hTrackCaloPOCAZ = new TH1F(
    "hTrackCaloPOCAZ",
    "Track-calo POCA z position;z [mm];Matched rank-0 electron tracks",
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

  // Dedicated two-electron track-momentum histograms.
  //
  // These mirror the energy plots above, but use the reconstructed track
  // momentum vector stored in trkcalohit at the calorimeter association.
  // The same event.tracks ordering defines electron 0 and electron 1 here.
  TH1F* hTwoElectronTrackCaloMomentumAll = new TH1F(
    "hTwoElectronTrackCaloMomentumAll",
    "Two-track events: track momentum at calo association;p_{track} [MeV/c];Selected reconstructed electrons",
    caloMomentumBins, caloMomentumMin, caloMomentumMax);
  TH1F* hTwoElectronTrackCaloMomentumElectron0 = new TH1F(
    "hTwoElectronTrackCaloMomentumElectron0",
    "Two-track events: electron 0 track momentum at calo association;p_{track} [MeV/c];Events",
    caloMomentumBins, caloMomentumMin, caloMomentumMax);
  TH1F* hTwoElectronTrackCaloMomentumElectron1 = new TH1F(
    "hTwoElectronTrackCaloMomentumElectron1",
    "Two-track events: electron 1 track momentum at calo association;p_{track} [MeV/c];Events",
    caloMomentumBins, caloMomentumMin, caloMomentumMax);
  TH1F* hTwoElectronTrackCaloMomentumSum = new TH1F(
    "hTwoElectronTrackCaloMomentumSum",
    "Two-track events: summed track momentum at calo association;p_{track,0}+p_{track,1} [MeV/c];Events",
    2 * caloMomentumBins, caloMomentumMin, twoElectronCaloMomentumSumMax);
  TH2F* hTwoElectronTrackCaloMomentumPair = new TH2F(
    "hTwoElectronTrackCaloMomentumPair",
    "Two-track events: track momentum pair at calo association;electron 0 p_{track} [MeV/c];electron 1 p_{track} [MeV/c]",
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

      // Momentum plots are filled with the matched reconstructed track momentum
      // at the calorimeter association, again one fill per selected electron.
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
           << " track_mom=" << trkcalohit->mom.R()
           << " track_mom_xyz=(" << trkcalohit->mom.x()
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
    // plots, because the pair sum is undefined if either matched momentum is
    // missing.
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

      // Convert this event's two selected electron impact points into the
      // generic CaloHitter overlay records.  CaloHitter only receives disk/x/y
      // marker data; all physics categorization stays in this analysis macro.
      auto makeSelectedHit = [](const SelectedElectronCaloPoint& point,
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

      auto saveFirstTopologyDisplay = [&](FirstDiskTopologyDisplay& display) {
        if (display.found)
        {
          return;
        }

        display.found = true;
        display.entry = i_event;
        display.run = run;
        display.subrun = subrun;
        display.event = eventNumber;
        display.selectedHits.clear();
        display.selectedHits.push_back(makeSelectedHit(selectedElectronCaloPoints[0], 0));
        display.selectedHits.push_back(makeSelectedHit(selectedElectronCaloPoints[1], 1));
      };

      string diskTopologyCategory;
      if (disk0 == disk1)
      {
        ++selectedEventsWithSameElectronCaloDisk;
        diskTopologyCategory = "same";

        if (disk0 == 0)
        {
          ++selectedEventsWithBothElectronsOnFrontDisk;
          saveFirstTopologyDisplay(firstFrontFrontDisplay);
        }
        else if (disk0 == 1)
        {
          ++selectedEventsWithBothElectronsOnBackDisk;
          saveFirstTopologyDisplay(firstBackBackDisplay);
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
          saveFirstTopologyDisplay(firstFrontBackDisplay);
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
        // Guard every index before dereferencing.  This keeps the macro useful
        // even if a file has an unusual cluster-hit reference: the bad index is
        // printed and the rest of the event can still be processed.
        for (const int hitIndex : cluster.hits_)
        {
          if (hitIndex < 0 || static_cast<size_t>(hitIndex) >= event.calohits->size())
          {
            ostringstream line;
            line << "    CRYSTAL_HIT invalid hit index " << hitIndex
                 << " for cluster_index=" << i_cluster;
            printSelectedLine(line.str());
            continue;
          }

          const auto& hit = event.calohits->at(hitIndex);
          hCrystalHitEnergy->Fill(hit.eDep_);
          if (hit.crystalId_ >= 0 && hit.crystalId_ < 1348)
          {
            // Fill by weight, not by count.  The bin content is summed deposited
            // energy for that crystal ID across all selected events.
            hCrystalEnergyById->Fill(hit.crystalId_, hit.eDep_);
          }

          if (printThisEvent)
          {
            ++printedCrystalHitCount;
          }

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
             << "# selected_events_without_calo_clusters " << selectedEventsWithoutCaloClusters << '\n'
             << "# two_electron_track_calo_energy_fills " << twoElectronTrackCaloEnergyFills << '\n'
             << "# selected_events_with_both_electron_calo_energies " << selectedEventsWithBothElectronCaloEnergies << '\n'
             << "# two_electron_track_calo_momentum_fills " << twoElectronTrackCaloMomentumFills << '\n'
             << "# selected_events_with_both_electron_calo_momenta " << selectedEventsWithBothElectronCaloMomenta << '\n'
             << "# selected_events_with_both_electron_calo_disks " << selectedEventsWithBothElectronCaloDisks << '\n'
             << "# selected_events_with_same_electron_calo_disk " << selectedEventsWithSameElectronCaloDisk << '\n'
             << "# selected_events_with_different_electron_calo_disks " << selectedEventsWithDifferentElectronCaloDisks << '\n'
             << "# selected_events_with_both_electrons_on_front_disk " << selectedEventsWithBothElectronsOnFrontDisk << '\n'
             << "# selected_events_with_both_electrons_on_back_disk " << selectedEventsWithBothElectronsOnBackDisk << '\n'
             << "# selected_events_with_front_back_electron_calo_disks " << selectedEventsWithFrontBackElectronCaloDisks << '\n'
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
             << "# histogram_track_calo_entries " << hTrackCaloEnergy->GetEntries() << '\n'
             << "# histogram_two_electron_track_calo_entries " << hTwoElectronTrackCaloEnergyAll->GetEntries() << '\n'
             << "# histogram_two_electron_track_calo_momentum_entries " << hTwoElectronTrackCaloMomentumAll->GetEntries() << '\n'
             << "# histogram_cluster_entries " << hClusterEnergy->GetEntries() << '\n'
             << "# histogram_crystal_hit_entries " << hCrystalHitEnergy->GetEntries() << '\n';

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
  cout << "  crystal-hit energy lines printed: " << printedCrystalHitCount << endl;
  cout << "  selected events without calo clusters: " << selectedEventsWithoutCaloClusters << endl;
  cout << "  two-electron track-calo energy histogram fills: " << twoElectronTrackCaloEnergyFills << endl;
  cout << "  selected events with both electron calo energies: " << selectedEventsWithBothElectronCaloEnergies << endl;
  cout << "  two-electron track-calo momentum histogram fills: " << twoElectronTrackCaloMomentumFills << endl;
  cout << "  selected events with both electron calo momenta: " << selectedEventsWithBothElectronCaloMomenta << endl;
  cout << "  selected events with both electron calo disk IDs: " << selectedEventsWithBothElectronCaloDisks << endl;
  cout << "  selected events with both electrons on the same disk: " << selectedEventsWithSameElectronCaloDisk << endl;
  cout << "  selected events with electrons on different disks: " << selectedEventsWithDifferentElectronCaloDisks << endl;
  cout << "  selected events with both electrons on the Front disk: " << selectedEventsWithBothElectronsOnFrontDisk << endl;
  cout << "  selected events with both electrons on the Back disk: " << selectedEventsWithBothElectronsOnBackDisk << endl;
  cout << "  selected events with one Front electron and one Back electron: " << selectedEventsWithFrontBackElectronCaloDisks << endl;
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

    // Electron-only matched-momentum spectra with the same selection and pair
    // logic as the energy histograms.
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
    histogramFile.Close();
    cout << "Wrote calorimeter histogram ROOT file to " << histogramRootFileName << endl;
  }
  else
  {
    cerr << "ERROR: could not create calorimeter histogram ROOT file: "
         << histogramRootFileName << endl;
  }

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

  // PDF group 2: linear-scale calorimeter-energy spectra for the two selected
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
  hTwoElectronTrackCaloEnergyAll->Draw("HIST");
  cTwoElectronTrackCalo->cd(2);
  hTwoElectronTrackCaloEnergyElectron0->Draw("HIST");
  cTwoElectronTrackCalo->cd(3);
  hTwoElectronTrackCaloEnergyElectron1->Draw("HIST");
  cTwoElectronTrackCalo->cd(4);
  hTwoElectronTrackCaloEnergyPair->Draw("COLZ");
  cTwoElectronTrackCalo->SaveAs(twoElectronTrackCaloPdfName.c_str());

  // PDF group 3: the same combined electron calorimeter-energy spectrum on a
  // log y-axis.  The histogram content is unchanged; only the drawing scale is
  // different so tails and low-count bins are easier to see.
  const string twoElectronTrackCaloLogPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_TwoElectronTrackCaloEnergyLogY.pdf";
  TCanvas* cTwoElectronTrackCaloLog = new TCanvas(
    "cTwoElectronTrackCaloLog",
    "Two selected reconstructed electrons: calorimeter energy log y",
    900, 700);
  cTwoElectronTrackCaloLog->SetLogy();
  hTwoElectronTrackCaloEnergyAll->SetMinimum(0.5);
  hTwoElectronTrackCaloEnergyAll->Draw("HIST");
  cTwoElectronTrackCaloLog->SaveAs(twoElectronTrackCaloLogPdfName.c_str());

  // PDF group 4: linear-scale matched-momentum spectra.  These plots are the
  // momentum counterparts of the energy plots above, using trkcalohit.mom.R().
  const string twoElectronTrackCaloMomentumPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_TwoElectronTrackCaloMomentum.pdf";
  TCanvas* cTwoElectronTrackCaloMomentum = new TCanvas(
    "cTwoElectronTrackCaloMomentum",
    "Two selected reconstructed electrons: track momentum at calo association",
    1400, 1000);
  cTwoElectronTrackCaloMomentum->Divide(2, 2);
  cTwoElectronTrackCaloMomentum->cd(1);
  hTwoElectronTrackCaloMomentumAll->Draw("HIST");
  cTwoElectronTrackCaloMomentum->cd(2);
  hTwoElectronTrackCaloMomentumElectron0->Draw("HIST");
  cTwoElectronTrackCaloMomentum->cd(3);
  hTwoElectronTrackCaloMomentumElectron1->Draw("HIST");
  cTwoElectronTrackCaloMomentum->cd(4);
  hTwoElectronTrackCaloMomentumPair->Draw("COLZ");
  cTwoElectronTrackCaloMomentum->SaveAs(twoElectronTrackCaloMomentumPdfName.c_str());

  // PDF group 5: log-scale version of the combined matched-momentum spectrum.
  // This is useful for checking small tails without changing the selected sample.
  const string twoElectronTrackCaloMomentumLogPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_TwoElectronTrackCaloMomentumLogY.pdf";
  TCanvas* cTwoElectronTrackCaloMomentumLog = new TCanvas(
    "cTwoElectronTrackCaloMomentumLog",
    "Two selected reconstructed electrons: track momentum at calo association log y",
    900, 700);
  cTwoElectronTrackCaloMomentumLog->SetLogy();
  hTwoElectronTrackCaloMomentumAll->SetMinimum(0.5);
  hTwoElectronTrackCaloMomentumAll->Draw("HIST");
  cTwoElectronTrackCaloMomentumLog->SaveAs(twoElectronTrackCaloMomentumLogPdfName.c_str());

  // PDF group 6: event-level reconstructed clusters.  These plots are not
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

  // PDF group 7: crystal-hit diagnostics.  The left plot counts individual
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

  cout << "Wrote calorimeter PDF plots to:" << endl;
  cout << "  " << trackCaloPdfName << endl;
  cout << "  " << twoElectronTrackCaloPdfName << endl;
  cout << "  " << twoElectronTrackCaloLogPdfName << endl;
  cout << "  " << twoElectronTrackCaloMomentumPdfName << endl;
  cout << "  " << twoElectronTrackCaloMomentumLogPdfName << endl;
  cout << "  " << clusterPdfName << endl;
  cout << "  " << crystalPdfName << endl;

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
         << " event=" << display.event << endl;
    calohitter::saveCalorimeterPdf(outputPdf, display.selectedHits, canvasName);
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

