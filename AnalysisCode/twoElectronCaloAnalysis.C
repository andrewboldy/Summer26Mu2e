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
//   4. Also draw and save a blank two-disk calorimeter map through CaloHitter.
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
                             const int maxSelectedEventsToPrint = -1)
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

  // Track-associated calorimeter histograms.
  //
  // These use trkcalohit, which is the reconstructed track-to-calo match.  The
  // energy is the matched calorimeter energy.  The momentum and POCA position
  // are the reconstructed track state at the calorimeter association.
  TH1F* hTrackCaloEnergy = new TH1F(
    "hTrackCaloEnergy",
    "Track-associated calorimeter energy;E_{calo} [MeV];Matched rank-0 electron tracks",
    150, 0.0, 150.0);
  TH1F* hTrackCaloMomentum = new TH1F(
    "hTrackCaloMomentum",
    "Track momentum at calorimeter association;p_{track} [MeV/c];Matched rank-0 electron tracks",
    150, 0.0, 150.0);
  TH2F* hTrackCaloPOCAXY = new TH2F(
    "hTrackCaloPOCAXY",
    "Track-calo POCA position;x [mm];y [mm]",
    200, -1000.0, 1000.0, 200, -1000.0, 1000.0);
  TH1F* hTrackCaloPOCAZ = new TH1F(
    "hTrackCaloPOCAZ",
    "Track-calo POCA z position;z [mm];Matched rank-0 electron tracks",
    240, -12000.0, 12000.0);

  // Event-level reconstructed calorimeter cluster histograms.
  //
  // caloclusters stores reconstructed energy and COG position.  EventNtuple's
  // COG is in the calorimeter disk front-face frame, not global Mu2e xyz.
  TH1F* hClusterEnergy = new TH1F(
    "hClusterEnergy",
    "Reconstructed calorimeter cluster energy;E_{cluster} [MeV];Clusters",
    150, 0.0, 150.0);
  TH1F* hClusterTime = new TH1F(
    "hClusterTime",
    "Reconstructed calorimeter cluster time;t_{cluster} [ns];Clusters",
    240, 0.0, 2400.0);
  TH1F* hClusterZ = new TH1F(
    "hClusterZ",
    "Reconstructed calorimeter cluster COG z;COG z [mm];Clusters",
    200, -1000.0, 1000.0);
  TH2F* hClusterXYByDisk[2] = {
    new TH2F(
      "hClusterXYDisk0",
      "Reconstructed cluster COG, Front disk;disk-local x [mm];disk-local y [mm]",
      200, -800.0, 800.0, 200, -800.0, 800.0),
    new TH2F(
      "hClusterXYDisk1",
      "Reconstructed cluster COG, Back disk;disk-local x [mm];disk-local y [mm]",
      200, -800.0, 800.0, 200, -800.0, 800.0)
  };

  // Crystal-hit histograms.  These use event.calohits, which provides crystal ID
  // and deposited energy.  The energy-by-ID histogram is the bridge to the
  // future CaloHitter overlay: each x bin corresponds to one calorimeter crystal.
  TH1F* hCrystalHitEnergy = new TH1F(
    "hCrystalHitEnergy",
    "Calorimeter crystal-hit energy;E_{hit} [MeV];Crystal hits",
    150, 0.0, 150.0);
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

  // Counters are separated into "selected" and "printed" because the optional
  // maxSelectedEventsToPrint limit can stop printing before the scan is done.
  long long selectedEventCount = 0;
  long long printedSelectedEventCount = 0;
  long long printedRank0ElectronCount = 0;
  long long printedTrackCaloCount = 0;
  long long printedClusterCount = 0;
  long long printedCrystalHitCount = 0;
  long long selectedEventsWithoutCaloClusters = 0;
  long long selectedRank0TracksWithoutCalo = 0;

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
      // association, not a per-crystal center position.
      hTrackCaloEnergy->Fill(trkcalohit->edep);
      hTrackCaloMomentum->Fill(trkcalohit->mom.R());
      hTrackCaloPOCAXY->Fill(trkcalohit->poca.x(), trkcalohit->poca.y());
      hTrackCaloPOCAZ->Fill(trkcalohit->poca.z());

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

    // The event-level caloclusters collection contains calorimeter clusters
    // independent of whether a particular track was matched to one.  A selected
    // two-electron event can still have no reconstructed calorimeter cluster.
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
        hClusterZ->Fill(cluster.cog_.z());
        if (cluster.diskID_ >= 0 && cluster.diskID_ < 2)
        {
          hClusterXYByDisk[cluster.diskID_]->Fill(cluster.cog_.x(), cluster.cog_.y());
        }

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
  outputFile << "\n# Summary\n"
             << "# selected_events " << selectedEventCount << '\n'
             << "# printed_selected_events " << printedSelectedEventCount << '\n'
             << "# printed_rank0_electrons " << printedRank0ElectronCount << '\n'
             << "# printed_track_calo_entries " << printedTrackCaloCount << '\n'
             << "# selected_rank0_tracks_without_track_calo " << selectedRank0TracksWithoutCalo << '\n'
             << "# printed_clusters " << printedClusterCount << '\n'
             << "# printed_crystal_hits " << printedCrystalHitCount << '\n'
             << "# selected_events_without_calo_clusters " << selectedEventsWithoutCaloClusters << '\n'
             << "# histogram_track_calo_entries " << hTrackCaloEnergy->GetEntries() << '\n'
             << "# histogram_cluster_entries " << hClusterEnergy->GetEntries() << '\n'
             << "# histogram_crystal_hit_entries " << hCrystalHitEnergy->GetEntries() << '\n';

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
  cout << "Wrote text output to " << outputFileName << endl;

  // Write the ROOT histograms to a file so the plotted distributions can be
  // reopened, rebinned, or overlaid later without rerunning the ntuple loop.
  const string histogramRootFileName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_RecoCaloHistograms.root";
  TFile histogramFile(histogramRootFileName.c_str(), "RECREATE");
  if (!histogramFile.IsZombie())
  {
    hTrackCaloEnergy->Write();
    hTrackCaloMomentum->Write();
    hTrackCaloPOCAXY->Write();
    hTrackCaloPOCAZ->Write();
    hClusterEnergy->Write();
    hClusterTime->Write();
    hClusterZ->Write();
    hClusterXYByDisk[0]->Write();
    hClusterXYByDisk[1]->Write();
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

  const string clusterPdfName =
    caloPlotsDirectory + "/twoElectronCaloAnalysis_" + generatorName + "_RecoCaloClusters.pdf";
  TCanvas* cClusters = new TCanvas(
    "cClusters",
    "Reconstructed calorimeter clusters",
    1600, 1000);
  cClusters->Divide(3, 2);
  cClusters->cd(1);
  hClusterEnergy->Draw("HIST");
  cClusters->cd(2);
  hClusterTime->Draw("HIST");
  cClusters->cd(3);
  hClusterZ->Draw("HIST");
  cClusters->cd(4);
  hClusterXYByDisk[0]->Draw("COLZ");
  cClusters->cd(5);
  hClusterXYByDisk[1]->Draw("COLZ");
  cClusters->SaveAs(clusterPdfName.c_str());

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
  cout << "  " << clusterPdfName << endl;
  cout << "  " << crystalPdfName << endl;

  // Draw the calorimeter after the event loop.  Right now CaloHitter only knows
  // the blank crystal geometry, but this location is intentional: later we can
  // use the accumulated event/crystal information from the analysis above to
  // alter crystal colors, labels, or hit markers before saving the diagram.
  const string blankCaloDiskPdfName =
    "Plots/CaloHitPlots/twoElectronCaloAnalysis_BlankCaloDisks_" + generatorName + ".pdf";
  cout << "Writing blank calorimeter disk PDF: " << blankCaloDiskPdfName << endl;
  calohitter::saveCalorimeterPdf(blankCaloDiskPdfName);

  timer.Stop();
  // Timing is useful when deciding whether to print all selected events or only
  // the first few during development.
  cout << "CPU time: " << timer.CpuTime() << " s, real time: " << timer.RealTime() << " s" << endl;

  TH1::AddDirectory(oldAddDirectoryStatus);
  gROOT->SetBatch(wasBatchMode);
}
