//----------------------------------------------------------------------------------
//
// twoElectronTrkSegAnalysisSelector.C
// Written by Andrew Boldy University of South Carolina, 2026
// Assisted by Codex
//
// Purpose:
//   Loop over an EventNtuple ROOT file or filelist and print reconstruction-only
//   information for tracks whose fitted particle hypothesis is a downstream
//   electron.
//
//   In EventNtuple, the reconstructed track summary is stored in a vector branch
//   called "trk" by default.  Each element of "trk" is one reconstructed track.
//   The corresponding track-state samples at detector surfaces are stored in
//   "trksegs", which is a vector of vectors:
//
//       trk[i_track]          -> reconstructed track summary
//       trksegs[i_track]      -> all stored surface intersections for that track
//       trksegs[i_track][j]   -> one fitted state at one surface
//
//   The "trksegs" objects contain the reconstructed momentum vector, position,
//   time, surface ID, and surface index.  This macro prints all of those segment
//   momenta for every reconstructed downstream electron track it finds.  The
//   segment printout is sorted by reconstructed segment time so the output is
//   easier to follow by hand as a track path.
//
// Important:
//   This macro intentionally does not use Monte Carlo truth information.  It turns
//   off all branches first, then enables only evtinfo, the requested track branch,
//   and the corresponding track-segment branch.
//
//   In particular, it never enables or reads:
//       trkmcsim
//       trkmc
//       trksegsmc
//       any other MC truth branch
//
//   The downstream-electron decision is based only on reconstructed quantities:
//
//       trk.pdg == 11
//       and
//       at least one reconstructed TT_Mid segment has p_z > 0
//
//   That means "downstream electron track" here means "the reconstruction fit
//   treated this track as an e-minus and the fitted track state is moving in the
//   +z direction at the tracker middle."  It does not mean "MC truth says this
//   was a downstream electron."
//
// Usage from ROOT:
//   .L CreatedCode/twoElectronTrkSegAnalysisSelector.C+
//   twoElectronTrkSegAnalysisSelector("path/to/nts.root")
//
// Optional arguments:
//   twoElectronTrkSegAnalysisSelector("filelist.txt", 100)          // first 100 events
//   twoElectronTrkSegAnalysisSelector("nts.root", -1, "de")         // use dedicated de/desegs
//   twoElectronTrkSegAnalysisSelector("nts.root", -1, "trk", -11)   // downstream e-plus fits
//   twoElectronTrkSegAnalysisSelector("nts.root", -1, "trk", 11, true)
//                                                           // full trkseg printout
//   twoElectronTrkSegAnalysisSelector("nts.root", -1, "trk", 11, false, "reduced.root")
//                                                           // write reduced candidate ntuple
//
// The default "trk" branch is useful for ntuples written with the common
// EventNtuple "All" branch.  If your file contains the dedicated downstream
// electron branch configured as "de", passing trackBranch="de" reads de/desegs
// directly and still applies the same reconstruction-only downstream check.
//
// By default, doFullTrkSegPrintout is false so the macro prints only the
// reconstructed-track and event-level summary counts.  Set it true when you want
// the full per-track, per-segment momentum/path printout.
//
// The final summary also counts events with more than one selected reconstructed
// downstream electron track and events with exactly two selected reconstructed
// downstream electron tracks.  These are reconstruction-only event counts: they
// do not prove, by MC truth, that distinct physical electrons were thrown.
//
// For exactly-two-track events, the summary also counts how many selected
// electron tracks have 0, 1, 2, 3, ... reconstructed trksegs on stopping-target
// foils, where ST_Foils is sid == 104.
//
// For each exactly-two-track event, the macro also prints which reconstructed
// surfaces are shared by both selected tracks when doFullTrkSegPrintout is true.
// A shared surface means the two selected tracks both have a reconstructed
// trkseg with the same surface id and surface index.
//
//----------------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include "EventNtuple/inc/EventInfo.hh"
#include "EventNtuple/inc/TrkCaloHitInfo.hh"
#include "EventNtuple/inc/TrkInfo.hh"
#include "EventNtuple/inc/TrkSegInfo.hh"

using namespace std;

namespace
{
  //============================================================================
  // Input Helpers
  //============================================================================

  // ROOT input can be either:
  //   1. one EventNtuple ROOT file, ending in ".root", or
  //   2. a text filelist, with one ROOT filename per line.
  //
  // This tiny helper decides which case we are in.
  bool hasRootSuffix(const string& path)
  {
      const string suffix = ".root";

      if (path.size() < suffix.size())
          return false;

      return path.substr(path.size() - suffix.size()) == suffix;
  }

  // EventNtuple branches are split ROOT branches.  Enabling only "trk" is not
  // always enough because leaves can appear as "trk.pdg", "trk.nhits", etc.
  // This helper enables both the top-level branch name and all split leaves
  // underneath it.
  void enableBranch(TChain& chain, const string& branchName)
  {
    chain.SetBranchStatus(branchName.c_str(), 1);
    chain.SetBranchStatus((branchName + ".*").c_str(), 1);
  }

  // Add either a single ROOT file or every ROOT file listed in a filelist to the
  // TChain.  A TChain behaves like one long TTree, even when it is built from
  // many files.
  bool addInput(TChain& chain, const string& inputName)
  {
    // If the input itself is a ROOT file, add it directly.
    if (hasRootSuffix(inputName))
    {
      chain.Add(inputName.c_str());
      return true;
    }

    // Otherwise treat the input as a plain text filelist.
    ifstream filelist(inputName);
    if (!filelist.is_open())
    {
      cerr << "ERROR: could not open input file or filelist: " << inputName << endl;
      return false;
    }

    string line;
    int nFiles = 0;
    while (getline(filelist, line))
    {
      // Allow blank lines and comments in the filelist.  This makes it easier
      // to temporarily disable an input file without editing the macro.
      if (line.empty() || line[0] == '#')
      {
        continue;
      }

      chain.Add(line.c_str());
      ++nFiles;
    }

    if (nFiles == 0)
    {
      cerr << "ERROR: filelist contains no ROOT files: " << inputName << endl;
      return false;
    }

    return true;
  }

  //============================================================================
  // Reconstruction-Only Track Selection Helpers
  //============================================================================

  // Reconstruction-only direction test.
  //
  // EventNtuple/RooUtil convention identifies downstream-going reconstructed
  // tracks by the sign of the reconstructed momentum's z component at the
  // tracker middle surface:
  //
  //   sid == mu2e::SurfaceIdDetail::TT_Mid
  //   pz  > 0
  //
  // This uses only trksegs.  It does not inspect trkmcsim, trkmc, trksegsmc, or
  // any other truth branch.
  bool hasDownstreamTrackerMiddleSegment(const vector<mu2e::TrkSegInfo>& segments)
  {
    for (const auto& segment : segments)
    {
      if (segment.sid == mu2e::SurfaceIdDetail::TT_Mid && segment.mom.z() > 0.0)
      {
        return true;
      }
    }

    return false;
  }

  // Full downstream-electron reconstructed-track selection.
  //
  // The PDG code is the fit hypothesis stored by reconstruction in TrkInfo.
  // The downstream requirement comes from the reconstructed segment direction.
  bool isSelectedDownstreamElectronTrack(const mu2e::TrkInfo& track,
                                         const vector<mu2e::TrkSegInfo>& segments,
                                         int electronPdg)
  {
    return track.pdg == electronPdg && hasDownstreamTrackerMiddleSegment(segments);
  }

  // Count how many reconstructed track-segment intersections for one track are
  // on stopping-target foils.  In the SurfaceId enum, ST_Foils is numeric sid
  // 104; sindex then identifies which foil was intersected.
  size_t countStoppingTargetFoilSegments(const vector<mu2e::TrkSegInfo>& segments)
  {
    size_t nFoilSegments = 0;

    for (const auto& segment : segments)
    {
      if (segment.sid == mu2e::SurfaceIdDetail::ST_Foils)
      {
        ++nFoilSegments;
      }
    }

    return nFoilSegments;
  }

  //============================================================================
  // Shared-Surface Helpers
  //============================================================================

  struct SharedSurfaceMatch
  {
    int sid = -1;
    int sindex = -1;
    size_t firstStoredSegmentIndex = 0;
    size_t secondStoredSegmentIndex = 0;
  };

  using SharedSurfaceMatches = vector<SharedSurfaceMatch>;

  SharedSurfaceMatches sharedSurfacesBetweenTracks(
    const vector<mu2e::TrkSegInfo>& firstTrackSegments,
    const vector<mu2e::TrkSegInfo>& secondTrackSegments)
  {
    map<pair<int, int>, size_t> firstSurfaceToSegmentIndex;
    map<pair<int, int>, size_t> secondSurfaceToSegmentIndex;

    for (size_t iFirstSegment = 0;
         iFirstSegment < firstTrackSegments.size();
         ++iFirstSegment)
    {
      const auto& firstSegment = firstTrackSegments.at(iFirstSegment);
      const auto key = make_pair(firstSegment.sid, firstSegment.sindex);
      if (firstSurfaceToSegmentIndex.find(key) == firstSurfaceToSegmentIndex.end())
      {
        firstSurfaceToSegmentIndex.emplace(key, iFirstSegment);
      }
    }

    for (size_t iSecondSegment = 0;
         iSecondSegment < secondTrackSegments.size();
         ++iSecondSegment)
    {
      const auto& secondSegment = secondTrackSegments.at(iSecondSegment);
      const auto key = make_pair(secondSegment.sid, secondSegment.sindex);
      if (secondSurfaceToSegmentIndex.find(key) == secondSurfaceToSegmentIndex.end())
      {
        secondSurfaceToSegmentIndex.emplace(key, iSecondSegment);
      }
    }

    SharedSurfaceMatches sharedSurfaces;
    for (const auto& firstEntry : firstSurfaceToSegmentIndex)
    {
      const auto secondIter = secondSurfaceToSegmentIndex.find(firstEntry.first);
      if (secondIter == secondSurfaceToSegmentIndex.end())
      {
        continue;
      }

      SharedSurfaceMatch sharedSurface;
      sharedSurface.sid = firstEntry.first.first;
      sharedSurface.sindex = firstEntry.first.second;
      sharedSurface.firstStoredSegmentIndex = firstEntry.second;
      sharedSurface.secondStoredSegmentIndex = secondIter->second;
      sharedSurfaces.push_back(sharedSurface);
    }

    return sharedSurfaces;
  }

  string formatSharedSurface(const SharedSurfaceMatch& sharedSurface)
  {
    ostringstream out;
    out << mu2e::SurfaceId(static_cast<mu2e::SurfaceIdDetail::enum_type>(sharedSurface.sid),
                           sharedSurface.sindex);
    return out.str();
  }

  //============================================================================
  // Printout Ordering Helpers
  //============================================================================

  // Return the original segment-vector indices sorted by reconstructed segment
  // time.  The EventNtuple producer stores trksegs in the KalSeed intersection
  // insertion order, not chronological order, so the macro sorts for printing
  // only.  Keeping indices instead of copying TrkSegInfo objects lets the
  // printout report the original stored segment index.
  vector<size_t> sortedSegmentIndicesByTime(const vector<mu2e::TrkSegInfo>& segments)
  {
    vector<size_t> sortedIndices;
    sortedIndices.reserve(segments.size());

    for (size_t iSegment = 0; iSegment < segments.size(); ++iSegment)
    {
      sortedIndices.push_back(iSegment);
    }

    stable_sort(sortedIndices.begin(),
                sortedIndices.end(),
                [&segments](size_t left, size_t right)
                {
                  return segments.at(left).time < segments.at(right).time;
                });

    return sortedIndices;
  }
}

void twoElectronTrkSegAnalysisSelector(const string& inputName,
                                       int maxEvents = -1,
                                       const string& trackBranch = "trk",
                                       int electronPdg = 11,
                                       bool doFullTrkSegPrintout = false,
                                       const string& reducedOutputName = "")
{
  //============================================================================
  // Main Macro Setup: Branch Names, Input Files, and Tree Validation
  //============================================================================

  // EventNtuple names the track-segment branch by appending "segs" to the track
  // branch prefix.  The common/all-track branch is "trk" -> "trksegs".  A
  // dedicated downstream-electron branch, when present in the ntuple, is
  // "de" -> "desegs".
  const string segmentBranch = trackBranch + "segs";

  // Build a TChain over the EventNtuple tree.  The tree path is fixed by the
  // EventNtuple maker: directory "EventNtuple", tree "ntuple".
  TChain ntuple("EventNtuple/ntuple");
  if (!addInput(ntuple, inputName))
  {
    return;
  }

  // A quick sanity check catches empty files or a filelist whose files do not
  // contain the expected EventNtuple tree.
  const Long64_t nEntries = ntuple.GetEntries();
  if (nEntries <= 0)
  {
    cerr << "ERROR: no entries found in EventNtuple/ntuple for input: " << inputName << endl;
    return;
  }

  // Confirm that the reconstruction branches requested by the user actually
  // exist before trying to set branch addresses.
  if (ntuple.GetBranch(trackBranch.c_str()) == nullptr)
  {
    cerr << "ERROR: missing requested track branch '" << trackBranch << "'." << endl;
    return;
  }

  if (ntuple.GetBranch(segmentBranch.c_str()) == nullptr)
  {
    cerr << "ERROR: missing requested track-segment branch '" << segmentBranch << "'." << endl;
    return;
  }

  //============================================================================
  // Branch Selection: Reconstruction Only, No Monte Carlo Truth
  //============================================================================

  // This is the key no-MC-truth safeguard:
  //   1. Turn off every branch.
  //   2. Turn on only evtinfo, the reconstructed track branch, and trksegs.
  //
  // After this, ROOT will not load trkmcsim or any other truth branch when
  // ntuple.GetEntry(...) is called.
  ntuple.SetBranchStatus("*", 0);
  enableBranch(ntuple, "evtinfo");
  enableBranch(ntuple, trackBranch);
  enableBranch(ntuple, segmentBranch);

  // Pointers that ROOT will fill for each TTree entry.  These point directly to
  // the branch data owned by ROOT; do not delete them.
  mu2e::EventInfo* evtinfo = nullptr;
  vector<mu2e::TrkInfo>* tracks = nullptr;
  vector<vector<mu2e::TrkSegInfo>>* trackSegments = nullptr;
  vector<mu2e::TrkCaloHitInfo>* trackCaloHits = nullptr;

  // evtinfo is helpful for printing run/subrun/event identifiers, but a few
  // minimal ntuples might omit it.  The analysis can still run without it.
  if (ntuple.GetBranch("evtinfo") != nullptr)
  {
    ntuple.SetBranchAddress("evtinfo", &evtinfo);
  }

  // Connect the requested reconstruction branches to the C++ pointers above.
  // After GetEntry(i), "tracks" and "trackSegments" point at entry i's data.
  ntuple.SetBranchAddress(trackBranch.c_str(), &tracks);
  ntuple.SetBranchAddress(segmentBranch.c_str(), &trackSegments);
  if (ntuple.GetBranch("trkcalohit") != nullptr)
  {
    enableBranch(ntuple, "trkcalohit");
    ntuple.SetBranchAddress("trkcalohit", &trackCaloHits);
  }

  //============================================================================
  // Run Configuration and Analysis Counters
  //============================================================================

  // maxEvents lets you test the macro on the first few events before printing a
  // very large file.  A negative maxEvents means scan everything.
  const Long64_t entriesToRead =
    (maxEvents >= 0 && static_cast<Long64_t>(maxEvents) < nEntries) ? maxEvents : nEntries;

  // Print a short configuration summary so the output log records exactly what
  // was scanned and which branch prefix/PDG criterion was used.
  cout << "Input: " << inputName << endl;
  cout << "Tree entries available: " << nEntries << endl;
  cout << "Tree entries being scanned: " << entriesToRead << endl;
  cout << "Track branch: " << trackBranch << endl;
  cout << "Track-segment branch: " << segmentBranch << endl;
  cout << "Downstream electron track criterion: " << trackBranch << ".pdg == " << electronPdg
       << " and " << segmentBranch << " has TT_Mid pz > 0"
       << " (reconstructed quantities only)" << endl;
  cout << "Full track-segment printout: "
       << (doFullTrkSegPrintout ? "enabled" : "disabled")
       << " (summary counts and percentages are always printed)" << endl;
  if (!reducedOutputName.empty())
  {
    cout << "Reduced output file: " << reducedOutputName << endl;
  }
  if (doFullTrkSegPrintout)
  {
    cout << "Segment print order: increasing reconstructed segment time; "
         << "stored_seg_index preserves original ntuple order." << endl;
  }
  if (ntuple.GetBranch("trkcalohit") != nullptr)
  {
    cout << "Track-calo association branch: trkcalohit (reduced output will record it)" << endl;
  }
  else if (!reducedOutputName.empty())
  {
    cout << "Track-calo association branch: trkcalohit not present; reduced output will mark those fields missing" << endl;
  }
  cout << "MC truth branches are not enabled or read by this macro." << endl;

  // First pass counter.  The user asked to begin by printing the number of
  // stored tracks that correspond to a downstream electron, so the macro counts
  // all selected reconstructed downstream electron tracks before printing any
  // per-track momenta.
  Long64_t totalDownstreamElectronTracks = 0;

  // Event-level counters for possible multi-electron reconstruction candidates.
  // These use the same reco-only downstream electron track selection used
  // everywhere else in this macro.
  Long64_t eventsWithMoreThanOneDownstreamElectronTrack = 0;
  Long64_t eventsWithExactlyTwoDownstreamElectronTracks = 0;

  // For exactly-two-track events, count how many selected electron tracks have
  // 0, 1, 2, 3, ... reconstructed trksegs at stopping-target foils.  This map's
  // key is the number of ST_Foils segments on one selected electron track, and
  // its value is the number of selected electron tracks with that multiplicity.
  map<size_t, Long64_t> foilSegmentMultiplicityForExactlyTwoTrackEvents;

  // Shared-surface counters for exactly-two-track events.  A shared surface is
  // one where both selected tracks have a reconstructed trkseg with the same sid
  // and sindex.
  Long64_t exactlyTwoTrackEventsWithSharedSurfaces = 0;
  Long64_t totalSharedSurfacesAcrossExactlyTwoTrackEvents = 0;
  Long64_t exactlyTwoTrackEventsWithSharedSTFoils = 0;
  Long64_t totalSharedSTFoilsAcrossExactlyTwoTrackEvents = 0;

  // Reduced-output tree.  This is the compact vertex-candidate schema.
  //
  // It stores one row per event, but only for exactly-two-track events with at
  // least one shared surface.  The row keeps:
  //   - event identity so the candidate can be traced back to the source ntuple
  //   - the two selected track indices and summary fit-quality fields
  //   - the track-calo association fields from trkcalohit, when available
  //   - the shared-surface state needed to extrapolate lines and vertex them
  Long64_t reducedEventsWritten = 0;
  bool writeReducedOutput = !reducedOutputName.empty();
  TFile* reducedFile = nullptr;
  TTree* reducedTree = nullptr;
  mu2e::EventInfo reducedEvtinfo;
  Long64_t reducedSourceEntry = -1;
  vector<int> reducedSelectedTrackIndices;
  vector<int> reducedSelectedTrackPdg;
  vector<int> reducedSelectedTrackNhits;
  vector<int> reducedSelectedTrackStatus;
  vector<int> reducedSelectedTrackGoodfit;
  vector<int> reducedSelectedTrackNdof;
  vector<int> reducedSelectedTrackNactive;
  vector<int> reducedSelectedTrackNseg;
  vector<int> reducedSelectedTrackNipadown;
  vector<int> reducedSelectedTrackNstdown;
  vector<int> reducedSelectedTrackFirstStInter;
  vector<int> reducedSelectedTrackNstup;
  vector<int> reducedSelectedTrackNipaup;
  vector<int> reducedSelectedTrackCaloActive;
  vector<int> reducedSelectedTrackCaloDid;
  vector<double> reducedSelectedTrackFitcon;
  vector<double> reducedSelectedTrackChisq;
  vector<double> reducedSelectedTrackCaloPOCAX;
  vector<double> reducedSelectedTrackCaloPOCAY;
  vector<double> reducedSelectedTrackCaloPOCAZ;
  vector<double> reducedSelectedTrackCaloMomX;
  vector<double> reducedSelectedTrackCaloMomY;
  vector<double> reducedSelectedTrackCaloMomZ;
  vector<double> reducedSelectedTrackCaloCDepth;
  vector<double> reducedSelectedTrackCaloTrkDepth;
  vector<double> reducedSelectedTrackCaloDphiDot;
  vector<double> reducedSelectedTrackCaloDoca;
  vector<double> reducedSelectedTrackCaloDt;
  vector<double> reducedSelectedTrackCaloPtoca;
  vector<double> reducedSelectedTrackCaloTocavar;
  vector<double> reducedSelectedTrackCaloTresid;
  vector<double> reducedSelectedTrackCaloTresidmvar;
  vector<double> reducedSelectedTrackCaloTresidpvar;
  vector<double> reducedSelectedTrackCaloCTime;
  vector<double> reducedSelectedTrackCaloCTimeErr;
  vector<double> reducedSelectedTrackCaloCSize;
  vector<double> reducedSelectedTrackCaloEdep;
  vector<double> reducedSelectedTrackCaloEdepErr;
  vector<int> reducedSelectedTrackTotalSegmentCount;
  vector<int> reducedSelectedTrackFoilSegmentCount;
  vector<int> reducedSharedSurfaceSid;
  vector<int> reducedSharedSurfaceSindex;
  vector<int> reducedSharedSurfaceFirstStoredSegmentIndex;
  vector<int> reducedSharedSurfaceSecondStoredSegmentIndex;
  vector<double> reducedSharedSurfaceFirstPosX;
  vector<double> reducedSharedSurfaceFirstPosY;
  vector<double> reducedSharedSurfaceFirstPosZ;
  vector<double> reducedSharedSurfaceFirstMomX;
  vector<double> reducedSharedSurfaceFirstMomY;
  vector<double> reducedSharedSurfaceFirstMomZ;
  vector<double> reducedSharedSurfaceFirstTime;
  vector<double> reducedSharedSurfaceSecondPosX;
  vector<double> reducedSharedSurfaceSecondPosY;
  vector<double> reducedSharedSurfaceSecondPosZ;
  vector<double> reducedSharedSurfaceSecondMomX;
  vector<double> reducedSharedSurfaceSecondMomY;
  vector<double> reducedSharedSurfaceSecondMomZ;
  vector<double> reducedSharedSurfaceSecondTime;
  int reducedSharedSurfaceCount = 0;
  int reducedSharedSTFoilCount = 0;

  if (writeReducedOutput)
  {
    reducedFile = TFile::Open(reducedOutputName.c_str(), "RECREATE");
    if (reducedFile == nullptr || reducedFile->IsZombie())
    {
      cerr << "ERROR: could not create reduced output file: " << reducedOutputName << endl;
      writeReducedOutput = false;
    }
    else
    {
      reducedFile->cd();
      reducedTree = new TTree("TwoElectronTrackVertexCandidates",
                              "Reduced two downstream electron track events");

      // Event identity and source bookkeeping.
      reducedTree->Branch("evtinfo", &reducedEvtinfo);
      reducedTree->Branch("sourceEntry", &reducedSourceEntry);

      // Track identity and fit-quality summary for the two selected tracks.
      reducedTree->Branch("selectedTrackIndices", &reducedSelectedTrackIndices);
      reducedTree->Branch("selectedTrackPdg", &reducedSelectedTrackPdg);
      reducedTree->Branch("selectedTrackNhits", &reducedSelectedTrackNhits);
      reducedTree->Branch("selectedTrackStatus", &reducedSelectedTrackStatus);
      reducedTree->Branch("selectedTrackGoodfit", &reducedSelectedTrackGoodfit);
      reducedTree->Branch("selectedTrackNdof", &reducedSelectedTrackNdof);
      reducedTree->Branch("selectedTrackNactive", &reducedSelectedTrackNactive);
      reducedTree->Branch("selectedTrackNseg", &reducedSelectedTrackNseg);
      reducedTree->Branch("selectedTrackNipadown", &reducedSelectedTrackNipadown);
      reducedTree->Branch("selectedTrackNstdown", &reducedSelectedTrackNstdown);
      reducedTree->Branch("selectedTrackFirstStInter", &reducedSelectedTrackFirstStInter);
      reducedTree->Branch("selectedTrackNstup", &reducedSelectedTrackNstup);
      reducedTree->Branch("selectedTrackNipaup", &reducedSelectedTrackNipaup);
      reducedTree->Branch("selectedTrackFitcon", &reducedSelectedTrackFitcon);
      reducedTree->Branch("selectedTrackChisq", &reducedSelectedTrackChisq);

      // Track-calo association fields.  These come from trkcalohit and capture
      // the reconstructed calorimeter impact point, direction, timing, and
      // deposited energy for each selected track.
      reducedTree->Branch("selectedTrackCaloActive", &reducedSelectedTrackCaloActive);
      reducedTree->Branch("selectedTrackCaloDid", &reducedSelectedTrackCaloDid);
      reducedTree->Branch("selectedTrackCaloPOCAX", &reducedSelectedTrackCaloPOCAX);
      reducedTree->Branch("selectedTrackCaloPOCAY", &reducedSelectedTrackCaloPOCAY);
      reducedTree->Branch("selectedTrackCaloPOCAZ", &reducedSelectedTrackCaloPOCAZ);
      reducedTree->Branch("selectedTrackCaloMomX", &reducedSelectedTrackCaloMomX);
      reducedTree->Branch("selectedTrackCaloMomY", &reducedSelectedTrackCaloMomY);
      reducedTree->Branch("selectedTrackCaloMomZ", &reducedSelectedTrackCaloMomZ);
      reducedTree->Branch("selectedTrackCaloCDepth", &reducedSelectedTrackCaloCDepth);
      reducedTree->Branch("selectedTrackCaloTrkDepth", &reducedSelectedTrackCaloTrkDepth);
      reducedTree->Branch("selectedTrackCaloDphiDot", &reducedSelectedTrackCaloDphiDot);
      reducedTree->Branch("selectedTrackCaloDoca", &reducedSelectedTrackCaloDoca);
      reducedTree->Branch("selectedTrackCaloDt", &reducedSelectedTrackCaloDt);
      reducedTree->Branch("selectedTrackCaloPtoca", &reducedSelectedTrackCaloPtoca);
      reducedTree->Branch("selectedTrackCaloTocavar", &reducedSelectedTrackCaloTocavar);
      reducedTree->Branch("selectedTrackCaloTresid", &reducedSelectedTrackCaloTresid);
      reducedTree->Branch("selectedTrackCaloTresidmvar", &reducedSelectedTrackCaloTresidmvar);
      reducedTree->Branch("selectedTrackCaloTresidpvar", &reducedSelectedTrackCaloTresidpvar);
      reducedTree->Branch("selectedTrackCaloCTime", &reducedSelectedTrackCaloCTime);
      reducedTree->Branch("selectedTrackCaloCTimeErr", &reducedSelectedTrackCaloCTimeErr);
      reducedTree->Branch("selectedTrackCaloCSize", &reducedSelectedTrackCaloCSize);
      reducedTree->Branch("selectedTrackCaloEdep", &reducedSelectedTrackCaloEdep);
      reducedTree->Branch("selectedTrackCaloEdepErr", &reducedSelectedTrackCaloEdepErr);

      // Segment-count bookkeeping for downstream-electron diagnostics.
      reducedTree->Branch("selectedTrackTotalSegmentCount", &reducedSelectedTrackTotalSegmentCount);
      reducedTree->Branch("selectedTrackFoilSegmentCount", &reducedSelectedTrackFoilSegmentCount);

      // Shared-surface geometry and per-track state at the shared surface.
      reducedTree->Branch("sharedSurfaceCount", &reducedSharedSurfaceCount);
      reducedTree->Branch("sharedSTFoilCount", &reducedSharedSTFoilCount);
      reducedTree->Branch("sharedSurfaceSid", &reducedSharedSurfaceSid);
      reducedTree->Branch("sharedSurfaceSindex", &reducedSharedSurfaceSindex);
      reducedTree->Branch("sharedSurfaceFirstStoredSegmentIndex",
                          &reducedSharedSurfaceFirstStoredSegmentIndex);
      reducedTree->Branch("sharedSurfaceSecondStoredSegmentIndex",
                          &reducedSharedSurfaceSecondStoredSegmentIndex);
      reducedTree->Branch("sharedSurfaceFirstPosX", &reducedSharedSurfaceFirstPosX);
      reducedTree->Branch("sharedSurfaceFirstPosY", &reducedSharedSurfaceFirstPosY);
      reducedTree->Branch("sharedSurfaceFirstPosZ", &reducedSharedSurfaceFirstPosZ);
      reducedTree->Branch("sharedSurfaceFirstMomX", &reducedSharedSurfaceFirstMomX);
      reducedTree->Branch("sharedSurfaceFirstMomY", &reducedSharedSurfaceFirstMomY);
      reducedTree->Branch("sharedSurfaceFirstMomZ", &reducedSharedSurfaceFirstMomZ);
      reducedTree->Branch("sharedSurfaceFirstTime", &reducedSharedSurfaceFirstTime);
      reducedTree->Branch("sharedSurfaceSecondPosX", &reducedSharedSurfaceSecondPosX);
      reducedTree->Branch("sharedSurfaceSecondPosY", &reducedSharedSurfaceSecondPosY);
      reducedTree->Branch("sharedSurfaceSecondPosZ", &reducedSharedSurfaceSecondPosZ);
      reducedTree->Branch("sharedSurfaceSecondMomX", &reducedSharedSurfaceSecondMomX);
      reducedTree->Branch("sharedSurfaceSecondMomY", &reducedSharedSurfaceSecondMomY);
      reducedTree->Branch("sharedSurfaceSecondMomZ", &reducedSharedSurfaceSecondMomZ);
      reducedTree->Branch("sharedSurfaceSecondTime", &reducedSharedSurfaceSecondTime);
    }
  }

  // The user wants a percentage relative to the number of generated/thrown
  // events.  Keep this explicit here so it is easy to change if the input sample
  // size changes later.
  const double nThrownEvents = 100000.0;

  // Use fixed-width decimal formatting for the printed momenta/positions so the
  // output is easier to scan by eye or redirect to a text file.
  cout << fixed << setprecision(3);

  //============================================================================
  // Pass 1: Count Tracks and Compute Exactly-Two-Track Event Diagnostics
  //============================================================================
  //
  // This pass reads the same reconstruction-only branches as the printout pass,
  // and inspects only track.pdg plus trksegs reconstructed momenta.  No MC truth
  // branch is enabled or touched.  It also computes all event-level quantities
  // needed for the final summary, so quiet mode can still print the summary even
  // when doFullTrkSegPrintout is false.
  for (Long64_t iEntry = 0; iEntry < entriesToRead; ++iEntry)
  {
    // Load this event's enabled branch data into evtinfo/tracks/trackSegments.
    ntuple.GetEntry(iEntry);

    // For the downstream count, both the reconstructed track vector and the
    // reconstructed track-segment vector are needed.
    if (tracks == nullptr || trackSegments == nullptr)
    {
      cerr << "WARNING: null track or track-segment data at entry " << iEntry
           << "; skipping." << endl;
      continue;
    }

    const size_t nTracks = tracks->size();
    Long64_t downstreamElectronTracksInThisEvent = 0;
    vector<size_t> selectedTrackIndices;
    vector<size_t> foilSegmentCountsForSelectedTracks;
    for (size_t iTrack = 0; iTrack < nTracks; ++iTrack)
    {
      if (iTrack >= trackSegments->size())
      {
        continue;
      }

      if (isSelectedDownstreamElectronTrack(tracks->at(iTrack),
                                            trackSegments->at(iTrack),
                                            electronPdg))
      {
        ++totalDownstreamElectronTracks;
        ++downstreamElectronTracksInThisEvent;
        selectedTrackIndices.push_back(iTrack);
        foilSegmentCountsForSelectedTracks.push_back(
          countStoppingTargetFoilSegments(trackSegments->at(iTrack)));
      }
    }

    if (downstreamElectronTracksInThisEvent > 1)
    {
      ++eventsWithMoreThanOneDownstreamElectronTrack;
    }

    if (downstreamElectronTracksInThisEvent == 2)
    {
      //----------------------------------------------------------------------
      // Exactly-Two-Track Event Bookkeeping
      //----------------------------------------------------------------------
      //
      // Only events with exactly two selected downstream electron tracks enter
      // the foil multiplicity and shared-surface diagnostics.  Events with three
      // or more selected tracks are counted above, but excluded here so there is
      // no ambiguity about which pair is being compared.
      ++eventsWithExactlyTwoDownstreamElectronTracks;

      for (const size_t nFoilSegments : foilSegmentCountsForSelectedTracks)
      {
        ++foilSegmentMultiplicityForExactlyTwoTrackEvents[nFoilSegments];
      }

      const size_t firstTrackIndex = selectedTrackIndices.at(0);
      const size_t secondTrackIndex = selectedTrackIndices.at(1);
      const auto& firstTrackSegments = trackSegments->at(firstTrackIndex);
      const auto& secondTrackSegments = trackSegments->at(secondTrackIndex);
      const SharedSurfaceMatches sharedSurfaces =
        sharedSurfacesBetweenTracks(firstTrackSegments, secondTrackSegments);

      if (!sharedSurfaces.empty())
      {
        ++exactlyTwoTrackEventsWithSharedSurfaces;
      }
      totalSharedSurfacesAcrossExactlyTwoTrackEvents +=
        static_cast<Long64_t>(sharedSurfaces.size());

      vector<SharedSurfaceMatch> sharedSTFoils;
      for (const auto& sharedSurface : sharedSurfaces)
      {
        if (sharedSurface.sid == mu2e::SurfaceIdDetail::ST_Foils)
        {
          sharedSTFoils.push_back(sharedSurface);
        }
      }
      if (!sharedSTFoils.empty())
      {
        ++exactlyTwoTrackEventsWithSharedSTFoils;
      }
      totalSharedSTFoilsAcrossExactlyTwoTrackEvents +=
        static_cast<Long64_t>(sharedSTFoils.size());

      if (writeReducedOutput && reducedTree != nullptr && !sharedSurfaces.empty())
      {
        reducedEvtinfo = evtinfo != nullptr ? *evtinfo : mu2e::EventInfo{};
        reducedSourceEntry = iEntry;

        reducedSelectedTrackIndices.clear();
        reducedSelectedTrackPdg.clear();
        reducedSelectedTrackNhits.clear();
        reducedSelectedTrackStatus.clear();
        reducedSelectedTrackGoodfit.clear();
        reducedSelectedTrackNdof.clear();
        reducedSelectedTrackNactive.clear();
        reducedSelectedTrackNseg.clear();
        reducedSelectedTrackNipadown.clear();
        reducedSelectedTrackNstdown.clear();
        reducedSelectedTrackFirstStInter.clear();
        reducedSelectedTrackNstup.clear();
        reducedSelectedTrackNipaup.clear();
        reducedSelectedTrackFitcon.clear();
        reducedSelectedTrackChisq.clear();
        reducedSelectedTrackCaloActive.clear();
        reducedSelectedTrackCaloDid.clear();
        reducedSelectedTrackCaloPOCAX.clear();
        reducedSelectedTrackCaloPOCAY.clear();
        reducedSelectedTrackCaloPOCAZ.clear();
        reducedSelectedTrackCaloMomX.clear();
        reducedSelectedTrackCaloMomY.clear();
        reducedSelectedTrackCaloMomZ.clear();
        reducedSelectedTrackCaloCDepth.clear();
        reducedSelectedTrackCaloTrkDepth.clear();
        reducedSelectedTrackCaloDphiDot.clear();
        reducedSelectedTrackCaloDoca.clear();
        reducedSelectedTrackCaloDt.clear();
        reducedSelectedTrackCaloPtoca.clear();
        reducedSelectedTrackCaloTocavar.clear();
        reducedSelectedTrackCaloTresid.clear();
        reducedSelectedTrackCaloTresidmvar.clear();
        reducedSelectedTrackCaloTresidpvar.clear();
        reducedSelectedTrackCaloCTime.clear();
        reducedSelectedTrackCaloCTimeErr.clear();
        reducedSelectedTrackCaloCSize.clear();
        reducedSelectedTrackCaloEdep.clear();
        reducedSelectedTrackCaloEdepErr.clear();
        reducedSelectedTrackTotalSegmentCount.clear();
        reducedSelectedTrackFoilSegmentCount.clear();
        reducedSharedSurfaceSid.clear();
        reducedSharedSurfaceSindex.clear();
        reducedSharedSurfaceFirstStoredSegmentIndex.clear();
        reducedSharedSurfaceSecondStoredSegmentIndex.clear();
        reducedSharedSurfaceFirstPosX.clear();
        reducedSharedSurfaceFirstPosY.clear();
        reducedSharedSurfaceFirstPosZ.clear();
        reducedSharedSurfaceFirstMomX.clear();
        reducedSharedSurfaceFirstMomY.clear();
        reducedSharedSurfaceFirstMomZ.clear();
        reducedSharedSurfaceFirstTime.clear();
        reducedSharedSurfaceSecondPosX.clear();
        reducedSharedSurfaceSecondPosY.clear();
        reducedSharedSurfaceSecondPosZ.clear();
        reducedSharedSurfaceSecondMomX.clear();
        reducedSharedSurfaceSecondMomY.clear();
        reducedSharedSurfaceSecondMomZ.clear();
        reducedSharedSurfaceSecondTime.clear();

        const auto& firstTrack = tracks->at(firstTrackIndex);
        const auto& secondTrack = tracks->at(secondTrackIndex);
        const bool hasTrackCaloHitData =
          trackCaloHits != nullptr &&
          firstTrackIndex < trackCaloHits->size() &&
          secondTrackIndex < trackCaloHits->size();

        // Helper used once per selected track.  It copies the compact
        // per-track summary fields and, if available, the track-calo match data.
        const auto appendTrackReducedData =
          [&](const mu2e::TrkInfo& track,
                 const vector<mu2e::TrkSegInfo>& segments,
                 size_t selectedTrackIndex,
                 size_t foilSegmentCount,
                 const mu2e::TrkCaloHitInfo* trkcalohit)
          {
            reducedSelectedTrackIndices.push_back(static_cast<int>(selectedTrackIndex));
            reducedSelectedTrackPdg.push_back(track.pdg);
            reducedSelectedTrackNhits.push_back(track.nhits);
            reducedSelectedTrackStatus.push_back(track.status);
            reducedSelectedTrackGoodfit.push_back(track.goodfit);
            reducedSelectedTrackNdof.push_back(track.ndof);
            reducedSelectedTrackNactive.push_back(track.nactive);
            reducedSelectedTrackNseg.push_back(track.nseg);
            reducedSelectedTrackNipadown.push_back(track.nipadown);
            reducedSelectedTrackNstdown.push_back(track.nstdown);
            reducedSelectedTrackFirstStInter.push_back(track.firststinter);
            reducedSelectedTrackNstup.push_back(track.nstup);
            reducedSelectedTrackNipaup.push_back(track.nipaup);
            reducedSelectedTrackFitcon.push_back(track.fitcon);
            reducedSelectedTrackChisq.push_back(track.chisq);
            reducedSelectedTrackTotalSegmentCount.push_back(static_cast<int>(segments.size()));
            reducedSelectedTrackFoilSegmentCount.push_back(static_cast<int>(foilSegmentCount));

            if (trkcalohit != nullptr)
            {
              // Track-associated calorimeter state.  This is the reconstruction
              // summary most useful for cross-checking vertexing against calo
              // timing, impact point, and deposited energy.
              reducedSelectedTrackCaloActive.push_back(trkcalohit->active ? 1 : 0);
              reducedSelectedTrackCaloDid.push_back(trkcalohit->did);
              reducedSelectedTrackCaloPOCAX.push_back(trkcalohit->poca.x());
              reducedSelectedTrackCaloPOCAY.push_back(trkcalohit->poca.y());
              reducedSelectedTrackCaloPOCAZ.push_back(trkcalohit->poca.z());
              reducedSelectedTrackCaloMomX.push_back(trkcalohit->mom.x());
              reducedSelectedTrackCaloMomY.push_back(trkcalohit->mom.y());
              reducedSelectedTrackCaloMomZ.push_back(trkcalohit->mom.z());
              reducedSelectedTrackCaloCDepth.push_back(trkcalohit->cdepth);
              reducedSelectedTrackCaloTrkDepth.push_back(trkcalohit->trkdepth);
              reducedSelectedTrackCaloDphiDot.push_back(trkcalohit->dphidot);
              reducedSelectedTrackCaloDoca.push_back(trkcalohit->doca);
              reducedSelectedTrackCaloDt.push_back(trkcalohit->dt);
              reducedSelectedTrackCaloPtoca.push_back(trkcalohit->ptoca);
              reducedSelectedTrackCaloTocavar.push_back(trkcalohit->tocavar);
              reducedSelectedTrackCaloTresid.push_back(trkcalohit->tresid);
              reducedSelectedTrackCaloTresidmvar.push_back(trkcalohit->tresidmvar);
              reducedSelectedTrackCaloTresidpvar.push_back(trkcalohit->tresidpvar);
              reducedSelectedTrackCaloCTime.push_back(trkcalohit->ctime);
              reducedSelectedTrackCaloCTimeErr.push_back(trkcalohit->ctimeerr);
              reducedSelectedTrackCaloCSize.push_back(trkcalohit->csize);
              reducedSelectedTrackCaloEdep.push_back(trkcalohit->edep);
              reducedSelectedTrackCaloEdepErr.push_back(trkcalohit->edeperr);
            }
            else
            {
              // If the input file does not carry trkcalohit, keep the schema
              // aligned and mark the calorimeter fields as missing.
              reducedSelectedTrackCaloActive.push_back(-1);
              reducedSelectedTrackCaloDid.push_back(-1);
              reducedSelectedTrackCaloPOCAX.push_back(-1000.0);
              reducedSelectedTrackCaloPOCAY.push_back(-1000.0);
              reducedSelectedTrackCaloPOCAZ.push_back(-1000.0);
              reducedSelectedTrackCaloMomX.push_back(-1000.0);
              reducedSelectedTrackCaloMomY.push_back(-1000.0);
              reducedSelectedTrackCaloMomZ.push_back(-1000.0);
              reducedSelectedTrackCaloCDepth.push_back(-1000.0);
              reducedSelectedTrackCaloTrkDepth.push_back(-1000.0);
              reducedSelectedTrackCaloDphiDot.push_back(-1000.0);
              reducedSelectedTrackCaloDoca.push_back(-1000.0);
              reducedSelectedTrackCaloDt.push_back(-1000.0);
              reducedSelectedTrackCaloPtoca.push_back(-1000.0);
              reducedSelectedTrackCaloTocavar.push_back(-1000.0);
              reducedSelectedTrackCaloTresid.push_back(-1000.0);
              reducedSelectedTrackCaloTresidmvar.push_back(-1000.0);
              reducedSelectedTrackCaloTresidpvar.push_back(-1000.0);
              reducedSelectedTrackCaloCTime.push_back(-1000.0);
              reducedSelectedTrackCaloCTimeErr.push_back(-1000.0);
              reducedSelectedTrackCaloCSize.push_back(-1000.0);
              reducedSelectedTrackCaloEdep.push_back(-1000.0);
              reducedSelectedTrackCaloEdepErr.push_back(-1000.0);
            }
          };

        appendTrackReducedData(firstTrack,
                               firstTrackSegments,
                               firstTrackIndex,
                               foilSegmentCountsForSelectedTracks.at(0),
                               hasTrackCaloHitData ? &trackCaloHits->at(firstTrackIndex) : nullptr);
        appendTrackReducedData(secondTrack,
                               secondTrackSegments,
                               secondTrackIndex,
                               foilSegmentCountsForSelectedTracks.at(1),
                               hasTrackCaloHitData ? &trackCaloHits->at(secondTrackIndex) : nullptr);

        reducedSharedSurfaceCount = static_cast<int>(sharedSurfaces.size());
        reducedSharedSTFoilCount = static_cast<int>(sharedSTFoils.size());
        for (const auto& sharedSurface : sharedSurfaces)
        {
          const auto& firstSegment =
            firstTrackSegments.at(sharedSurface.firstStoredSegmentIndex);
          const auto& secondSegment =
            secondTrackSegments.at(sharedSurface.secondStoredSegmentIndex);

          reducedSharedSurfaceSid.push_back(sharedSurface.sid);
          reducedSharedSurfaceSindex.push_back(sharedSurface.sindex);
          reducedSharedSurfaceFirstStoredSegmentIndex.push_back(
            static_cast<int>(sharedSurface.firstStoredSegmentIndex));
          reducedSharedSurfaceSecondStoredSegmentIndex.push_back(
            static_cast<int>(sharedSurface.secondStoredSegmentIndex));
          reducedSharedSurfaceFirstPosX.push_back(firstSegment.pos.x());
          reducedSharedSurfaceFirstPosY.push_back(firstSegment.pos.y());
          reducedSharedSurfaceFirstPosZ.push_back(firstSegment.pos.z());
          reducedSharedSurfaceFirstMomX.push_back(firstSegment.mom.x());
          reducedSharedSurfaceFirstMomY.push_back(firstSegment.mom.y());
          reducedSharedSurfaceFirstMomZ.push_back(firstSegment.mom.z());
          reducedSharedSurfaceFirstTime.push_back(firstSegment.time);
          reducedSharedSurfaceSecondPosX.push_back(secondSegment.pos.x());
          reducedSharedSurfaceSecondPosY.push_back(secondSegment.pos.y());
          reducedSharedSurfaceSecondPosZ.push_back(secondSegment.pos.z());
          reducedSharedSurfaceSecondMomX.push_back(secondSegment.mom.x());
          reducedSharedSurfaceSecondMomY.push_back(secondSegment.mom.y());
          reducedSharedSurfaceSecondMomZ.push_back(secondSegment.mom.z());
          reducedSharedSurfaceSecondTime.push_back(secondSegment.time);
        }

        reducedTree->Fill();
        ++reducedEventsWritten;
      }

      if (doFullTrkSegPrintout)
      {
        //--------------------------------------------------------------------
        // Optional Event-Level Shared-Surface Printout
        //--------------------------------------------------------------------
        //
        // This block is intentionally inside doFullTrkSegPrintout.  The user
        // gets quiet summary-only output by default, but can enable this when
        // hand-checking individual events and shared surfaces.
        cout << "\nExactly-two-track shared-surface summary"
             << " | entry=" << iEntry;
        if (evtinfo != nullptr)
        {
          cout << " run=" << evtinfo->run
               << " subrun=" << evtinfo->subrun
               << " event=" << evtinfo->event;
        }
        cout << " track_indices=(" << firstTrackIndex
             << ", " << secondTrackIndex << ")"
             << " foil_trksegs=(" << foilSegmentCountsForSelectedTracks.at(0)
             << ", " << foilSegmentCountsForSelectedTracks.at(1) << ")";

        cout << " shared_surface_count=" << sharedSurfaces.size();
        if (sharedSurfaces.empty())
        {
          cout << " shared_surfaces=(none)";
        }
        else
        {
          cout << " shared_surfaces=[";
          for (size_t iSharedSurface = 0; iSharedSurface < sharedSurfaces.size(); ++iSharedSurface)
          {
            if (iSharedSurface > 0)
            {
              cout << ", ";
            }
            cout << formatSharedSurface(sharedSurfaces.at(iSharedSurface));
          }
          cout << "]";
        }

        cout << " shared_st_foils_count=" << sharedSTFoils.size();
        if (sharedSTFoils.empty())
        {
          cout << " shared_st_foils=(none)";
        }
        else
        {
          cout << " shared_st_foils=[";
          for (size_t iSharedSurface = 0; iSharedSurface < sharedSTFoils.size(); ++iSharedSurface)
          {
            if (iSharedSurface > 0)
            {
              cout << ", ";
            }
            cout << formatSharedSurface(sharedSTFoils.at(iSharedSurface));
          }
          cout << "]";
        }
        cout << endl;
      }
    }
  }

  cout << "\nTotal reconstructed downstream electron tracks found: "
       << totalDownstreamElectronTracks << endl;

  if (doFullTrkSegPrintout)
  {
    cout << "Beginning reconstructed track-segment momentum printout." << endl;

    // Separate print counter.  This should end at the same value as
    // totalDownstreamElectronTracks unless a malformed entry prevents safe printing.
    Long64_t printedDownstreamElectronTracks = 0;

    //==========================================================================
    // Pass 2: Optional Detailed Track and Segment Printout
    //==========================================================================
    //
    // Each TTree entry corresponds to one Mu2e event stored in the EventNtuple.
    // This pass is skipped in the default quiet mode.  The summary numbers are
    // already computed in Pass 1, so disabling this pass only suppresses the
    // event-by-event and track-by-track text.
    for (Long64_t iEntry = 0; iEntry < entriesToRead; ++iEntry)
    {
      // Load this event's enabled branch data into evtinfo/tracks/trackSegments.
      ntuple.GetEntry(iEntry);

      // If ROOT could not populate the requested branches for this entry, skip
      // it rather than dereferencing a null pointer.
      if (tracks == nullptr || trackSegments == nullptr)
      {
        cerr << "WARNING: null track or track-segment data at entry " << iEntry
             << "; skipping momentum printout for this entry." << endl;
        continue;
      }

      // In a well-formed EventNtuple, the outer vector sizes should match:
      //   tracks->at(i)       belongs with
      //   trackSegments->at(i)
      //
      // If they do not match, we still process the tracks that are safe to
      // access, but we print a warning because the ntuple layout is suspicious.
      const size_t nTracks = tracks->size();
      if (trackSegments->size() != nTracks)
      {
        cerr << "WARNING: entry " << iEntry << " has " << nTracks << " tracks but "
             << trackSegments->size() << " track-segment vectors." << endl;
      }

      //----------------------------------------------------------------------
      // Track Loop: Print Only Selected Downstream Electron Fits
      //----------------------------------------------------------------------
      //
      // Each element in "tracks" is one reconstructed track, usually
      // corresponding to a particular fit hypothesis.  The same reco-only
      // selection used in Pass 1 is repeated here so the detailed printout and
      // the summary counts refer to the same selected tracks.
      for (size_t iTrack = 0; iTrack < nTracks; ++iTrack)
      {
        const auto& track = tracks->at(iTrack);

        // Protect against a mismatched outer-vector size before asking for this
        // track's segment vector.  The downstream selection needs trksegs, so a
        // track with no matching segment vector cannot pass this reconstructed
        // direction requirement.
        if (iTrack >= trackSegments->size())
        {
          cout << "\nSkipping track at entry=" << iEntry
               << " track_index=" << iTrack
               << " because no " << segmentBranch
               << " vector is available for the downstream-direction test." << endl;
          continue;
        }

        // This vector contains all stored reconstructed surface intersections
        // for the current track.  Examples of surfaces include tracker
        // front/mid/back, stopping-target surfaces, and absorber surfaces,
        // depending on how the upstream reconstruction FCL sampled/extrapolated
        // the track.
        const auto& segments = trackSegments->at(iTrack);

        // Reconstruction-only downstream electron selection.
        //
        // PDG 11  = e-minus hypothesis.
        // PDG -11 = e-plus hypothesis, if requested through the electronPdg arg.
        //
        // The downstream requirement is made from the reconstructed trksegs
        // momentum at TT_Mid, not from MC truth.
        if (!isSelectedDownstreamElectronTrack(track, segments, electronPdg))
        {
          continue;
        }

        ++printedDownstreamElectronTracks;

        cout << "\nDownstream electron track " << printedDownstreamElectronTracks
             << " | entry=" << iEntry;
        if (evtinfo != nullptr)
        {
          cout << " run=" << evtinfo->run
               << " subrun=" << evtinfo->subrun
               << " event=" << evtinfo->event;
        }
        cout << " track_index=" << iTrack
             << " nhits=" << track.nhits
             << " nactive=" << track.nactive
             << " fitcon=" << track.fitcon
             << endl;

        if (segments.empty())
        {
          cout << "  This downstream electron track has no stored "
               << segmentBranch << " entries." << endl;
          continue;
        }

        //--------------------------------------------------------------------
        // Segment Loop: Chronological Reconstructed Track-State Printout
        //--------------------------------------------------------------------
        //
        // The momentum is printed both as a vector (px, py, pz) and as a scalar
        // magnitude |p|.  The surface ID and index are included so you can later
        // select specific locations, for example:
        //   sid == 0    -> TT_Front
        //   sid == 1    -> TT_Mid
        //   sid == 2    -> TT_Back
        //   sid == 104  -> ST_Foils
        //
        // EventNtuple does not store these segment vectors in chronological
        // order.  Sort by reconstructed segment time for printing, while
        // preserving the original stored segment index in the "stored_seg_index"
        // field.
        const vector<size_t> sortedSegmentIndices = sortedSegmentIndicesByTime(segments);
        for (size_t iTimeOrder = 0;
             iTimeOrder < sortedSegmentIndices.size();
             ++iTimeOrder)
        {
          const size_t storedSegmentIndex = sortedSegmentIndices.at(iTimeOrder);
          const auto& segment = segments.at(storedSegmentIndex);
          cout << "  time_order=" << setw(3) << iTimeOrder
               << " stored_seg_index=" << setw(3) << storedSegmentIndex
               << " sid=" << setw(4) << segment.sid
               << " sindex=" << setw(3) << segment.sindex
               << " |p|=" << setw(9) << segment.mom.R() << " MeV/c"
               << " px=" << setw(9) << segment.mom.x()
               << " py=" << setw(9) << segment.mom.y()
               << " pz=" << setw(9) << segment.mom.z()
               << " time=" << setw(9) << segment.time << " ns"
               << " pos=(" << segment.pos.x()
               << ", " << segment.pos.y()
               << ", " << segment.pos.z() << ") mm"
               << " inbounds=" << segment.inbounds
               << " gap=" << segment.gap
               << endl;
        }
      }
    }

    // Final cross-check.  This repeats the count after printing so the end of the
    // log tells you whether all counted downstream electron tracks were
    // successfully printed.
    cout << "\nPrinted reconstructed downstream electron tracks: "
         << printedDownstreamElectronTracks
         << " of " << totalDownstreamElectronTracks << endl;
  }
  else
  {
    cout << "Full reconstructed track-segment momentum printout disabled; "
         << "summary counts and percentages will still be printed." << endl;
  }

  //============================================================================
  // Summary Calculations: Convert Counts into Percentages and Averages
  //============================================================================

  const double percentOfScannedEvents =
    entriesToRead > 0
      ? 100.0 * static_cast<double>(eventsWithExactlyTwoDownstreamElectronTracks) /
          static_cast<double>(entriesToRead)
      : 0.0;
  const double percentOfThrownEvents =
    100.0 * static_cast<double>(eventsWithExactlyTwoDownstreamElectronTracks) /
      nThrownEvents;
  const double percentWithSharedSurfacesOfScannedEvents =
    entriesToRead > 0
      ? 100.0 * static_cast<double>(exactlyTwoTrackEventsWithSharedSurfaces) /
          static_cast<double>(entriesToRead)
      : 0.0;
  const double percentWithSharedSurfacesOfThrownEvents =
    100.0 * static_cast<double>(exactlyTwoTrackEventsWithSharedSurfaces) /
      nThrownEvents;
  const double averageSharedSurfacesPerExactlyTwoTrackEvent =
    eventsWithExactlyTwoDownstreamElectronTracks > 0
      ? static_cast<double>(totalSharedSurfacesAcrossExactlyTwoTrackEvents) /
          static_cast<double>(eventsWithExactlyTwoDownstreamElectronTracks)
      : 0.0;

  //============================================================================
  // Final Summary Printout
  //============================================================================

  cout << "\nEvents with more than one reconstructed downstream electron track: "
       << eventsWithMoreThanOneDownstreamElectronTrack << endl;
  cout << "Events with exactly two reconstructed downstream electron tracks: "
       << eventsWithExactlyTwoDownstreamElectronTracks << endl;

  cout << "Electron tracks in exactly-two-track events: "
       << 2 * eventsWithExactlyTwoDownstreamElectronTracks << endl;
  cout << "ST_Foils trkseg multiplicity for electron tracks in exactly-two-track events:"
       << endl;
  if (foilSegmentMultiplicityForExactlyTwoTrackEvents.empty())
  {
    cout << "  No electron tracks from exactly-two-track events were found." << endl;
  }
  else
  {
    const size_t maxFoilSegments =
      foilSegmentMultiplicityForExactlyTwoTrackEvents.rbegin()->first;
    for (size_t nFoilSegments = 0; nFoilSegments <= maxFoilSegments; ++nFoilSegments)
    {
      const auto countIter =
        foilSegmentMultiplicityForExactlyTwoTrackEvents.find(nFoilSegments);
      const Long64_t nTracksWithThisMultiplicity =
        countIter != foilSegmentMultiplicityForExactlyTwoTrackEvents.end()
          ? countIter->second
          : 0;

      cout << "  " << setw(3) << nFoilSegments
           << " ST_Foils trksegs: "
           << nTracksWithThisMultiplicity
           << " electron tracks" << endl;
    }
  }

  cout << "Shared-surface diagnostic:" << endl;
  cout << "  Exactly-two-track events with at least one shared surface: "
       << exactlyTwoTrackEventsWithSharedSurfaces << endl;
  cout << "  Total shared surfaces across exactly-two-track events: "
       << totalSharedSurfacesAcrossExactlyTwoTrackEvents << endl;
  cout << "  Average shared surfaces per exactly-two-track event: "
       << averageSharedSurfacesPerExactlyTwoTrackEvent << endl;
  cout << "  Events with at least one shared surface: "
       << exactlyTwoTrackEventsWithSharedSurfaces
       << " events | " << percentWithSharedSurfacesOfScannedEvents
       << "% of scanned ntuple events | "
       << percentWithSharedSurfacesOfThrownEvents
       << "% of thrown events" << endl;
  cout << "  Exactly-two-track events with at least one shared ST_Foils surface: "
       << exactlyTwoTrackEventsWithSharedSTFoils << endl;
  cout << "  Total shared ST_Foils surfaces across exactly-two-track events: "
       << totalSharedSTFoilsAcrossExactlyTwoTrackEvents << endl;

  cout << "Percent of scanned ntuple events with exactly two reconstructed downstream electron tracks ("
       << entriesToRead << "): "
       << percentOfScannedEvents << "%" << endl;
  cout << "Percent of thrown events with exactly two reconstructed downstream electron tracks ("
       << static_cast<Long64_t>(nThrownEvents) << "): "
       << percentOfThrownEvents << "%" << endl;

  if (writeReducedOutput && reducedFile != nullptr && reducedTree != nullptr)
  {
    reducedFile->cd();
    reducedTree->Write();
    reducedFile->Write();
    reducedFile->Close();
    cout << "Reduced vertex-candidate ntuple written with " << reducedEventsWritten
         << " events to " << reducedOutputName << endl;
  }
}

