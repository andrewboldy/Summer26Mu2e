//----------------------------------------------------------------------------------
//
// TrkSegInvestigator.C
// Written by Andrew Boldy University of South Carolina, 2026
// Assisted by Codex
//
// Purpose:
//   Loop over an EventNtuple ROOT file or filelist with RooUtil and inspect the
//   reconstructed trksegs branch.
//
//   In the common EventNtuple branch layout:
//
//       trk[i_track]          -> reconstructed track summary
//       trksegs[i_track]      -> all stored surface intersections for that track
//       trksegs[i_track][j]   -> one fitted track state at one surface
//
//   This first version is deliberately printout-focused.  The hard-coded manual
//   switches below let you turn pieces of the trkseg printout on and off while
//   keeping the event/track/segment loop unchanged.
//
// Usage from ROOT:
//   .L CreatedCode/TrkSegInvestigator.C+
//   TrkSegInvestigator("path/to/nts.root")
//   TrkSegInvestigator("filelist.txt", 25)   // print only the first 25 entries
//
// Notes:
//   - This macro leaves all EventNtuple branches enabled for RooUtil.
//     RooUtil::GetEvent() calls Event::Update(), which can touch companion
//     branches while it builds full Track wrappers.
//   - The printout below only uses reconstructed trk/trksegs information.
//   - EventNtuple stores trksegs in the KalSeed intersection order.  This macro
//     sorts the printout by reconstructed trkseg time so the path is easier to
//     follow by eye, while still showing the original stored segment index.
//
//----------------------------------------------------------------------------------

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "EventNtuple/inc/TrkInfo.hh"
#include "EventNtuple/inc/TrkSegInfo.hh"
#include "EventNtuple/rooutil/inc/RooUtil.hh"

using namespace std;
using namespace rooutil;

namespace
{
  //============================================================================
  // Manual Switches
  //============================================================================

  // Master switch for all text printout.  The detailed PRINT_* switches below
  // are only checked when DO_PRINTOUT is true.
  const bool DO_PRINTOUT = true;

  // Individual trkseg printout switches.  These all start as true so the first
  // run prints the complete reconstructed trkseg payload.  Manually change these
  // to false when you want a quieter output log.
  const bool PRINT_MOM_MOMERR_DMOM = true;
  const bool PRINT_POS             = true;
  const bool PRINT_TIME            = true;
  const bool PRINT_INBOUNDS        = true;
  const bool PRINT_GAP             = true;
  const bool PRINT_SID_SINDEX      = true;

  // Placeholder for the next stage.  Leave this false until we add specific
  // histogram definitions and fill logic.
  const bool DO_HISTOGRAMS = false;

  //============================================================================
  // Small Formatting Helpers
  //============================================================================

  string surfaceName(int sid)
  {
    switch (sid)
    {
      case mu2e::SurfaceIdDetail::TT_Front: return "TT_Front";
      case mu2e::SurfaceIdDetail::TT_Mid:   return "TT_Mid";
      case mu2e::SurfaceIdDetail::TT_Back:  return "TT_Back";
      case mu2e::SurfaceIdDetail::TT_Inner: return "TT_Inner";
      case mu2e::SurfaceIdDetail::TT_Outer: return "TT_Outer";
      case mu2e::SurfaceIdDetail::DS_Front: return "DS_Front";
      case mu2e::SurfaceIdDetail::DS_Back:  return "DS_Back";
      case mu2e::SurfaceIdDetail::DS_Inner: return "DS_Inner";
      case mu2e::SurfaceIdDetail::DS_Outer: return "DS_Outer";
      case mu2e::SurfaceIdDetail::IPA:      return "IPA";
      case mu2e::SurfaceIdDetail::IPA_Front:return "IPA_Front";
      case mu2e::SurfaceIdDetail::IPA_Back: return "IPA_Back";
      case mu2e::SurfaceIdDetail::OPA:      return "OPA";
      case mu2e::SurfaceIdDetail::TSDA:     return "TSDA";
      case mu2e::SurfaceIdDetail::ST_Front: return "ST_Front";
      case mu2e::SurfaceIdDetail::ST_Back:  return "ST_Back";
      case mu2e::SurfaceIdDetail::ST_Inner: return "ST_Inner";
      case mu2e::SurfaceIdDetail::ST_Outer: return "ST_Outer";
      case mu2e::SurfaceIdDetail::ST_Foils: return "ST_Foils";
      case mu2e::SurfaceIdDetail::ST_Wires: return "ST_Wires";
      case mu2e::SurfaceIdDetail::TCRV:     return "TCRV";
      default:                              return "unknown";
    }
  }

  void printEventSeparator()
  {
    cout << "----------------------------------------------------------------------" << endl;
  }

  void printTrackHeader(size_t iTrack,
                        const mu2e::TrkInfo* track,
                        size_t nSegments)
  {
    cout << "Track " << iTrack
         << ": n_trksegs = " << nSegments;

    if (track != nullptr)
    {
      cout << ", pdg = " << track->pdg
           << ", status = " << track->status
           << ", goodfit = " << track->goodfit
           << ", nactive = " << track->nactive;
    }

    cout << endl;
  }

  bool selectionCriteria(const mu2e::TrkInfo* track,
                         const vector<mu2e::TrkSegInfo>& trackSegments)
  {
    // Reconstructed downstream-electron selection.
    //
    // This intentionally uses only reconstructed information:
    //   1. The fitted track hypothesis must be an electron: trk.pdg == 11.
    //   2. The track must have a tracker-middle segment moving downstream:
    //      trkseg.sid == TT_Mid and trkseg.mom.z() > 0.
    //
    // MC-truth branches such as trkmcsim/trkmc/trksegsmc are not used here.
    if (track == nullptr)
    {
      return false;
    }

    if (track->pdg != 11)
    {
      return false;
    }

    for (const auto& segment : trackSegments)
    {
      if (segment.sid == mu2e::SurfaceIdDetail::TT_Mid &&
          segment.mom.z() > 0.0)
      {
        return true;
      }
    }

    return false;
  }

  vector<size_t> sortedSegmentIndicesByTime(const vector<mu2e::TrkSegInfo>& trackSegments)
  {
    // EventNtuple does not guarantee that trksegs are stored in time order.
    // Keep indices rather than copying segments so we can still report the
    // original stored index from trksegs[i_track][storedIndex].
    vector<size_t> sortedIndices;
    sortedIndices.reserve(trackSegments.size());

    for (size_t iSegment = 0; iSegment < trackSegments.size(); ++iSegment)
    {
      sortedIndices.push_back(iSegment);
    }

    stable_sort(sortedIndices.begin(),
                sortedIndices.end(),
                [&trackSegments](size_t lhs, size_t rhs)
                {
                  return trackSegments.at(lhs).time < trackSegments.at(rhs).time;
                });

    return sortedIndices;
  }

  void printSegment(const mu2e::TrkSegInfo& segment,
                    size_t sortedIndex,
                    size_t storedIndex)
  {
    cout << "  trkseg " << setw(3) << sortedIndex
         << " | stored_index = " << storedIndex;

    if (PRINT_MOM_MOMERR_DMOM)
    {
      cout << " | px = " << segment.mom.x() << " MeV/c"
           << ", py = " << segment.mom.y() << " MeV/c"
           << ", pz = " << segment.mom.z() << " MeV/c"
           << ", |p| = " << segment.mom.R() << " MeV/c"
           << ", momerr = " << segment.momerr << " MeV/c"
           << ", dmom = " << segment.dmom << " MeV/c";
    }

    if (PRINT_POS)
    {
      cout << " | pos = ("
           << segment.pos.x() << ", "
           << segment.pos.y() << ", "
           << segment.pos.z() << ") mm";
    }

    if (PRINT_TIME)
    {
      cout << " | time = " << segment.time << " ns";
    }

    if (PRINT_INBOUNDS)
    {
      cout << " | inbounds = " << boolalpha << segment.inbounds << noboolalpha;
    }

    if (PRINT_GAP)
    {
      cout << " | gap = " << boolalpha << segment.gap << noboolalpha;
    }

    if (PRINT_SID_SINDEX)
    {
      cout << " | sid = " << segment.sid
           << " (" << surfaceName(segment.sid) << ")"
           << ", sindex = " << segment.sindex;
    }

    cout << endl;
  }

  void runHistogramPlaceholder()
  {
    // This is intentionally empty for now.  The DO_HISTOGRAMS switch is already
    // in place so we can add one histogram at a time after the printout is
    // confirmed to be reading the intended trkseg content.
  }
}

void TrkSegInvestigator(const string& inputName, int maxEvents = -1)
{
  //----------------------------------------------------------------------------
  // Open EventNtuple With RooUtil
  //----------------------------------------------------------------------------

  RooUtil util(inputName);
  const int nEntries = util.GetNEvents();

  if (nEntries <= 0)
  {
    cerr << "ERROR: no entries found in EventNtuple/ntuple for input: "
         << inputName << endl;
    return;
  }

  // Leave all discovered branches enabled.  RooUtil::GetEvent() calls
  // Event::Update(), and Event::Update() can touch companion branches while it
  // builds the convenient Track wrappers.

  //ternary operator to determine the number of entries to read
  const int entriesToRead = (maxEvents >= 0) ? min(maxEvents, nEntries) : nEntries;

  cout << "Input: " << inputName << endl;
  cout << "Tree entries available: " << nEntries << endl;
  cout << "Tree entries being scanned: " << entriesToRead << endl;
  cout << "RooUtil branch mode: all discovered branches left enabled." << endl;
  cout << "DO_PRINTOUT = " << boolalpha << DO_PRINTOUT << endl;
  if (DO_PRINTOUT)
  {
    cout << "  PRINT_MOM_MOMERR_DMOM = " << PRINT_MOM_MOMERR_DMOM << endl;
    cout << "  PRINT_POS             = " << PRINT_POS << endl;
    cout << "  PRINT_TIME            = " << PRINT_TIME << endl;
    cout << "  PRINT_INBOUNDS        = " << PRINT_INBOUNDS << endl;
    cout << "  PRINT_GAP             = " << PRINT_GAP << endl;
    cout << "  PRINT_SID_SINDEX      = " << PRINT_SID_SINDEX << endl;
  }
  cout << "DO_HISTOGRAMS = " << DO_HISTOGRAMS << noboolalpha << endl;

  // Counters are useful both for sanity checks and for deciding what histogram
  // to add next.
  size_t nEventsRead = 0;
  size_t nEventsWithSelectedTracks = 0;
  size_t nTracksExamined = 0;
  size_t nTracksSelected = 0;
  size_t nSegmentsSelected = 0;

  //----------------------------------------------------------------------------
  // Main Event Loop
  //----------------------------------------------------------------------------

  for (int iEntry = 0; iEntry < entriesToRead; ++iEntry)
  {
    const auto& event = util.GetEvent(iEntry);
    ++nEventsRead;

    if (event.trksegs == nullptr)
    {
      cerr << "ERROR: event.trksegs is null.  This file may not contain the "
           << "standard trksegs branch expected by RooUtil." << endl;
      return;
    }

    bool eventHasSelectedTrack = false;
    bool printedEventHeader = false;

    for (size_t iTrack = 0; iTrack < event.trksegs->size(); ++iTrack)
    {
      const auto& trackSegments = event.trksegs->at(iTrack);
      const mu2e::TrkInfo* track = nullptr;

      if (event.trk != nullptr && iTrack < event.trk->size())
      {
        track = &(event.trk->at(iTrack));
      }

      ++nTracksExamined;

      if (!selectionCriteria(track, trackSegments))
      {
        continue;
      }

      if (!eventHasSelectedTrack)
      {
        ++nEventsWithSelectedTracks;
        eventHasSelectedTrack = true;
      }

      ++nTracksSelected;
      nSegmentsSelected += trackSegments.size();

      if (DO_PRINTOUT)
      {
        if (!printedEventHeader)
        {
          printEventSeparator();
          cout << "Entry " << iEntry;

          if (event.evtinfo != nullptr)
          {
            cout << " | run = " << event.evtinfo->run
                 << ", subrun = " << event.evtinfo->subrun
                 << ", event = " << event.evtinfo->event;
          }

          cout << " | n_tracks = " << event.trksegs->size() << endl;
          printedEventHeader = true;
        }

        printTrackHeader(iTrack, track, trackSegments.size());

        const vector<size_t> sortedSegmentIndices =
          sortedSegmentIndicesByTime(trackSegments);

        for (size_t iSortedSegment = 0;
             iSortedSegment < sortedSegmentIndices.size();
             ++iSortedSegment)
        {
          const size_t storedSegmentIndex = sortedSegmentIndices.at(iSortedSegment);
          printSegment(trackSegments.at(storedSegmentIndex),
                       iSortedSegment,
                       storedSegmentIndex);
        }
      }

      if (DO_HISTOGRAMS)
      {
        runHistogramPlaceholder();
      }
    }
  }

  if (DO_PRINTOUT)
  {
    printEventSeparator();
  }

  cout << "Summary:" << endl;
  cout << "  events read: " << nEventsRead << endl;
  cout << "  events with selected tracks: " << nEventsWithSelectedTracks << endl;
  cout << "  tracks examined: " << nTracksExamined << endl;
  cout << "  tracks selected: " << nTracksSelected << endl;
  cout << "  selected-track trksegs read: " << nSegmentsSelected << endl;
} // end of TrkSegInvestigator function

