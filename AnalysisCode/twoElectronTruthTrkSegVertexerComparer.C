//----------------------------------------------------------------------------------
//
// twoElectronTruthTrkSegVertexerComparer.C
// Written by Andrew Boldy University of South Carolina, 2026
// Assisted by Codex
//
// Purpose:
//   Read an EventNtuple `nts.*.root` file directly and compare the first useful
//   pieces of MC truth information against reconstructed track/calo information.
//
//   This macro is the new ntuple-based replacement for the art-module attempt.
//   It does not call `mu2e`, does not use art::RootInput, and does not require
//   an original art event file.  It reads the ROOT tree:
//
//       EventNtuple/ntuple
//
//   from an `nts.*.root` file or from a text filelist of such files.
//
// First diagnostic pass:
//   For each scanned event, print:
//     - event identity, when evtinfo is available
//     - the MC truth `trkmcsim` origin for each stored truth particle:
//         origin = (t, x, y, z)
//     - whether reconstructed tracks exist in the requested track branch
//     - whether each reconstructed track has an associated calorimeter object
//     - whether each reconstructed momentum lies between 50 and 53 MeV/c
//     - whether the event contains two reconstructed tracks with those same
//       requested properties
//
// Histogram pass:
//   The macro also writes a ROOT file containing the first comparison
//   histograms.  Histogram entries are filled only for events with exactly two
//   selected reconstructed downstream electron tracks.  In the default
//   selection, each selected reconstructed track must:
//
//       trk.pdg == 11
//       have TT_Mid pz > 0
//       have an associated trkcalohit
//       have reconstructed momentum in 50-53 MeV/c
//
//   The selected reconstructed and truth momentum histograms are still drawn
//   over 30-55 MeV/c.  That plotting range is intentionally wider than the
//   selection cut.
//
//   The macro also checks the reconstructed momentum z signs on the actual
//   shared ST_Foils segment pair selected for the vertex:
//
//       one selected segment has mom.z() > 0
//       the other selected segment has mom.z() < 0
//
//   This is deliberately less restrictive than classifying a full track by all
//   of its stored foil crossings.
//
//   Within those selected events, the reconstructed vertex is built from the
//   shared ST_Foils surface pair with the smallest closest-line distance.  A
//   parallel TEST histogram set is also filled for the shared ST_Foils surface
//   pair with the smallest absolute time difference.  MC truth histograms are
//   filled only from trkmcsim particles that are rank 0 downstream electrons:
//   The summary also compares each constructed vertex z position with the
//   nearest configured stopping-target foil center and reports whether that
//   geometry-nearest foil index matches the shared ST_Foils sindex used to
//   build the vertex.
//
//       sim.valid
//       sim.rank == 0
//       sim.pdg == 11
//       sim.mom.z() > 0
//
//   The output file name is controlled by HISTOGRAM_OUTPUT_FILE near the top of
//   the implementation.  The PDF plots are written under
//   Plots/TruthVsRecoPlots/.
//
// Branch conventions:
//   The default branch prefix is "trk", which uses:
//
//       trk        -> vector<mu2e::TrkInfo>
//       trksegs    -> vector<vector<mu2e::TrkSegInfo>>
//       trkmcsim   -> vector<vector<mu2e::SimInfo>>
//       trkcalohit -> vector<mu2e::TrkCaloHitInfo>
//
//   If a different trackBranch is supplied, branch names are derived as:
//
//       <trackBranch>segs
//       <trackBranch>mcsim
//       <trackBranch>calohit
//
//   The derived names can be overridden by passing explicit simBranch and
//   caloBranch values.
//
// Momentum convention for this first pass:
//   The reconstructed momentum test uses the best diagnostic momentum available
//   for each track:
//
//     1. trkcalohit[i].mom.R(), when an active calo association exists
//     2. the reconstructed TT_Mid segment momentum, when present
//     3. the first stored trkseg momentum, as a last diagnostic fallback
//
//   The printout labels which source was used.  This keeps the first pass
//   robust across ntuples while making the choice explicit in the terminal log.
//
// Usage from ROOT:
//   .L CreatedCode/HistogramMakers/twoElectronTruthTrkSegVertexerComparer.C+
//   twoElectronTruthTrkSegVertexerComparer("path/to/nts.root")
//
// Optional examples:
//   twoElectronTruthTrkSegVertexerComparer("filelist.txt", 25)
//   twoElectronTruthTrkSegVertexerComparer("nts.root", -1, "trk", 50.0, 53.0)
//   twoElectronTruthTrkSegVertexerComparer("nts.root", -1, "trk", 50.0, 53.0,
//                                          11, true, 3)
//
//   The fourth and fifth arguments are the reconstructed momentum selection cut.
//
//   The sixth argument is the required reconstructed-track PDG hypothesis.
//   Use requiredRecoPdg = 0 to disable that requirement.
//
//   The seventh argument requires downstream TT_Mid pz > 0 when true.
//
//   The eighth argument limits how many `trkmcsim` entries are printed per
//   reconstructed track.  Use -1 to print all.
//
//----------------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <TChain.h>

#include "EventNtuple/inc/EventInfo.hh"
#include "EventNtuple/inc/SimInfo.hh"
#include "EventNtuple/inc/TrkCaloHitInfo.hh"
#include "EventNtuple/inc/TrkInfo.hh"
#include "EventNtuple/inc/TrkSegInfo.hh"
#include "Helpers/twoElectronSelectedParticleHistograms.hh"
#include "Helpers/twoElectronSelectedParticlePlots.hh"
#include "Helpers/twoParticleVertexer.hh"

using namespace std;

namespace
{
  //============================================================================
  // Manual Printout Switch
  //============================================================================

  // Set this to false by hand when you want the macro to scan the input and
  // print only the final summary counters.  This avoids having to pass a long
  // optional argument list from the ROOT command line while the macro is still
  // in its first diagnostic-printout stage.
  //
  // Example manual edit:
  //
  //   static const bool DO_FULL_PRINT = false;
  //
  static const bool DO_FULL_PRINT = false;

  // This prints the comparison between the selected reconstructed vertex and
  // the representative rank-0 downstream-electron MC truth origin.
  static const bool DO_VERTEX_ORIGIN_RESIDUAL_PRINT = true;

  // Histogram output is controlled here for the same reason as DO_FULL_PRINT:
  // the macro is still evolving, and editing one clearly named constant is less
  // error-prone than passing a long positional argument list from ROOT.
  static const bool WRITE_HISTOGRAM_FILE = true;
  static const char* HISTOGRAM_OUTPUT_FILE =
    "twoElectronTruthTrkSegVertexerComparer_histograms.root";
  static const char* HISTOGRAM_PLOT_OUTPUT_DIR = "Plots/TruthVsRecoPlots";

  // Current 37-foil stopping-target geometry used by the recent run-1/phase-1
  // geometry files:
  //   Offline/Mu2eG4/geom/stoppingTargetHoles_v02.txt
  //   Offline/Mu2eG4/geom/stoppingTargetHoles_DOE_review_2017.txt
  //
  // EventNtuple TrkSegInfo positions are in detector coordinates, so the
  // Mu2e-coordinate foil centers are shifted by the detector-system z origin.
  static const int CURRENT_ST_GEOMETRY_N_FOILS = 37;
  static const double CURRENT_ST_GEOMETRY_Z0_IN_MU2E = 5871.0;
  static const double CURRENT_ST_GEOMETRY_DELTA_Z = 22.222222;
  static const double CURRENT_ST_GEOMETRY_DETECTOR_SYSTEM_Z0_IN_MU2E = 10171.0;

  //============================================================================
  // Input Helpers
  //============================================================================

  bool hasRootSuffix(const string& path)
  {
    const string suffix = ".root";
    if (path.size() < suffix.size())
    {
      return false;
    }
    return path.substr(path.size() - suffix.size()) == suffix;
  }

  void enableBranch(TChain& chain, const string& branchName)
  {
    // EventNtuple branches are usually split.  Enabling both the top-level
    // branch and its split leaves avoids missing fields such as trk.pdg or
    // trkcalohit.mom.fCoordinates.fX.
    chain.SetBranchStatus(branchName.c_str(), 1);
    chain.SetBranchStatus((branchName + ".*").c_str(), 1);
  }

  bool addInput(TChain& chain, const string& inputName)
  {
    if (hasRootSuffix(inputName))
    {
      chain.Add(inputName.c_str());
      return true;
    }

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

  string deriveSegmentBranchName(const string& trackBranch)
  {
    return trackBranch + "segs";
  }

  string deriveSimBranchName(const string& trackBranch)
  {
    return trackBranch + "mcsim";
  }

  string deriveCaloBranchName(const string& trackBranch)
  {
    return trackBranch + "calohit";
  }

  string formatEventLabel(const mu2e::EventInfo* evtinfo, Long64_t entry)
  {
    ostringstream out;
    if (evtinfo != nullptr)
    {
      out << "run=" << evtinfo->run
          << " subrun=" << evtinfo->subrun
          << " event=" << evtinfo->event
          << " entry=" << entry;
    }
    else
    {
      out << "entry=" << entry;
    }
    return out.str();
  }

  string formatVector3(const XYZVectorF& vector)
  {
    ostringstream out;
    out << "(" << vector.x() << ", " << vector.y() << ", " << vector.z() << ")";
    return out.str();
  }

  struct GeometryFoilMatch
  {
    bool valid = false;
    int foilIndex = -1;
    double foilCenterZ = numeric_limits<double>::quiet_NaN();
    double deltaZ = numeric_limits<double>::quiet_NaN();
  };

  struct VertexFoilIndexCheck
  {
    bool valid = false;
    bool matchesSharedFoilIndex = false;
    int vertexNearestFoilIndex = -1;
    int sharedFoilIndex = -1;
    double vertexZ = numeric_limits<double>::quiet_NaN();
    double nearestFoilCenterZ = numeric_limits<double>::quiet_NaN();
    double vertexMinusFoilCenterZ = numeric_limits<double>::quiet_NaN();
    string failureReason;
  };

  const vector<double>& configuredStoppingTargetFoilCentersZ()
  {
    static const vector<double> centers = []() {
      vector<double> result;
      result.reserve(CURRENT_ST_GEOMETRY_N_FOILS);

      const int n0 = CURRENT_ST_GEOMETRY_N_FOILS / 2;
      const double offset =
        (CURRENT_ST_GEOMETRY_N_FOILS % 2 == 1)
          ? CURRENT_ST_GEOMETRY_Z0_IN_MU2E
          : CURRENT_ST_GEOMETRY_Z0_IN_MU2E +
              0.5 * CURRENT_ST_GEOMETRY_DELTA_Z;

      for (int foilIndex = 0;
           foilIndex < CURRENT_ST_GEOMETRY_N_FOILS;
           ++foilIndex)
      {
        const double foilCenterZInMu2e =
          offset +
          static_cast<double>(foilIndex - n0) * CURRENT_ST_GEOMETRY_DELTA_Z;
        result.push_back(
          foilCenterZInMu2e -
          CURRENT_ST_GEOMETRY_DETECTOR_SYSTEM_Z0_IN_MU2E);
      }

      return result;
    }();

    return centers;
  }

  GeometryFoilMatch nearestConfiguredStoppingTargetFoil(double detectorZ)
  {
    GeometryFoilMatch match;
    if (!std::isfinite(detectorZ))
    {
      return match;
    }

    const vector<double>& foilCentersZ = configuredStoppingTargetFoilCentersZ();
    double bestAbsDeltaZ = numeric_limits<double>::infinity();
    for (size_t foilIndex = 0; foilIndex < foilCentersZ.size(); ++foilIndex)
    {
      const double deltaZ = detectorZ - foilCentersZ.at(foilIndex);
      const double absDeltaZ = fabs(deltaZ);
      if (absDeltaZ < bestAbsDeltaZ)
      {
        bestAbsDeltaZ = absDeltaZ;
        match.valid = true;
        match.foilIndex = static_cast<int>(foilIndex);
        match.foilCenterZ = foilCentersZ.at(foilIndex);
        match.deltaZ = deltaZ;
      }
    }

    return match;
  }

  VertexFoilIndexCheck compareVertexToSharedFoilIndex(
    const twoparticlevertexer::VertexResult& vertex,
    const mu2e::TrkSegInfo* firstSegment,
    const mu2e::TrkSegInfo* secondSegment)
  {
    VertexFoilIndexCheck check;
    if (!vertex.valid)
    {
      check.failureReason = "vertex is invalid";
      return check;
    }

    if (firstSegment == nullptr || secondSegment == nullptr)
    {
      check.failureReason = "missing shared ST_Foils segment pointer";
      return check;
    }

    if (firstSegment->sid != mu2e::SurfaceIdDetail::ST_Foils ||
        secondSegment->sid != mu2e::SurfaceIdDetail::ST_Foils)
    {
      check.failureReason = "vertex seed segments are not both ST_Foils";
      return check;
    }

    if (firstSegment->sindex != secondSegment->sindex)
    {
      check.failureReason = "vertex seed segments do not share the same foil sindex";
      return check;
    }

    const GeometryFoilMatch nearestFoil =
      nearestConfiguredStoppingTargetFoil(vertex.vertex.z());
    if (!nearestFoil.valid)
    {
      check.failureReason = "could not match vertex z to configured foil geometry";
      return check;
    }

    check.valid = true;
    check.sharedFoilIndex = firstSegment->sindex;
    check.vertexNearestFoilIndex = nearestFoil.foilIndex;
    check.vertexZ = vertex.vertex.z();
    check.nearestFoilCenterZ = nearestFoil.foilCenterZ;
    check.vertexMinusFoilCenterZ = nearestFoil.deltaZ;
    check.matchesSharedFoilIndex =
      check.vertexNearestFoilIndex == check.sharedFoilIndex;
    return check;
  }

  void countVertexFoilIndexCheck(const VertexFoilIndexCheck& check,
                                 Long64_t& checkedVertices,
                                 Long64_t& matchingVertices)
  {
    if (!check.valid)
    {
      return;
    }

    ++checkedVertices;
    if (check.matchesSharedFoilIndex)
    {
      ++matchingVertices;
    }
  }

  void printVertexFoilIndexCheck(const string& label,
                                 const VertexFoilIndexCheck& check)
  {
    cout << "    " << label << " nearest configured foil check: ";
    if (!check.valid)
    {
      cout << "unavailable (" << check.failureReason << ")" << endl;
      return;
    }

    cout << "vertex_z=" << check.vertexZ << " mm"
         << " nearest_geometry_foil_sindex=" << check.vertexNearestFoilIndex
         << " foil_center_z=" << check.nearestFoilCenterZ << " mm"
         << " vertex_minus_foil_center_z=" << check.vertexMinusFoilCenterZ
         << " mm"
         << " shared_seed_sindex=" << check.sharedFoilIndex
         << " match="
         << (check.matchesSharedFoilIndex ? "yes" : "no")
         << endl;
  }

  void printVertexFoilIndexSummary(const string& label,
                                   Long64_t matchingVertices,
                                   Long64_t checkedVertices)
  {
    cout << "  " << label << ": "
         << matchingVertices << " / " << checkedVertices;
    if (checkedVertices > 0)
    {
      cout << " ("
           << 100.0 * static_cast<double>(matchingVertices) /
                static_cast<double>(checkedVertices)
           << "%)";
    }
    cout << endl;
  }

  struct SharedSurfaceMatch
  {
    int sid = -1;
    int sindex = -1;
    size_t firstStoredSegmentIndex = 0;
    size_t secondStoredSegmentIndex = 0;
  };

  using SharedSurfaceMatches = vector<SharedSurfaceMatch>;

  enum class SharedFoilDeltaZScatterCategory
  {
    kAllCandidatePairs,
    kSpaceSelectedPair,
    kMinTimeSelectedPair
  };

  SharedSurfaceMatches sharedSurfacesBetweenTracks(
    const vector<mu2e::TrkSegInfo>& firstTrackSegments,
    const vector<mu2e::TrkSegInfo>& secondTrackSegments)
  {
    // Keep every stored segment for a surface key.  A track can intersect the
    // same ST_Foils sindex more than once, and the vertexer must be able to
    // choose the smallest-distance pair among those repeated crossings.
    map<pair<int, int>, vector<size_t>> firstSurfaceToSegmentIndices;
    map<pair<int, int>, vector<size_t>> secondSurfaceToSegmentIndices;

    for (size_t iFirstSegment = 0; iFirstSegment < firstTrackSegments.size(); ++iFirstSegment)
    {
      const auto& firstSegment = firstTrackSegments.at(iFirstSegment);
      const auto key = make_pair(firstSegment.sid, firstSegment.sindex);
      firstSurfaceToSegmentIndices[key].push_back(iFirstSegment);
    }

    for (size_t iSecondSegment = 0; iSecondSegment < secondTrackSegments.size(); ++iSecondSegment)
    {
      const auto& secondSegment = secondTrackSegments.at(iSecondSegment);
      const auto key = make_pair(secondSegment.sid, secondSegment.sindex);
      secondSurfaceToSegmentIndices[key].push_back(iSecondSegment);
    }

    SharedSurfaceMatches sharedSurfaces;
    for (const auto& firstEntry : firstSurfaceToSegmentIndices)
    {
      const auto secondIter = secondSurfaceToSegmentIndices.find(firstEntry.first);
      if (secondIter == secondSurfaceToSegmentIndices.end())
      {
        continue;
      }

      for (const size_t firstSegmentIndex : firstEntry.second)
      {
        for (const size_t secondSegmentIndex : secondIter->second)
        {
          SharedSurfaceMatch sharedSurface;
          sharedSurface.sid = firstEntry.first.first;
          sharedSurface.sindex = firstEntry.first.second;
          sharedSurface.firstStoredSegmentIndex = firstSegmentIndex;
          sharedSurface.secondStoredSegmentIndex = secondSegmentIndex;
          sharedSurfaces.push_back(sharedSurface);
        }
      }
    }

    return sharedSurfaces;
  }

  void fillSharedFoilDeltaZScatterDiagnostics(
    twoelectronhist::HistogramBook& histograms,
    const twoparticlevertexer::VertexResult& vertex,
    const mu2e::TrkSegInfo* firstSegment,
    const mu2e::TrkSegInfo* secondSegment,
    const mu2e::SimInfo& truthOrigin,
    SharedFoilDeltaZScatterCategory category)
  {
    if (!vertex.valid || firstSegment == nullptr || secondSegment == nullptr)
    {
      return;
    }

    if (firstSegment->sid != mu2e::SurfaceIdDetail::ST_Foils ||
        secondSegment->sid != mu2e::SurfaceIdDetail::ST_Foils ||
        firstSegment->sindex != secondSegment->sindex)
    {
      return;
    }

    const double l1 = (firstSegment->pos - truthOrigin.pos).R();
    const double l2 = (secondSegment->pos - truthOrigin.pos).R();
    const double maxPointTruthDistance = std::max(l1, l2);
    const double averageAbsLineParameter =
      0.5 * (std::fabs(vertex.firstLineParameter) +
             std::fabs(vertex.secondLineParameter));
    const double recoMinusTruthZ = vertex.vertex.z() - truthOrigin.pos.z();

    switch (category)
    {
      case SharedFoilDeltaZScatterCategory::kAllCandidatePairs:
        twoelectronhist::fillRecoAllSharedFoilCandidateMaxLvsDeltaZ(
          histograms,
          maxPointTruthDistance,
          recoMinusTruthZ);
        twoelectronhist::fillRecoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZ(
          histograms,
          averageAbsLineParameter,
          recoMinusTruthZ);
        break;
      case SharedFoilDeltaZScatterCategory::kSpaceSelectedPair:
        twoelectronhist::fillRecoSpaceSelectedSharedFoilMaxLvsDeltaZ(
          histograms,
          maxPointTruthDistance,
          recoMinusTruthZ);
        twoelectronhist::fillRecoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZ(
          histograms,
          averageAbsLineParameter,
          recoMinusTruthZ);
        break;
      case SharedFoilDeltaZScatterCategory::kMinTimeSelectedPair:
        twoelectronhist::fillRecoMinTimeSharedFoilMaxLvsDeltaZTest(
          histograms,
          maxPointTruthDistance,
          recoMinusTruthZ);
        twoelectronhist::fillRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZTest(
          histograms,
          averageAbsLineParameter,
          recoMinusTruthZ);
        break;
    }
  }

  void fillAllSharedFoilCandidateDeltaZScatterDiagnostics(
    twoelectronhist::HistogramBook& histograms,
    const vector<vector<mu2e::TrkSegInfo>>* trackSegments,
    size_t firstTrackIndex,
    size_t secondTrackIndex,
    const mu2e::SimInfo& truthOrigin)
  {
    if (trackSegments == nullptr ||
        firstTrackIndex >= trackSegments->size() ||
        secondTrackIndex >= trackSegments->size())
    {
      return;
    }

    const auto& firstTrackSegments = trackSegments->at(firstTrackIndex);
    const auto& secondTrackSegments = trackSegments->at(secondTrackIndex);
    const SharedSurfaceMatches sharedSurfaces =
      sharedSurfacesBetweenTracks(firstTrackSegments, secondTrackSegments);

    for (const auto& sharedSurface : sharedSurfaces)
    {
      if (sharedSurface.sid != mu2e::SurfaceIdDetail::ST_Foils)
      {
        continue;
      }

      const auto& firstSharedSegment =
        firstTrackSegments.at(sharedSurface.firstStoredSegmentIndex);
      const auto& secondSharedSegment =
        secondTrackSegments.at(sharedSurface.secondStoredSegmentIndex);

      const auto firstState =
        twoparticlevertexer::makeParticleStateFromTrackSegment(
          firstSharedSegment,
          static_cast<int>(firstTrackIndex),
          "first selected reco track shared-foil candidate");
      const auto secondState =
        twoparticlevertexer::makeParticleStateFromTrackSegment(
          secondSharedSegment,
          static_cast<int>(secondTrackIndex),
          "second selected reco track shared-foil candidate");

      const auto candidateVertex =
        twoparticlevertexer::vertexFromParticleStates(firstState, secondState);
      fillSharedFoilDeltaZScatterDiagnostics(
        histograms,
        candidateVertex,
        &firstSharedSegment,
        &secondSharedSegment,
        truthOrigin,
        SharedFoilDeltaZScatterCategory::kAllCandidatePairs);
    }
  }

  int countUniqueStoppingTargetFoils(
    const vector<vector<mu2e::TrkSegInfo>>* trackSegments,
    size_t trackIndex);

  int countSharedStoppingTargetFoils(
    const vector<vector<mu2e::TrkSegInfo>>* trackSegments,
    size_t firstTrackIndex,
    size_t secondTrackIndex);

  void fillSpaceSelectedSharedFoilDeltaZRelationDiagnostics(
    twoelectronhist::HistogramBook& histograms,
    const vector<vector<mu2e::TrkSegInfo>>* trackSegments,
    size_t firstTrackIndex,
    size_t secondTrackIndex,
    const twoparticlevertexer::VertexResult& vertex,
    const mu2e::TrkSegInfo* firstSegment,
    const mu2e::TrkSegInfo* secondSegment,
    const mu2e::SimInfo& truthOrigin)
  {
    if (!vertex.valid || firstSegment == nullptr || secondSegment == nullptr)
    {
      return;
    }

    if (firstSegment->sid != mu2e::SurfaceIdDetail::ST_Foils ||
        secondSegment->sid != mu2e::SurfaceIdDetail::ST_Foils ||
        firstSegment->sindex != secondSegment->sindex)
    {
      return;
    }

    const int firstTrackFoilsHit =
      countUniqueStoppingTargetFoils(trackSegments, firstTrackIndex);
    const int secondTrackFoilsHit =
      countUniqueStoppingTargetFoils(trackSegments, secondTrackIndex);
    const int sharedFoilsHit =
      countSharedStoppingTargetFoils(trackSegments,
                                     firstTrackIndex,
                                     secondTrackIndex);
    const int maxFoilsHit = std::max(firstTrackFoilsHit, secondTrackFoilsHit);
    const int absDeltaFoilsHit =
      firstTrackFoilsHit >= secondTrackFoilsHit
        ? firstTrackFoilsHit - secondTrackFoilsHit
        : secondTrackFoilsHit - firstTrackFoilsHit;
    const double recoMinusTruthZ = vertex.vertex.z() - truthOrigin.pos.z();
    const double openingAngleDegrees =
      twoelectronhist::momentumOpeningAngle(vertex);

    twoelectronhist::fillRecoSpaceSelectedSharedFoilNumberVsDeltaZ(
      histograms,
      firstSegment->sindex,
      recoMinusTruthZ);
    twoelectronhist::fillRecoSpaceSelectedSharedFoilCountVsDeltaZ(
      histograms,
      sharedFoilsHit,
      recoMinusTruthZ);
    twoelectronhist::fillRecoSpaceSelectedMaxFoilsHitVsDeltaZ(
      histograms,
      maxFoilsHit,
      recoMinusTruthZ);
    twoelectronhist::fillRecoSpaceSelectedAbsDeltaFoilsHitVsDeltaZ(
      histograms,
      absDeltaFoilsHit,
      recoMinusTruthZ);
    twoelectronhist::fillRecoSpaceSelectedOpeningAngleVsDeltaZ(
      histograms,
      openingAngleDegrees,
      recoMinusTruthZ);
  }

  //============================================================================
  // Reconstruction Helpers
  //============================================================================

  struct MomentumChoice
  {
    bool valid = false;
    double momentum = -1.0;
    string source;
  };

  bool hasAssociatedCaloHit(const vector<mu2e::TrkCaloHitInfo>* trackCaloHits,
                            size_t trackIndex)
  {
    if (trackCaloHits == nullptr || trackIndex >= trackCaloHits->size())
    {
      return false;
    }

    const auto& caloHit = trackCaloHits->at(trackIndex);
    return caloHit.active && caloHit.did >= 0;
  }

  const mu2e::TrkSegInfo* findTrackerMiddleSegment(
    const vector<mu2e::TrkSegInfo>& segments)
  {
    for (const auto& segment : segments)
    {
      if (segment.sid == mu2e::SurfaceIdDetail::TT_Mid)
      {
        return &segment;
      }
    }
    return nullptr;
  }

  bool hasDownstreamTrackerMiddleSegment(const vector<mu2e::TrkSegInfo>& segments)
  {
    const mu2e::TrkSegInfo* trackerMiddleSegment = findTrackerMiddleSegment(segments);
    return trackerMiddleSegment != nullptr && trackerMiddleSegment->mom.z() > 0.0;
  }

  int countUniqueStoppingTargetFoils(const vector<mu2e::TrkSegInfo>& segments)
  {
    set<int> foilIndices;
    for (const auto& segment : segments)
    {
      if (segment.sid == mu2e::SurfaceIdDetail::ST_Foils &&
          segment.sindex >= 0)
      {
        foilIndices.insert(segment.sindex);
      }
    }
    return static_cast<int>(foilIndices.size());
  }

  int countUniqueStoppingTargetFoils(
    const vector<vector<mu2e::TrkSegInfo>>* trackSegments,
    size_t trackIndex)
  {
    if (trackSegments == nullptr || trackIndex >= trackSegments->size())
    {
      return 0;
    }
    return countUniqueStoppingTargetFoils(trackSegments->at(trackIndex));
  }

  int countSharedStoppingTargetFoils(
    const vector<vector<mu2e::TrkSegInfo>>* trackSegments,
    size_t firstTrackIndex,
    size_t secondTrackIndex)
  {
    if (trackSegments == nullptr ||
        firstTrackIndex >= trackSegments->size() ||
        secondTrackIndex >= trackSegments->size())
    {
      return 0;
    }

    set<int> firstFoilIndices;
    for (const auto& segment : trackSegments->at(firstTrackIndex))
    {
      if (segment.sid == mu2e::SurfaceIdDetail::ST_Foils &&
          segment.sindex >= 0)
      {
        firstFoilIndices.insert(segment.sindex);
      }
    }

    set<int> sharedFoilIndices;
    for (const auto& segment : trackSegments->at(secondTrackIndex))
    {
      if (segment.sid == mu2e::SurfaceIdDetail::ST_Foils &&
          segment.sindex >= 0 &&
          firstFoilIndices.find(segment.sindex) != firstFoilIndices.end())
      {
        sharedFoilIndices.insert(segment.sindex);
      }
    }

    return static_cast<int>(sharedFoilIndices.size());
  }

  bool selectedSharedFoilSegmentsHaveOppositePz(
    const mu2e::TrkSegInfo* firstSegment,
    const mu2e::TrkSegInfo* secondSegment)
  {
    if (firstSegment == nullptr || secondSegment == nullptr)
    {
      return false;
    }

    if (firstSegment->sid != mu2e::SurfaceIdDetail::ST_Foils ||
        secondSegment->sid != mu2e::SurfaceIdDetail::ST_Foils ||
        firstSegment->sindex != secondSegment->sindex)
    {
      return false;
    }

    const double firstPz = firstSegment->mom.z();
    const double secondPz = secondSegment->mom.z();
    return (firstPz > 0.0 && secondPz < 0.0) ||
           (firstPz < 0.0 && secondPz > 0.0);
  }

  bool isRankZeroDownstreamElectron(const mu2e::SimInfo& sim, int electronPdg = 11)
  {
    return sim.valid &&
           sim.rank == 0 &&
           sim.pdg == electronPdg &&
           sim.mom.z() > 0.0;
  }

  MomentumChoice chooseRecoMomentum(const vector<mu2e::TrkCaloHitInfo>* trackCaloHits,
                                    const vector<vector<mu2e::TrkSegInfo>>* trackSegments,
                                    size_t trackIndex)
  {
    MomentumChoice choice;

    // Prefer the momentum at the track-calo point of closest approach when a
    // usable track-calo association exists.  This ties the momentum cut directly
    // to the same reconstructed object used by the "has calo" decision.
    if (hasAssociatedCaloHit(trackCaloHits, trackIndex))
    {
      const auto& caloHit = trackCaloHits->at(trackIndex);
      const double caloMomentum = caloHit.mom.R();
      if (caloMomentum > 0.0)
      {
        choice.valid = true;
        choice.momentum = caloMomentum;
        choice.source = "trkcalohit.mom";
        return choice;
      }
    }

    if (trackSegments == nullptr || trackIndex >= trackSegments->size())
    {
      return choice;
    }

    const auto& segments = trackSegments->at(trackIndex);
    const mu2e::TrkSegInfo* trackerMiddleSegment = findTrackerMiddleSegment(segments);
    if (trackerMiddleSegment != nullptr && trackerMiddleSegment->mom.R() > 0.0)
    {
      choice.valid = true;
      choice.momentum = trackerMiddleSegment->mom.R();
      choice.source = "trksegs TT_Mid mom";
      return choice;
    }

    if (!segments.empty() && segments.front().mom.R() > 0.0)
    {
      choice.valid = true;
      choice.momentum = segments.front().mom.R();
      choice.source = "first trkseg mom";
    }

    return choice;
  }

  bool momentumInWindow(const MomentumChoice& choice,
                        double momentumMin,
                        double momentumMax)
  {
    return choice.valid &&
           choice.momentum >= momentumMin &&
           choice.momentum <= momentumMax;
  }

  struct RecoTrackDecision
  {
    size_t trackIndex = 0;
    bool pdgPass = true;
    bool downstreamPass = true;
    bool hasCalo = false;
    bool momentumPass = false;
    bool candidatePass = false;
    MomentumChoice momentum;
  };

  RecoTrackDecision evaluateRecoTrack(
    const vector<mu2e::TrkInfo>* tracks,
    const vector<vector<mu2e::TrkSegInfo>>* trackSegments,
    const vector<mu2e::TrkCaloHitInfo>* trackCaloHits,
    size_t trackIndex,
    double momentumMin,
    double momentumMax,
    int requiredRecoPdg,
    bool requireDownstream)
  {
    RecoTrackDecision decision;
    decision.trackIndex = trackIndex;

    if (tracks != nullptr && trackIndex < tracks->size() && requiredRecoPdg != 0)
    {
      decision.pdgPass = tracks->at(trackIndex).pdg == requiredRecoPdg;
    }

    if (requireDownstream)
    {
      decision.downstreamPass =
        trackSegments != nullptr &&
        trackIndex < trackSegments->size() &&
        hasDownstreamTrackerMiddleSegment(trackSegments->at(trackIndex));
    }

    decision.hasCalo = hasAssociatedCaloHit(trackCaloHits, trackIndex);
    decision.momentum = chooseRecoMomentum(trackCaloHits, trackSegments, trackIndex);
    decision.momentumPass = momentumInWindow(decision.momentum, momentumMin, momentumMax);
    decision.candidatePass =
      decision.pdgPass &&
      decision.downstreamPass &&
      decision.hasCalo &&
      decision.momentumPass;

    return decision;
  }

  void printTruthForTrack(const vector<vector<mu2e::SimInfo>>* truthSimByTrack,
                          size_t trackIndex,
                          int maxTruthPerTrack)
  {
    if (truthSimByTrack == nullptr || trackIndex >= truthSimByTrack->size())
    {
      cout << "      truth trkmcsim: unavailable for this track index" << endl;
      return;
    }

    const auto& truthSims = truthSimByTrack->at(trackIndex);
    cout << "      truth trkmcsim entries: " << truthSims.size() << endl;

    const size_t nToPrint =
      maxTruthPerTrack < 0
        ? truthSims.size()
        : min(truthSims.size(), static_cast<size_t>(maxTruthPerTrack));

    for (size_t iTruth = 0; iTruth < nToPrint; ++iTruth)
    {
      const auto& sim = truthSims.at(iTruth);
      cout << "        truth[" << iTruth << "]"
           << " valid=" << sim.valid
           << " rank=" << sim.rank
           << " id=" << sim.id
           << " pdg=" << sim.pdg
           << " nhits=" << sim.nhits
           << " nactive=" << sim.nactive
           << " origin(t,x,y,z)=("
           << sim.time << ", "
           << sim.pos.x() << ", "
           << sim.pos.y() << ", "
           << sim.pos.z() << ")"
           << " origin_mom=" << formatVector3(sim.mom)
           << endl;
    }

    if (nToPrint < truthSims.size())
    {
      cout << "        ... " << (truthSims.size() - nToPrint)
           << " additional trkmcsim entries suppressed by maxTruthPerTrack"
           << endl;
    }
  }

  const mu2e::SimInfo* findRepresentativeTruthOrigin(
    const vector<vector<mu2e::SimInfo>>* truthSimByTrack,
    const vector<size_t>& selectedTrackIndices,
    int electronPdg = 11)
  {
    if (truthSimByTrack == nullptr)
    {
      return nullptr;
    }

    for (const size_t trackIndex : selectedTrackIndices)
    {
      if (trackIndex >= truthSimByTrack->size())
      {
        continue;
      }

      const auto& truthSims = truthSimByTrack->at(trackIndex);
      for (const auto& sim : truthSims)
      {
        if (isRankZeroDownstreamElectron(sim, electronPdg))
        {
          return &sim;
        }
      }
    }

    return nullptr;
  }

  enum class SharedSurfaceSelectionMetric
  {
    kMinDistance,
    kMinAbsTimeDifference
  };

  // Build the reconstructed vertex from shared ST_Foils only, and pick the
  // shared foil with the minimum score under the requested metric.
  twoparticlevertexer::VertexResult buildVertexForSelectedTrackPairByMetric(
    const vector<mu2e::TrkInfo>* tracks,
    const vector<vector<mu2e::TrkSegInfo>>* trackSegments,
    size_t firstTrackIndex,
    size_t secondTrackIndex,
    const mu2e::TrkSegInfo*& firstSegment,
    const mu2e::TrkSegInfo*& secondSegment,
    SharedSurfaceSelectionMetric metric)
  {
    twoparticlevertexer::VertexResult vertex;
    firstSegment = nullptr;
    secondSegment = nullptr;

    if (trackSegments == nullptr ||
        firstTrackIndex >= trackSegments->size() ||
        secondTrackIndex >= trackSegments->size())
    {
      vertex.failureReason = "missing track-segment data";
      return vertex;
    }

    if (tracks == nullptr ||
        firstTrackIndex >= tracks->size() ||
        secondTrackIndex >= tracks->size())
    {
      vertex.failureReason = "missing reconstructed track summary";
      return vertex;
    }

    const auto& firstTrackSegments = trackSegments->at(firstTrackIndex);
    const auto& secondTrackSegments = trackSegments->at(secondTrackIndex);
    const SharedSurfaceMatches sharedSurfaces =
      sharedSurfacesBetweenTracks(firstTrackSegments, secondTrackSegments);

    bool sawSharedSTFoil = false;
    bool foundValidSharedSTFoil = false;
    twoparticlevertexer::VertexResult bestVertex;
    double bestScore = numeric_limits<double>::infinity();

    for (const auto& sharedSurface : sharedSurfaces)
    {
      if (sharedSurface.sid != mu2e::SurfaceIdDetail::ST_Foils)
      {
        continue;
      }

      sawSharedSTFoil = true;
      const auto& firstSharedSegment =
        firstTrackSegments.at(sharedSurface.firstStoredSegmentIndex);
      const auto& secondSharedSegment =
        secondTrackSegments.at(sharedSurface.secondStoredSegmentIndex);

      const auto firstState =
        twoparticlevertexer::makeParticleStateFromTrackSegment(
          firstSharedSegment,
          static_cast<int>(firstTrackIndex),
          "first selected reco track");
      const auto secondState =
        twoparticlevertexer::makeParticleStateFromTrackSegment(
          secondSharedSegment,
          static_cast<int>(secondTrackIndex),
          "second selected reco track");

      const auto candidateVertex =
        twoparticlevertexer::vertexFromParticleStates(firstState, secondState);
      if (!candidateVertex.valid)
      {
        continue;
      }

      foundValidSharedSTFoil = true;
      const double candidateScore =
        metric == SharedSurfaceSelectionMetric::kMinDistance
          ? candidateVertex.distance
          : fabs(candidateVertex.deltaInputTime);

      if (candidateScore < bestScore)
      {
        bestScore = candidateScore;
        bestVertex = candidateVertex;
        firstSegment = &firstSharedSegment;
        secondSegment = &secondSharedSegment;
      }
    }

    if (!sawSharedSTFoil)
    {
      vertex.failureReason = "no shared ST_Foils surface";
      return vertex;
    }

    if (!foundValidSharedSTFoil)
    {
      vertex.failureReason = "shared ST_Foils lines are parallel or invalid";
      return vertex;
    }

    return bestVertex;
  }

  twoparticlevertexer::VertexResult buildVertexForSelectedTrackPair(
    const vector<mu2e::TrkInfo>* tracks,
    const vector<vector<mu2e::TrkSegInfo>>* trackSegments,
    size_t firstTrackIndex,
    size_t secondTrackIndex,
    const mu2e::TrkSegInfo*& firstSegment,
    const mu2e::TrkSegInfo*& secondSegment)
  {
    return buildVertexForSelectedTrackPairByMetric(
      tracks,
      trackSegments,
      firstTrackIndex,
      secondTrackIndex,
      firstSegment,
      secondSegment,
      SharedSurfaceSelectionMetric::kMinDistance);
  }

  twoparticlevertexer::VertexResult buildVertexForSelectedTrackPairMinTimeDifference(
    const vector<mu2e::TrkInfo>* tracks,
    const vector<vector<mu2e::TrkSegInfo>>* trackSegments,
    size_t firstTrackIndex,
    size_t secondTrackIndex,
    const mu2e::TrkSegInfo*& firstSegment,
    const mu2e::TrkSegInfo*& secondSegment)
  {
    return buildVertexForSelectedTrackPairByMetric(
      tracks,
      trackSegments,
      firstTrackIndex,
      secondTrackIndex,
      firstSegment,
      secondSegment,
      SharedSurfaceSelectionMetric::kMinAbsTimeDifference);
  }
}

void twoElectronTruthTrkSegVertexerComparer(const string& inputName,
                                            int maxEvents = -1,
                                            const string& trackBranch = "trk",
                                            double momentumCutMin = 50.0,
                                            double momentumCutMax = 53.0,
                                            int requiredRecoPdg = 11,
                                            bool requireDownstream = true,
                                            int maxTruthPerTrack = -1,
                                            const string& explicitSimBranch = "",
                                            const string& explicitCaloBranch = "")
{
  //============================================================================
  // Main Macro Setup
  //============================================================================

  const string segmentBranch = deriveSegmentBranchName(trackBranch);
  const string simBranch =
    explicitSimBranch.empty() ? deriveSimBranchName(trackBranch) : explicitSimBranch;
  const string caloBranch =
    explicitCaloBranch.empty() ? deriveCaloBranchName(trackBranch) : explicitCaloBranch;

  TChain ntuple("EventNtuple/ntuple");
  if (!addInput(ntuple, inputName))
  {
    return;
  }

  const Long64_t nEntries = ntuple.GetEntries();
  if (nEntries <= 0)
  {
    cerr << "ERROR: no entries found in EventNtuple/ntuple for input: "
         << inputName << endl;
    return;
  }

  if (ntuple.GetBranch(trackBranch.c_str()) == nullptr)
  {
    cerr << "ERROR: missing requested reconstructed track branch '"
         << trackBranch << "'." << endl;
    return;
  }

  if (ntuple.GetBranch(segmentBranch.c_str()) == nullptr)
  {
    cerr << "ERROR: missing requested reconstructed track-segment branch '"
         << segmentBranch << "'." << endl;
    return;
  }

  if (ntuple.GetBranch(simBranch.c_str()) == nullptr)
  {
    cerr << "ERROR: missing requested truth branch '" << simBranch << "'." << endl;
    return;
  }

  const bool hasCaloBranch = ntuple.GetBranch(caloBranch.c_str()) != nullptr;
  if (!hasCaloBranch)
  {
    cout << "WARNING: missing requested track-calo branch '" << caloBranch
         << "'. Reco-calo decisions will be false." << endl;
  }

  ntuple.SetBranchStatus("*", 0);
  enableBranch(ntuple, "evtinfo");
  enableBranch(ntuple, trackBranch);
  enableBranch(ntuple, segmentBranch);
  enableBranch(ntuple, simBranch);
  if (hasCaloBranch)
  {
    enableBranch(ntuple, caloBranch);
  }

  mu2e::EventInfo* evtinfo = nullptr;
  vector<mu2e::TrkInfo>* tracks = nullptr;
  vector<vector<mu2e::TrkSegInfo>>* trackSegments = nullptr;
  vector<vector<mu2e::SimInfo>>* truthSimByTrack = nullptr;
  vector<mu2e::TrkCaloHitInfo>* trackCaloHits = nullptr;

  if (ntuple.GetBranch("evtinfo") != nullptr)
  {
    ntuple.SetBranchAddress("evtinfo", &evtinfo);
  }
  ntuple.SetBranchAddress(trackBranch.c_str(), &tracks);
  ntuple.SetBranchAddress(segmentBranch.c_str(), &trackSegments);
  ntuple.SetBranchAddress(simBranch.c_str(), &truthSimByTrack);
  if (hasCaloBranch)
  {
    ntuple.SetBranchAddress(caloBranch.c_str(), &trackCaloHits);
  }

  const Long64_t entriesToRead =
    (maxEvents >= 0 && static_cast<Long64_t>(maxEvents) < nEntries) ? maxEvents : nEntries;

  if (DO_FULL_PRINT)
  {
    cout << "Input: " << inputName << endl;
    cout << "Tree entries available: " << nEntries << endl;
    cout << "Tree entries being scanned: " << entriesToRead << endl;
    cout << "Track branch: " << trackBranch << endl;
    cout << "Track-segment branch: " << segmentBranch << endl;
    cout << "Truth branch: " << simBranch << endl;
    cout << "Track-calo branch: " << (hasCaloBranch ? caloBranch : string("missing")) << endl;
    cout << "Reco momentum selection cut: ["
         << momentumCutMin << ", " << momentumCutMax << "] MeV/c" << endl;
    cout << "Required reco PDG: "
         << (requiredRecoPdg == 0 ? string("disabled") : to_string(requiredRecoPdg)) << endl;
    cout << "Require downstream TT_Mid pz > 0: "
         << (requireDownstream ? "yes" : "no") << endl;
    cout << "maxTruthPerTrack: "
         << (maxTruthPerTrack < 0 ? string("all") : to_string(maxTruthPerTrack)) << endl;
  }

  //============================================================================
  // Event Loop
  //============================================================================

  twoelectronhist::HistogramBook histograms = twoelectronhist::bookHistograms();

  Long64_t eventsWithRecoTracks = 0;
  Long64_t eventsWithAtLeastOneRecoTrackInMomentumWindow = 0;
  Long64_t eventsWithExactlyTwoRecoTracksInMomentumWindow = 0;
  Long64_t eventsWithExactlyTwoRecoTracksInMomentumWindowAndGoodCalo = 0;
  Long64_t eventsWithAtLeastOneCandidateTrack = 0;
  Long64_t eventsWithExactlyTwoCandidateTracks = 0;
  Long64_t eventsWithAtLeastTwoCandidateTracks = 0;
  Long64_t totalRecoTracks = 0;
  Long64_t totalCandidateTracks = 0;
  Long64_t histogramSelectedSpaceSharedFoilOppositePzPairs = 0;
  Long64_t histogramSelectedSpaceSharedFoilOppositePzTracks = 0;
  Long64_t histogramSelectedTimeSharedFoilOppositePzPairs = 0;
  Long64_t histogramSelectedTimeSharedFoilOppositePzTracks = 0;
  Long64_t spaceSelectedVertexFoilIndexChecks = 0;
  Long64_t spaceSelectedVertexFoilIndexMatches = 0;
  Long64_t timeSelectedVertexFoilIndexChecks = 0;
  Long64_t timeSelectedVertexFoilIndexMatches = 0;
  Long64_t histogramSelectedSpaceVertexFoilIndexChecks = 0;
  Long64_t histogramSelectedSpaceVertexFoilIndexMatches = 0;
  Long64_t histogramSelectedTimeVertexFoilIndexChecks = 0;
  Long64_t histogramSelectedTimeVertexFoilIndexMatches = 0;

  for (Long64_t entry = 0; entry < entriesToRead; ++entry)
  {
    ntuple.GetEntry(entry);

    const size_t nTracks = tracks != nullptr ? tracks->size() : 0;
    const bool eventHasRecoTracks = nTracks > 0;
    if (eventHasRecoTracks)
    {
      ++eventsWithRecoTracks;
    }
    totalRecoTracks += static_cast<Long64_t>(nTracks);

    if (DO_FULL_PRINT)
    {
      cout << "-------------------------------------------------------------------------------" << endl;
      cout << "Event " << formatEventLabel(evtinfo, entry) << endl;
      cout << "  reconstructed track branch present: yes" << endl;
      cout << "  reconstructed tracks stored: " << nTracks << endl;
      if (truthSimByTrack != nullptr)
      {
        cout << "  trkmcsim outer-vector size: " << truthSimByTrack->size() << endl;
      }
      else
      {
        cout << "  trkmcsim outer-vector size: unavailable" << endl;
      }
    }

    vector<RecoTrackDecision> candidateTrackDecisions;
    size_t recoTracksInMomentumWindow = 0;
    size_t recoTracksInMomentumWindowAndGoodCalo = 0;

    for (size_t iTrack = 0; iTrack < nTracks; ++iTrack)
    {
      const auto& track = tracks->at(iTrack);
      const RecoTrackDecision decision =
        evaluateRecoTrack(tracks,
                          trackSegments,
                          trackCaloHits,
                          iTrack,
                          momentumCutMin,
                          momentumCutMax,
                          requiredRecoPdg,
                          requireDownstream);

      if (decision.momentumPass)
      {
        ++recoTracksInMomentumWindow;
        if (decision.hasCalo)
        {
          ++recoTracksInMomentumWindowAndGoodCalo;
        }
      }

      if (DO_FULL_PRINT)
      {
        cout << "    reco track[" << iTrack << "]"
             << " pdg=" << track.pdg
             << " status=" << track.status
             << " goodfit=" << track.goodfit
             << " nhits=" << track.nhits
             << " nactive=" << track.nactive
             << " fitcon=" << track.fitcon
             << endl;

        printTruthForTrack(truthSimByTrack, iTrack, maxTruthPerTrack);

        cout << "      reco branch for track: yes" << endl;
        cout << "      associated calo hit: " << (decision.hasCalo ? "yes" : "no") << endl;
        cout << "      reco momentum: ";
        if (decision.momentum.valid)
        {
          cout << decision.momentum.momentum
               << " MeV/c from " << decision.momentum.source;
        }
        else
        {
          cout << "unavailable";
        }
        cout << endl;
        cout << "      reco momentum in [" << momentumCutMin << ", " << momentumCutMax
             << "] MeV/c: " << (decision.momentumPass ? "yes" : "no") << endl;
        cout << "      reco PDG requirement: " << (decision.pdgPass ? "pass" : "fail") << endl;
        cout << "      downstream requirement: "
             << (decision.downstreamPass ? "pass" : "fail") << endl;
        cout << "      selected candidate track: "
             << (decision.candidatePass ? "yes" : "no") << endl;
      }

      if (decision.candidatePass)
      {
        candidateTrackDecisions.push_back(decision);
      }
    }

    if (recoTracksInMomentumWindow >= 1)
    {
      ++eventsWithAtLeastOneRecoTrackInMomentumWindow;
    }
    if (recoTracksInMomentumWindow == 2)
    {
      ++eventsWithExactlyTwoRecoTracksInMomentumWindow;
    }
    if (recoTracksInMomentumWindowAndGoodCalo == 2)
    {
      ++eventsWithExactlyTwoRecoTracksInMomentumWindowAndGoodCalo;
    }

    totalCandidateTracks += static_cast<Long64_t>(candidateTrackDecisions.size());
    if (!candidateTrackDecisions.empty())
    {
      ++eventsWithAtLeastOneCandidateTrack;
    }
    if (candidateTrackDecisions.size() == 2)
    {
      ++eventsWithExactlyTwoCandidateTracks;
    }
    if (candidateTrackDecisions.size() >= 2)
    {
      ++eventsWithAtLeastTwoCandidateTracks;
    }

    if (DO_FULL_PRINT)
    {
      cout << "  event candidate-track count: " << candidateTrackDecisions.size() << endl;
      cout << "  event has two reconstructed tracks with requested properties: "
           << (candidateTrackDecisions.size() >= 2 ? "yes" : "no") << endl;
    }

    // Build vertices for the first two selected tracks when at least two pass.
    // The printout uses these as diagnostic previews.  The histogram fill below
    // is stricter and requires exactly two selected tracks in the event.
    bool vertexWasAttempted = false;
    const mu2e::TrkSegInfo* firstVertexSegment = nullptr;
    const mu2e::TrkSegInfo* secondVertexSegment = nullptr;
    const mu2e::TrkSegInfo* firstTimeVertexSegment = nullptr;
    const mu2e::TrkSegInfo* secondTimeVertexSegment = nullptr;
    twoparticlevertexer::VertexResult selectedPairVertex;
    twoparticlevertexer::VertexResult timeSelectedPairVertex;
    VertexFoilIndexCheck selectedPairVertexFoilCheck;
    VertexFoilIndexCheck timeSelectedPairVertexFoilCheck;
    size_t firstVertexTrackIndex = 0;
    size_t secondVertexTrackIndex = 0;

    if (candidateTrackDecisions.size() >= 2)
    {
      firstVertexTrackIndex = candidateTrackDecisions.at(0).trackIndex;
      secondVertexTrackIndex = candidateTrackDecisions.at(1).trackIndex;
      selectedPairVertex =
        buildVertexForSelectedTrackPair(tracks,
                                        trackSegments,
                                        firstVertexTrackIndex,
                                        secondVertexTrackIndex,
                                        firstVertexSegment,
                                        secondVertexSegment);
      timeSelectedPairVertex =
        buildVertexForSelectedTrackPairMinTimeDifference(
          tracks,
          trackSegments,
          firstVertexTrackIndex,
          secondVertexTrackIndex,
          firstTimeVertexSegment,
          secondTimeVertexSegment);
      vertexWasAttempted = true;

      selectedPairVertexFoilCheck =
        compareVertexToSharedFoilIndex(selectedPairVertex,
                                       firstVertexSegment,
                                       secondVertexSegment);
      timeSelectedPairVertexFoilCheck =
        compareVertexToSharedFoilIndex(timeSelectedPairVertex,
                                       firstTimeVertexSegment,
                                       secondTimeVertexSegment);
      countVertexFoilIndexCheck(selectedPairVertexFoilCheck,
                                spaceSelectedVertexFoilIndexChecks,
                                spaceSelectedVertexFoilIndexMatches);
      countVertexFoilIndexCheck(timeSelectedPairVertexFoilCheck,
                                timeSelectedVertexFoilIndexChecks,
                                timeSelectedVertexFoilIndexMatches);
    }

    // Histogram selection:
    //   - exactly two selected reconstructed tracks in the event
    //   - each selected track is a downstream electron by the reco selection
    //   - each selected track has an associated calo hit
    //   - each selected track has reconstructed momentum in the 50-53 MeV/c
    //     selection cut configured above
    //
    // Only after this event-level reco selection passes do we fill any truth,
    // reco-momentum, or vertex histogram.
    if (candidateTrackDecisions.size() == 2)
    {
      ++histograms.selectedEvents;

      countVertexFoilIndexCheck(selectedPairVertexFoilCheck,
                                histogramSelectedSpaceVertexFoilIndexChecks,
                                histogramSelectedSpaceVertexFoilIndexMatches);
      countVertexFoilIndexCheck(timeSelectedPairVertexFoilCheck,
                                histogramSelectedTimeVertexFoilIndexChecks,
                                histogramSelectedTimeVertexFoilIndexMatches);

      if (vertexWasAttempted &&
          selectedPairVertex.valid &&
          selectedSharedFoilSegmentsHaveOppositePz(firstVertexSegment,
                                                  secondVertexSegment))
      {
        ++histogramSelectedSpaceSharedFoilOppositePzPairs;
        histogramSelectedSpaceSharedFoilOppositePzTracks += 2;
      }

      if (vertexWasAttempted &&
          timeSelectedPairVertex.valid &&
          selectedSharedFoilSegmentsHaveOppositePz(firstTimeVertexSegment,
                                                  secondTimeVertexSegment))
      {
        ++histogramSelectedTimeSharedFoilOppositePzPairs;
        histogramSelectedTimeSharedFoilOppositePzTracks += 2;
      }

      vector<size_t> selectedTrackIndices;
      selectedTrackIndices.push_back(candidateTrackDecisions.at(0).trackIndex);
      selectedTrackIndices.push_back(candidateTrackDecisions.at(1).trackIndex);

      const mu2e::SimInfo* truthOrigin =
        findRepresentativeTruthOrigin(truthSimByTrack, selectedTrackIndices, 11);

      if (truthOrigin != nullptr && vertexWasAttempted)
      {
        fillAllSharedFoilCandidateDeltaZScatterDiagnostics(
          histograms,
          trackSegments,
          firstVertexTrackIndex,
          secondVertexTrackIndex,
          *truthOrigin);
      }

      twoelectronhist::fillRankZeroDownstreamElectronTruthFromSelectedTracks(
        histograms,
        truthSimByTrack,
        selectedTrackIndices,
        11);

      for (const auto& decision : candidateTrackDecisions)
      {
        twoelectronhist::fillRecoMomentum(histograms, decision.momentum.momentum);
      }

      if (vertexWasAttempted && selectedPairVertex.valid)
      {
        twoelectronhist::fillRecoVertex(histograms, selectedPairVertex);
        twoelectronhist::fillRecoVertexSelectedSegmentTimeDifference(
          histograms,
          selectedPairVertex);

        if (selectedPairVertexFoilCheck.valid &&
            selectedPairVertexFoilCheck.matchesSharedFoilIndex)
        {
          twoelectronhist::fillRecoVertexFoilIndexMatchedMaps(
            histograms,
            selectedPairVertex);
        }

        if (truthOrigin != nullptr)
        {
          const XYZVectorF truthToReco = selectedPairVertex.vertex - truthOrigin->pos;
          twoelectronhist::fillRecoVertexTruthResidual(histograms, truthToReco);
          fillSharedFoilDeltaZScatterDiagnostics(
            histograms,
            selectedPairVertex,
            firstVertexSegment,
            secondVertexSegment,
            *truthOrigin,
            SharedFoilDeltaZScatterCategory::kSpaceSelectedPair);
          fillSpaceSelectedSharedFoilDeltaZRelationDiagnostics(
            histograms,
            trackSegments,
            firstVertexTrackIndex,
            secondVertexTrackIndex,
            selectedPairVertex,
            firstVertexSegment,
            secondVertexSegment,
            *truthOrigin);

          if (selectedPairVertexFoilCheck.valid &&
              selectedPairVertexFoilCheck.matchesSharedFoilIndex)
          {
            twoelectronhist::fillRecoVertexFoilIndexMatchedTruthResidual(
              histograms,
              truthToReco);
          }
        }
      }

      if (vertexWasAttempted && timeSelectedPairVertex.valid)
      {
        twoelectronhist::fillRecoVertexMinTimeChoiceTest(
          histograms,
          timeSelectedPairVertex);
        twoelectronhist::fillRecoVertexMinTimeDifferenceTest(
          histograms,
          timeSelectedPairVertex);

        if (timeSelectedPairVertexFoilCheck.valid &&
            timeSelectedPairVertexFoilCheck.matchesSharedFoilIndex)
        {
          twoelectronhist::fillRecoVertexMinTimeFoilIndexMatchedMapsTest(
            histograms,
            timeSelectedPairVertex);
        }

        if (truthOrigin != nullptr)
        {
          const XYZVectorF truthToReco =
            timeSelectedPairVertex.vertex - truthOrigin->pos;
          twoelectronhist::fillRecoVertexMinTimeTruthResidualTest(
            histograms,
            truthToReco);
          fillSharedFoilDeltaZScatterDiagnostics(
            histograms,
            timeSelectedPairVertex,
            firstTimeVertexSegment,
            secondTimeVertexSegment,
            *truthOrigin,
            SharedFoilDeltaZScatterCategory::kMinTimeSelectedPair);

          if (timeSelectedPairVertexFoilCheck.valid &&
              timeSelectedPairVertexFoilCheck.matchesSharedFoilIndex)
          {
            twoelectronhist::fillRecoVertexMinTimeFoilIndexMatchedTruthResidualTest(
              histograms,
              truthToReco);
          }
        }
      }
    }

    // Demonstrate the generic helper on the first two passing tracks.  This is
    // a textual preview only; histograms fill from the stricter exactly-two
    // event selection above.
    if (DO_FULL_PRINT && vertexWasAttempted)
    {
        cout << "  twoParticleVertexer preview using selected tracks "
             << firstVertexTrackIndex << " and " << secondVertexTrackIndex << ":" << endl;
        if (selectedPairVertex.valid)
        {
          cout << "    seed surfaces: first sid=" << firstVertexSegment->sid
               << " sindex=" << firstVertexSegment->sindex
               << ", second sid=" << secondVertexSegment->sid
               << " sindex=" << secondVertexSegment->sindex << endl;
          cout << "    vertex midpoint(x,y,z)=" << formatVector3(selectedPairVertex.vertex)
               << " closest-line distance=" << selectedPairVertex.distance << " mm"
               << " delta_input_time=" << selectedPairVertex.deltaInputTime << " ns"
               << endl;
          printVertexFoilIndexCheck("space-selected vertex",
                                    selectedPairVertexFoilCheck);
        }
        else
        {
          cout << "    vertex unavailable: " << selectedPairVertex.failureReason << endl;
          printVertexFoilIndexCheck("space-selected vertex",
                                    selectedPairVertexFoilCheck);
        }

        cout << "  twoParticleVertexer preview minimizing |delta_input_time| using selected tracks "
             << firstVertexTrackIndex << " and " << secondVertexTrackIndex << ":" << endl;
        if (timeSelectedPairVertex.valid)
        {
          cout << "    seed surfaces: first sid=" << firstTimeVertexSegment->sid
               << " sindex=" << firstTimeVertexSegment->sindex
               << ", second sid=" << secondTimeVertexSegment->sid
               << " sindex=" << secondTimeVertexSegment->sindex << endl;
          cout << "    vertex midpoint(x,y,z)=" << formatVector3(timeSelectedPairVertex.vertex)
               << " closest-line distance=" << timeSelectedPairVertex.distance << " mm"
               << " delta_input_time=" << timeSelectedPairVertex.deltaInputTime << " ns"
               << endl;
          printVertexFoilIndexCheck("TEST minimum-time vertex",
                                    timeSelectedPairVertexFoilCheck);
        }
        else
        {
          cout << "    vertex unavailable: " << timeSelectedPairVertex.failureReason << endl;
          printVertexFoilIndexCheck("TEST minimum-time vertex",
                                    timeSelectedPairVertexFoilCheck);
        }
    }

    if (DO_VERTEX_ORIGIN_RESIDUAL_PRINT &&
        vertexWasAttempted &&
        selectedPairVertex.valid)
    {
      const vector<size_t> selectedTrackIndicesForTruth = {
        firstVertexTrackIndex,
        secondVertexTrackIndex
      };
      const mu2e::SimInfo* truthOrigin =
        findRepresentativeTruthOrigin(truthSimByTrack, selectedTrackIndicesForTruth, 11);

      if (truthOrigin != nullptr)
      {
        const XYZVectorF delta = selectedPairVertex.vertex - truthOrigin->pos;
        cout << "  truth-origin vs reconstructed-vertex residual:" << endl;
        cout << "    truth origin(t,x,y,z)=("
             << truthOrigin->time << ", "
             << truthOrigin->pos.x() << ", "
             << truthOrigin->pos.y() << ", "
             << truthOrigin->pos.z() << ")" << endl;
        cout << "    reco vertex(x,y,z)="
             << formatVector3(selectedPairVertex.vertex) << endl;
        cout << "    delta(reco - truth)(dx,dy,dz)="
             << formatVector3(delta) << endl;
        cout << "    abs(delta)(|dx|,|dy|,|dz|)=("
             << fabs(delta.x()) << ", "
             << fabs(delta.y()) << ", "
             << fabs(delta.z()) << ")" << endl;
        cout << "    delta_z=" << delta.z()
             << " mm, |delta_z|=" << fabs(delta.z()) << " mm" << endl;
      }
      else
      {
        cout << "  truth-origin vs reconstructed-vertex residual: unavailable"
             << endl;
      }
    }
  }

  bool histogramFileWasWritten = false;
  if (WRITE_HISTOGRAM_FILE)
  {
    histogramFileWasWritten =
      twoelectronhist::writeHistograms(histograms, HISTOGRAM_OUTPUT_FILE);
    if (histogramFileWasWritten)
    {
      twoelectronplots::saveTruthVsRecoPlotsFromFile(HISTOGRAM_OUTPUT_FILE,
                                                     HISTOGRAM_PLOT_OUTPUT_DIR);
    }
  }

  //============================================================================
  // Summary
  //============================================================================

  cout << "===============================================================================" << endl;
  cout << "twoElectronTruthTrkSegVertexerComparer summary" << endl;
  cout << "  events scanned: " << entriesToRead << endl;
  cout << "  events with at least one reconstructed track: "
       << eventsWithRecoTracks << endl;
  cout << "  total reconstructed tracks: " << totalRecoTracks << endl;
  cout << "  events with at least one reconstructed track with momentum in ["
       << momentumCutMin << ", " << momentumCutMax << "] MeV/c: "
       << eventsWithAtLeastOneRecoTrackInMomentumWindow << endl;
  cout << "  events with exactly two reconstructed tracks with momentum in ["
       << momentumCutMin << ", " << momentumCutMax << "] MeV/c: "
       << eventsWithExactlyTwoRecoTracksInMomentumWindow << endl;
  cout << "  events with exactly two reconstructed tracks with momentum in ["
       << momentumCutMin << ", " << momentumCutMax
       << "] MeV/c and a good calo hit: "
       << eventsWithExactlyTwoRecoTracksInMomentumWindowAndGoodCalo << endl;
  cout << "  total candidate tracks: " << totalCandidateTracks << endl;
  cout << "  events with at least one candidate track: "
       << eventsWithAtLeastOneCandidateTrack << endl;
  cout << "  events with exactly two candidate tracks: "
       << eventsWithExactlyTwoCandidateTracks << endl;
  cout << "  events with at least two candidate tracks: "
       << eventsWithAtLeastTwoCandidateTracks << endl;
  cout << "  histogram-selected space vertices with opposite-pz selected shared ST_Foils pair: "
       << histogramSelectedSpaceSharedFoilOppositePzPairs << endl;
  cout << "  histogram-selected tracks in opposite-pz space selected shared ST_Foils pairs: "
       << histogramSelectedSpaceSharedFoilOppositePzTracks << endl;
  cout << "  histogram-selected TEST minimum-time vertices with opposite-pz selected shared ST_Foils pair: "
       << histogramSelectedTimeSharedFoilOppositePzPairs << endl;
  cout << "  histogram-selected tracks in opposite-pz TEST minimum-time selected shared ST_Foils pairs: "
       << histogramSelectedTimeSharedFoilOppositePzTracks << endl;
  const vector<double>& configuredFoilCentersZ =
    configuredStoppingTargetFoilCentersZ();
  if (!configuredFoilCentersZ.empty())
  {
    cout << "  configured stopping-target foil geometry for vertex-z check: "
         << configuredFoilCentersZ.size()
         << " foils, detector z centers ["
         << configuredFoilCentersZ.front() << ", "
         << configuredFoilCentersZ.back() << "] mm, spacing "
         << CURRENT_ST_GEOMETRY_DELTA_Z << " mm" << endl;
  }
  printVertexFoilIndexSummary(
    "constructed first-two-candidate space vertex foil-index matches",
    spaceSelectedVertexFoilIndexMatches,
    spaceSelectedVertexFoilIndexChecks);
  printVertexFoilIndexSummary(
    "constructed first-two-candidate TEST minimum-time vertex foil-index matches",
    timeSelectedVertexFoilIndexMatches,
    timeSelectedVertexFoilIndexChecks);
  printVertexFoilIndexSummary(
    "histogram-selected exactly-two-candidate space vertex foil-index matches",
    histogramSelectedSpaceVertexFoilIndexMatches,
    histogramSelectedSpaceVertexFoilIndexChecks);
  printVertexFoilIndexSummary(
    "histogram-selected exactly-two-candidate TEST minimum-time vertex foil-index matches",
    histogramSelectedTimeVertexFoilIndexMatches,
    histogramSelectedTimeVertexFoilIndexChecks);
  cout << "  histogram-selected events, exactly two candidate tracks: "
       << histograms.selectedEvents << endl;
  cout << "  histogram MC truth entries, rank-0 downstream e-: "
       << histograms.mcTruthEntries << endl;
  cout << "  histogram reconstructed momentum entries: "
       << histograms.recoMomentumEntries << endl;
  cout << "  histogram reconstructed vertex entries: "
       << histograms.recoVertexEntries << endl;
  cout << "  histogram reconstructed vertex momentum theta-pair entries: "
       << histograms.recoVertexMomentumThetaEntries << endl;
  cout << "  histogram reconstructed vertex momentum opening-angle entries: "
       << histograms.recoVertexMomentumOpeningAngleEntries << endl;
  cout << "  histogram reconstructed vertex line-parameter entries: "
       << histograms.recoVertexLineParameterEntries << endl;
  cout << "  histogram foil-index matched reconstructed vertex map entries: "
       << histograms.recoVertexFoilIndexMatchedEntries << endl;
  cout << "  histogram reconstructed-vs-truth vertex residual entries: "
       << histograms.recoVertexTruthResidualEntries << endl;
  cout << "  histogram foil-index matched reconstructed-vs-truth vertex residual entries: "
       << histograms.recoVertexFoilIndexMatchedTruthResidualEntries << endl;
  cout << "  scatter all shared same-foil candidate pair max(L1,L2) vs delta-z points: "
       << histograms.recoAllSharedFoilCandidateMaxLvsDeltaZEntries << endl;
  cout << "  scatter space-selected shared same-foil max(L1,L2) vs delta-z points: "
       << histograms.recoSpaceSelectedSharedFoilMaxLvsDeltaZEntries << endl;
  cout << "  scatter all shared same-foil candidate pair (|s|+|t|)/2 vs delta-z points: "
       << histograms.recoAllSharedFoilCandidateAvgAbsLineParameterVsDeltaZEntries << endl;
  cout << "  scatter space-selected shared same-foil (|s|+|t|)/2 vs delta-z points: "
       << histograms.recoSpaceSelectedSharedFoilAvgAbsLineParameterVsDeltaZEntries << endl;
  cout << "  scatter space-selected shared foil sindex vs delta-z points: "
       << histograms.recoSpaceSelectedSharedFoilNumberVsDeltaZEntries << endl;
  cout << "  scatter space-selected number of shared foil indices vs delta-z points: "
       << histograms.recoSpaceSelectedSharedFoilCountVsDeltaZEntries << endl;
  cout << "  scatter space-selected max unique foils hit vs delta-z points: "
       << histograms.recoSpaceSelectedMaxFoilsHitVsDeltaZEntries << endl;
  cout << "  scatter space-selected absolute foil-count difference vs delta-z points: "
       << histograms.recoSpaceSelectedAbsDeltaFoilsHitVsDeltaZEntries << endl;
  cout << "  scatter space-selected opening angle vs delta-z points: "
       << histograms.recoSpaceSelectedOpeningAngleVsDeltaZEntries << endl;
  cout << "  histogram TEST minimum-time reconstructed vertex entries: "
       << histograms.testRecoVertexMinTimeEntries << endl;
  cout << "  histogram TEST minimum-time momentum theta-pair entries: "
       << histograms.testRecoVertexMinTimeMomentumThetaEntries << endl;
  cout << "  histogram TEST minimum-time momentum opening-angle entries: "
       << histograms.testRecoVertexMinTimeMomentumOpeningAngleEntries << endl;
  cout << "  histogram TEST foil-index matched minimum-time reconstructed vertex map entries: "
       << histograms.testRecoVertexMinTimeFoilIndexMatchedEntries << endl;
  cout << "  histogram TEST minimum-time vertex delta-t entries: "
       << histograms.recoVertexMinTimeDifferenceEntries << endl;
  cout << "  histogram TEST minimum-time reconstructed-vs-truth vertex residual entries: "
       << histograms.testRecoVertexMinTimeTruthResidualEntries << endl;
  cout << "  histogram TEST foil-index matched minimum-time reconstructed-vs-truth vertex residual entries: "
       << histograms.testRecoVertexMinTimeFoilIndexMatchedTruthResidualEntries << endl;
  cout << "  scatter TEST minimum-time shared same-foil max(L1,L2) vs delta-z points: "
       << histograms.testRecoMinTimeSharedFoilMaxLvsDeltaZEntries << endl;
  cout << "  scatter TEST minimum-time shared same-foil (|s|+|t|)/2 vs delta-z points: "
       << histograms.testRecoMinTimeSharedFoilAvgAbsLineParameterVsDeltaZEntries << endl;
  if (WRITE_HISTOGRAM_FILE)
  {
    cout << "  histogram output file: "
         << (histogramFileWasWritten ? HISTOGRAM_OUTPUT_FILE : "not written") << endl;
    cout << "  histogram PDF output directory: "
         << (histogramFileWasWritten ? HISTOGRAM_PLOT_OUTPUT_DIR : "not written")
         << endl;
  }
}

