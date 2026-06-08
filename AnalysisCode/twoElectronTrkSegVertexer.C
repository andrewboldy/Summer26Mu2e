//----------------------------------------------------------------------------------
//
// twoElectronTrkSegVertexer.C
// Written by Andrew Boldy University of South Carolina, 2026
// Assisted by Codex
//
// Purpose:
//   Read the reduced candidate ntuple produced by
//   twoElectronTrkSegAnalysisSelector.C and print the shared-surface state for
//   each selected event in time order.
//
//   The reduced tree is expected to contain only exactly-two-downstream-electron
//   candidate events with at least one shared surface.  For each event, this
//   macro prints:
//     - event identity and source entry
//     - the two selected track indices
//     - the reconstructed calorimeter energy associated with each track
//     - each shared surface, ordered by reconstructed time
//     - the reconstructed position, momentum, and time for both tracks at that
//       shared surface
//     - ST foil labels when the shared surface is on ST_Foils
//
//   In addition to the printout, this macro makes a tentative diagnostic
//   correlation plot of event-level closest ST_Foils line distance versus the
//   time difference between the two track segments at the selected shared
//   foil.  This is an exploratory setup and is labeled as TEST in the output
//   PDF name.
//
// Usage from ROOT:
//   .L CreatedCode/twoElectronTrkSegVertexer.C+
//   twoElectronTrkSegVertexer("reduced.root")
//   twoElectronTrkSegVertexer("reduced.root", -1, false, false)
//
//   The optional fourth argument controls whether ROOT displays canvases
//   interactively while saving the PDF files.  Set it to false for quieter,
//   lower-lag batch-style plotting.
//
//   The input can also be a filelist containing reduced ROOT files.
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
#include <sstream>
#include <string>
#include <vector>

#include <TChain.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMarker.h>
#include <TROOT.h>
#include <TVirtualPad.h>

#include "EventNtuple/inc/EventInfo.hh"
#include "EventNtuple/inc/TrkSegInfo.hh"

using namespace std;

namespace
{
  // Stopping-target overlay geometry, expressed in the same tracker/detector
  // coordinates as the vertex maps:
  //   - EventNtuple/inc/TrkSegInfo.hh documents TrkSegInfo::pos as positions
  //     in mm WRT the tracker center.
  //   - EventNtuple/src/InfoStructHelper.cc stores kinter.position3()
  //     directly into TrkSegInfo::pos in fillTrkSegInfo().
  //   - Offline/GeometryService/src/StoppingTargetMaker.cc says stopping
  //     target positions are in detector coordinates and builds TargetFoil
  //     centers from the geometry config plus the DetectorSystem origin.
  //   - Offline/Mu2eG4/geom/geom_run2.txt sets mu2e.detectorSystemZ0=10171
  //     and includes stoppingTargetHoles_v02.txt.
  //   - stoppingTargetHoles_v02.txt includes
  //     stoppingTargetHoles_DOE_review_2017.txt, which sets
  //     stoppingTarget.holeRadius=21.5 mm, stoppingTarget.rOut=75 mm for the
  //     foils, and stoppingTarget.deltaZ=22.222222 mm.  The base target file
  //     stoppingTarget_CD3C_34foils.txt sets stoppingTarget.z0InMu2e=5871 mm,
  //     while stoppingTargetHoles_v02.txt sets the foil half-thickness to
  //     0.0528 mm.
  //
  // For geom_run2, detector/tracker z = Mu2e z - 10171 mm.  With 37 foils,
  // these config values give foil centers from about -4700 to -3900 mm and
  // physical foil faces from -4700.053 to -3899.947 mm.  The XY annulus is
  // centered at (0,0) because the same geometry config leaves the foil x/y
  // offsets at 0.
  const double kStoppingTargetCenterX = 0.0;
  const double kStoppingTargetCenterY = 0.0;
  const double kStoppingTargetHoleRadius = 21.5;
  const double kStoppingTargetOuterRadius = 75.0;
  const double kStoppingTargetXYMin = -kStoppingTargetOuterRadius;
  const double kStoppingTargetXYMax = kStoppingTargetOuterRadius;
  const double kStoppingTargetZMin = -4700.053;
  const double kStoppingTargetZMax = -3899.947;

  enum VertexProjection
  {
    kVertexXY,
    kVertexXZ,
    kVertexYZ
  };

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
    // The reduced vertex tree stores these as direct branches, not split
    // branch wildcards.  Enabling the exact branch name avoids noisy ROOT
    // warnings for branchName.* patterns that do not exist in this tree.
    if (chain.GetBranch(branchName.c_str()) != nullptr)
    {
      chain.SetBranchStatus(branchName.c_str(), 1);
    }
  }

  bool addInput(TChain& chain, const string& inputName)
  {
    // Support both a single reduced ROOT file and a plain text filelist.
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

  //============================================================================
  // Shared-Surface Printing Helpers
  //============================================================================

  struct SharedSurfacePrintOrder
  {
    size_t index = 0;
    double averageTime = 0.0;
  };

  // Build a stable print order from the mean time of the two track states at
  // the shared surface.  This keeps the output chronological even when the
  // stored branches are not ordered that way.
  vector<SharedSurfacePrintOrder> sortedSharedSurfaceIndicesByTime(
    const vector<double>& firstTimes,
    const vector<double>& secondTimes)
  {
    vector<SharedSurfacePrintOrder> ordered;
    const size_t nSharedSurfaces = min(firstTimes.size(), secondTimes.size());
    ordered.reserve(nSharedSurfaces);

    for (size_t i = 0; i < nSharedSurfaces; ++i)
    {
      SharedSurfacePrintOrder entry;
      entry.index = i;
      entry.averageTime = 0.5 * (firstTimes.at(i) + secondTimes.at(i));
      ordered.push_back(entry);
    }

    stable_sort(ordered.begin(),
                ordered.end(),
                [](const SharedSurfacePrintOrder& left,
                   const SharedSurfacePrintOrder& right)
                {
                  return left.averageTime < right.averageTime;
                });

    return ordered;
  }

  string formatSurfaceLabel(int sid, int sindex)
  {
    ostringstream out;
    // ST foils get a special label because the foil index matters downstream.
    if (sid == mu2e::SurfaceIdDetail::ST_Foils)
    {
      out << "ST_Foils[" << sindex << "]";
    }
    else
    {
      out << "sid=" << sid << ", sindex=" << sindex;
    }
    return out.str();
  }

  string formatVector3(double x, double y, double z)
  {
    // Keep vector formatting in one place so the printout is easy to change.
    ostringstream out;
    out << "(" << x << ", " << y << ", " << z << ")";
    return out.str();
  }

  // One reconstructed line in 3D.  The label "3D line" is represented here as
  // Line3D because C++ identifiers cannot begin with a digit.
  struct Line3D
  {
    XYZVectorF posOnSurface = XYZVectorF();
    XYZVectorF unitMom = XYZVectorF();
  };

  // Build a line from a surface position and a momentum vector.  The position
  // is kept as-is; the momentum is normalized so downstream vertexing can work
  // from a direction vector instead of a full momentum scale.
  Line3D makeLine3D(const XYZVectorF& posOnSurface, const XYZVectorF& mom)
  {
    Line3D line;
    line.posOnSurface = posOnSurface;
    const double momMag = mom.R();
    if (momMag > 0.0)
    {
      line.unitMom = mom / momMag;
    }
    else
    {
      line.unitMom = XYZVectorF(0.0, 0.0, 0.0);
    }
    return line;
  }

  // Distance between two 3D lines using the standard skew-line formula:
  //
  //   d = |(r2 - r1) · (u1 x u2)| / |u1 x u2|
  //
  // where r1/r2 are points on the lines and u1/u2 are unit direction vectors.
  // If the lines are parallel or nearly parallel, the denominator vanishes and
  // the distance is not numerically stable.  In that case we return -1.0 as a
  // sentinel and let the printout label it as undefined.
  double lineLineDistance(const Line3D& line1, const Line3D& line2)
  {
    const XYZVectorF delta = line2.posOnSurface - line1.posOnSurface;
    const XYZVectorF cross = line1.unitMom.Cross(line2.unitMom);
    const double crossMag = cross.R();

    if (crossMag <= 0.0)
    {
      return -1.0;
    }

    return fabs(delta.Dot(cross)) / crossMag;
  }

  // Full closest-point separation between two lines.  This gives the signed
  // x/y/z offset from the closest point on line1 to the closest point on line2.
  // The scalar distance is the magnitude of that separation vector.
  struct LineSeparation3D
  {
    bool valid = false;
    XYZVectorF closestPointOnFirst = XYZVectorF();
    XYZVectorF closestPointOnSecond = XYZVectorF();
    XYZVectorF separation = XYZVectorF();
    double distance = -1.0;
  };

  LineSeparation3D lineLineSeparation(const Line3D& line1, const Line3D& line2)
  {
    LineSeparation3D result;

    const XYZVectorF w0 = line1.posOnSurface - line2.posOnSurface;
    const double a = line1.unitMom.Dot(line1.unitMom);
    const double b = line1.unitMom.Dot(line2.unitMom);
    const double c = line2.unitMom.Dot(line2.unitMom);
    const double d = line1.unitMom.Dot(w0);
    const double e = line2.unitMom.Dot(w0);
    const double denom = a * c - b * b;

    if (fabs(denom) <= 1e-12)
    {
      return result;
    }

    const double sc = (b * e - c * d) / denom;
    const double tc = (a * e - b * d) / denom;
    result.closestPointOnFirst = line1.posOnSurface + sc * line1.unitMom;
    result.closestPointOnSecond = line2.posOnSurface + tc * line2.unitMom;
    result.separation = result.closestPointOnSecond - result.closestPointOnFirst;
    result.distance = result.separation.R();
    result.valid = true;
    return result;
  }

  XYZVectorF midpointBetweenClosestApproachPoints(const LineSeparation3D& lineSeparation)
  {
    // The closest points are built from TrkSegInfo surface positions, which are
    // already stored in tracker/detector coordinates.  This midpoint therefore
    // stays in the same coordinate frame as the stopping-target overlay below.
    return 0.5 * (lineSeparation.closestPointOnFirst +
                  lineSeparation.closestPointOnSecond);
  }

  void drawStoppingTargetBox(VertexProjection projection)
  {
    if (projection == kVertexXY)
    {
      TEllipse outerFoilEdge(kStoppingTargetCenterX,
                             kStoppingTargetCenterY,
                             kStoppingTargetOuterRadius,
                             kStoppingTargetOuterRadius);
      outerFoilEdge.SetFillStyle(0);
      outerFoilEdge.SetLineColor(kRed + 1);
      outerFoilEdge.SetLineWidth(3);
      outerFoilEdge.DrawClone("same");

      TEllipse innerFoilHole(kStoppingTargetCenterX,
                             kStoppingTargetCenterY,
                             kStoppingTargetHoleRadius,
                             kStoppingTargetHoleRadius);
      innerFoilHole.SetFillStyle(0);
      innerFoilHole.SetLineColor(kRed + 1);
      innerFoilHole.SetLineWidth(3);
      innerFoilHole.SetLineStyle(2);
      innerFoilHole.DrawClone("same");

      TMarker targetCenter(kStoppingTargetCenterX, kStoppingTargetCenterY, 2);
      targetCenter.SetMarkerColor(kRed + 1);
      targetCenter.SetMarkerSize(1.2);
      targetCenter.DrawClone("same");
      return;
    }

    double boxXMin = kStoppingTargetXYMin;
    double boxXMax = kStoppingTargetXYMax;
    double boxYMin = kStoppingTargetXYMin;
    double boxYMax = kStoppingTargetXYMax;

    // For XZ and YZ maps, the red box spans the foil's outer transverse
    // radius in the horizontal coordinate and the physical foil z-face bounds
    // in tracker/detector z.  Those numbers come from the geometry files
    // listed in the constants block above.
    if (projection == kVertexXZ || projection == kVertexYZ)
    {
      boxYMin = kStoppingTargetZMin;
      boxYMax = kStoppingTargetZMax;
    }

    TBox targetBox(boxXMin, boxYMin, boxXMax, boxYMax);
    targetBox.SetFillStyle(0);
    targetBox.SetLineColor(kRed + 1);
    targetBox.SetLineWidth(3);
    targetBox.DrawClone("same");
  }

  void draw2DHistogram(const string& canvasName,
                       const string& canvasTitle,
                       TH2F* histogram,
                       VertexProjection targetProjection,
                       const string& outputPath)
  {
    TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 900, 700);
    histogram->Draw("COLZ");
    drawStoppingTargetBox(targetProjection);
    canvas->SaveAs(outputPath.c_str());
  }

  void drawHistogramPad(TCanvas* canvas, int padNumber, TH1F* histogram, bool useLogY)
  {
    TVirtualPad* pad = canvas->cd(padNumber);
    if (pad != nullptr && useLogY)
    {
      pad->SetLogy();
      if (histogram->GetMaximum() > 0.0)
      {
        histogram->SetMinimum(0.5);
      }
    }

    histogram->SetLineWidth(2);
    histogram->Draw("HIST");
  }

  void drawSeparationComponentCanvas(const string& canvasName,
                                     const string& canvasTitle,
                                     TH1F* hDx,
                                     TH1F* hDy,
                                     TH1F* hDz,
                                     TH1F* hDistance,
                                     bool useLogY,
                                     const string& outputPath)
  {
    TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 1200, 900);
    canvas->Divide(2, 2);

    drawHistogramPad(canvas, 1, hDx, useLogY);
    drawHistogramPad(canvas, 2, hDy, useLogY);
    drawHistogramPad(canvas, 3, hDz, useLogY);
    drawHistogramPad(canvas, 4, hDistance, useLogY);

    canvas->SaveAs(outputPath.c_str());
  }

  bool passesTrackMomentumCut(double firstTrackMomentum,
                              double secondTrackMomentum,
                              double cutMin,
                              double cutMax)
  {
    return firstTrackMomentum >= cutMin &&
           firstTrackMomentum <= cutMax &&
           secondTrackMomentum >= cutMin &&
           secondTrackMomentum <= cutMax;
  }
}

void twoElectronTrkSegVertexer(const string& inputName,
                               int maxEvents = -1,
                               bool doFullPrintout = true,
                               bool displayCanvases = true)
{
  // The reduced tree created by the selector stores one row per selected event.
  const string treeName = "TwoElectronTrackVertexCandidates";

  // Read the reduced tree with a TChain so the vertexer can accept either a
  // single file or a filelist.
  TChain ntuple(treeName.c_str());
  if (!addInput(ntuple, inputName))
  {
    return;
  }

  // Quick sanity check before wiring branches.  A zero-entry tree usually
  // means the input path is wrong or the selector did not write candidates.
  const Long64_t nEntries = ntuple.GetEntries();
  if (nEntries <= 0)
  {
    cerr << "ERROR: no entries found in " << treeName << " for input: " << inputName << endl;
    return;
  }

  // Disable everything first, then turn on only the reduced candidate fields
  // that this vertexer actually reads.
  ntuple.SetBranchStatus("*", 0);
  enableBranch(ntuple, "evtinfo");
  enableBranch(ntuple, "sourceEntry");
  enableBranch(ntuple, "selectedTrackIndices");
  enableBranch(ntuple, "selectedTrackPdg");
  enableBranch(ntuple, "selectedTrackNhits");
  enableBranch(ntuple, "selectedTrackStatus");
  enableBranch(ntuple, "selectedTrackGoodfit");
  enableBranch(ntuple, "selectedTrackNdof");
  enableBranch(ntuple, "selectedTrackNactive");
  enableBranch(ntuple, "selectedTrackNseg");
  enableBranch(ntuple, "selectedTrackNipadown");
  enableBranch(ntuple, "selectedTrackNstdown");
  enableBranch(ntuple, "selectedTrackFirstStInter");
  enableBranch(ntuple, "selectedTrackNstup");
  enableBranch(ntuple, "selectedTrackNipaup");
  enableBranch(ntuple, "selectedTrackFitcon");
  enableBranch(ntuple, "selectedTrackChisq");
  enableBranch(ntuple, "selectedTrackCaloActive");
  enableBranch(ntuple, "selectedTrackCaloDid");
  enableBranch(ntuple, "selectedTrackCaloPOCAX");
  enableBranch(ntuple, "selectedTrackCaloPOCAY");
  enableBranch(ntuple, "selectedTrackCaloPOCAZ");
  enableBranch(ntuple, "selectedTrackCaloMomX");
  enableBranch(ntuple, "selectedTrackCaloMomY");
  enableBranch(ntuple, "selectedTrackCaloMomZ");
  enableBranch(ntuple, "selectedTrackCaloCDepth");
  enableBranch(ntuple, "selectedTrackCaloTrkDepth");
  enableBranch(ntuple, "selectedTrackCaloDphiDot");
  enableBranch(ntuple, "selectedTrackCaloDoca");
  enableBranch(ntuple, "selectedTrackCaloDt");
  enableBranch(ntuple, "selectedTrackCaloPtoca");
  enableBranch(ntuple, "selectedTrackCaloTocavar");
  enableBranch(ntuple, "selectedTrackCaloTresid");
  enableBranch(ntuple, "selectedTrackCaloTresidmvar");
  enableBranch(ntuple, "selectedTrackCaloTresidpvar");
  enableBranch(ntuple, "selectedTrackCaloCSize");
  enableBranch(ntuple, "selectedTrackCaloEdep");
  enableBranch(ntuple, "selectedTrackCaloEdepErr");
  enableBranch(ntuple, "selectedTrackTotalSegmentCount");
  enableBranch(ntuple, "selectedTrackFoilSegmentCount");
  enableBranch(ntuple, "sharedSurfaceCount");
  enableBranch(ntuple, "sharedSTFoilCount");
  enableBranch(ntuple, "sharedSurfaceSid");
  enableBranch(ntuple, "sharedSurfaceSindex");
  enableBranch(ntuple, "sharedSurfaceFirstStoredSegmentIndex");
  enableBranch(ntuple, "sharedSurfaceSecondStoredSegmentIndex");
  enableBranch(ntuple, "sharedSurfaceFirstPosX");
  enableBranch(ntuple, "sharedSurfaceFirstPosY");
  enableBranch(ntuple, "sharedSurfaceFirstPosZ");
  enableBranch(ntuple, "sharedSurfaceFirstMomX");
  enableBranch(ntuple, "sharedSurfaceFirstMomY");
  enableBranch(ntuple, "sharedSurfaceFirstMomZ");
  enableBranch(ntuple, "sharedSurfaceFirstTime");
  enableBranch(ntuple, "sharedSurfaceSecondPosX");
  enableBranch(ntuple, "sharedSurfaceSecondPosY");
  enableBranch(ntuple, "sharedSurfaceSecondPosZ");
  enableBranch(ntuple, "sharedSurfaceSecondMomX");
  enableBranch(ntuple, "sharedSurfaceSecondMomY");
  enableBranch(ntuple, "sharedSurfaceSecondMomZ");
  enableBranch(ntuple, "sharedSurfaceSecondTime");

  // Pointers for the branches that ROOT will populate on each GetEntry call.
  // The reduced tree is deliberately flat, so each field is read back into a
  // direct pointer rather than a nested event object.
  mu2e::EventInfo* evtinfo = nullptr;
  Long64_t sourceEntry = -1;
  vector<int>* selectedTrackIndices = nullptr;
  vector<int>* selectedTrackPdg = nullptr;
  vector<int>* selectedTrackNhits = nullptr;
  vector<int>* selectedTrackStatus = nullptr;
  vector<int>* selectedTrackGoodfit = nullptr;
  vector<int>* selectedTrackNdof = nullptr;
  vector<int>* selectedTrackNactive = nullptr;
  vector<int>* selectedTrackNseg = nullptr;
  vector<int>* selectedTrackNipadown = nullptr;
  vector<int>* selectedTrackNstdown = nullptr;
  vector<int>* selectedTrackFirstStInter = nullptr;
  vector<int>* selectedTrackNstup = nullptr;
  vector<int>* selectedTrackNipaup = nullptr;
  vector<double>* selectedTrackFitcon = nullptr;
  vector<double>* selectedTrackChisq = nullptr;
  vector<int>* selectedTrackCaloActive = nullptr;
  vector<int>* selectedTrackCaloDid = nullptr;
  vector<double>* selectedTrackCaloPOCAX = nullptr;
  vector<double>* selectedTrackCaloPOCAY = nullptr;
  vector<double>* selectedTrackCaloPOCAZ = nullptr;
  vector<double>* selectedTrackCaloMomX = nullptr;
  vector<double>* selectedTrackCaloMomY = nullptr;
  vector<double>* selectedTrackCaloMomZ = nullptr;
  vector<double>* selectedTrackCaloCDepth = nullptr;
  vector<double>* selectedTrackCaloTrkDepth = nullptr;
  vector<double>* selectedTrackCaloDphiDot = nullptr;
  vector<double>* selectedTrackCaloDoca = nullptr;
  vector<double>* selectedTrackCaloDt = nullptr;
  vector<double>* selectedTrackCaloPtoca = nullptr;
  vector<double>* selectedTrackCaloTocavar = nullptr;
  vector<double>* selectedTrackCaloTresid = nullptr;
  vector<double>* selectedTrackCaloTresidmvar = nullptr;
  vector<double>* selectedTrackCaloTresidpvar = nullptr;
  vector<double>* selectedTrackCaloCSize = nullptr;
  vector<double>* selectedTrackCaloEdep = nullptr;
  vector<double>* selectedTrackCaloEdepErr = nullptr;
  vector<int>* selectedTrackTotalSegmentCount = nullptr;
  vector<int>* selectedTrackFoilSegmentCount = nullptr;
  int sharedSurfaceCount = 0;
  int sharedSTFoilCount = 0;
  vector<int>* sharedSurfaceSid = nullptr;
  vector<int>* sharedSurfaceSindex = nullptr;
  vector<int>* sharedSurfaceFirstStoredSegmentIndex = nullptr;
  vector<int>* sharedSurfaceSecondStoredSegmentIndex = nullptr;
  vector<double>* sharedSurfaceFirstPosX = nullptr;
  vector<double>* sharedSurfaceFirstPosY = nullptr;
  vector<double>* sharedSurfaceFirstPosZ = nullptr;
  vector<double>* sharedSurfaceFirstMomX = nullptr;
  vector<double>* sharedSurfaceFirstMomY = nullptr;
  vector<double>* sharedSurfaceFirstMomZ = nullptr;
  vector<double>* sharedSurfaceFirstTime = nullptr;
  vector<double>* sharedSurfaceSecondPosX = nullptr;
  vector<double>* sharedSurfaceSecondPosY = nullptr;
  vector<double>* sharedSurfaceSecondPosZ = nullptr;
  vector<double>* sharedSurfaceSecondMomX = nullptr;
  vector<double>* sharedSurfaceSecondMomY = nullptr;
  vector<double>* sharedSurfaceSecondMomZ = nullptr;
  vector<double>* sharedSurfaceSecondTime = nullptr;

  // Bind the local pointers to the reduced-tree branches.
  ntuple.SetBranchAddress("evtinfo", &evtinfo);
  ntuple.SetBranchAddress("sourceEntry", &sourceEntry);
  ntuple.SetBranchAddress("selectedTrackIndices", &selectedTrackIndices);
  ntuple.SetBranchAddress("selectedTrackPdg", &selectedTrackPdg);
  ntuple.SetBranchAddress("selectedTrackNhits", &selectedTrackNhits);
  ntuple.SetBranchAddress("selectedTrackStatus", &selectedTrackStatus);
  ntuple.SetBranchAddress("selectedTrackGoodfit", &selectedTrackGoodfit);
  ntuple.SetBranchAddress("selectedTrackNdof", &selectedTrackNdof);
  ntuple.SetBranchAddress("selectedTrackNactive", &selectedTrackNactive);
  ntuple.SetBranchAddress("selectedTrackNseg", &selectedTrackNseg);
  ntuple.SetBranchAddress("selectedTrackNipadown", &selectedTrackNipadown);
  ntuple.SetBranchAddress("selectedTrackNstdown", &selectedTrackNstdown);
  ntuple.SetBranchAddress("selectedTrackFirstStInter", &selectedTrackFirstStInter);
  ntuple.SetBranchAddress("selectedTrackNstup", &selectedTrackNstup);
  ntuple.SetBranchAddress("selectedTrackNipaup", &selectedTrackNipaup);
  ntuple.SetBranchAddress("selectedTrackFitcon", &selectedTrackFitcon);
  ntuple.SetBranchAddress("selectedTrackChisq", &selectedTrackChisq);
  ntuple.SetBranchAddress("selectedTrackCaloActive", &selectedTrackCaloActive);
  ntuple.SetBranchAddress("selectedTrackCaloDid", &selectedTrackCaloDid);
  ntuple.SetBranchAddress("selectedTrackCaloPOCAX", &selectedTrackCaloPOCAX);
  ntuple.SetBranchAddress("selectedTrackCaloPOCAY", &selectedTrackCaloPOCAY);
  ntuple.SetBranchAddress("selectedTrackCaloPOCAZ", &selectedTrackCaloPOCAZ);
  ntuple.SetBranchAddress("selectedTrackCaloMomX", &selectedTrackCaloMomX);
  ntuple.SetBranchAddress("selectedTrackCaloMomY", &selectedTrackCaloMomY);
  ntuple.SetBranchAddress("selectedTrackCaloMomZ", &selectedTrackCaloMomZ);
  ntuple.SetBranchAddress("selectedTrackCaloCDepth", &selectedTrackCaloCDepth);
  ntuple.SetBranchAddress("selectedTrackCaloTrkDepth", &selectedTrackCaloTrkDepth);
  ntuple.SetBranchAddress("selectedTrackCaloDphiDot", &selectedTrackCaloDphiDot);
  ntuple.SetBranchAddress("selectedTrackCaloDoca", &selectedTrackCaloDoca);
  ntuple.SetBranchAddress("selectedTrackCaloDt", &selectedTrackCaloDt);
  ntuple.SetBranchAddress("selectedTrackCaloPtoca", &selectedTrackCaloPtoca);
  ntuple.SetBranchAddress("selectedTrackCaloTocavar", &selectedTrackCaloTocavar);
  ntuple.SetBranchAddress("selectedTrackCaloTresid", &selectedTrackCaloTresid);
  ntuple.SetBranchAddress("selectedTrackCaloTresidmvar", &selectedTrackCaloTresidmvar);
  ntuple.SetBranchAddress("selectedTrackCaloTresidpvar", &selectedTrackCaloTresidpvar);
  ntuple.SetBranchAddress("selectedTrackCaloCSize", &selectedTrackCaloCSize);
  ntuple.SetBranchAddress("selectedTrackCaloEdep", &selectedTrackCaloEdep);
  ntuple.SetBranchAddress("selectedTrackCaloEdepErr", &selectedTrackCaloEdepErr);
  ntuple.SetBranchAddress("selectedTrackTotalSegmentCount", &selectedTrackTotalSegmentCount);
  ntuple.SetBranchAddress("selectedTrackFoilSegmentCount", &selectedTrackFoilSegmentCount);
  ntuple.SetBranchAddress("sharedSurfaceCount", &sharedSurfaceCount);
  ntuple.SetBranchAddress("sharedSTFoilCount", &sharedSTFoilCount);
  ntuple.SetBranchAddress("sharedSurfaceSid", &sharedSurfaceSid);
  ntuple.SetBranchAddress("sharedSurfaceSindex", &sharedSurfaceSindex);
  ntuple.SetBranchAddress("sharedSurfaceFirstStoredSegmentIndex",
                          &sharedSurfaceFirstStoredSegmentIndex);
  ntuple.SetBranchAddress("sharedSurfaceSecondStoredSegmentIndex",
                          &sharedSurfaceSecondStoredSegmentIndex);
  ntuple.SetBranchAddress("sharedSurfaceFirstPosX", &sharedSurfaceFirstPosX);
  ntuple.SetBranchAddress("sharedSurfaceFirstPosY", &sharedSurfaceFirstPosY);
  ntuple.SetBranchAddress("sharedSurfaceFirstPosZ", &sharedSurfaceFirstPosZ);
  ntuple.SetBranchAddress("sharedSurfaceFirstMomX", &sharedSurfaceFirstMomX);
  ntuple.SetBranchAddress("sharedSurfaceFirstMomY", &sharedSurfaceFirstMomY);
  ntuple.SetBranchAddress("sharedSurfaceFirstMomZ", &sharedSurfaceFirstMomZ);
  ntuple.SetBranchAddress("sharedSurfaceFirstTime", &sharedSurfaceFirstTime);
  ntuple.SetBranchAddress("sharedSurfaceSecondPosX", &sharedSurfaceSecondPosX);
  ntuple.SetBranchAddress("sharedSurfaceSecondPosY", &sharedSurfaceSecondPosY);
  ntuple.SetBranchAddress("sharedSurfaceSecondPosZ", &sharedSurfaceSecondPosZ);
  ntuple.SetBranchAddress("sharedSurfaceSecondMomX", &sharedSurfaceSecondMomX);
  ntuple.SetBranchAddress("sharedSurfaceSecondMomY", &sharedSurfaceSecondMomY);
  ntuple.SetBranchAddress("sharedSurfaceSecondMomZ", &sharedSurfaceSecondMomZ);
  ntuple.SetBranchAddress("sharedSurfaceSecondTime", &sharedSurfaceSecondTime);

  // Respect maxEvents if the caller wants a short test run.
  const Long64_t entriesToRead =
    (maxEvents >= 0 && static_cast<Long64_t>(maxEvents) < nEntries) ? maxEvents : nEntries;

  // Keep floating-point output compact and consistent across the whole macro.
  cout << fixed << setprecision(3);
  cout << "Input: " << inputName << endl;
  cout << "Tree: " << treeName << endl;
  cout << "Tree entries available: " << nEntries << endl;
  cout << "Tree entries being scanned: " << entriesToRead << endl;
  cout << "Output mode: "
       << (doFullPrintout ? "shared-surface vertex printout" : "summary only")
       << endl;
  const bool originalRootBatchMode = gROOT->IsBatch();
  gROOT->SetBatch(!displayCanvases);
  cout << "Canvas display: "
       << (displayCanvases ? "enabled" : "disabled; PDFs are still saved")
       << endl;

  // Summary counters for the final footer.
  Long64_t eventsPrinted = 0;
  Long64_t sharedSurfacesPrinted = 0;
  Long64_t sharedSTFoilsPrinted = 0;
  // Event-level ST-foil minima are what we summarize at the end of the scan.
  Long64_t eventsWithSharedSTFoils = 0;
  Long64_t eventsWithSharedSTFoilsMomentumCut = 0;
  Long64_t eventsWithMinDistanceUnder1mm = 0;
  Long64_t eventsWithMinDistanceUnder0p5mm = 0;
  Long64_t eventsWithMinDistanceUnder0p1mm = 0;
  map<int, Long64_t> eventsBySharedFoilMultiplicity;
  double smallestEventMinSTFoilDistance = numeric_limits<double>::infinity();
  double largestEventMinSTFoilDistance = -numeric_limits<double>::infinity();
  const double trackMomentumCutMin = 50.0;
  const double trackMomentumCutMax = 53.0;
  const double momentumCutReferenceEventCount = 100000.0;
  const double vertexXYMin = -600.0;
  const double vertexXYMax = 600.0;
  const double vertexZMin = -5000.0;
  const double vertexZMax = -3000.0;
  const int vertexXYBins = 200;
  const int vertexZBins = 405;
  TH1F* hEventMinSTFoilDistance = new TH1F(
    "hEventMinSTFoilDistance",
    "Event-level closest shared ST_Foils line distance;closest shared ST_Foils line distance [mm];Events",
    200, 0.0, 20.0);
  TH1F* hAllSTFoilLineDistance = new TH1F(
    "hAllSTFoilLineDistance",
    "All shared ST_Foils line distances;shared ST_Foils line distance [mm];Pairs",
    200, 0.0, 20.0);
  // Exploratory component-level separation histograms.  These are kept
  // alongside the scalar line-distance plots so we can inspect the geometry in
  // x, y, z, and transverse radius without changing the existing resolution
  // study.
  TH1F* hEventMinSTFoilLineDx = new TH1F(
    "hEventMinSTFoilLineDx",
    "Event-level closest shared ST_Foils line dx;line_dx [mm];Events",
    200, -20.0, 20.0);
  TH1F* hEventMinSTFoilLineDy = new TH1F(
    "hEventMinSTFoilLineDy",
    "Event-level closest shared ST_Foils line dy;line_dy [mm];Events",
    200, -20.0, 20.0);
  TH1F* hEventMinSTFoilLineDz = new TH1F(
    "hEventMinSTFoilLineDz",
    "Event-level closest shared ST_Foils line dz;line_dz [mm];Events",
    200, -20.0, 20.0);
  TH1F* hEventMinSTFoilLineRxy = new TH1F(
    "hEventMinSTFoilLineRxy",
    "Event-level closest shared ST_Foils transverse separation;#sqrt{line_dx^{2}+line_dy^{2}} [mm];Events",
    200, 0.0, 20.0);
  TH1F* hAllSTFoilLineDx = new TH1F(
    "hAllSTFoilLineDx",
    "All shared ST_Foils line dx;line_dx [mm];Pairs",
    200, -20.0, 20.0);
  TH1F* hAllSTFoilLineDy = new TH1F(
    "hAllSTFoilLineDy",
    "All shared ST_Foils line dy;line_dy [mm];Pairs",
    200, -20.0, 20.0);
  TH1F* hAllSTFoilLineDz = new TH1F(
    "hAllSTFoilLineDz",
    "All shared ST_Foils line dz;line_dz [mm];Pairs",
    200, -20.0, 20.0);
  TH1F* hAllSTFoilLineRxy = new TH1F(
    "hAllSTFoilLineRxy",
    "All shared ST_Foils transverse separation;#sqrt{line_dx^{2}+line_dy^{2}} [mm];Pairs",
    200, 0.0, 20.0);
  TH2F* hEventMinSTFoilDistanceVsElectronTimeDiff = new TH2F(
    "hEventMinSTFoilDistanceVsElectronTimeDiff_TEST",
    "TEST: event-level closest shared ST_Foils line distance vs absolute shared-foil track-segment time difference;|track-segment time 1 - track-segment time 2| [ns];closest shared ST_Foils line distance [mm]",
    200, 0.0, 600.0,
    200, 0.0, 20.0);
  TH2F* hEventMinSTFoilVertexXY = new TH2F(
    "hEventMinSTFoilVertexXY",
    "Event-level closest shared ST_Foils midpoint vertex XY;vertex x [mm];vertex y [mm]",
    vertexXYBins, vertexXYMin, vertexXYMax,
    vertexXYBins, vertexXYMin, vertexXYMax);
  TH2F* hEventMinSTFoilVertexXZ = new TH2F(
    "hEventMinSTFoilVertexXZ",
    "Event-level closest shared ST_Foils midpoint vertex XZ;vertex x [mm];vertex z [mm]",
    vertexXYBins, vertexXYMin, vertexXYMax,
    vertexZBins, vertexZMin, vertexZMax);
  TH2F* hEventMinSTFoilVertexYZ = new TH2F(
    "hEventMinSTFoilVertexYZ",
    "Event-level closest shared ST_Foils midpoint vertex YZ;vertex y [mm];vertex z [mm]",
    vertexXYBins, vertexXYMin, vertexXYMax,
    vertexZBins, vertexZMin, vertexZMax);
  TH1F* hEventMinSTFoilDistanceMomentumCut = new TH1F(
    "hEventMinSTFoilDistanceMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils line distance;closest shared ST_Foils line distance [mm];Events",
    200, 0.0, 20.0);
  TH1F* hAllSTFoilLineDistanceMomentumCut = new TH1F(
    "hAllSTFoilLineDistanceMomentumCut",
    "Momentum cut: all shared ST_Foils line distances;shared ST_Foils line distance [mm];Pairs",
    200, 0.0, 20.0);
  TH1F* hEventMinSTFoilLineDxMomentumCut = new TH1F(
    "hEventMinSTFoilLineDxMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils line dx;line_dx [mm];Events",
    200, -20.0, 20.0);
  TH1F* hEventMinSTFoilLineDyMomentumCut = new TH1F(
    "hEventMinSTFoilLineDyMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils line dy;line_dy [mm];Events",
    200, -20.0, 20.0);
  TH1F* hEventMinSTFoilLineDzMomentumCut = new TH1F(
    "hEventMinSTFoilLineDzMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils line dz;line_dz [mm];Events",
    200, -20.0, 20.0);
  TH1F* hEventMinSTFoilLineRxyMomentumCut = new TH1F(
    "hEventMinSTFoilLineRxyMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils transverse separation;#sqrt{line_dx^{2}+line_dy^{2}} [mm];Events",
    200, 0.0, 20.0);
  TH1F* hAllSTFoilLineDxMomentumCut = new TH1F(
    "hAllSTFoilLineDxMomentumCut",
    "Momentum cut: all shared ST_Foils line dx;line_dx [mm];Pairs",
    200, -20.0, 20.0);
  TH1F* hAllSTFoilLineDyMomentumCut = new TH1F(
    "hAllSTFoilLineDyMomentumCut",
    "Momentum cut: all shared ST_Foils line dy;line_dy [mm];Pairs",
    200, -20.0, 20.0);
  TH1F* hAllSTFoilLineDzMomentumCut = new TH1F(
    "hAllSTFoilLineDzMomentumCut",
    "Momentum cut: all shared ST_Foils line dz;line_dz [mm];Pairs",
    200, -20.0, 20.0);
  TH1F* hAllSTFoilLineRxyMomentumCut = new TH1F(
    "hAllSTFoilLineRxyMomentumCut",
    "Momentum cut: all shared ST_Foils transverse separation;#sqrt{line_dx^{2}+line_dy^{2}} [mm];Pairs",
    200, 0.0, 20.0);
  TH2F* hEventMinSTFoilDistanceVsElectronTimeDiffMomentumCut = new TH2F(
    "hEventMinSTFoilDistanceVsElectronTimeDiffMomentumCut_TEST",
    "TEST: Momentum cut: event-level closest shared ST_Foils line distance vs absolute shared-foil track-segment time difference;|track-segment time 1 - track-segment time 2| [ns];closest shared ST_Foils line distance [mm]",
    200, 0.0, 600.0,
    200, 0.0, 20.0);
  TH2F* hEventMinSTFoilVertexXYMomentumCut = new TH2F(
    "hEventMinSTFoilVertexXYMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils midpoint vertex XY;vertex x [mm];vertex y [mm]",
    vertexXYBins, vertexXYMin, vertexXYMax,
    vertexXYBins, vertexXYMin, vertexXYMax);
  TH2F* hEventMinSTFoilVertexXZMomentumCut = new TH2F(
    "hEventMinSTFoilVertexXZMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils midpoint vertex XZ;vertex x [mm];vertex z [mm]",
    vertexXYBins, vertexXYMin, vertexXYMax,
    vertexZBins, vertexZMin, vertexZMax);
  TH2F* hEventMinSTFoilVertexYZMomentumCut = new TH2F(
    "hEventMinSTFoilVertexYZMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils midpoint vertex YZ;vertex y [mm];vertex z [mm]",
    vertexXYBins, vertexXYMin, vertexXYMax,
    vertexZBins, vertexZMin, vertexZMax);

  for (Long64_t iEntry = 0; iEntry < entriesToRead; ++iEntry)
  {
    // Load one reduced-candidate event from the chain.
    ntuple.GetEntry(iEntry);

    // Guard the readout: if a required branch is missing, skip this row rather
    // than dereferencing a null pointer or walking malformed vectors.
    if (evtinfo == nullptr ||
        selectedTrackIndices == nullptr ||
        selectedTrackPdg == nullptr ||
        selectedTrackCaloEdep == nullptr ||
        selectedTrackCaloMomX == nullptr ||
        selectedTrackCaloMomY == nullptr ||
        selectedTrackCaloMomZ == nullptr ||
        sharedSurfaceSid == nullptr ||
        sharedSurfaceSindex == nullptr ||
        sharedSurfaceFirstPosX == nullptr ||
        sharedSurfaceFirstPosY == nullptr ||
        sharedSurfaceFirstPosZ == nullptr ||
        sharedSurfaceFirstMomX == nullptr ||
        sharedSurfaceFirstMomY == nullptr ||
        sharedSurfaceFirstMomZ == nullptr ||
        sharedSurfaceFirstTime == nullptr ||
        sharedSurfaceSecondPosX == nullptr ||
        sharedSurfaceSecondPosY == nullptr ||
        sharedSurfaceSecondPosZ == nullptr ||
        sharedSurfaceSecondMomX == nullptr ||
        sharedSurfaceSecondMomY == nullptr ||
        sharedSurfaceSecondMomZ == nullptr ||
        sharedSurfaceSecondTime == nullptr)
    {
      cerr << "WARNING: missing reduced-ntuple data at entry " << iEntry << "; skipping." << endl;
      continue;
    }

    // The selector is supposed to have already filtered to two-track
    // candidates.  These checks are just a defensive guard against malformed
    // input or stale files.
    if (selectedTrackIndices->size() < 2 ||
        selectedTrackCaloEdep->size() < 2 ||
        sharedSurfaceSid->size() != sharedSurfaceFirstTime->size() ||
        sharedSurfaceSid->size() != sharedSurfaceSecondTime->size())
    {
      cerr << "WARNING: malformed reduced candidate at entry " << iEntry << "; skipping." << endl;
      continue;
    }

    ++eventsPrinted;

    // Store the counts locally so the event header stays easy to read.
    const int nSharedSurfaces = sharedSurfaceCount;
    const int nSharedFoils = sharedSTFoilCount;

    sharedSurfacesPrinted += nSharedSurfaces;
    if (nSharedFoils > 0)
    {
      ++eventsBySharedFoilMultiplicity[nSharedFoils];
    }

    if (doFullPrintout)
    {
      // Event header.  This stays off in summary-only mode.
      cout << "\nVertex candidate event " << eventsPrinted
           << " | entry=" << iEntry;
      if (evtinfo != nullptr)
      {
        cout << " run=" << evtinfo->run
             << " subrun=" << evtinfo->subrun
             << " event=" << evtinfo->event;
      }
      cout << " sourceEntry=" << sourceEntry
           << " track_indices=(" << selectedTrackIndices->at(0)
           << ", " << selectedTrackIndices->at(1) << ")"
           << " shared_surface_count=" << nSharedSurfaces
           << " shared_st_foils_count=" << nSharedFoils
           << endl;

      // The track-level calo momentum is a convenient single-vector summary for
      // the line construction step and for debugging against calorimeter data.
      const double track0MomMag = sqrt(
        selectedTrackCaloMomX->at(0) * selectedTrackCaloMomX->at(0) +
        selectedTrackCaloMomY->at(0) * selectedTrackCaloMomY->at(0) +
        selectedTrackCaloMomZ->at(0) * selectedTrackCaloMomZ->at(0));
      const double track1MomMag = sqrt(
        selectedTrackCaloMomX->at(1) * selectedTrackCaloMomX->at(1) +
        selectedTrackCaloMomY->at(1) * selectedTrackCaloMomY->at(1) +
        selectedTrackCaloMomZ->at(1) * selectedTrackCaloMomZ->at(1));

      cout << "  track0: pdg=" << selectedTrackPdg->at(0)
           << " caloEdep=" << selectedTrackCaloEdep->at(0)
           << " caloMom=" << formatVector3(selectedTrackCaloMomX->at(0),
                                           selectedTrackCaloMomY->at(0),
                                           selectedTrackCaloMomZ->at(0))
           << " |p|=" << track0MomMag << " MeV/c"
           << endl;
      cout << "  track1: pdg=" << selectedTrackPdg->at(1)
           << " caloEdep=" << selectedTrackCaloEdep->at(1)
           << " caloMom=" << formatVector3(selectedTrackCaloMomX->at(1),
                                           selectedTrackCaloMomY->at(1),
                                           selectedTrackCaloMomZ->at(1))
           << " |p|=" << track1MomMag << " MeV/c"
           << endl;
    }

    // Sort the shared surfaces chronologically by the mean of the two track
    // times at that surface.  This is done even in summary-only mode because
    // the same ordered loop fills the distance and vertex histograms.
    const vector<SharedSurfacePrintOrder> orderedSharedSurfaces =
      sortedSharedSurfaceIndicesByTime(*sharedSurfaceFirstTime, *sharedSurfaceSecondTime);

    // Shared ST foil lines are stored explicitly so the event can later be fed
    // into a vertexing calculation without re-parsing the full printout.  One
    // event can have multiple shared foil surfaces, so we cache them locally as
    // line objects while scanning the event.
    vector<Line3D> firstTrackFoilLines;
    vector<Line3D> secondTrackFoilLines;
    // Track the closest ST-foil pair for this event.  The vertex study uses
    // the minimum separation across all shared foils, not every raw value.
    double eventMinSTFoilDistance = numeric_limits<double>::infinity();
    XYZVectorF eventMinSeparation = XYZVectorF();
    XYZVectorF eventMinVertex = XYZVectorF();
    int eventMinFoilSindex = -1;
    int eventMinFoilSid = -1;
    double eventMinFirstTime = numeric_limits<double>::quiet_NaN();
    double eventMinSecondTime = numeric_limits<double>::quiet_NaN();
    double eventMinSTFoilDistanceMomentumCut = numeric_limits<double>::infinity();
    XYZVectorF eventMinSeparationMomentumCut = XYZVectorF();
    XYZVectorF eventMinVertexMomentumCut = XYZVectorF();
    double eventMinFirstTimeMomentumCut = numeric_limits<double>::quiet_NaN();
    double eventMinSecondTimeMomentumCut = numeric_limits<double>::quiet_NaN();

    for (size_t iOrder = 0; iOrder < orderedSharedSurfaces.size(); ++iOrder)
    {
      // Pick out the stored row for this shared surface.
      const size_t iShared = orderedSharedSurfaces.at(iOrder).index;
      const int sid = sharedSurfaceSid->at(iShared);
      const int sindex = sharedSurfaceSindex->at(iShared);

      const string surfaceLabel = formatSurfaceLabel(sid, sindex);
      const double firstTime = sharedSurfaceFirstTime->at(iShared);
      const double secondTime = sharedSurfaceSecondTime->at(iShared);
      const double firstMomX = sharedSurfaceFirstMomX->at(iShared);
      const double firstMomY = sharedSurfaceFirstMomY->at(iShared);
      const double firstMomZ = sharedSurfaceFirstMomZ->at(iShared);
      const double secondMomX = sharedSurfaceSecondMomX->at(iShared);
      const double secondMomY = sharedSurfaceSecondMomY->at(iShared);
      const double secondMomZ = sharedSurfaceSecondMomZ->at(iShared);
      const double firstPosX = sharedSurfaceFirstPosX->at(iShared);
      const double firstPosY = sharedSurfaceFirstPosY->at(iShared);
      const double firstPosZ = sharedSurfaceFirstPosZ->at(iShared);
      const double secondPosX = sharedSurfaceSecondPosX->at(iShared);
      const double secondPosY = sharedSurfaceSecondPosY->at(iShared);
      const double secondPosZ = sharedSurfaceSecondPosZ->at(iShared);
      if (sid == mu2e::SurfaceIdDetail::ST_Foils)
      {
        // For foil surfaces, the line geometry is what we care about: the
        // position where each track crosses the foil and the direction of
        // that track at the crossing.
        const double firstMomMag = sqrt(firstMomX * firstMomX +
                                        firstMomY * firstMomY +
                                        firstMomZ * firstMomZ);
        const double secondMomMag = sqrt(secondMomX * secondMomX +
                                         secondMomY * secondMomY +
                                         secondMomZ * secondMomZ);
        const bool passesMomentumCut =
          passesTrackMomentumCut(firstMomMag,
                                 secondMomMag,
                                 trackMomentumCutMin,
                                 trackMomentumCutMax);

        // For shared foil surfaces, construct the 3D line objects from the
        // surface position and the unit momentum vector at that surface.
        const Line3D firstLine = makeLine3D(
          XYZVectorF(firstPosX, firstPosY, firstPosZ),
          XYZVectorF(firstMomX, firstMomY, firstMomZ));
        const Line3D secondLine = makeLine3D(
          XYZVectorF(secondPosX, secondPosY, secondPosZ),
          XYZVectorF(secondMomX, secondMomY, secondMomZ));
        const LineSeparation3D lineSeparation = lineLineSeparation(firstLine, secondLine);
        firstTrackFoilLines.push_back(firstLine);
        secondTrackFoilLines.push_back(secondLine);

        if (doFullPrintout)
        {
          // Print the signed component-wise closest-point separation together
          // with the scalar distance and the line geometry used to compute it.
          cout << "  shared_surface_time_order=" << setw(3) << iOrder
               << " surface=" << surfaceLabel
               << " avg_time=" << orderedSharedSurfaces.at(iOrder).averageTime << " ns"
               << " first_time=" << firstTime << " ns"
               << " first_pos=" << formatVector3(firstPosX, firstPosY, firstPosZ) << " mm"
               << " first_mom=" << formatVector3(firstMomX, firstMomY, firstMomZ)
               << " MeV/c first_|p|=" << firstMomMag << " MeV/c"
               << " first_energy=" << selectedTrackCaloEdep->at(0) << " MeV"
               << " second_time=" << secondTime << " ns"
               << " second_pos=" << formatVector3(secondPosX, secondPosY, secondPosZ) << " mm"
               << " second_mom=" << formatVector3(secondMomX, secondMomY, secondMomZ)
               << " MeV/c second_|p|=" << secondMomMag << " MeV/c"
               << " second_energy=" << selectedTrackCaloEdep->at(1) << " MeV"
               << " st_foil_index=" << sindex
               << " first_line_pos=" << formatVector3(firstLine.posOnSurface.x(),
                                                       firstLine.posOnSurface.y(),
                                                       firstLine.posOnSurface.z()) << " mm"
               << " first_line_unit_mom=" << formatVector3(firstLine.unitMom.x(),
                                                           firstLine.unitMom.y(),
                                                           firstLine.unitMom.z())
               << " second_line_pos=" << formatVector3(secondLine.posOnSurface.x(),
                                                       secondLine.posOnSurface.y(),
                                                       secondLine.posOnSurface.z()) << " mm"
               << " second_line_unit_mom=" << formatVector3(secondLine.unitMom.x(),
                                                            secondLine.unitMom.y(),
                                                            secondLine.unitMom.z());
        }
        if (lineSeparation.valid)
        {
          const XYZVectorF vertexMidpoint =
            midpointBetweenClosestApproachPoints(lineSeparation);
          if (doFullPrintout)
          {
            cout << " line_dx=" << lineSeparation.separation.x() << " mm"
                 << " line_dy=" << lineSeparation.separation.y() << " mm"
                 << " line_dz=" << lineSeparation.separation.z() << " mm"
                 << " line_distance=" << lineSeparation.distance << " mm";
          }

          // This histogram is the full population of ST-foil line separations,
          // not just the event-level minimum.  It is useful for understanding
          // the raw spread of the geometry before the minimum-pick step.
          hAllSTFoilLineDx->Fill(lineSeparation.separation.x());
          hAllSTFoilLineDy->Fill(lineSeparation.separation.y());
          hAllSTFoilLineDz->Fill(lineSeparation.separation.z());
          hAllSTFoilLineRxy->Fill(
            sqrt(lineSeparation.separation.x() * lineSeparation.separation.x() +
                 lineSeparation.separation.y() * lineSeparation.separation.y()));
          hAllSTFoilLineDistance->Fill(lineSeparation.distance);
          if (passesMomentumCut)
          {
            hAllSTFoilLineDxMomentumCut->Fill(lineSeparation.separation.x());
            hAllSTFoilLineDyMomentumCut->Fill(lineSeparation.separation.y());
            hAllSTFoilLineDzMomentumCut->Fill(lineSeparation.separation.z());
            hAllSTFoilLineRxyMomentumCut->Fill(
              sqrt(lineSeparation.separation.x() * lineSeparation.separation.x() +
                   lineSeparation.separation.y() * lineSeparation.separation.y()));
            hAllSTFoilLineDistanceMomentumCut->Fill(lineSeparation.distance);

            if (lineSeparation.distance < eventMinSTFoilDistanceMomentumCut)
            {
              eventMinSTFoilDistanceMomentumCut = lineSeparation.distance;
              eventMinSeparationMomentumCut = lineSeparation.separation;
              eventMinVertexMomentumCut = vertexMidpoint;
              eventMinFirstTimeMomentumCut = firstTime;
              eventMinSecondTimeMomentumCut = secondTime;
            }
          }

          if (lineSeparation.distance < eventMinSTFoilDistance)
          {
            eventMinSTFoilDistance = lineSeparation.distance;
            eventMinSeparation = lineSeparation.separation;
            eventMinVertex = vertexMidpoint;
            eventMinFoilSindex = sindex;
            eventMinFoilSid = sid;
            eventMinFirstTime = firstTime;
            eventMinSecondTime = secondTime;
          }
        }
        else
        {
          if (doFullPrintout)
          {
            cout << " line_dx=undefined"
                 << " line_dy=undefined"
                 << " line_dz=undefined"
                 << " line_distance=undefined_parallel";
          }
        }
        if (doFullPrintout)
        {
          cout << endl;
        }

        // Shared-ST-foil counters are only for the final summary, so they are
        // updated even when the verbose printout is suppressed.
        ++sharedSTFoilsPrinted;
      }
    }

    if (eventMinSTFoilDistanceMomentumCut < numeric_limits<double>::infinity())
    {
      ++eventsWithSharedSTFoilsMomentumCut;
      hEventMinSTFoilLineDxMomentumCut->Fill(eventMinSeparationMomentumCut.x());
      hEventMinSTFoilLineDyMomentumCut->Fill(eventMinSeparationMomentumCut.y());
      hEventMinSTFoilLineDzMomentumCut->Fill(eventMinSeparationMomentumCut.z());
      hEventMinSTFoilLineRxyMomentumCut->Fill(
        sqrt(eventMinSeparationMomentumCut.x() * eventMinSeparationMomentumCut.x() +
             eventMinSeparationMomentumCut.y() * eventMinSeparationMomentumCut.y()));
      hEventMinSTFoilDistanceMomentumCut->Fill(eventMinSTFoilDistanceMomentumCut);
      hEventMinSTFoilDistanceVsElectronTimeDiffMomentumCut->Fill(
        fabs(eventMinFirstTimeMomentumCut - eventMinSecondTimeMomentumCut),
        eventMinSTFoilDistanceMomentumCut);
      hEventMinSTFoilVertexXYMomentumCut->Fill(eventMinVertexMomentumCut.x(),
                                               eventMinVertexMomentumCut.y());
      hEventMinSTFoilVertexXZMomentumCut->Fill(eventMinVertexMomentumCut.x(),
                                               eventMinVertexMomentumCut.z());
      hEventMinSTFoilVertexYZMomentumCut->Fill(eventMinVertexMomentumCut.y(),
                                               eventMinVertexMomentumCut.z());
    }

    // Keep the foil-line caches alive through this block so it is obvious
    // they are intentionally gathered as part of the vertexing workflow.
    // They are not consumed yet, but they define the exact data the next
    // stage will use.
    (void)firstTrackFoilLines;
    (void)secondTrackFoilLines;

    // Update the scan-wide min/max from the event-level minimum separation.
    // That keeps the summary focused on the best ST-foil candidate per event.
    if (eventMinSTFoilDistance < numeric_limits<double>::infinity())
    {
      ++eventsWithSharedSTFoils;
      hEventMinSTFoilLineDx->Fill(eventMinSeparation.x());
      hEventMinSTFoilLineDy->Fill(eventMinSeparation.y());
      hEventMinSTFoilLineDz->Fill(eventMinSeparation.z());
      hEventMinSTFoilLineRxy->Fill(
        sqrt(eventMinSeparation.x() * eventMinSeparation.x() +
             eventMinSeparation.y() * eventMinSeparation.y()));
      hEventMinSTFoilDistance->Fill(eventMinSTFoilDistance);
      hEventMinSTFoilDistanceVsElectronTimeDiff->Fill(
        fabs(eventMinFirstTime - eventMinSecondTime),
        eventMinSTFoilDistance);
      hEventMinSTFoilVertexXY->Fill(eventMinVertex.x(), eventMinVertex.y());
      hEventMinSTFoilVertexXZ->Fill(eventMinVertex.x(), eventMinVertex.z());
      hEventMinSTFoilVertexYZ->Fill(eventMinVertex.y(), eventMinVertex.z());
      if (eventMinSTFoilDistance < 1.0)
      {
        ++eventsWithMinDistanceUnder1mm;
      }
      if (eventMinSTFoilDistance < 0.5)
      {
        ++eventsWithMinDistanceUnder0p5mm;
      }
      if (eventMinSTFoilDistance < 0.1)
      {
        ++eventsWithMinDistanceUnder0p1mm;
      }
      if (eventMinSTFoilDistance < smallestEventMinSTFoilDistance)
      {
        smallestEventMinSTFoilDistance = eventMinSTFoilDistance;
      }
      if (eventMinSTFoilDistance > largestEventMinSTFoilDistance)
      {
        largestEventMinSTFoilDistance = eventMinSTFoilDistance;
      }

      if (doFullPrintout)
      {
        cout << "  event_smallest_shared_ST_Foils_line_distance="
             << eventMinSTFoilDistance << " mm"
             << " event_smallest_line_dx=" << eventMinSeparation.x() << " mm"
             << " event_smallest_line_dy=" << eventMinSeparation.y() << " mm"
             << " event_smallest_line_dz=" << eventMinSeparation.z() << " mm"
             << " surface=" << formatSurfaceLabel(eventMinFoilSid, eventMinFoilSindex)
             << endl;
      }
    }
  }

  // Draw and save the event-level minimum-distance histogram.  This acts as a
  // direct resolution gauge for the vertexing separation study.
  TCanvas* cEventMinSTFoilDistance = new TCanvas(
    "cEventMinSTFoilDistance",
    "Event-level closest shared ST_Foils line distance",
    900, 700);
  hEventMinSTFoilDistance->SetLineWidth(2);
  hEventMinSTFoilDistance->Draw("HIST");
  cEventMinSTFoilDistance->SaveAs(
    "Plots/DistancePlots/twoElectronTrkSegVertexer_EventMinSTFoilDistance.pdf");

  TCanvas* cAllSTFoilLineDistance = new TCanvas(
    "cAllSTFoilLineDistance",
    "All shared ST_Foils line distances",
    900, 700);
  hAllSTFoilLineDistance->SetLineWidth(2);
  hAllSTFoilLineDistance->Draw("HIST");
  cAllSTFoilLineDistance->SaveAs(
    "Plots/DistancePlots/twoElectronTrkSegVertexer_AllSTFoilLineDistance.pdf");

  // Exploratory component plots.  These are diagnostic additions to the scalar
  // separation study.
  TCanvas* cEventMinSTFoilLineRxy = new TCanvas(
    "cEventMinSTFoilLineRxy",
    "Event-level closest shared ST_Foils transverse separation",
    900, 700);
  hEventMinSTFoilLineRxy->Draw("HIST");
  cEventMinSTFoilLineRxy->SaveAs(
    "Plots/DistancePlots/twoElectronTrkSegVertexer_EventMinSTFoilLineRxy.pdf");

  drawSeparationComponentCanvas(
    "cEventMinSTFoilSeparationComponents",
    "Event-level closest shared ST_Foils line separation components",
    hEventMinSTFoilLineDx,
    hEventMinSTFoilLineDy,
    hEventMinSTFoilLineDz,
    hEventMinSTFoilDistance,
    false,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_EventMinSTFoilSeparationComponents.pdf");

  drawSeparationComponentCanvas(
    "cEventMinSTFoilSeparationComponentsLogY",
    "Event-level closest shared ST_Foils line separation components (log-y)",
    hEventMinSTFoilLineDx,
    hEventMinSTFoilLineDy,
    hEventMinSTFoilLineDz,
    hEventMinSTFoilDistance,
    true,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_EventMinSTFoilSeparationComponents_LogY.pdf");

  TCanvas* cAllSTFoilLineRxy = new TCanvas(
    "cAllSTFoilLineRxy",
    "All shared ST_Foils transverse separation",
    900, 700);
  hAllSTFoilLineRxy->Draw("HIST");
  cAllSTFoilLineRxy->SaveAs(
    "Plots/DistancePlots/twoElectronTrkSegVertexer_AllSTFoilLineRxy.pdf");

  drawSeparationComponentCanvas(
    "cAllSTFoilLineSeparationComponents",
    "All shared ST_Foils line separation components",
    hAllSTFoilLineDx,
    hAllSTFoilLineDy,
    hAllSTFoilLineDz,
    hAllSTFoilLineDistance,
    false,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_AllSTFoilLineSeparationComponents.pdf");

  drawSeparationComponentCanvas(
    "cAllSTFoilLineSeparationComponentsLogY",
    "All shared ST_Foils line separation components (log-y)",
    hAllSTFoilLineDx,
    hAllSTFoilLineDy,
    hAllSTFoilLineDz,
    hAllSTFoilLineDistance,
    true,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_AllSTFoilLineSeparationComponents_LogY.pdf");

  TCanvas* cEventMinSTFoilDistanceVsElectronTimeDiff = new TCanvas(
    "cEventMinSTFoilDistanceVsElectronTimeDiff",
    "TEST: event-level closest shared ST_Foils line distance vs absolute shared-foil track-segment time difference",
    900, 700);
  hEventMinSTFoilDistanceVsElectronTimeDiff->Draw("COLZ");
  cEventMinSTFoilDistanceVsElectronTimeDiff->SaveAs(
    "Plots/DistancePlots/twoElectronTrkSegVertexer_TEST_EventMinSTFoilDistanceVsElectronTimeDiff.pdf");

  draw2DHistogram(
    "cEventMinSTFoilVertexXY",
    "Event-level closest shared ST_Foils midpoint vertex XY",
    hEventMinSTFoilVertexXY,
    kVertexXY,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_EventMinSTFoilVertexXY.pdf");

  draw2DHistogram(
    "cEventMinSTFoilVertexXZ",
    "Event-level closest shared ST_Foils midpoint vertex XZ",
    hEventMinSTFoilVertexXZ,
    kVertexXZ,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_EventMinSTFoilVertexXZ.pdf");

  draw2DHistogram(
    "cEventMinSTFoilVertexYZ",
    "Event-level closest shared ST_Foils midpoint vertex YZ",
    hEventMinSTFoilVertexYZ,
    kVertexYZ,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_EventMinSTFoilVertexYZ.pdf");

  TCanvas* cEventMinSTFoilDistanceMomentumCut = new TCanvas(
    "cEventMinSTFoilDistanceMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils line distance",
    900, 700);
  hEventMinSTFoilDistanceMomentumCut->SetLineWidth(2);
  hEventMinSTFoilDistanceMomentumCut->Draw("HIST");
  cEventMinSTFoilDistanceMomentumCut->SaveAs(
    "Plots/DistancePlots/twoElectronTrkSegVertexer_MomentumCut_EventMinSTFoilDistance.pdf");

  TCanvas* cAllSTFoilLineDistanceMomentumCut = new TCanvas(
    "cAllSTFoilLineDistanceMomentumCut",
    "Momentum cut: all shared ST_Foils line distances",
    900, 700);
  hAllSTFoilLineDistanceMomentumCut->SetLineWidth(2);
  hAllSTFoilLineDistanceMomentumCut->Draw("HIST");
  cAllSTFoilLineDistanceMomentumCut->SaveAs(
    "Plots/DistancePlots/twoElectronTrkSegVertexer_MomentumCut_AllSTFoilLineDistance.pdf");

  TCanvas* cEventMinSTFoilLineRxyMomentumCut = new TCanvas(
    "cEventMinSTFoilLineRxyMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils transverse separation",
    900, 700);
  hEventMinSTFoilLineRxyMomentumCut->Draw("HIST");
  cEventMinSTFoilLineRxyMomentumCut->SaveAs(
    "Plots/DistancePlots/twoElectronTrkSegVertexer_MomentumCut_EventMinSTFoilLineRxy.pdf");

  drawSeparationComponentCanvas(
    "cEventMinSTFoilSeparationComponentsMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils line separation components",
    hEventMinSTFoilLineDxMomentumCut,
    hEventMinSTFoilLineDyMomentumCut,
    hEventMinSTFoilLineDzMomentumCut,
    hEventMinSTFoilDistanceMomentumCut,
    false,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_MomentumCut_EventMinSTFoilSeparationComponents.pdf");

  drawSeparationComponentCanvas(
    "cEventMinSTFoilSeparationComponentsMomentumCutLogY",
    "Momentum cut: event-level closest shared ST_Foils line separation components (log-y)",
    hEventMinSTFoilLineDxMomentumCut,
    hEventMinSTFoilLineDyMomentumCut,
    hEventMinSTFoilLineDzMomentumCut,
    hEventMinSTFoilDistanceMomentumCut,
    true,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_MomentumCut_EventMinSTFoilSeparationComponents_LogY.pdf");

  TCanvas* cAllSTFoilLineRxyMomentumCut = new TCanvas(
    "cAllSTFoilLineRxyMomentumCut",
    "Momentum cut: all shared ST_Foils transverse separation",
    900, 700);
  hAllSTFoilLineRxyMomentumCut->Draw("HIST");
  cAllSTFoilLineRxyMomentumCut->SaveAs(
    "Plots/DistancePlots/twoElectronTrkSegVertexer_MomentumCut_AllSTFoilLineRxy.pdf");

  drawSeparationComponentCanvas(
    "cAllSTFoilLineSeparationComponentsMomentumCut",
    "Momentum cut: all shared ST_Foils line separation components",
    hAllSTFoilLineDxMomentumCut,
    hAllSTFoilLineDyMomentumCut,
    hAllSTFoilLineDzMomentumCut,
    hAllSTFoilLineDistanceMomentumCut,
    false,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_MomentumCut_AllSTFoilLineSeparationComponents.pdf");

  drawSeparationComponentCanvas(
    "cAllSTFoilLineSeparationComponentsMomentumCutLogY",
    "Momentum cut: all shared ST_Foils line separation components (log-y)",
    hAllSTFoilLineDxMomentumCut,
    hAllSTFoilLineDyMomentumCut,
    hAllSTFoilLineDzMomentumCut,
    hAllSTFoilLineDistanceMomentumCut,
    true,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_MomentumCut_AllSTFoilLineSeparationComponents_LogY.pdf");

  TCanvas* cEventMinSTFoilDistanceVsElectronTimeDiffMomentumCut = new TCanvas(
    "cEventMinSTFoilDistanceVsElectronTimeDiffMomentumCut",
    "TEST: Momentum cut: event-level closest shared ST_Foils line distance vs absolute shared-foil track-segment time difference",
    900, 700);
  hEventMinSTFoilDistanceVsElectronTimeDiffMomentumCut->Draw("COLZ");
  cEventMinSTFoilDistanceVsElectronTimeDiffMomentumCut->SaveAs(
    "Plots/DistancePlots/twoElectronTrkSegVertexer_TEST_MomentumCut_EventMinSTFoilDistanceVsElectronTimeDiff.pdf");

  draw2DHistogram(
    "cEventMinSTFoilVertexXYMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils midpoint vertex XY",
    hEventMinSTFoilVertexXYMomentumCut,
    kVertexXY,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_MomentumCut_EventMinSTFoilVertexXY.pdf");

  draw2DHistogram(
    "cEventMinSTFoilVertexXZMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils midpoint vertex XZ",
    hEventMinSTFoilVertexXZMomentumCut,
    kVertexXZ,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_MomentumCut_EventMinSTFoilVertexXZ.pdf");

  draw2DHistogram(
    "cEventMinSTFoilVertexYZMomentumCut",
    "Momentum cut: event-level closest shared ST_Foils midpoint vertex YZ",
    hEventMinSTFoilVertexYZMomentumCut,
    kVertexYZ,
    "Plots/DistancePlots/twoElectronTrkSegVertexer_MomentumCut_EventMinSTFoilVertexYZ.pdf");

  // Final footer: keep this visible even in summary-only mode so the scan
  // result is still useful when the verbose printout is disabled.
  cout << "\nSummary" << endl;
  cout << "  events processed: " << eventsPrinted << endl;
  cout << "  shared surfaces processed: " << sharedSurfacesPrinted << endl;
  cout << "  shared ST_Foils surfaces processed: " << sharedSTFoilsPrinted << endl;
  cout << "  events with at least one shared ST_Foils surface: "
       << eventsWithSharedSTFoils << endl;
  cout << "  events with at least one shared ST_Foils surface passing "
       << trackMomentumCutMin << "-" << trackMomentumCutMax
       << " MeV/c momentum cut: "
       << eventsWithSharedSTFoilsMomentumCut << endl;
  const Long64_t eventsCutByMomentumRange =
    eventsWithSharedSTFoils - eventsWithSharedSTFoilsMomentumCut;
  cout << "  momentum cut range: "
       << trackMomentumCutMin << "-" << trackMomentumCutMax
       << " MeV/c" << endl;
  cout << "  events cut by momentum range: "
       << eventsCutByMomentumRange << endl;
  cout << "  events remaining after momentum cut: "
       << eventsWithSharedSTFoilsMomentumCut << endl;
  if (eventsWithSharedSTFoils > 0)
  {
    cout << "  percent remaining after momentum cut among shared ST_Foils events: "
         << 100.0 * static_cast<double>(eventsWithSharedSTFoilsMomentumCut) /
              static_cast<double>(eventsWithSharedSTFoils)
         << "%" << endl;
  }
  else
  {
    cout << "  percent remaining after momentum cut among shared ST_Foils events: undefined" << endl;
  }
  cout << "  percent of 100,000 remaining after momentum cut: "
       << 100.0 * static_cast<double>(eventsWithSharedSTFoilsMomentumCut) /
            momentumCutReferenceEventCount
       << "%" << endl;
  cout << "  events by shared ST_Foils multiplicity:" << endl;
  if (eventsBySharedFoilMultiplicity.empty())
  {
    cout << "    none" << endl;
  }
  else
  {
    const int maxSharedFoils = eventsBySharedFoilMultiplicity.rbegin()->first;
    for (int nSharedFoils = 1; nSharedFoils <= maxSharedFoils; ++nSharedFoils)
    {
      const auto iter = eventsBySharedFoilMultiplicity.find(nSharedFoils);
      const Long64_t count =
        iter != eventsBySharedFoilMultiplicity.end() ? iter->second : 0;
      cout << "    " << nSharedFoils << " shared foils: "
           << count << " events" << endl;
    }
  }
  cout << "  events with closest ST_Foils line distance < 1 mm: "
       << eventsWithMinDistanceUnder1mm << endl;
  cout << "  events with closest ST_Foils line distance < 0.5 mm: "
       << eventsWithMinDistanceUnder0p5mm << endl;
  cout << "  events with closest ST_Foils line distance < 0.1 mm: "
       << eventsWithMinDistanceUnder0p1mm << endl;
  // The footer reports the extrema of the event-level minimum ST-foil line
  // distances, not the extrema of every individual shared foil line pair.
  if (eventsWithSharedSTFoils > 0)
  {
    cout << "  smallest event-level ST_Foils line distance: "
         << smallestEventMinSTFoilDistance << " mm" << endl;
    cout << "  largest event-level ST_Foils line distance: "
         << largestEventMinSTFoilDistance << " mm" << endl;
  }
  else
  {
    cout << "  smallest event-level ST_Foils line distance: undefined" << endl;
    cout << "  largest event-level ST_Foils line distance: undefined" << endl;
  }

  gROOT->SetBatch(originalRootBatchMode);
}

