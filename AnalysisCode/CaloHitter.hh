#ifndef CREATEDCODE_HISTOGRAMMAKERS_CALOHITTER_HH
#define CREATEDCODE_HISTOGRAMMAKERS_CALOHITTER_HH

//----------------------------------------------------------------------------------
//
// CaloHitter.hh
//
// Lightweight ROOT helper for drawing the Mu2e CsI calorimeter disk geometry.
//
// The important analysis-facing structures are:
//   - CaloCrystal: one crystal object; for now it stores only crystalNumber.
//   - Calorimeter: one disk object; it stores disk number and its crystals.
//
// The drawing code below reproduces the ideal Mu2e calorimeter arrangement used
// by Offline/CalorimeterGeom:
//   - two disks
//   - 674 crystals per disk
//   - shifted square lattice
//   - annular inner/outer crystal acceptance
//   - four manually removed edge crystals
//
// This is intentionally header-only so it can be included directly by ROOT
// macros without adding a build-system target.
//
//----------------------------------------------------------------------------------

// Standard C++ utilities for vectors, math, strings, and diagnostic printing.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// ROOT drawing primitives.  TBox draws each square crystal, TEllipse draws the
// inner/outer annulus boundaries, and TCanvas owns the final PDF page.
#include <TBox.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TObject.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

namespace calohitter {

// Public crystal-count constant used by both analysis-facing helpers and the
// internal drawing code.  Keeping this outside detail lets SelectedHit infer a
// disk number from a global crystal ID without depending on private internals.
static const int kCaloHitterCrystalsPerDisk = 674;

// One calorimeter crystal.  This deliberately starts small: the crystal number
// is enough to map EventNtuple calohits onto a drawn location.  Future fields
// can include deposited energy, hit time, MC truth labels, or display colors.
struct CaloCrystal {
  // Global EventNtuple-style crystal number.  For the current two-disk geometry
  // this is 0-673 on the front disk and 674-1347 on the back disk.
  int crystalNumber = -1;

  CaloCrystal() = default;
  explicit CaloCrystal(const int number) : crystalNumber(number) {}
};

using CaloCrystals = std::vector<CaloCrystal>;

// One calorimeter disk.  The user's requested "Calorimeter" object is modeled as
// a disk because Mu2e's calorimeter is a pair of disks, each with its own crystal
// numbering offset.
struct Calorimeter {
  // Raw Mu2e/EventNtuple disk number.  User-facing plots translate this to
  // "Front" and "Back", but storing the raw ID keeps the object compatible
  // with trkcalohit.did and calocluster.diskID_.
  int number = -1;          // disk number: 0 or 1

  // Collection of crystals assigned to this disk.  The crystal objects are
  // intentionally lightweight right now; future passes can add energy, time,
  // display color, MC truth labels, or selected-track associations here.
  CaloCrystals crystals;    // crystals belonging to this disk

  Calorimeter() = default;
  explicit Calorimeter(const int diskNumber) : number(diskNumber) {}
};

// Optional overlay marker supplied by analysis code.
//
// CaloHitter deliberately does not decide which physics objects are important.
// Analysis macros should build SelectedHit records for the objects they want to
// display, then pass the vector to saveCalorimeterPdf().  The marker can be
// placed either by disk-local x/y or, for future crystal-level overlays, by a
// global crystal number.
struct SelectedHit {
  // Raw disk number: 0 = Front, 1 = Back.  If a crystal number is supplied and
  // diskNumber is left negative, the disk can be inferred from the crystal ID.
  int diskNumber = -1;

  // Optional global crystal number.  Set useCrystalCenter=true to draw the
  // marker at this crystal's ideal center instead of using x/y directly.
  int crystalNumber = -1;
  bool useCrystalCenter = false;

  // Disk-local marker position in millimeters.  These coordinates are used when
  // useCrystalCenter is false.
  double x = 0.0;
  double y = 0.0;

  // Optional display metadata controlled by the caller.
  std::string label;
  Color_t markerColor = kRed + 1;
  Style_t markerStyle = 20;
  double markerSize = 1.35;

  SelectedHit() = default;

  static SelectedHit fromXY(const int diskNumberValue,
                            const double xValue,
                            const double yValue,
                            const std::string& labelValue = "",
                            const Color_t markerColorValue = kRed + 1,
                            const Style_t markerStyleValue = 20,
                            const double markerSizeValue = 1.35)
  {
    SelectedHit hit;
    hit.diskNumber = diskNumberValue;
    hit.x = xValue;
    hit.y = yValue;
    hit.label = labelValue;
    hit.markerColor = markerColorValue;
    hit.markerStyle = markerStyleValue;
    hit.markerSize = markerSizeValue;
    return hit;
  }

  static SelectedHit fromCrystal(const int crystalNumberValue,
                                 const std::string& labelValue = "",
                                 const Color_t markerColorValue = kRed + 1,
                                 const Style_t markerStyleValue = 20,
                                 const double markerSizeValue = 1.35)
  {
    SelectedHit hit;
    hit.crystalNumber = crystalNumberValue;
    hit.diskNumber = crystalNumberValue >= 0
                       ? crystalNumberValue / kCaloHitterCrystalsPerDisk
                       : -1;
    hit.useCrystalCenter = true;
    hit.label = labelValue;
    hit.markerColor = markerColorValue;
    hit.markerStyle = markerStyleValue;
    hit.markerSize = markerSizeValue;
    return hit;
  }
};

// Human-facing disk labels.  Keep the raw Mu2e disk numbers internally because
// crystal IDs and EventNtuple disk IDs use 0 and 1.  Only plot/text labels are
// translated to Front and Back for readability.
inline std::string diskShortLabel(const int diskNumber)
{
  if (diskNumber == 0) return "Front";
  if (diskNumber == 1) return "Back";
  return "Unknown";
}

inline std::string diskDisplayLabel(const int diskNumber)
{
  return diskShortLabel(diskNumber) + " disk (raw ID " + std::to_string(diskNumber) + ")";
}

// Internal drawing/geometry record.  CaloCrystal only stores crystalNumber for
// now, but the drawing needs x/y placement, local crystal number, and row.
struct CrystalPlacement {
  int crystalNumber = -1;       // global crystal number: disk * 674 + local
  int diskNumber = -1;          // disk number used to build this placement
  int localCrystalNumber = -1;  // 0-673 within a single disk
  int mapIndex = -1;            // shifted-square mapper index before annulus cuts
  int row = 0;                  // shifted-square row number
  double x = 0.0;               // disk-local x center in mm
  double y = 0.0;               // disk-local y center in mm
};

namespace detail {

// The Offline geometry mapper uses two integer coordinates, l and k, to label
// cells on a shifted square lattice.  ShiftIndex is a minimal local copy of that
// concept so this helper does not depend on CLHEP or compiled Offline libraries.
struct ShiftIndex {
  int l = 0;
  int k = 0;

  ShiftIndex() = default;
  ShiftIndex(const int lValue, const int kValue) : l(lValue), k(kValue) {}

  void add(const ShiftIndex& step)
  {
    l += step.l;
    k += step.k;
  }
};

// Geometry constants from the CsI calorimeter geometry currently used in the
// Offline checkout.  Units are millimeters where applicable.
static const int kNumberOfDisks = 2;
static const int kCrystalsPerDisk = kCaloHitterCrystalsPerDisk;
static const double kInnerCrystalRadius = 374.0;  // mm
static const double kOuterCrystalRadius = 660.0;  // mm
static const double kCrystalXYLength = 34.0;      // mm
static const double kWrapperThickness = 0.150;    // mm
static const double kCellSize = kCrystalXYLength + 2.0 * kWrapperThickness;
static const double kDrawRange = 725.0;           // mm
static const char* kDefaultPlotDirectory = "Plots/CaloHitPlots";

// Step directions used to walk around each ring of the shifted-square lattice.
// These mirror SquareShiftMapper.cc in Offline/CalorimeterGeom.
inline std::vector<ShiftIndex> shiftedSquareSteps()
{
  // The six step directions walk around one complete ring of the shifted-square
  // lattice.  Starting from the bottom of a ring and applying these directions
  // reproduces the ordering used by the Offline square-shift mapper.
  std::vector<ShiftIndex> steps;
  steps.push_back(ShiftIndex(1, 1));
  steps.push_back(ShiftIndex(0, 1));
  steps.push_back(ShiftIndex(-1, 0));
  steps.push_back(ShiftIndex(-1, -1));
  steps.push_back(ShiftIndex(0, -1));
  steps.push_back(ShiftIndex(1, 0));
  return steps;
}

// Number of mapper cells inside all rings up to maxRing.  This is larger than
// the number of real crystals because the annulus cut is applied afterward.
inline int nCrystalMax(const int maxRing)
{
  return 3 * maxRing * (maxRing + 1) + 1;
}

// Convert a mapper index into the shifted-square l,k coordinates.  The indexing
// starts at the center and then walks ring-by-ring around the disk.
inline ShiftIndex lkFromIndex(const int thisIndex)
{
  // Mapper index 0 is the central lattice cell.  The Mu2e calorimeter has a
  // hole there, so this cell will be rejected later by the annulus check, but
  // keeping the mapping complete makes the helper match the Offline convention.
  if (thisIndex < 1) {
    return ShiftIndex(0, 0);
  }

  const std::vector<ShiftIndex> steps = shiftedSquareSteps();
  const int nRing = static_cast<int>(0.5 + std::sqrt(0.25 + static_cast<double>(thisIndex - 1) / 3.0));
  const int nSeg = (thisIndex - 3 * nRing * (nRing - 1) - 1) / nRing;
  const int nPos = (thisIndex - 3 * nRing * (nRing - 1) - 1) % nRing;

  if (nSeg < 0) {
    return ShiftIndex(0, 0);
  }

  int l = nPos * steps.at(nSeg).l;
  int k = -nRing + nPos * steps.at(nSeg).k;

  for (int i = 0; i < nSeg; ++i) {
    l += steps.at(i).l * nRing;
    k += steps.at(i).k * nRing;
  }

  return ShiftIndex(l, k);
}

// Convert mapper l,k coordinates into a disk-local x coordinate.  Adjacent rows
// are shifted by half a cell because x uses (l+k)/2.
inline double mapXFromIndex(const int mapIndex)
{
  const ShiftIndex lk = lkFromIndex(mapIndex);
  return kCellSize * static_cast<double>(lk.l + lk.k) / 2.0;
}

// Convert mapper l,k coordinates into a disk-local y coordinate.  y uses l-k,
// which separates neighboring shifted rows by one full cell pitch.
inline double mapYFromIndex(const int mapIndex)
{
  const ShiftIndex lk = lkFromIndex(mapIndex);
  return kCellSize * static_cast<double>(lk.l - lk.k);
}

// Row number in the shifted-square grid.  This is useful for debugging the disk
// arrangement and later for row-based drawing or summaries.
inline int rowFromIndex(const int mapIndex)
{
  const ShiftIndex lk = lkFromIndex(mapIndex);
  return lk.l - lk.k;
}

// Test whether the full square crystal cell lies inside the allowed annulus.
// The check follows Disk::isInsideDisk: for each crystal edge, compute the
// closest and farthest distance from the disk center and reject edge crossings.
inline bool isInsideDisk(const double x, const double y, const double widthX, const double widthY)
{
  // Check all four edges of the square cell.  A crystal is accepted only if its
  // full square footprint is outside the inner hole and inside the outer radius.
  // This avoids drawing partial crystals that would cross either boundary.
  const double apexX[5] = {-0.5, 0.5, 0.5, -0.5, -0.5};
  const double apexY[5] = {-0.5, -0.5, 0.5, 0.5, -0.5};

  for (int i = 1; i < 5; ++i) {
    const double p0x = x + widthX * apexX[i - 1];
    const double p0y = y + widthY * apexY[i - 1];
    const double p1x = x + widthX * apexX[i];
    const double p1y = y + widthY * apexY[i];

    const double vx = p1x - p0x;
    const double vy = p1y - p0y;
    const double vv = vx * vx + vy * vy;
    double t = -(p0x * vx + p0y * vy) / vv;

    double minDist = 0.0;
    if (t < 0.0) {
      minDist = std::sqrt(p0x * p0x + p0y * p0y);
    } else if (t > 1.0) {
      minDist = std::sqrt(p1x * p1x + p1y * p1y);
    } else {
      const double cx = p0x + t * vx;
      const double cy = p0y + t * vy;
      minDist = std::sqrt(cx * cx + cy * cy);
    }

    const double p0Dist = std::sqrt(p0x * p0x + p0y * p0y);
    const double p1Dist = std::sqrt(p1x * p1x + p1y * p1y);
    const double maxDist = std::max(p0Dist, p1Dist);

    if (minDist < kInnerCrystalRadius || maxDist > kOuterCrystalRadius) {
      return false;
    }
  }

  return true;
}

// The Offline ideal geometry removes these four edge crystals by hand after the
// annulus cut.  Keeping that rule here makes the blank drawing count match Mu2e:
// 674 crystals per disk instead of the pre-removal count.
inline bool isManuallyRemovedCrystal(const double x, const double y)
{
  if (std::abs(x - 257.25) < 1.0 && std::abs(y - 583.1) < 1.0) return true;
  if (std::abs(x - 257.25) < 1.0 && std::abs(y + 583.1) < 1.0) return true;
  if (std::abs(x + 257.25) < 1.0 && std::abs(y - 583.1) < 1.0) return true;
  if (std::abs(x + 257.25) < 1.0 && std::abs(y + 583.1) < 1.0) return true;
  return false;
}

// Allocate and draw a ROOT box owned by the current pad.  ROOT macros often use
// raw pointers for drawn objects; TObject::kCanDelete lets the pad clean them up.
inline void drawOwnedBox(const double xMin, const double yMin, const double xMax, const double yMax,
                         const Color_t lineColor, const Style_t fillStyle, const Color_t fillColor)
{
  // ROOT pads own drawn primitives only when requested explicitly.  kCanDelete
  // lets the pad clean this TBox up when the canvas is destroyed or cleared.
  TBox* box = new TBox(xMin, yMin, xMax, yMax);
  box->SetLineColor(lineColor);
  box->SetLineWidth(1);
  box->SetFillStyle(fillStyle);
  box->SetFillColor(fillColor);
  box->SetBit(TObject::kCanDelete);
  box->Draw("l");
}

// Draw a circular boundary.  Both calorimeter radii are shown as unfilled
// TEllipse objects with equal x/y radii.
inline void drawOwnedCircle(const double radius, const Color_t color, const Width_t width)
{
  // TEllipse is used as a circle by giving equal x and y radii.  The fill is
  // disabled so the crystal grid remains visible underneath the boundaries.
  TEllipse* circle = new TEllipse(0.0, 0.0, radius, radius);
  circle->SetFillStyle(0);
  circle->SetLineColor(color);
  circle->SetLineWidth(width);
  circle->SetBit(TObject::kCanDelete);
  circle->Draw("same");
}

// Draw a light reference axis through the disk center.  This helps distinguish
// the disk-local x and y directions on the otherwise symmetric annulus.
inline void drawOwnedLine(const double x1, const double y1, const double x2, const double y2,
                          const Color_t color)
{
  TLine* line = new TLine(x1, y1, x2, y2);
  line->SetLineColor(color);
  line->SetLineStyle(3);
  line->SetBit(TObject::kCanDelete);
  line->Draw("same");
}

// Draw a small text label inside the pad.  It reports disk number and crystal
// count on the saved PDF.
inline void drawOwnedLatex(const double x, const double y, const std::string& text,
                           const double textSize, const int align)
{
  TLatex* latex = new TLatex(x, y, text.c_str());
  latex->SetTextAlign(align);
  latex->SetTextSize(textSize);
  latex->SetTextFont(42);
  latex->SetBit(TObject::kCanDelete);
  latex->Draw("same");
}

// Draw a caller-supplied selected hit marker.  The marker is kept independent
// from the crystal boxes so analyses can overlay tracks, clusters, crystal
// centers, or any other selected object without changing the geometry code.
inline void drawOwnedMarker(const double x, const double y, const SelectedHit& hit)
{
  TMarker* marker = new TMarker(x, y, hit.markerStyle);
  marker->SetMarkerColor(hit.markerColor);
  marker->SetMarkerSize(hit.markerSize);
  marker->SetBit(TObject::kCanDelete);
  marker->Draw("same");
}

// ROOT's canvas printing does not create missing directories.  Before writing a
// PDF, extract the parent directory from the output path and create it.
inline void ensureOutputDirectory(const std::string& outputPath)
{
  // If the caller passes just "file.pdf", there is no directory component to
  // create.  Otherwise create the parent directory tree before ROOT prints.
  const std::string::size_type slashPosition = outputPath.find_last_of("/\\");
  if (slashPosition == std::string::npos) {
    return;
  }

  const std::string outputDirectory = outputPath.substr(0, slashPosition);
  if (!outputDirectory.empty() && gSystem != nullptr) {
    gSystem->mkdir(outputDirectory.c_str(), true);
  }
}

// Resolve which disk should receive a SelectedHit.  For crystal-center overlays,
// the disk can be inferred from the global crystal number if the caller did not
// set diskNumber explicitly.
inline int selectedHitDiskNumber(const SelectedHit& hit)
{
  if (hit.diskNumber >= 0) {
    return hit.diskNumber;
  }
  if (hit.crystalNumber >= 0) {
    return hit.crystalNumber / kCrystalsPerDisk;
  }
  return -1;
}

// Convert a SelectedHit into a drawable x/y position for this disk.  Most track
// overlays will provide x/y directly.  Crystal overlays can instead request the
// ideal crystal center using the global crystal number.
inline bool selectedHitPosition(const SelectedHit& hit,
                                const std::vector<CrystalPlacement>& placements,
                                double& x,
                                double& y)
{
  if (!hit.useCrystalCenter) {
    x = hit.x;
    y = hit.y;
    return true;
  }

  for (const CrystalPlacement& placement : placements) {
    if (placement.crystalNumber == hit.crystalNumber) {
      x = placement.x;
      y = placement.y;
      return true;
    }
  }

  return false;
}

// Draw all SelectedHit overlays that belong to one disk.  The physics selection
// remains entirely in the caller; this routine only filters by disk and draws
// the provided markers in the same coordinate frame as the blank geometry.
inline void drawSelectedHitsForDisk(const int diskNumber,
                                    const std::vector<SelectedHit>& selectedHits,
                                    const std::vector<CrystalPlacement>& placements)
{
  for (const SelectedHit& hit : selectedHits) {
    if (selectedHitDiskNumber(hit) != diskNumber) {
      continue;
    }

    double x = 0.0;
    double y = 0.0;
    if (!selectedHitPosition(hit, placements, x, y)) {
      std::cerr << "WARNING: could not place selected hit";
      if (hit.crystalNumber >= 0) {
        std::cerr << " with crystal number " << hit.crystalNumber;
      }
      std::cerr << " on disk " << diskNumber << std::endl;
      continue;
    }

    drawOwnedMarker(x, y, hit);
    if (!hit.label.empty()) {
      drawOwnedLatex(x + 18.0, y + 18.0, hit.label, 0.025, 11);
    }
  }
}

}  // namespace detail

// Build the full list of drawable crystal placements for one disk.  This is the
// core geometry routine: scan the shifted-square map, keep cells inside the
// annulus, remove the four hand-removed cells, then assign local/global numbers.
inline std::vector<CrystalPlacement> crystalPlacementsForDisk(const int diskNumber)
{
  // This vector is the detailed geometry representation.  It has more
  // information than Calorimeter::crystals because drawing needs x/y locations
  // and local numbering, while the public analysis object currently stores only
  // crystal number.
  std::vector<CrystalPlacement> placements;
  placements.reserve(detail::kCrystalsPerDisk);

  // Match Disk::fillCrystalsIdeal: choose a mapper range large enough to cover
  // the outer radius, then reject any cells that do not fit.
  const int nRingsMax = static_cast<int>(1.5 * detail::kOuterCrystalRadius / detail::kCellSize);
  const int nCrystalMap = detail::nCrystalMax(nRingsMax);

  int localCrystalNumber = 0;
  for (int mapIndex = 0; mapIndex < nCrystalMap; ++mapIndex) {
    // Convert this mapper index into an ideal crystal center in the disk plane.
    const double x = detail::mapXFromIndex(mapIndex);
    const double y = detail::mapYFromIndex(mapIndex);

    // Reject cells whose square footprint crosses the inner or outer boundary.
    if (!detail::isInsideDisk(x, y, detail::kCellSize, detail::kCellSize)) {
      continue;
    }
    // Reject the four explicitly removed edge cells.
    if (detail::isManuallyRemovedCrystal(x, y)) {
      continue;
    }

    // Local numbering is sequential after all geometry cuts.  The global
    // crystal number follows Mu2e convention: disk * 674 + localCrystalNumber.
    CrystalPlacement placement;
    placement.diskNumber = diskNumber;
    placement.localCrystalNumber = localCrystalNumber;
    placement.crystalNumber = diskNumber * detail::kCrystalsPerDisk + localCrystalNumber;
    placement.mapIndex = mapIndex;
    placement.row = detail::rowFromIndex(mapIndex);
    placement.x = x;
    placement.y = y;
    placements.push_back(placement);

    ++localCrystalNumber;
  }

  // This warning protects against accidental geometry drift if constants or the
  // mapper reproduction are edited later.
  if (static_cast<int>(placements.size()) != detail::kCrystalsPerDisk) {
    std::cerr << "WARNING: CaloHitter built " << placements.size()
              << " crystals for disk " << diskNumber
              << ", expected " << detail::kCrystalsPerDisk << std::endl;
  }

  return placements;
}

// Build the user-facing Calorimeter object for one disk.  The drawable placement
// list contains geometry details, while Calorimeter currently stores only the
// crystal objects requested for the analysis interface.
inline Calorimeter buildCalorimeter(const int diskNumber)
{
  // Convert detailed placements into the compact analysis-facing object.  This
  // keeps the public Calorimeter/CaloCrystal structures simple while allowing
  // the drawing code to keep using full geometry details internally.
  Calorimeter calorimeter(diskNumber);
  const std::vector<CrystalPlacement> placements = crystalPlacementsForDisk(diskNumber);
  calorimeter.crystals.reserve(placements.size());

  for (const CrystalPlacement& placement : placements) {
    calorimeter.crystals.push_back(CaloCrystal(placement.crystalNumber));
  }

  return calorimeter;
}

// Build both Mu2e calorimeter disks.
inline std::vector<Calorimeter> buildCalorimeters()
{
  // Build disks in raw-ID order.  The drawing labels convert these to Front and
  // Back, but the vector order remains disk 0 followed by disk 1.
  std::vector<Calorimeter> calorimeters;
  calorimeters.reserve(detail::kNumberOfDisks);

  for (int diskNumber = 0; diskNumber < detail::kNumberOfDisks; ++diskNumber) {
    calorimeters.push_back(buildCalorimeter(diskNumber));
  }

  return calorimeters;
}

// Draw one calorimeter disk into the current ROOT pad.  The pad must already
// exist; drawCalorimeterDisks creates and divides the canvas before calling this.
// selectedHits is an optional caller-provided overlay; the helper only draws
// records assigned to this disk.
inline void drawCalorimeterDisk(const Calorimeter& calorimeter,
                                const std::vector<SelectedHit>& selectedHits)
{
  // Basic pad styling.  The margins leave room for axis titles.
  gPad->SetTicks(1, 1);
  gPad->SetLeftMargin(0.10);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.10);

  const std::string title = "Calorimeter " + diskDisplayLabel(calorimeter.number) +
                            ";disk-local x [mm];disk-local y [mm]";
  // DrawFrame creates the axes and establishes the coordinate system for the
  // rest of the disk.  All subsequent boxes, circles, and labels are drawn in
  // this same disk-local millimeter frame.
  gPad->DrawFrame(-detail::kDrawRange, -detail::kDrawRange,
                  detail::kDrawRange, detail::kDrawRange, title.c_str());

  // Draw reference axes first so the crystal grid stays visible above them.
  detail::drawOwnedLine(-detail::kDrawRange, 0.0, detail::kDrawRange, 0.0, kGray + 2);
  detail::drawOwnedLine(0.0, -detail::kDrawRange, 0.0, detail::kDrawRange, kGray + 2);

  // Draw every crystal as one square cell.  At this stage the cells are blank;
  // later we can color them by hit energy or mark selected crystals.
  const std::vector<CrystalPlacement> placements = crystalPlacementsForDisk(calorimeter.number);
  for (const CrystalPlacement& placement : placements) {
    const double halfSize = 0.5 * detail::kCellSize;
    // Blank white fill for now.  Later, this is the natural place to translate
    // a crystal's stored energy or hit state into a color before drawing.
    detail::drawOwnedBox(placement.x - halfSize, placement.y - halfSize,
                         placement.x + halfSize, placement.y + halfSize,
                         kGray + 1, 0, kWhite);
  }

  // Draw the annulus boundaries on top of the cells so the physical aperture is
  // clear in the saved PDF.
  detail::drawOwnedCircle(detail::kOuterCrystalRadius, kBlack, 2);
  detail::drawOwnedCircle(detail::kInnerCrystalRadius, kBlack, 2);

  // Overlay caller-selected objects after the blank disk is drawn so the
  // markers and labels sit on top of the crystal grid.
  detail::drawSelectedHitsForDisk(calorimeter.number, selectedHits, placements);

  // Small label confirming which disk is shown and how many crystals were built.
  detail::drawOwnedLatex(-690.0, 665.0,
                         diskDisplayLabel(calorimeter.number) +
                           "  crystals: " + std::to_string(calorimeter.crystals.size()),
                         0.032, 11);
}

// Backward-compatible blank-disk draw entry point.
inline void drawCalorimeterDisk(const Calorimeter& calorimeter)
{
  const std::vector<SelectedHit> selectedHits;
  drawCalorimeterDisk(calorimeter, selectedHits);
}

// Draw all supplied Calorimeter objects on one canvas.  For the current Mu2e
// geometry this makes a two-panel canvas: Front disk on the left, Back disk on the right.
inline TCanvas* drawCalorimeterDisks(const std::vector<Calorimeter>& calorimeters,
                                     const std::vector<SelectedHit>& selectedHits,
                                     const std::string& canvasName = "cCaloHitter")
{
  // Disable ROOT's automatic statistics box.  It is useful for histograms but
  // distracting for a geometry drawing where there is no statistical sample.
  gStyle->SetOptStat(0);

  TCanvas* canvas = new TCanvas(canvasName.c_str(), "CaloHitter calorimeter disks", 1500, 750);
  canvas->Divide(static_cast<int>(calorimeters.size()), 1);

  // ROOT pads are numbered from 1, not 0.
  for (size_t i = 0; i < calorimeters.size(); ++i) {
    canvas->cd(static_cast<int>(i + 1));
    drawCalorimeterDisk(calorimeters.at(i), selectedHits);
  }

  canvas->cd();
  canvas->Modified();
  canvas->Update();
  return canvas;
}

// Backward-compatible blank two-disk draw entry point.
inline TCanvas* drawCalorimeterDisks(const std::vector<Calorimeter>& calorimeters,
                                     const std::string& canvasName = "cCaloHitter")
{
  const std::vector<SelectedHit> selectedHits;
  return drawCalorimeterDisks(calorimeters, selectedHits, canvasName);
}

// Convenience entry point for scripts and macros.  It builds the blank
// calorimeter, draws it, saves the canvas as a PDF, and prints a count summary.
inline TCanvas* saveCalorimeterPdf(const std::string& outputPdf =
                                     std::string(detail::kDefaultPlotDirectory) +
                                     "/CaloHitter_CalorimeterDisks.pdf",
                                   const std::vector<SelectedHit>& selectedHits =
                                     std::vector<SelectedHit>(),
                                   const std::string& canvasName = "cCaloHitter")
{
  // Batch mode prevents ROOT from opening GUI windows when the macro is run from
  // a terminal or batch job.  Restore the user's previous setting before return.
  const bool wasBatchMode = gROOT->IsBatch();
  gROOT->SetBatch(true);

  const std::vector<Calorimeter> calorimeters = buildCalorimeters();
  TCanvas* canvas = drawCalorimeterDisks(calorimeters, selectedHits, canvasName);

  // ROOT will not create missing directories during Print(), so prepare the
  // output directory first and then write the PDF.
  detail::ensureOutputDirectory(outputPdf);
  canvas->Print(outputPdf.c_str());

  // Print disk and total counts so the terminal output confirms the geometry.
  int totalCrystals = 0;
  for (const Calorimeter& calorimeter : calorimeters) {
    totalCrystals += static_cast<int>(calorimeter.crystals.size());
    std::cout << "CaloHitter " << diskDisplayLabel(calorimeter.number)
              << ": " << calorimeter.crystals.size() << " crystals" << std::endl;
  }
  std::cout << "CaloHitter total crystals: " << totalCrystals << std::endl;
  std::cout << "CaloHitter wrote " << outputPdf << std::endl;

  gROOT->SetBatch(wasBatchMode);
  return canvas;
}

}  // namespace calohitter

#endif

