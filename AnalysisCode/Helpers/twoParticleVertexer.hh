#ifndef CREATEDCODE_HISTOGRAMMAKERS_TWOPARTICLEVERTEXER_HH
#define CREATEDCODE_HISTOGRAMMAKERS_TWOPARTICLEVERTEXER_HH

//----------------------------------------------------------------------------------
//
// twoParticleVertexer.hh
// Written by Andrew Boldy University of South Carolina, 2026
// Assisted by Codex
//
// Purpose:
//   Provide a small, reusable two-particle vertex helper for ROOT/EventNtuple
//   analysis macros.
//
//   The helper treats each reconstructed particle state as a 3D line:
//
//       line(s) = position + s * unit(momentum)
//
//   For two such lines, it calculates the pair of closest approach points, the
//   midpoint between those points, and the separation distance between the two
//   lines.  That midpoint is the first-pass reconstructed two-particle vertex.
//
//   This is deliberately generic.  It does not assume the particles are
//   electrons, does not assume a particular EventNtuple branch prefix, and does
//   not decide which track segment should be used.  The calling analysis chooses
//   the two states and passes only position, momentum, time, and labels here.
//
// Units and coordinates:
//   - position is assumed to be in millimeters
//   - momentum is assumed to be in MeV/c
//   - time is assumed to be in nanoseconds
//   - the coordinate frame is whatever frame the caller supplies
//
// EventNtuple convention:
//   TrkSegInfo::pos and TrkSegInfo::mom are already in the tracker/detector
//   coordinate frame used by the existing vertexer macro, so they can be passed
//   directly into makeParticleStateFromTrackSegment(...).
//
//----------------------------------------------------------------------------------

#include <cmath>
#include <string>

#include "EventNtuple/inc/RootVectors.hh"
#include "EventNtuple/inc/TrkSegInfo.hh"

namespace twoparticlevertexer
{
  // Numerical guard for line-line closest approach.  When two momentum
  // directions are parallel or nearly parallel, the closest-approach equations
  // become singular.  In that case VertexResult::valid is false.
  static const double kDefaultParallelTolerance = 1e-12;

  // One particle state at one chosen reference point.
  //
  // The reference point can be a reconstructed track segment, a calorimeter
  // POCA, a virtual-detector step, or any other point/direction pair.  Extra
  // metadata fields are optional and exist only to make printouts and debugging
  // easier for caller code.
  struct ParticleState
  {
    bool valid = false;
    int objectIndex = -1;
    int pdg = 0;
    int surfaceId = -1;
    int surfaceIndex = -1;
    std::string label;
    XYZVectorF position = XYZVectorF();
    XYZVectorF momentum = XYZVectorF();
    double time = 0.0;
  };

  // Internal line representation used by the vertex calculation.
  //
  // The point is copied from ParticleState::position.  The direction is the
  // unit momentum vector.  The original scalar momentum is kept because it is
  // useful for diagnostics and future histogram labels.
  struct Line3D
  {
    bool valid = false;
    int objectIndex = -1;
    int pdg = 0;
    int surfaceId = -1;
    int surfaceIndex = -1;
    std::string label;
    XYZVectorF point = XYZVectorF();
    XYZVectorF unitDirection = XYZVectorF();
    double momentum = -1.0;
    double time = 0.0;
  };

  // Result of the two-line closest-approach calculation.
  //
  // closestPointOnFirst and closestPointOnSecond are the two points, one on
  // each reconstructed particle line, that minimize the line-line distance.
  // vertex is the midpoint between them.  separation points from the first
  // closest point to the second closest point.
  struct VertexResult
  {
    bool valid = false;
    std::string failureReason;
    Line3D firstLine;
    Line3D secondLine;
    XYZVectorF closestPointOnFirst = XYZVectorF();
    XYZVectorF closestPointOnSecond = XYZVectorF();
    XYZVectorF vertex = XYZVectorF();
    XYZVectorF separation = XYZVectorF();
    double distance = -1.0;
    double firstLineParameter = 0.0;
    double secondLineParameter = 0.0;
    double averageInputTime = 0.0;
    double deltaInputTime = 0.0;
  };

  inline double magnitude(const XYZVectorF& vector)
  {
    return vector.R();
  }

  inline bool hasNonZeroDirection(const XYZVectorF& momentum)
  {
    return magnitude(momentum) > 0.0;
  }

  inline XYZVectorF unitVector(const XYZVectorF& vector)
  {
    const double mag = magnitude(vector);
    if (mag <= 0.0)
    {
      return XYZVectorF();
    }
    return vector / mag;
  }

  inline ParticleState makeParticleState(const XYZVectorF& position,
                                         const XYZVectorF& momentum,
                                         const double time,
                                         const int objectIndex = -1,
                                         const int pdg = 0,
                                         const int surfaceId = -1,
                                         const int surfaceIndex = -1,
                                         const std::string& label = "")
  {
    ParticleState state;
    state.valid = hasNonZeroDirection(momentum);
    state.objectIndex = objectIndex;
    state.pdg = pdg;
    state.surfaceId = surfaceId;
    state.surfaceIndex = surfaceIndex;
    state.label = label;
    state.position = position;
    state.momentum = momentum;
    state.time = time;
    return state;
  }

  inline ParticleState makeParticleStateFromTrackSegment(const mu2e::TrkSegInfo& segment,
                                                         const int trackIndex = -1,
                                                         const int pdg = 11,
                                                         const std::string& label = "")
  {
    return makeParticleState(segment.pos,
                             segment.mom,
                             segment.time,
                             trackIndex,
                             pdg,
                             segment.sid,
                             segment.sindex,
                             label);
  }

  inline ParticleState makeParticleStateFromTrackSegment(const mu2e::TrkSegInfo& segment,
                                                         const int trackIndex,
                                                         const std::string& label)
  {
    return makeParticleStateFromTrackSegment(segment, trackIndex, 11, label);
  }

  inline Line3D makeLine(const ParticleState& state)
  {
    Line3D line;
    line.valid = state.valid && hasNonZeroDirection(state.momentum);
    line.objectIndex = state.objectIndex;
    line.pdg = state.pdg;
    line.surfaceId = state.surfaceId;
    line.surfaceIndex = state.surfaceIndex;
    line.label = state.label;
    line.point = state.position;
    line.unitDirection = unitVector(state.momentum);
    line.momentum = magnitude(state.momentum);
    line.time = state.time;
    return line;
  }

  inline VertexResult vertexFromLines(const Line3D& firstLine,
                                      const Line3D& secondLine,
                                      const double parallelTolerance = kDefaultParallelTolerance)
  {
    VertexResult result;
    result.firstLine = firstLine;
    result.secondLine = secondLine;
    result.averageInputTime = 0.5 * (firstLine.time + secondLine.time);
    result.deltaInputTime = firstLine.time - secondLine.time;

    if (!firstLine.valid || !secondLine.valid)
    {
      result.failureReason = "one or both input lines are invalid";
      return result;
    }

    // Closest approach between two lines:
    //
    //   p1(s) = r1 + s * u1
    //   p2(t) = r2 + t * u2
    //
    // The closest pair satisfies:
    //
    //   (p1 - p2) dot u1 = 0
    //   (p1 - p2) dot u2 = 0
    //
    // Solving those two linear equations gives the parameters s and t below.
    const XYZVectorF w0 = firstLine.point - secondLine.point;
    const double a = firstLine.unitDirection.Dot(firstLine.unitDirection);
    const double b = firstLine.unitDirection.Dot(secondLine.unitDirection);
    const double c = secondLine.unitDirection.Dot(secondLine.unitDirection);
    const double d = firstLine.unitDirection.Dot(w0);
    const double e = secondLine.unitDirection.Dot(w0);
    const double denominator = a * c - b * b;

    if (std::fabs(denominator) <= parallelTolerance)
    {
      result.failureReason = "input lines are parallel or nearly parallel";
      return result;
    }

    result.firstLineParameter = (b * e - c * d) / denominator;
    result.secondLineParameter = (a * e - b * d) / denominator;
    result.closestPointOnFirst =
      firstLine.point + result.firstLineParameter * firstLine.unitDirection;
    result.closestPointOnSecond =
      secondLine.point + result.secondLineParameter * secondLine.unitDirection;
    result.separation = result.closestPointOnSecond - result.closestPointOnFirst;
    result.vertex = 0.5 * (result.closestPointOnFirst + result.closestPointOnSecond);
    result.distance = result.separation.R();
    result.valid = true;
    return result;
  }

  inline VertexResult vertexFromParticleStates(const ParticleState& firstState,
                                               const ParticleState& secondState,
                                               const double parallelTolerance = kDefaultParallelTolerance)
  {
    return vertexFromLines(makeLine(firstState),
                           makeLine(secondState),
                           parallelTolerance);
  }
}

#endif

