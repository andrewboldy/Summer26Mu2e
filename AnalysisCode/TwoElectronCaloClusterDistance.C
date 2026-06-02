//----------------------------------------------------------------------------------
//
// TwoElectronCaloClusterDistance.C
//
// Implementation for the two-electron calorimeter-cluster COG distance helper.
// This file is included by twoElectronCaloAnalysis.C so ROOT macro execution has
// the method definitions without requiring a separate build-system target.
//
//----------------------------------------------------------------------------------

#include "TwoElectronCaloClusterDistance.hh"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace twoelectroncalocog {

namespace {

bool isUsableTrackCaloHit(const mu2e::TrkCaloHitInfo* trkcalohit)
{
  return trkcalohit != nullptr && trkcalohit->did >= 0 && trkcalohit->edep >= 0.0;
}

double absoluteDifference(const double lhs, const double rhs)
{
  return std::fabs(lhs - rhs);
}

double clusterMatchScore(const mu2e::TrkCaloHitInfo& trkcalohit,
                         const mu2e::CaloClusterInfo& cluster)
{
  const double deltaEnergy = absoluteDifference(cluster.energyDep_, trkcalohit.edep);
  const double deltaTime = absoluteDifference(cluster.time_, trkcalohit.ctime);
  const double deltaSize = absoluteDifference(
    static_cast<double>(cluster.size_), static_cast<double>(trkcalohit.csize));

  // Energy, time, and size are all copied from the associated cluster into
  // trkcalohit.  In normal files the correct cluster should be nearly exact.
  // The weighted score keeps the helper useful if minor float formatting or
  // branch-production differences prevent exact equality.
  return deltaEnergy + 0.01 * deltaTime + 0.001 * deltaSize;
}

ClusterCogMatch makeClusterCogMatch(const mu2e::TrkCaloHitInfo& trkcalohit,
                                    const mu2e::CaloClusterInfo& cluster,
                                    const int clusterIndex)
{
  ClusterCogMatch match;
  match.valid = true;
  match.clusterIndex = clusterIndex;
  match.disk = cluster.diskID_;
  match.energy = cluster.energyDep_;
  match.time = cluster.time_;
  match.size = cluster.size_;
  match.x = cluster.cog_.x();
  match.y = cluster.cog_.y();
  match.z = cluster.cog_.z();
  match.deltaEnergy = absoluteDifference(cluster.energyDep_, trkcalohit.edep);
  match.deltaTime = absoluteDifference(cluster.time_, trkcalohit.ctime);
  match.deltaSize = absoluteDifference(
    static_cast<double>(cluster.size_), static_cast<double>(trkcalohit.csize));
  match.matchScore = clusterMatchScore(trkcalohit, cluster);
  return match;
}

void appendClusterMatchText(std::ostringstream& line,
                            const std::string& prefix,
                            const ClusterCogMatch& match)
{
  line << " " << prefix << "_cluster_index=" << match.clusterIndex
       << " " << prefix << "_disk=" << match.disk
       << " " << prefix << "_energy=" << match.energy
       << " " << prefix << "_time=" << match.time
       << " " << prefix << "_size=" << match.size
       << " " << prefix << "_cog_xyz=("
       << match.x << ", " << match.y << ", " << match.z << ")"
       << " " << prefix << "_match_delta_energy=" << match.deltaEnergy
       << " " << prefix << "_match_delta_time=" << match.deltaTime
       << " " << prefix << "_match_delta_size=" << match.deltaSize
       << " " << prefix << "_match_score=" << match.matchScore;
}

}  // namespace

ClusterCogMatch findBestClusterCogMatch(
  const mu2e::TrkCaloHitInfo* trkcalohit,
  const std::vector<mu2e::CaloClusterInfo>* caloclusters)
{
  ClusterCogMatch bestMatch;
  if (!isUsableTrackCaloHit(trkcalohit) || caloclusters == nullptr || caloclusters->empty())
  {
    return bestMatch;
  }

  double bestScore = std::numeric_limits<double>::max();
  for (size_t iCluster = 0; iCluster < caloclusters->size(); ++iCluster)
  {
    const auto& cluster = caloclusters->at(iCluster);
    if (cluster.diskID_ != trkcalohit->did)
    {
      continue;
    }

    const double score = clusterMatchScore(*trkcalohit, cluster);
    if (score < bestScore)
    {
      bestScore = score;
      bestMatch = makeClusterCogMatch(
        *trkcalohit, cluster, static_cast<int>(iCluster));
    }
  }

  return bestMatch;
}

CogDistanceResult calculateTwoElectronClusterCogDistance(
  const mu2e::TrkCaloHitInfo* electron0TrkCaloHit,
  const mu2e::TrkCaloHitInfo* electron1TrkCaloHit,
  const std::vector<mu2e::CaloClusterInfo>* caloclusters)
{
  CogDistanceResult result;

  if (!isUsableTrackCaloHit(electron0TrkCaloHit) ||
      !isUsableTrackCaloHit(electron1TrkCaloHit))
  {
    result.failureReason = "one_or_both_selected_tracks_lack_usable_trkcalohit";
    return result;
  }

  if (caloclusters == nullptr)
  {
    result.failureReason = "caloclusters_branch_missing";
    return result;
  }
  if (caloclusters->empty())
  {
    result.failureReason = "caloclusters_branch_empty";
    return result;
  }

  result.electron0 = findBestClusterCogMatch(electron0TrkCaloHit, caloclusters);
  result.electron1 = findBestClusterCogMatch(electron1TrkCaloHit, caloclusters);

  if (!result.electron0.valid || !result.electron1.valid)
  {
    result.failureReason = "could_not_match_both_trkcalohits_to_caloclusters";
    return result;
  }

  result.valid = true;
  result.sameDisk = result.electron0.disk == result.electron1.disk;
  result.dx = result.electron1.x - result.electron0.x;
  result.dy = result.electron1.y - result.electron0.y;
  result.dz = result.electron1.z - result.electron0.z;
  result.distance3D = std::sqrt(
    result.dx * result.dx + result.dy * result.dy + result.dz * result.dz);
  result.distanceXY = std::sqrt(result.dx * result.dx + result.dy * result.dy);
  return result;
}

TH1F* makeClusterCogDistance3DHistogram()
{
  return new TH1F(
    "hTwoElectronTrackMatchedCaloClusterCogDistance3D",
    "Two-electron events: distance between matched calo cluster COGs;COG distance in disk front-face coordinates [mm];Events",
    300, 0.0, 1500.0);
}

TH1F* makeClusterCogDistanceXYHistogram()
{
  return new TH1F(
    "hTwoElectronTrackMatchedCaloClusterCogDistanceXY",
    "Two-electron events: transverse distance between matched calo cluster COGs;COG #DeltaR_{xy} in disk front-face coordinates [mm];Events",
    300, 0.0, 1500.0);
}

void fillClusterCogDistanceHistograms(const CogDistanceResult& result,
                                      TH1F* distance3DHistogram,
                                      TH1F* distanceXYHistogram)
{
  if (!result.valid)
  {
    return;
  }
  if (distance3DHistogram != nullptr)
  {
    distance3DHistogram->Fill(result.distance3D);
  }
  if (distanceXYHistogram != nullptr)
  {
    distanceXYHistogram->Fill(result.distanceXY);
  }
}

std::string formatClusterCogDistanceLine(const CogDistanceResult& result)
{
  std::ostringstream line;
  line << "  CALO_CLUSTER_COG_DISTANCE";

  if (!result.valid)
  {
    line << " valid=0 reason=" << result.failureReason;
    return line.str();
  }

  line << " valid=1"
       << " same_disk=" << (result.sameDisk ? 1 : 0)
       << std::fixed << std::setprecision(6)
       << " dx=" << result.dx
       << " dy=" << result.dy
       << " dz=" << result.dz
       << " distance_3d_diskff=" << result.distance3D
       << " distance_xy_diskff=" << result.distanceXY;
  appendClusterMatchText(line, "electron0", result.electron0);
  appendClusterMatchText(line, "electron1", result.electron1);
  return line.str();
}

}  // namespace twoelectroncalocog


