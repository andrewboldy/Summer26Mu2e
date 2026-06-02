#ifndef CREATEDCODE_HISTOGRAMMAKERS_TWOELECTRONCALOCLUSTERDISTANCE_HH
#define CREATEDCODE_HISTOGRAMMAKERS_TWOELECTRONCALOCLUSTERDISTANCE_HH

//----------------------------------------------------------------------------------
//
// TwoElectronCaloClusterDistance.hh
//
// Helper interface for matching the two selected track-associated calorimeter
// objects back to event-level caloclusters and calculating the distance between
// their reconstructed cluster COGs.
//
// EventNtuple calohits do not store positions.  The COG used here is therefore
// caloclusters.cog_, matched to each selected track's trkcalohit by disk, energy,
// time, and cluster size.
//
//----------------------------------------------------------------------------------

#include <string>
#include <vector>

#include <TH1F.h>

#include "EventNtuple/inc/CaloClusterInfo.hh"
#include "EventNtuple/inc/TrkCaloHitInfo.hh"

namespace twoelectroncalocog {

struct ClusterCogMatch {
  bool valid = false;
  int clusterIndex = -1;
  int disk = -1;
  double energy = -1.0;
  double time = -1.0;
  unsigned size = 0;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double deltaEnergy = -1.0;
  double deltaTime = -1.0;
  double deltaSize = -1.0;
  double matchScore = -1.0;
};

struct CogDistanceResult {
  bool valid = false;
  bool sameDisk = false;
  std::string failureReason;
  ClusterCogMatch electron0;
  ClusterCogMatch electron1;
  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
  double distance3D = -1.0;
  double distanceXY = -1.0;
};

ClusterCogMatch findBestClusterCogMatch(
  const mu2e::TrkCaloHitInfo* trkcalohit,
  const std::vector<mu2e::CaloClusterInfo>* caloclusters);

CogDistanceResult calculateTwoElectronClusterCogDistance(
  const mu2e::TrkCaloHitInfo* electron0TrkCaloHit,
  const mu2e::TrkCaloHitInfo* electron1TrkCaloHit,
  const std::vector<mu2e::CaloClusterInfo>* caloclusters);

TH1F* makeClusterCogDistance3DHistogram();

TH1F* makeClusterCogDistanceXYHistogram();

void fillClusterCogDistanceHistograms(const CogDistanceResult& result,
                                      TH1F* distance3DHistogram,
                                      TH1F* distanceXYHistogram);

std::string formatClusterCogDistanceLine(const CogDistanceResult& result);

}  // namespace twoelectroncalocog

#endif

