// momentumHists.C
//
// Plot reconstructed TT_Mid trkseg momentum through the same progressive cuts
// used by twoElectronTruthTrkSegVertexerComparer.C, overlaid with rank-0
// downstream-electron trkmcsim truth.  The left pad is linear and the right pad
// is logarithmic.
//
// ROOT usage:
//   .L CreatedCode/momentumHists.C+
//   momentumHists("path/to/nts.root")
//
// A text file containing one ROOT file per line may also be supplied.  Optional
// arguments select another track branch and PDF output file name:
//   momentumHists("files.txt", "trk", "other.pdf")

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TChain.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "EventNtuple/inc/SimInfo.hh"
#include "EventNtuple/inc/TrkCaloHitInfo.hh"
#include "EventNtuple/inc/TrkInfo.hh"
#include "EventNtuple/inc/TrkSegInfo.hh"

using namespace std;

namespace momentum_hists_detail
{
  bool hasRootSuffix(const string& name)
  {
    const string suffix = ".root";
    return name.size() >= suffix.size() &&
           name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0;
  }

  bool addInput(TChain& chain, const string& inputName)
  {
    if (hasRootSuffix(inputName))
    {
      return chain.Add(inputName.c_str()) > 0;
    }

    ifstream fileList(inputName);
    if (!fileList.is_open())
    {
      cerr << "ERROR: could not open input file or filelist: " << inputName << endl;
      return false;
    }

    string line;
    int filesAdded = 0;
    while (getline(fileList, line))
    {
      const size_t first = line.find_first_not_of(" \t\r");
      if (first == string::npos || line[first] == '#')
      {
        continue;
      }
      const size_t last = line.find_last_not_of(" \t\r");
      filesAdded += chain.Add(line.substr(first, last - first + 1).c_str()) > 0;
    }
    return filesAdded > 0;
  }

  void enableBranch(TChain& chain, const string& name)
  {
    chain.SetBranchStatus(name.c_str(), 1);
    chain.SetBranchStatus((name + ".*").c_str(), 1);
  }

  const mu2e::TrkSegInfo* downstreamTrackerMiddleSegment(
    const vector<mu2e::TrkSegInfo>& segments)
  {
    for (const auto& segment : segments)
    {
      if (segment.sid == mu2e::SurfaceIdDetail::TT_Mid && segment.mom.z() > 0.0)
      {
        return &segment;
      }
    }
    return nullptr;
  }

  bool hasGoodCaloHit(const vector<mu2e::TrkCaloHitInfo>* caloHits,
                      size_t trackIndex)
  {
    if (caloHits == nullptr || trackIndex >= caloHits->size())
    {
      return false;
    }
    const auto& hit = caloHits->at(trackIndex);
    return hit.active && hit.did >= 0;
  }

  bool isTruthElectron(const mu2e::SimInfo& sim)
  {
    return sim.valid && sim.rank == 0 && sim.pdg == 11 && sim.mom.z() > 0.0;
  }

  void styleHistogram(TH1D& histogram, Color_t color, Style_t lineStyle,
                      Width_t lineWidth = 3)
  {
    histogram.SetLineColor(color);
    histogram.SetLineStyle(lineStyle);
    histogram.SetLineWidth(lineWidth);
    histogram.SetFillStyle(0);
    histogram.SetStats(false);
  }

  void drawPad(TH1D& truth, TH1D& allReco, TH1D& twoTracks,
               TH1D& goodCalo, TH1D& momentumCut, bool logarithmic)
  {
    gPad->SetTicks(1, 1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.12);
    gPad->SetLogy(logarithmic);

    const double maximum = max({truth.GetMaximum(), allReco.GetMaximum(),
                                twoTracks.GetMaximum(), goodCalo.GetMaximum(),
                                momentumCut.GetMaximum()});
    const double yMinimum = logarithmic ? 0.5 : 0.0;
    const double yMaximum = maximum > 0.0
                              ? maximum * (logarithmic ? 8.0 : 1.25)
                              : 1.0;
    const char* title = logarithmic
                          ? "Momentum cut progression (log scale);Momentum [MeV/c];Tracks"
                          : "Momentum cut progression;Momentum [MeV/c];Tracks";
    TH1F* frame = gPad->DrawFrame(30.0, yMinimum, 60.0, yMaximum, title);
    frame->SetStats(false);

    // Draw truth first as the common reference, then all reconstructed stages.
    truth.Draw("HIST SAME");
    allReco.Draw("HIST SAME");
    twoTracks.Draw("HIST SAME");
    goodCalo.Draw("HIST SAME");
    momentumCut.Draw("HIST SAME");

    // Upper-left avoids the expected signal peak near 52 MeV/c.
    TLegend* legend = new TLegend(0.14, 0.64, 0.55, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.031);
    legend->AddEntry(&truth, "trkmcsim truth (rank 0 e^{-})", "l");
    legend->AddEntry(&allReco, "All reco downstream e^{-} tracks", "l");
    legend->AddEntry(&twoTracks, "Exactly two tracks", "l");
    legend->AddEntry(&goodCalo, "+ good calo hit for each", "l");
    legend->AddEntry(&momentumCut, "+ both 50 #leq p #leq 53 MeV/c", "l");
    legend->Draw();
    gPad->RedrawAxis();
  }
}

void momentumHists(const string& inputName,
                   const string& trackBranch = "trk",
                   const string& plotOutputName =
                     "Plots/GoodTimingData/2ndRun/BasicAnalysisPlots/momentumHists.pdf")
{
  using namespace momentum_hists_detail;

  const string defaultOutputDirectory =
    "Plots/GoodTimingData/2ndRun/BasicAnalysisPlots";
  if (gSystem->mkdir(defaultOutputDirectory.c_str(), true) != 0 &&
      gSystem->AccessPathName(defaultOutputDirectory.c_str()))
  {
    cerr << "ERROR: could not create output directory: "
         << defaultOutputDirectory << endl;
    return;
  }

  const string segmentBranch = trackBranch + "segs";
  const string simBranch = trackBranch + "mcsim";
  const string caloBranch = trackBranch + "calohit";

  TChain ntuple("EventNtuple/ntuple");
  if (!addInput(ntuple, inputName) || ntuple.GetEntries() <= 0)
  {
    cerr << "ERROR: no EventNtuple/ntuple entries found for: " << inputName << endl;
    return;
  }

  for (const string& branch : {trackBranch, segmentBranch, simBranch, caloBranch})
  {
    if (ntuple.GetBranch(branch.c_str()) == nullptr)
    {
      cerr << "ERROR: required branch '" << branch << "' is missing." << endl;
      return;
    }
  }

  ntuple.SetBranchStatus("*", 0);
  enableBranch(ntuple, trackBranch);
  enableBranch(ntuple, segmentBranch);
  enableBranch(ntuple, simBranch);
  enableBranch(ntuple, caloBranch);

  vector<mu2e::TrkInfo>* tracks = nullptr;
  vector<vector<mu2e::TrkSegInfo>>* segmentsByTrack = nullptr;
  vector<vector<mu2e::SimInfo>>* truthByTrack = nullptr;
  vector<mu2e::TrkCaloHitInfo>* caloHits = nullptr;
  ntuple.SetBranchAddress(trackBranch.c_str(), &tracks);
  ntuple.SetBranchAddress(segmentBranch.c_str(), &segmentsByTrack);
  ntuple.SetBranchAddress(simBranch.c_str(), &truthByTrack);
  ntuple.SetBranchAddress(caloBranch.c_str(), &caloHits);

  const int nBins = 120; // 0.25 MeV/c per bin from 30 to 60 MeV/c.
  TH1D hTruth("hTruthMomentum", "", nBins, 30.0, 60.0);
  TH1D hAllReco("hAllRecoMomentum", "", nBins, 30.0, 60.0);
  TH1D hExactlyTwo("hExactlyTwoMomentum", "", nBins, 30.0, 60.0);
  TH1D hGoodCalo("hGoodCaloMomentum", "", nBins, 30.0, 60.0);
  TH1D hMomentumCut("hMomentumCut", "", nBins, 30.0, 60.0);

  styleHistogram(hTruth, kBlack, 1, 4);
  styleHistogram(hAllReco, kGray + 2, 2);
  styleHistogram(hExactlyTwo, kBlue + 1, 1);
  styleHistogram(hGoodCalo, kOrange + 7, 1);
  styleHistogram(hMomentumCut, kRed + 1, 1, 4);

  Long64_t eventsExactlyTwo = 0;
  Long64_t eventsGoodCalo = 0;
  Long64_t eventsMomentumCut = 0;

  for (Long64_t entry = 0; entry < ntuple.GetEntries(); ++entry)
  {
    ntuple.GetEntry(entry);
    if (tracks == nullptr || segmentsByTrack == nullptr)
    {
      continue;
    }

    // Truth is an all-event reference.  Count each stored rank-0 downstream
    // electron match once per reconstructed track.
    if (truthByTrack != nullptr)
    {
      for (const auto& trackTruth : *truthByTrack)
      {
        for (const auto& sim : trackTruth)
        {
          if (isTruthElectron(sim))
          {
            hTruth.Fill(sim.mom.R());
          }
        }
      }
    }

    vector<size_t> selectedIndices;
    vector<double> selectedMomenta;
    const size_t nTracks = min(tracks->size(), segmentsByTrack->size());
    for (size_t iTrack = 0; iTrack < nTracks; ++iTrack)
    {
      if (tracks->at(iTrack).pdg != 11)
      {
        continue;
      }
      const mu2e::TrkSegInfo* middle =
        downstreamTrackerMiddleSegment(segmentsByTrack->at(iTrack));
      if (middle == nullptr)
      {
        continue;
      }

      const double momentum = middle->mom.R();
      hAllReco.Fill(momentum);
      selectedIndices.push_back(iTrack);
      selectedMomenta.push_back(momentum);
    }

    if (selectedIndices.size() != 2)
    {
      continue;
    }
    ++eventsExactlyTwo;
    hExactlyTwo.Fill(selectedMomenta[0]);
    hExactlyTwo.Fill(selectedMomenta[1]);

    if (!hasGoodCaloHit(caloHits, selectedIndices[0]) ||
        !hasGoodCaloHit(caloHits, selectedIndices[1]))
    {
      continue;
    }
    ++eventsGoodCalo;
    hGoodCalo.Fill(selectedMomenta[0]);
    hGoodCalo.Fill(selectedMomenta[1]);

    if (selectedMomenta[0] < 50.0 || selectedMomenta[0] > 53.0 ||
        selectedMomenta[1] < 50.0 || selectedMomenta[1] > 53.0)
    {
      continue;
    }
    ++eventsMomentumCut;
    hMomentumCut.Fill(selectedMomenta[0]);
    hMomentumCut.Fill(selectedMomenta[1]);
  }

  const bool oldBatch = gROOT->IsBatch();
  const int oldOptStat = gStyle->GetOptStat();
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  TCanvas canvas("cMomentumHists", "Momentum histograms", 1500, 650);
  canvas.Divide(2, 1, 0.008, 0.008);
  canvas.cd(1);
  drawPad(hTruth, hAllReco, hExactlyTwo, hGoodCalo, hMomentumCut, false);
  canvas.cd(2);
  drawPad(hTruth, hAllReco, hExactlyTwo, hGoodCalo, hMomentumCut, true);
  canvas.SaveAs(plotOutputName.c_str());
  gStyle->SetOptStat(oldOptStat);
  gROOT->SetBatch(oldBatch);

  cout << "Scanned " << ntuple.GetEntries() << " events." << endl;
  cout << "Exactly-two-track events: " << eventsExactlyTwo << endl;
  cout << "After good-calo cut: " << eventsGoodCalo << endl;
  cout << "After 50-53 MeV/c cut: " << eventsMomentumCut << endl;
  cout << "Wrote " << plotOutputName << endl;
}

