// Generates two electrons at half endpoint energy that will be attached to a mu- in
// the input SimParticleCollection.
// This module throws an exception if no suitable muon is found.
//
// Andrew Boldy Spring 2026

#include <iostream>
#include <string>
#include <cmath>
#include <memory>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"

namespace mu2e {

  //================================================================
  class B2BCeEndpointDiffEnergy : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"),Comment("A SimParticleCollection with input stopped muons.")};
        fhicl::Atom<std::string> stoppingTargetMaterial{
        Name("stoppingTargetMaterial"),Comment("Material determines endpoint energy and muon life time.  Material must be known to the GlobalConstantsService."),"Al" };
        fhicl::Atom<unsigned> verbosity{Name("verbosity"),0};
        fhicl::Atom<int> pdgId{Name("pdgId"),Comment("pdg id of daughter particle")};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit B2BCeEndpointDiffEnergy(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:
    const PDGCode::type electronId_ = PDGCode::e_minus; // for mass only
    double electronMass_;
    double endPointEnergy_;
    double endPointMomentum_;
    double newEnergy_;
    double newMomentum_;
    double muonLifeTime_;

    art::ProductToken<SimParticleCollection> const simsToken_;

    unsigned verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandExponential randExp_;
    RandomUnitSphere   randomUnitSphere_;
    ProcessCode process;
    int pdgId_;
    PDGCode::type pid;
  };

  //================================================================
  B2BCeEndpointDiffEnergy::B2BCeEndpointDiffEnergy(const Parameters& conf)
    : EDProducer{conf}
    , electronMass_(GlobalConstantsHandle<ParticleDataList>()->particle(electronId_).mass())
    , endPointEnergy_()
    , endPointMomentum_ ()
    , newEnergy_(0.0)
    , newMomentum_ (0.0)
    , muonLifeTime_{GlobalConstantsHandle<PhysicsParams>()->getDecayTime(conf().stoppingTargetMaterial())}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
    , randomUnitSphere_{eng_}
    , pdgId_(conf().pdgId())
  {
    produces<mu2e::StageParticleCollection>();
    pid = static_cast<PDGCode::type>(pdgId_);

    if (pid == PDGCode::e_minus) {
      process = ProcessCode::mu2eCeEndpointDiffEnergyB2B;
      endPointEnergy_ = GlobalConstantsHandle<PhysicsParams>()->getEndpointEnergy(conf().stoppingTargetMaterial());
      newEnergy_ = endPointEnergy_/2;
    }
    else if (pid == PDGCode::e_plus) {
      process = ProcessCode::mu2eCePlusEndpoint;
      endPointEnergy_ = GlobalConstantsHandle<PhysicsParams>()->getePlusEndpointEnergy(conf().stoppingTargetMaterial());
      newEnergy_ = endPointEnergy_/2;
    }
    else {
      throw   cet::exception("BADINPUT")
        <<"CeEndpointDiffEnergyGenerator::produce(): No process associated with chosen PDG id\n";
    }
    endPointMomentum_ = endPointEnergy_*sqrt(1 - std::pow(electronMass_/endPointEnergy_,2));
    newMomentum_ = newEnergy_*sqrt(1-std::pow(electronMass_/newEnergy_,2));
    if(verbosity_ > 0) {
      mf::LogInfo log("CeEndpointDiffEnergy");
      log<<"stoppingTargetMaterial = "<<conf().stoppingTargetMaterial()
         <<", electron energy = "<<newEnergy_
         <<", muon lifetime = "<<muonLifeTime_
         <<std::endl;
    }
  }

  //================================================================
  void B2BCeEndpointDiffEnergy::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto mus = stoppedMuMinusList(simh);

    if(mus.empty()) {
      throw   cet::exception("BADINPUT")
        <<"CeEndpoint::produce(): no suitable stopped muon in the input SimParticleCollection\n";

    }

    // Normally we have exactly one mu stop, but it is not impossible to get more.
    // Pick one of them; we don't want more than one CE per event.
    // Note that Rmue normalization is per muon, not per primary proton.

    const auto mustop = mus.at(eng_.operator unsigned int() % mus.size());

    const auto p3 = randomUnitSphere_.fire(newMomentum_);
    const auto mom1 = CLHEP::HepLorentzVector{p3,newEnergy_};
    const auto mom2 = CLHEP::HepLorentzVector{-1*p3,newEnergy_};
    const auto daughterTime = mustop->endGlobalTime() + randExp_.fire(muonLifeTime_);
    output->emplace_back(mustop,
                         process,
                         pid,
                         mustop->endPosition(),
                         //CLHEP::HepLorentzVector{randomUnitSphere_.fire(endPointMomentum_), endPointEnergy_},
                         mom1,
			 //CLHEP::HepLorentzVector{p3,newEnergy_},
                         daughterTime
                         );
    output->emplace_back(mustop,
                        process,
                        pid,
                        mustop->endPosition(),
                        //CLHEP::HepLorentzVector{randomUnitSphere_.fire(endPointMomentum_), endPointEnergy_},
                        mom2,
			//CLHEP::HepLorentzVector{-1*p3,newEnergy_},
                        daughterTime
                                              );
    event.put(std::move(output));
      //Make Outputs that can be read into a txt file. 
//cout << "#$%#$%" << endl;
//cout << "Printing the following information for each pair of Lorentz vectors being added to the stack: " 
//<< "momNumber (1 or 2) || threeMom vector Cartesian (x,y,z) || threeMomSpherical (r,theta,phi) || magnitude (of 4 vector) || Energy (MeV/c)" << endl;

    std::cout << "MomNumber: 1 || " 
     << "Cartesian vector: (" << mom1.x() << ", " << mom1.y() << ", " << mom1.z() << ")" 
     << " || Spherical Vector: (" << mom1.rho() << ", " << mom1.theta() << ", " << mom1.phi() << ")"
     << " || Magnitude of 4 Vector: " << mom1.mag() 
     << " || Energy: " << mom1.e() << std::endl; 

std::cout << "MomNumber: 2 || " 
     << "Cartesian vector: (" << mom2.x() << ", " << mom2.y() << ", " << mom2.z() << ")" 
     << " || Spherical Vector: (" << mom2.rho() << ", " << mom2.theta() << ", " << mom2.phi() << ")"
     << " || Magnitude of 4 Vector: " << mom2.mag() 
     << " || Energy: " << mom2.e() << std::endl; 

//cout << "%$#%$#" << endl;
  }

  //================================================================
} // namespace mu2e

//DEFINE_ART_MODULE(mu2e::CeEndpoint)
DEFINE_ART_MODULE(mu2e::B2BCeEndpointDiffEnergy)
