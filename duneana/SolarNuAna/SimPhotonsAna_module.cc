/**
 * @author      : Pablo Barham Alzás 
 * @file        : SimPhotonsAna_module.cc
 * @created     : Thursday Feb 22, 2024 08:46:59 CST
 * @brief       : This analyzer writes out a TTree containing the properties of
 *                each reconstructed flash
 */

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

// C++ includes
#include "math.h"
#include <cstring>

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "dunecore/DuneInterface/Service/RawDigitExtractService.h"
#include "dunecore/DuneInterface/Service/PedestalEvaluationService.h"

#include "lardataobj/Simulation/SimPhotons.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"

namespace sim {

  class SimPhotonsAna : public art::EDAnalyzer {
  public:
    // Standard constructor and destructor for an ART module.
    SimPhotonsAna(const fhicl::ParameterSet&);

    // The analyzer routine, called once per event.
    void analyze(const art::Event&);

  private:
    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fSimPhotonsLiteModuleLabel; // Input tag for SimPhotonsLite object

    // Define the tree
    TTree* fPhotonsTree;

    // Define the variables that will go into the tree
    Int_t fEventID;
    Int_t fOpChannel;
    
    std::vector<int> fTickTime; // Tick time in ns (I think? pretty irrelevant)
    std::vector<int> fDetectedPhotonsCountPerTick; // Number of detected photons per tick
    Int_t fDetectedPhotonsCount; // Total number of detected photons

    // Object that we'll actually read from the artroot file
    std::map<int, int> fDetectedPhotons;
  };

}

namespace sim {

  //-----------------------------------------------------------------------
  // Constructor
  SimPhotonsAna::SimPhotonsAna(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
  {
    // Indicate that the Input Module comes from .fcl
    fSimPhotonsLiteModuleLabel = pset.get<std::string>("SimPhotonsLiteModuleLabel");

    art::ServiceHandle<art::TFileService const> tfs;

    // Make and add the branches to the tree
    fPhotonsTree = tfs->make<TTree>("PhotonsTree", "PhotonsTree");
    fPhotonsTree->Branch("EventID", &fEventID, "EventID/I");
    fPhotonsTree->Branch("OpChannel", &fOpChannel, "OpChannel/I");
    fPhotonsTree->Branch("TickTime", &fTickTime);
    fPhotonsTree->Branch("DetectedPhotonsCountPerTick", &fDetectedPhotonsCountPerTick);
    fPhotonsTree->Branch("DetectedPhotonsCount", &fDetectedPhotonsCount, "DetectedPhotonsCount/I");
  }

  //-----------------------------------------------------------------------
  void
  SimPhotonsAna::analyze(const art::Event& evt)
  {
   
    // auto const* geo = lar::providerFrom<geo::Geometry>();

    // Access ART's TFileService, which will handle creating and writing
    // histograms for us.
    art::ServiceHandle<art::TFileService const> tfs;

    fEventID = evt.id().event();

    // Get the SimPhotonsLite object
    art::Handle<std::vector<sim::SimPhotonsLite>> SimPhotonsLiteHandle;
    evt.getByLabel(fSimPhotonsLiteModuleLabel, SimPhotonsLiteHandle);

    // Loop over all the SimPhotonsLite objects
    for (unsigned int i = 0; i < SimPhotonsLiteHandle->size(); i++) {
      // const sim::SimPhotonsLite& SimPhotonsLite = SimPhotonsLiteHandle->at(i);

      // fOpChannel = SimPhotonsLite.OpChannel();
      // fDetectedPhotons = SimPhotonsLite.DetectedPhotons();

      fOpChannel = SimPhotonsLiteHandle->at(i).OpChannel;
      fDetectedPhotons = SimPhotonsLiteHandle->at(i).DetectedPhotons;

      fTickTime.clear();
      fDetectedPhotonsCountPerTick.clear();
      fDetectedPhotonsCount = 0;

      // Separate the DetectedPhotons map<int,int> into two vectors
      for (std::map<int, int>::iterator it = fDetectedPhotons.begin(); it != fDetectedPhotons.end(); ++it) {
        fTickTime.push_back(it->first);
        fDetectedPhotonsCountPerTick.push_back(it->second);
        fDetectedPhotonsCount += it->second;
      }

      fPhotonsTree->Fill();
    }

  }

} // namespace sim

namespace sim {
  DEFINE_ART_MODULE(SimPhotonsAna)
}
