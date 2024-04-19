/**
 * @author      : Pablo Barham Alzás, Daniele Guffanti
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
    struct Config_t {
      fhicl::Atom<std::string>     module_type{ fhicl::Name("module_type"), "SimPhotonsAna" };
      fhicl::Sequence<std::string> simphotons_labels{ fhicl::Name("SimPhotonsLiteModuleLabels") };
    };

    using Parameters = art::EDAnalyzer::Table<Config_t>;

    // Standard constructor and destructor for an ART module.
    explicit SimPhotonsAna(Parameters const&);

    // The analyzer routine, called once per event.
    void analyze(const art::Event&);

  private:
    const std::map<std::string, UShort_t> fSimPhotonsLabelIDMap = {
      {"PDFastSimAr", 1}, 
      {"PDFastSimXe", 2}, 
      {"PDFastSimArExternal", 3}, 
      {"PDFastSimXeExternal", 4}
    };

    inline UShort_t GetSimPhotonsLabelID(const std::string& label) const {
      if (fSimPhotonsLabelIDMap.find(label) != fSimPhotonsLabelIDMap.end()) {
        return fSimPhotonsLabelIDMap.find(label)->second; 
      }
      return 0;
    };

    // The parameters we'll read from the .fcl file.
    std::vector<std::string> fSimPhotonsLiteModuleLabels = {}; // Input tags for SimPhotonsLite object

    // Define the tree
    TTree* fPhotonsTree = {};

    // Define the variables that will go into the tree
    Int_t fEventID = {};
    Int_t fOpChannel = {};
    UShort_t fSimPhotonsLabelID = {};
    
    std::vector<int> fTickTime = {}; // Tick time in ns (I think? pretty irrelevant)
    std::vector<int> fDetectedPhotonsCountPerTick = {}; // Number of detected photons per tick
    Int_t fDetectedPhotonsCount = {}; // Total number of detected photons

    // Object that we'll actually read from the artroot file
    std::map<int, int> fDetectedPhotons;
  };

}

namespace sim {

  //-----------------------------------------------------------------------
  // Constructor
  SimPhotonsAna::SimPhotonsAna(Parameters const& pset) : art::EDAnalyzer{pset}, 
    fSimPhotonsLiteModuleLabels{ pset().simphotons_labels() }
  {
    art::ServiceHandle<art::TFileService const> tfs;

    // Make and add the branches to the tree
    fPhotonsTree = tfs->make<TTree>("PhotonsTree", "PhotonsTree");
    fPhotonsTree->Branch("EventID", &fEventID, "EventID/I");
    fPhotonsTree->Branch("OpChannel", &fOpChannel, "OpChannel/I");
    fPhotonsTree->Branch("TickTime", &fTickTime);
    fPhotonsTree->Branch("DetectedPhotonsCountPerTick", &fDetectedPhotonsCountPerTick);
    fPhotonsTree->Branch("DetectedPhotonsCount", &fDetectedPhotonsCount, "DetectedPhotonsCount/I");
    fPhotonsTree->Branch("SimPhotonsLabel", &fSimPhotonsLabelID, "SimPhotonsLabel/s");
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

    for (const auto& simphotons_label : fSimPhotonsLiteModuleLabels) {
      evt.getByLabel(simphotons_label, SimPhotonsLiteHandle);

      fSimPhotonsLabelID = GetSimPhotonsLabelID( simphotons_label ); 
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
  }

} // namespace sim

namespace sim {
  DEFINE_ART_MODULE(SimPhotonsAna)
}
