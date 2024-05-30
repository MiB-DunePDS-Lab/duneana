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
#include "larcorealg/CoreUtils/counter.h"
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
      public: 
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

    struct SimPhotonsLabel_t {
      std::string fLabel = {}; 
      std::string fInstance = {};

      SimPhotonsLabel_t() {}
      SimPhotonsLabel_t(const std::string label, const std::string instance) :
        fLabel(label), fInstance(instance) {}
      SimPhotonsLabel_t(const std::string label) {
        size_t pos = label.find(':');
        if (pos != std::string::npos) {
            fLabel = label.substr(0, pos);
            fInstance = label.substr(pos + 1);
        } else {
          fLabel = label;
        }
      }
      inline std::string to_string() const {
        if (fInstance.empty()) return fLabel;

        std::string str_id = fLabel + ":" + fInstance;
        return str_id;
      }
    }; 

    // The parameters we'll read from the .fcl file.
    std::vector<SimPhotonsLabel_t> fSimPhotonsLiteModuleLabels = {}; // Input tags for SimPhotonsLite object
    std::vector<std::string> fSimPhotonsLiteLabelString = {};


    // Define the tree
    TTree* fPhotonsTree = {};
    TTree* fOpDetTree = {}; 

    // Define the variables that will go into the tree
    Int_t event_counter = 0;
    Float_t opDetH = 0.0; 
    Float_t opDetL = 0.0; 
    Float_t opDetW = 0.0; 
    Float_t opDetPos[3] = {0.0, 0.0, 0.0};
    UInt_t  opDetOrientation = 9;
    std::vector<size_t> opChannel; 

    Int_t fEventID = {};
    Int_t fOpChannel = {};
    std::string fSimPhotonsLabelID = {};
    
    std::vector<int> fTickTime = {}; // Tick time in ns (I think? pretty irrelevant)
    std::vector<int> fDetectedPhotonsCountPerTick = {}; // Number of detected photons per tick
    Int_t fDetectedPhotonsCount = {}; // Total number of detected photons

    // Object that we'll actually read from the artroot file
    std::map<int, int> fDetectedPhotons;

    void fill_opdet_tree(); 
    UInt_t get_opdet_orientation(const Float_t height, const Float_t length, const Float_t width) const;
  };

}

namespace sim {

  //-----------------------------------------------------------------------
  // Constructor
  SimPhotonsAna::SimPhotonsAna(Parameters const& pset) : art::EDAnalyzer{pset}, 
    fSimPhotonsLiteLabelString( pset().simphotons_labels() )
  {
    for (const auto& label : fSimPhotonsLiteLabelString ) {
      fSimPhotonsLiteModuleLabels.push_back( SimPhotonsLabel_t(label) ); 
    }

    art::ServiceHandle<art::TFileService const> tfs;

    // Make and add the branches to the tree
    fPhotonsTree = tfs->make<TTree>("PhotonsTree", "PhotonsTree");
    fPhotonsTree->Branch("EventID", &fEventID, "EventID/I");
    fPhotonsTree->Branch("OpChannel", &fOpChannel, "OpChannel/I");
    fPhotonsTree->Branch("TickTime", &fTickTime);
    fPhotonsTree->Branch("DetectedPhotonsCountPerTick", &fDetectedPhotonsCountPerTick);
    fPhotonsTree->Branch("DetectedPhotonsCount", &fDetectedPhotonsCount, "DetectedPhotonsCount/I");
    fPhotonsTree->Branch("SimPhotonsLabel", &fSimPhotonsLabelID);

    fOpDetTree = tfs->make<TTree>("opDetMap", "opDetMap");
    fOpDetTree->Branch("opDetH", &opDetH); 
    fOpDetTree->Branch("opDetL", &opDetL); 
    fOpDetTree->Branch("opDetW", &opDetW);
    fOpDetTree->Branch("opDetPos", &opDetPos, "opDetPos[3]/F");
    fOpDetTree->Branch("opDetOrientation", &opDetOrientation); 
    fOpDetTree->Branch("opDetCh", &opChannel); 

    return;
  }

  //-----------------------------------------------------------------------
    
  UInt_t SimPhotonsAna::get_opdet_orientation(Float_t height, Float_t length, Float_t width) const 
  {
    UInt_t orientation = 9;
    //Float_t local_height = 0; 
    //Float_t local_length = 0; 
    if (width > length) { // laterals along x-y plane, Z dimension smallest
      orientation = 2;
      //local_length = width;
    }
    else if (width > height) { // laterals along x-z plane, Y dimension smallest
      orientation = 1;
      //local_height = width;
    }
    else { // anode/cathode (default), X dimension smallest
      orientation = 0;
    }

    return orientation;
  }


  void SimPhotonsAna::fill_opdet_tree()
    {
      auto const* geom = lar::providerFrom<geo::Geometry>();
      UInt_t nOpDets = geom->NOpDets();
      for (size_t i : util::counter(nOpDets)) {
        opChannel.clear(); 
        geo::OpDetGeo const& opDet = geom->OpDetGeoFromOpDet(i);
        auto center = opDet.GetCenter();
        center.GetCoordinates( opDetPos );
        opDetH = opDet.Height();
        opDetW = opDet.Width();
        opDetL = opDet.Length();

        opDetOrientation = get_opdet_orientation(opDetH, opDetL, opDetW);

        size_t n_ch = geom->NOpHardwareChannels(i); 
        opChannel.resize(n_ch, 0); 
        for (size_t ich = 0; ich < n_ch; ich++) {
          opChannel.at(ich) = geom->OpChannel(i, ich);  
        }

        fOpDetTree->Fill();
      }

      return;
    }

  //-----------------------------------------------------------------------
  void SimPhotonsAna::analyze(const art::Event& evt)
  {
   
     //auto const* geo = lar::providerFrom<geo::Geometry>();

    if (event_counter == 0) {
      fill_opdet_tree(); 
    }

    fEventID = evt.id().event();

    // Get the SimPhotonsLite object
    art::Handle<std::vector<sim::SimPhotonsLite>> SimPhotonsLiteHandle;

    for (const auto& simphotons_label : fSimPhotonsLiteModuleLabels) {
      if ( simphotons_label.fInstance.empty() ) {
        evt.getByLabel(simphotons_label.fLabel, SimPhotonsLiteHandle);
      }
      else {
        evt.getByLabel(simphotons_label.fLabel, simphotons_label.fInstance, SimPhotonsLiteHandle);
      }
      fSimPhotonsLabelID = simphotons_label.to_string(); 
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

    event_counter++;
  }

} // namespace sim

namespace sim {
  DEFINE_ART_MODULE(SimPhotonsAna)
}
