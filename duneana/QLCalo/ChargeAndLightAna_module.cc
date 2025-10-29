// -*- mode: c++; c-basic-offset: 2; -*-
// This analyzer writes out a TTree for evaluating deposited energy using both Charge and Light
// M.Torti 

#ifndef ChargeAndLightAna_H
#define ChargeAndLightAna_H 1

// ROOT includes
#include "TH1.h"
#include "TEfficiency.h"
#include "TTree.h"
#include "TH3.h"
#include "TFile.h"

// C++ includes
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/SpacePointSolver/Solver.h"
#include "larsim/IonizationScintillation/ISCalcCorrelated.h"


// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
//DUNE includes
#include "duneopdet/OpticalDetector/OpFlashSort.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaUtilsBase.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"

namespace opdet {

  class ChargeAndLightAna : public art::EDAnalyzer{
  public:

    // Standard constructor and destructor for an ART module.
    ChargeAndLightAna(const fhicl::ParameterSet&);
    virtual ~ChargeAndLightAna();

    // This method is called once, at the start of the job. In this
    // example, it will define the histogram we'll write.
    void beginJob();

    // The analyzer routine, called once per event.
    void analyze (const art::Event&);

    void endJob();
    
    float GetFvisFromHisto(float x, float y, float z, float q);

  private:

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fEdepLabel;                // Input tag for Energy deposit collection
    std::string felecDriftLabel;           // Input tag for electron Drift collection
    std::string fOpFlashModuleLabel;       // Input tag for OpFlash collection
    std::string fOpHitModuleLabel;         // Input tag for OpHit collection
    std::string fSignalLabel;              // Input tag for the signal generator label
    std::string fGeantLabel;               // Input tag for GEANT
    std::string fHitsLabel;                //Input tag for Charge hits
    std::string fParticleModuleLabel;       //Input tag for reco Particles
    std::string fTrackLabel;                //Input tag for reco Tracks
    std::string fShowerLabel;                //Input tag for reco Showers
    std::string fHitToSpacePointLabel;      //Input tag for SpacePoints
    bool	fBeam;                      // Simulated events are beam neutrinos
    bool        fIsVD;                     // Is it a FD2-VD sample?

    TTree * fChargeLightTree;
    
    //LightMap variables
    TFile *fLightMap;
    TH3D *h3LightMap;
    Float_t  fFvis;
    
    // Parameters from the fhicl
    int   fNBinsE;
    float fLowE;
    float fHighE;
    int   fNBinsX;
    float fLowX;
    float fHighX;
    float fDistanceCut;

    Int_t fEventID;

    Float_t fTrueX;
    Float_t fTrueY;
    Float_t fTrueZ;
    Float_t fTrueT;
    Float_t fDetectedT;
    Float_t fTrueE;
    Int_t   fTruePDG;
    Int_t   fTrueCCNC;
    Float_t fRecoX;
    
//********Ion and scint
    Float_t fEdep;
    Float_t fTotEdep;
    Float_t fLdepSim; //track length from simulation
    Float_t fLdepGeom; // track lengt from geometrical calculation
    Int_t fTotEion; //number of ionizing electrons
    Int_t fTotGammaScint; //number of scintillation gammas
    Float_t fEiondE; //number of ionizing electrons divided by deposited energy
    Float_t fGammaScintdE; //number of scintillation gammas divided by deposited energy
   
    
//*********Charge from reco hits
    Float_t fTotalCharge; //Integral under the calibrated signal waveform of the hit, in tick x ADC units
    Float_t fTotalChargeCorr;//Total charge corrected for the electron lifetime
    Float_t fADCSum; //The sum of calibrated ADC counts of the hit (0. by default)
    Float_t fPeakAmplitude; //The estimated amplitude of the hit at its peak, in ADC units
    // Float_t fFirstHitTime; //First Time recorded among hits, Time of the signal peak,in us
    //Float_t fLastHitTime; //Last Time recorded among hits, Time of the signal peak, in us
    Float_t fMeanHitTime; //Last Time recorded among hits, Time of the signal peak, in us
    Float_t fHitDist; //First hit time * 1.6 mm/us --> theor distance travelled by charge (cm)
    
    std::vector< Float_t > fTruePxallpart;
    std::vector< Float_t > fTruePyallpart;
    std::vector< Float_t > fTruePzallpart;
    std::vector< Float_t > fTrueEallpart;
    std::vector< Int_t >   fTrueAllPDG;

    Int_t fNFlashes;
    /*std::vector< Int_t >   fFlashIDVector;
    std::vector< Float_t > fYCenterVector;
    std::vector< Float_t > fZCenterVector;
    std::vector< Float_t > fYWidthVector;
    std::vector< Float_t > fZWidthVector;
    std::vector< Float_t > fTimeVector;
    std::vector< Float_t > fRecoXVector;
    std::vector< Float_t > fTimeWidthVector;
    std::vector< Float_t > fTimeDiffVector;
    std::vector< Float_t > fTotalPEVector;
    std::vector< Float_t > fPurityVector;
    std::vector< Float_t > fDistanceVector; */
    Int_t fNOpDets;
    //std::vector<Int_t> fNHitOpDetVector;
    std::vector< Float_t > fOpHitPeakTime; //Time of OpHit in us (?)

    //*******Ion and Scint  
    std::vector< Float_t > fEnergyDepositionVector;
    std::vector< Float_t > fPointX;
    std::vector< Float_t > fPointY;
    std::vector< Float_t > fPointZ;
    std::vector< Int_t >   fGammaScint; //scint gamma each step
    std::vector< Int_t >   fPEperOpDet; //scint gamma each step
    std::vector< Float_t > fStepLength;
    std::vector< Float_t > fStepLCumVector;
    std::vector< Float_t > fStepEdepCumVector;

    std::vector< Float_t > fHitCharge;
    std::vector< Int_t > fHitMultiplicity; //How many hits could this one be shared with. Index of this hit among the Multiplicity() hits in the signal window
    std::vector< Float_t > fHitPeakTime; //Time of the signal peak, converted in us
    //std::vector< Float_t > fHitPeakTimeTicks; //Time of the signal peak, in tick units.
    
    Int_t    fFlashID;
    Float_t  fYCenter;
    Float_t  fZCenter;
    Float_t  fYWidth;
    Float_t  fZWidth;
    Float_t  fTime;
    Float_t  fTimeWidth;
    Float_t  fTimeDiff;
    Float_t  fTotalPE;    
    Float_t  fSumPE;//Sum of PE in all
    Float_t  fOpHitArea; //Sum of OpHit area
    Float_t  fPurity;
    Float_t  fDistance;
    Int_t    fNHitOpDets;
    std::vector< Float_t > fPEsPerOpDetVector;
    
        
    //Reco variables
    Int_t fNTrack;
    Float_t fPandoraVtxX,fPandoraVtxY,fPandoraVtxZ;
    //std::vector< Float_t > fTrkEnVector;
    //std::vector< Float_t > fTrkdEdxVector;
    std::vector< Float_t > fTrkLengthVector;
    
    //Float_t fSelTrkEn;
    Float_t fSelTrkLength;
    std::vector< Float_t > fSelTrkPointX;
    std::vector< Float_t > fSelTrkPointY;
    std::vector< Float_t > fSelTrkPointZ;
    
    std::vector< Float_t >  fHitToXVector;
    std::vector< Float_t >  fHitToYVector;
    std::vector< Float_t >  fHitToZVector;
    
    // SpacePoint data
    Short_t nspacepoints;
    std::vector<Float_t> fSpacePointX;   // X position of this SpacePoint
    std::vector<Float_t> fSpacePointY;   // Y position of this SpacePoint
    std::vector<Float_t> fSpacePointZ;   // Z position of this SpacePoint
    std::vector<Float_t> fHitSPCharge;   // charge of this SpacePoint

    // For counting waveforms
    std::string fOpDetWaveformLabel;
    float fBaseline;
    float fPE;
    TTree * fCountTree;
    Int_t fnwaveforms1pe;
    Int_t fnwaveforms2pe;
    Int_t fnwaveforms3pe;
    
  };

}

#endif // ChargeAndLightAna_H

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  ChargeAndLightAna::ChargeAndLightAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {

    // Indicate that the Input Module comes from .fcl
    fEdepLabel          = pset.get<std::string>("EdepLabel","IonAndScint");
    felecDriftLabel     = pset.get<std::string>("elecDriftLabel","elecDrift");
    fOpFlashModuleLabel = pset.get<std::string>("OpFlashModuleLabel");
    fOpHitModuleLabel   = pset.get<std::string>("OpHitModuleLabel");
    fHitsLabel          = pset.get<std::string>("HitsLabel");
    fParticleModuleLabel= pset.get<std::string>("ParticleModuleLabel");
    fShowerLabel        = pset.get<std::string>("ShowerLabel");
    fTrackLabel         = pset.get<std::string>("TrackLabel");
    fSignalLabel        = pset.get<std::string>("SignalLabel");
    fGeantLabel         = pset.get<std::string>("GeantLabel");
    fHitToSpacePointLabel  = pset.get<std::string>("HitToSpacePointLabel");
    fIsVD               = pset.get<bool>("IsVD");
    fBeam               = pset.get<bool>("Beam");
    fNBinsE             = pset.get<int>("NBinsE");
    fLowE               = pset.get<float>("LowE");
    fHighE              = pset.get<float>("HighE");
    fNBinsX             = pset.get<int>("NBinsX");
    fLowX               = pset.get<float>("LowX");
    fHighX              = pset.get<float>("HighX");
    fDistanceCut        = pset.get<float>("DistanceCut");

    fOpDetWaveformLabel = pset.get<std::string>("OpDetWaveformLabel","");
    fBaseline           = pset.get<float>("Baseline", 1500.);
    fPE                 = pset.get<float>("PE", 18.);

    art::ServiceHandle< art::TFileService > tfs;

    fChargeLightTree = tfs->make<TTree>("ChargeLightTree","ChargeLightTree");
    fChargeLightTree->Branch("EventID",                     &fEventID,   "EventID/I");
    fChargeLightTree->Branch("TrueX",                       &fTrueX,     "TrueX/F");
    fChargeLightTree->Branch("TrueY",                       &fTrueY,     "TrueY/F");
    fChargeLightTree->Branch("TrueZ",                       &fTrueZ,     "TrueZ/F");
    fChargeLightTree->Branch("TrueT",                       &fTrueT,     "TrueT/F");
    fChargeLightTree->Branch("DetectedT",                   &fDetectedT, "DetectedT/F");
    fChargeLightTree->Branch("TrueE",                       &fTrueE,     "TrueE/F");
    fChargeLightTree->Branch("TruePDG",                     &fTruePDG,   "TruePDG/I");
    fChargeLightTree->Branch("TrueCCNC",                    &fTrueCCNC,  "TrueCCNC/I");
    fChargeLightTree->Branch("NFlashes",                    &fNFlashes,  "NFlashes/I");
    /*fChargeLightTree->Branch("FlashIDVector",               &fFlashIDVector);
    fChargeLightTree->Branch("YCenterVector",               &fYCenterVector);
    fChargeLightTree->Branch("ZCenterVector",               &fZCenterVector);
    fChargeLightTree->Branch("YWidthVector",                &fYWidthVector);
    fChargeLightTree->Branch("ZWidthVector",                &fZWidthVector);
    fChargeLightTree->Branch("TimeVector",                  &fTimeVector);
    fChargeLightTree->Branch("TimeWidthVector",             &fTimeWidthVector);
    fChargeLightTree->Branch("TimeDiffVector",              &fTimeDiffVector);
    fChargeLightTree->Branch("TotalPEVector",               &fTotalPEVector); */
    fChargeLightTree->Branch("SumPE",                       &fSumPE,      "SumPE/F");
    fChargeLightTree->Branch("Fvis",                        &fFvis,      "Fvis/F");
    fChargeLightTree->Branch("NOpDets",                     &fNOpDets, "NOpDets/I");
   // fChargeLightTree->Branch("NHitOpDetVector",             &fNHitOpDetVector);
    fChargeLightTree->Branch("OpHitPeakTime",               &fOpHitPeakTime);
    fChargeLightTree->Branch("OpHitArea",    		      &fOpHitArea ,      "OpHitArea/F");
    //fChargeLightTree->Branch("Purity",                      &fPurityVector);
    //fChargeLightTree->Branch("Distance",                    &fDistanceVector);
    //fChargeLightTree->Branch("RecoXVector",                 &fRecoXVector);
    fChargeLightTree->Branch("TruePxallpart",               &fTruePxallpart);
    fChargeLightTree->Branch("TruePyallpart",               &fTruePyallpart);
    fChargeLightTree->Branch("TruePzallpart",               &fTruePzallpart);
    fChargeLightTree->Branch("TrueEallpart",                &fTrueEallpart);
    fChargeLightTree->Branch("TrueAllPDG",                  &fTrueAllPDG);
    fChargeLightTree->Branch("PointX",     		     &fPointX);
    fChargeLightTree->Branch("PointY",     		     &fPointY);
    fChargeLightTree->Branch("PointZ",     		     &fPointZ);
    fChargeLightTree->Branch("GammaScint",     	     &fGammaScint);
    fChargeLightTree->Branch("PEperOpDet",     	     &fPEperOpDet);
    fChargeLightTree->Branch("EnergyDepositionVector",      &fEnergyDepositionVector);
    fChargeLightTree->Branch("fStepLCumVector",             &fStepLCumVector);
    //fChargeLightTree->Branch("fStepEdepCumVector",          &fStepEdepCumVector);
    fChargeLightTree->Branch("TotEdep",                     &fTotEdep, "TotEdep/F");
    //fChargeLightTree->Branch("StepLengthEdep",              &fStepLength);
    fChargeLightTree->Branch("LdepSim",                     &fLdepSim,   "LdepSim/F");
    //fChargeLightTree->Branch("LdepGeom",                    &fLdepGeom,"LdepGeom/F");
    fChargeLightTree->Branch("TotEion",                     &fTotEion, "TotEion/I");
    fChargeLightTree->Branch("TotGammaScint",               &fTotGammaScint,"TotGammaScint/I");
    fChargeLightTree->Branch("HitCharge",                   &fHitCharge);    
    fChargeLightTree->Branch("TotalCharge",                 &fTotalCharge, "TotalCharge/F");
    fChargeLightTree->Branch("TotalChargeCorr",            &fTotalChargeCorr, "TotalChargeCorr/F");
    fChargeLightTree->Branch("ADCSum",                      &fADCSum,      "ADCSum/F");
    fChargeLightTree->Branch("PeakAmplitude",               &fPeakAmplitude,"PeakAmplitude/F");
    fChargeLightTree->Branch("HitMultiplicity",             &fHitMultiplicity);
    //fChargeLightTree->Branch("FirstHitTime",                &fFirstHitTime,  "FirstHitTime/F");
    //fChargeLightTree->Branch("LastHitTime",                 &fLastHitTime,  "LastHitTime/F");
    fChargeLightTree->Branch("MeanHitTime",                 &fMeanHitTime, "MeanHitTime/F");
    fChargeLightTree->Branch("HitPeakTime",                 &fHitPeakTime);
    //fChargeLightTree->Branch("HitPeakTimeTicks",            &fHitPeakTimeTicks);
    fChargeLightTree->Branch("HitDist",                     &fHitDist,  "HitDist/F");
    fChargeLightTree->Branch("NTrack",                      &fNTrack, "NTrack/I");
    fChargeLightTree->Branch("TrkLengthVector",             &fTrkLengthVector);
    //fChargeLightTree->Branch("TrkEnVector",                 &fTrkEnVector);
    //fChargeLightTree->Branch("TrkdEdxVector",               &fTrkdEdxVector);
    //fChargeLightTree->Branch("SelTrkEn",                    &fSelTrkEn, "SelTrkEn/F");
    fChargeLightTree->Branch("SelTrkLength",                &fSelTrkLength, "SelTrkLength/F");
    fChargeLightTree->Branch("SelTrkPointX",     	      &fSelTrkPointX);
    fChargeLightTree->Branch("SelTrkPointY",     	      &fSelTrkPointY);
    fChargeLightTree->Branch("SelTrkPointZ",     	      &fSelTrkPointZ);        
    fChargeLightTree->Branch("PandoraVtxX",                 &fPandoraVtxX,"PandoraVtxX/F");
    fChargeLightTree->Branch("PandoraVtxY",                 &fPandoraVtxY,"PandoraVtxY/F");
    fChargeLightTree->Branch("PandoraVtxZ",                 &fPandoraVtxZ,"PandoraVtxZ/F");
    fChargeLightTree->Branch("HitToXVector",		      &fHitToXVector);
    fChargeLightTree->Branch("HitToYVector",		      &fHitToYVector);
    fChargeLightTree->Branch("HitToZVector",		      &fHitToZVector);
    fChargeLightTree->Branch("SpacePointX",		      &fSpacePointX);
    fChargeLightTree->Branch("SpacePointY",		      &fSpacePointY);
    fChargeLightTree->Branch("SpacePointZ",		      &fSpacePointZ);
    fChargeLightTree->Branch("HitSPCharge",                 &fHitSPCharge);


    if (!fOpDetWaveformLabel.empty()) {
      fCountTree = tfs->make<TTree>("CountWaveforms","CountWaveforms");
      fCountTree->Branch("EventID",       &fEventID,      "EventID/I");
      fCountTree->Branch("nwaveforms1pe", &fnwaveforms1pe, "nwaveforms1pe/I");
      fCountTree->Branch("nwaveforms2pe", &fnwaveforms2pe, "nwaveforms2pe/I");
      fCountTree->Branch("nwaveforms3pe", &fnwaveforms3pe, "nwaveforms3pe/I");
    }

  }

  //-----------------------------------------------------------------------
  // Destructor
  ChargeAndLightAna::~ChargeAndLightAna()
  {}

  //-----------------------------------------------------------------------
  void ChargeAndLightAna::beginJob()
  {
  
   std::cout << "Opening the light map file" << std::endl;
   //Open the light map rootfile   
   if (fIsVD) {
    fLightMap = new TFile ("/dune/data2/users/dguffant/QLcalo/combined_visibility_map.root", "READ");     
      if (fLightMap->IsOpen()) h3LightMap = (TH3D *) fLightMap->Get("h3VisMap_Ar_Xe10ppm");
      else std::cout << "Light Map VD not found!!" << std::endl; 
      std::cout<< "Light Map VD found!" << std::endl;
     } //open the light map for the VD
   
    else { 
      fLightMap = new TFile ("/dune/app/users/dguffant/test001/lightmap_h3.root","READ");     
       if (fLightMap->IsOpen()) h3LightMap = (TH3D *) fLightMap->Get("h3VisMap");
       else std::cout << "Light Map HD not found!!" << std::endl;
       std::cout<< "Light Map HD found!" << std::endl;
      } //light map for the HD
         
  }

  //-----------------------------------------------------------------------
  void ChargeAndLightAna::analyze(const art::Event& evt)
  {
    // Get the required services
    art::ServiceHandle< geo::Geometry > geom;
    art::ServiceHandle< cheat::PhotonBackTrackerService > pbt;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle< cheat::ParticleInventoryService > pinv;
    art::ServiceHandle< art::TFileService > tfs;
    //bool isMC = !evt.isRealData();

    //pbt->Rebuild(evt);


    // Record the event ID
    fEventID = evt.id().event();
    
  std::cout << "=============== EVENT ID " << fEventID << " ================" << std::endl;

    ///////////////////////////////////////////////////
    // Count waveforms if a waveform label was given //
    ///////////////////////////////////////////////////

    if (!fOpDetWaveformLabel.empty()) {
      fnwaveforms1pe = 0;
      fnwaveforms2pe = 0;
      fnwaveforms3pe = 0;
      art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
      if (evt.getByLabel(fOpDetWaveformLabel, wfHandle)) {
        fnwaveforms1pe = wfHandle->size();

        for (auto wf: *wfHandle) {
          auto it = max_element(std::begin(wf), std::end(wf));
          double peak = *it - fBaseline;
          if ( peak > (1.5*fPE)) {
            ++fnwaveforms2pe;

            if ( peak > (2.5*fPE) )
              ++fnwaveforms3pe;
          }
        }
        fCountTree->Fill();
      }
    }


    //////////////////////////////////////
    // Access all the Flash Information //
    //////////////////////////////////////

    // Get flashes from event
    art::Handle< std::vector< recob::OpFlash > > FlashHandle;
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    if (evt.getByLabel(fOpFlashModuleLabel, FlashHandle)) {
      art::fill_ptr_vector(flashlist, FlashHandle);
      std::sort(flashlist.begin(), flashlist.end(), recob::OpFlashPtrSortByPE);
    }
    else {
      mf::LogWarning("ChargeAndLightAna") << "Cannot load any flashes. Failing";
      return;
    }
    
 
     //////////////////////////////////////
    // Access all the OpHit Information //
    //////////////////////////////////////
       
    art::Handle< std::vector< recob::OpHit > > HitHandle;
    std::vector<art::Ptr<recob::OpHit> > hitlist;
    if (evt.getByLabel(fOpHitModuleLabel, HitHandle)) {
      art::fill_ptr_vector(hitlist, HitHandle);
    }

    // Get total PE in all Ophits
    fSumPE = 0;
    fOpHitArea = 0;
    for (auto hit: hitlist) {
        fSumPE += hit->PE();
        fOpHitPeakTime.emplace_back(hit->PeakTime());  
        fOpHitArea += hit->Area();                 
        }
        
      fPEperOpDet.clear();
      for(unsigned int iOD = 0; iOD < geom->NOpDets(); ++iOD){
        fPEperOpDet.emplace_back(0);
      }
      
     unsigned int iC = 0;   
     for (auto hit:hitlist) {     
      iC = hit->OpChannel();
      unsigned int iOD = geom->OpDetFromOpChannel(iC);
      fPEperOpDet[iOD] += hit->PE();
     }   

    // Get assosciations between flashes and hits
    //art::FindManyP< recob::OpHit > Assns(flashlist, evt, fOpFlashModuleLabel);

    /////////////////////////
    // G4 Deposited Energy
    //////////////////////
    fTotEdep = 0;
    
    ///// From IonAndScint
    fTotGammaScint = 0;
    fTotEion = 0;
    fLdepSim=0;
    fLdepGeom=0;
    fGammaScintdE = 0;
    fEiondE = 0;
    
    std::cout<<"Looking for deposited energy of IonAndScint"<<std::endl;
    int nSimEnergyDeposits = 0;
    double dedxsteps=0;
    double dEdx=0;
    	
    art::Handle< std::vector<sim::SimEnergyDeposit> > energyDepositHandle;
    std::vector<art::Ptr<sim::SimEnergyDeposit> > energyDepositlist;
      if(evt.getByLabel(fEdepLabel, energyDepositHandle)){
	art::fill_ptr_vector(energyDepositlist, energyDepositHandle);
	nSimEnergyDeposits = energyDepositlist.size();

std::cout<< "nSimEnergyDeposits " << nSimEnergyDeposits << std::endl;

if (nSimEnergyDeposits == 0) return;

	//look for start point and end point of total deposition, make a scan in Z (beam dir) and get corresponding x,y
	double edepstartx = 0;
	double edepstarty = 0;
	double edepstartz = energyDepositlist[0]->StartZ();
	double edependx = 0;
	double edependy = 0;
	double edependz =  energyDepositlist[0]->EndZ();
	double thisposition = 0;
	int istart = 0;
	int iend = 0;


       for(int i = 0; i < nSimEnergyDeposits; i++){
       
          //points of the simulated track
          fPointX.emplace_back(energyDepositlist[i]->StartX());
          fPointY.emplace_back(energyDepositlist[i]->StartY());
          fPointZ.emplace_back(energyDepositlist[i]->StartZ());
          
          fGammaScint.emplace_back(energyDepositlist[i]->NumPhotons());
              
	  fEdep = energyDepositlist[i]->E(); 
	  fEnergyDepositionVector.emplace_back(fEdep);
	  fTotEdep += fEdep;
          fStepEdepCumVector.emplace_back(fTotEdep);
	  fTotGammaScint += energyDepositlist[i]->NumPhotons();
	  fTotEion += energyDepositlist[i]->NumElectrons();
	  //find starting point, scan along z
	  thisposition = energyDepositlist[i]->StartZ();
	  if(thisposition < edepstartz ){
	    edepstartz=thisposition;
	    istart=i;
	  }
	  //find end point, scan along z
	  thisposition = energyDepositlist[i]->EndZ();
	  if(thisposition > edependz ){
	    edependz=thisposition;
	    iend=i;
	  }
	  //save single steps of edep length
	  fStepLength.emplace_back( energyDepositlist[i]->StepLength()); //in cm
	 //sum of single steps of edep length
	  fLdepSim += energyDepositlist[i]->StepLength(); //in cm
	  fStepLCumVector.emplace_back(fLdepSim); //cumulative of energy in each step
	  
	  dedxsteps += fEdep/fLdepSim; //in MeV/cm
	}
	
	std::cout<<"Number of deposits: "<<nSimEnergyDeposits<<" Single step L(cm): "<<energyDepositlist[2]->StepLength()<<" Total number of emittend photons: "<<fTotGammaScint<<" and ionization electrons: "<<fTotEion<<" ---> Check w.r.t path length: "<<std::endl;

	edepstartx=energyDepositlist[istart]->StartX();
	edepstarty=energyDepositlist[istart]->StartY();
	edependx=energyDepositlist[iend]->EndX();
	edependy=energyDepositlist[iend]->EndY();
	double diffx=TMath::Power((edependx-edepstartx),2);
	double diffy=TMath::Power((edependy-edepstarty),2);
	double diffz=TMath::Power((edependz-edepstartz),2);
	fLdepGeom=TMath::Sqrt(diffx+diffy+diffz);

	//std::cout<<"Scan in Z, Start of energy deposition (x,y,z): "<<edepstartx<<" "<<edepstarty<<" "<<edepstartz<<" End (x,y,z): "<<edependx<<" "<<edependy<<" "<<edependz<<" cm --> length: "<<fLdepGeom<<" cm"<<std::endl;
	//std::cout<<"Gamma/path length: "<<fTotGammaScint/fLdepSim<<" Ion. electrons/path length: "<<fTotEion/fLdepSim<<std::endl;
	fGammaScintdE=fTotGammaScint/fTotEdep;
	fEiondE=fTotEion/fTotEdep;
	//std::cout<<"Gamma/TotEdep: "<<fGammaScintdE<<" Ion. electrons/TotEdep: "<<fEiondE<<std::endl;
	//std::cout<<"(Gamma/TotEdep)/dx: "<<fGammaScintdE/fLdepSim<<" (Ion.electrons/TotEdep)/dx: "<<fEiondE/fLdepSim<<std::endl;
	dEdx=fTotEdep/fLdepSim;
	//std::cout<<"Sum of every step de/dx: "<<dedxsteps<<" MeV/cm, ...divinding TotEdep by Total L--> "<<std::endl;  
	std::cout<<"Total dE/dX: "<<dEdx<<" (dx is sum of steps L). MeV/cm"<<std::endl;
	
      }else{
	mf::LogWarning("ChargeAndLightAna") << "Cannot Find Deposited Energy. Failing";
	return;
      }

    std::cout<<"TotEdep is: "<<fTotEdep<<" (MeV) ";
              
     
    //////////////////////////////////////
    // Access all the truth information //
    //////////////////////////////////////

    std::set<int> signal_trackids;
    geo::PlaneID planeid;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    try {
      auto MClistHandle = evt.getValidHandle<std::vector<simb::MCTruth> >(fSignalLabel);

      art::Ptr<simb::MCTruth> mctruth(MClistHandle, 0);
      if (mctruth->NParticles() == 0) {
        mf::LogError("ChargeAndLightAna") << "No MCTruth Particles";
      }

      // Get all the track ids associated with the signal event.
      art::FindManyP<simb::MCParticle> SignalGeantAssns(MClistHandle,evt,fGeantLabel);
      for ( size_t i = 0; i < SignalGeantAssns.size(); i++) {
        auto parts = SignalGeantAssns.at(i);
        for (auto part = parts.begin(); part != parts.end(); part++) {
          signal_trackids.emplace((*part)->TrackId());
        }
      }

      // Get just the neutrino, entry 0 from the list, and record its properties
      const simb::MCParticle& part(mctruth->GetParticle(0));
      fTrueX     = part.Vx();
      fTrueY     = part.Vy();
      fTrueZ     = part.Vz();
      fTrueT     = part.T()*1000; // ns -> us
      fTrueE     = part.E();
      fTruePDG   = part.PdgCode();
      
      std::cout << "Vertex X " << fTrueX << " Y " << fTrueY << " Z " << fTrueZ << std::endl;
      std::cout << "True Energy " << fTrueE << " PDG code " << fTruePDG <<  std::endl;
      
      //CCNC only for the neutrino otherwise set to 1000.
      // 0=CC 1=NC
      int intType=-1000.;
      if (mctruth->Origin() == simb::kBeamNeutrino){
	fTrueCCNC  = mctruth->GetNeutrino().CCNC();
	intType=mctruth->GetNeutrino().InteractionType();
	std::cout<<"Interaction Type is: "<<intType;
	if(intType==1091) std::cout<<" --> CCDIS"<<std::endl;
	if(intType==1001) std::cout<<" --> CCQE"<<std::endl;
	if(intType==1003) std::cout<<" --> ResCCNuProtonPiPlus"<<std::endl;
	if(intType==1004) std::cout<<" --> ResCCNuNeutronPi0"<<std::endl;
	if(intType==1005) std::cout<<" --> ResCCNuNeutronPiPlus"<<std::endl;
	std::cout<<std::endl;
      }else{
	fTrueCCNC=1000;
      }

      if (mctruth->Origin() == simb::kBeamNeutrino && fTruePDG==14 &&  fTrueCCNC==0) { std::cout<<" *********** numu CC ********* "<<std::endl;
      std::cout<<"TrueE is: "<<fTrueE<<" (GeV)"<<std::endl;}
      
      if (mctruth->Origin() == simb::kBeamNeutrino && fTruePDG==12 &&  fTrueCCNC==0) { std::cout<<" *********** nue CC ********* "<<std::endl;
      std::cout<<"TrueE is: "<<fTrueE<<" (GeV)"<<std::endl;}
      
      // Get all the paricle including neutrino, and record its properties
      unsigned int const nParticles = mctruth->NParticles();
      std::cout<<"There are: "<<nParticles<<" secondary particles"<<std::endl;

      for (unsigned int i = 0; i < nParticles; ++i) {
	simb::MCParticle const& particle = mctruth->GetParticle(i);
        fTruePxallpart    .emplace_back(particle.Px());
        fTruePyallpart    .emplace_back(particle.Py());
        fTruePzallpart    .emplace_back(particle.Pz());
        fTrueEallpart     .emplace_back(particle.E());
        fTrueAllPDG       .emplace_back(particle.PdgCode());
      }

      
      // Get the PlaneID which describes the location of the true vertex
      int plane = 0;
      geo::TPCID tpc = geom->FindTPCAtPosition(geo::Point_t{part.Vx(), part.Vy(), part.Vz()});
      if (! geom->HasTPC(tpc) ) {
        if(fIsVD){ //allowing Flashes to be recorded outside active volume for VD but with no match
          std::cout << "No valid TPC" << std::endl;
          fDetectedT=-1;
        } else{ //Flash match only if inside active volume for HD                  
           mf::LogWarning("ChargeAndLightAna") << "No valid TPC for " << tpc;
           std::cout << "No find TPC at Position vertex" << std::endl;
           return;
      }
      } else {
      geo::PlaneID tempid(tpc, plane);
      planeid = tempid;
      std::cout<<" Plane ID from true Vertex: "<<planeid<<" TPC: "<<tpc<<'\n';
      
      // Convert true X to would-be charge arrival time, and convert from ticks to us, add to MC time
      double deltaTicks = detProp.ConvertXToTicks(part.Vx(), planeid);
      double deltaT = clockData.TPCTick2Time(deltaTicks);
      fDetectedT = fTrueT + deltaT;
      std::cout<<"True Vertex X to get the deltaTicks: "<<part.Vx()<<'\n';
      std::cout<<"MC True T: "<< fTrueT<<" ticks and coversion to us: "<<deltaTicks<<" "<<deltaT<< '\n';
      std::cout<<"Detected time (would-be charge arrival time) (us): "<< fDetectedT << '\n';
    }
    }
    catch (art::Exception const& err) 
    {
      // If the error isn't that a product wasn't found, throw it back up.
      if ( err.categoryCode() != art::errors::ProductNotFound ) throw;

      // Otherwise, this event just doesn't have signal. Fill default values
      fTrueX = 0;
      fTrueY = 0;
      fTrueZ = 0;
      fTrueT = 0;
      fTrueE = 0;
      fTruePDG = 0;
      fTrueCCNC = 0;
      fDetectedT = 0;

      mf::LogError("ChargeAndLightAna") << "Event doesn't have signal ";
    }

    // Get the maximum possible time difference by getting number of ticks corresponding to
    // one full drift distance, and converting to time.
    //double maxT=0;
    //if(fBeam) maxT = clockData.TPCTick2Time(detProp.NumberTimeSamples());


    //////////////////////////////////////
    // CHARGE collected                 //
    ////////////////////////////////////// 
//    unsigned short wirePlane = geo::kZ;

    // Total charge collected in the event and time info of the hits
    fTotalCharge = 0.0;
    fTotalChargeCorr = 0.0;
    fADCSum = 0.0;
    fPeakAmplitude = 0.0;
    //fFirstHitTime = -1000;
    fHitDist = -1000;
    fMeanHitTime = -1000;
    //double firstHit=100000;
    double meanHit=0;
    //double lastHit=-100000;
    int collpl=0;
    auto hitListHandle = evt.getValidHandle<std::vector<recob::Hit>>(fHitsLabel);

    for(size_t iHit = 0; iHit < hitListHandle->size(); ++iHit){
      //thishit=static_cast<int>(iHit);  
      art::Ptr<recob::Hit> hitPtr(hitListHandle, iHit);
	if(geom->SignalType(hitPtr->Channel()) == geo::kCollection){
     // if(hitPtr->View() == wirePlane){
	collpl++;
	fHitCharge.emplace_back(hitPtr->Integral());
	fTotalCharge += hitPtr->Integral();        
	fADCSum += hitPtr->SummedADC();
	fPeakAmplitude += hitPtr->PeakAmplitude();
	fHitMultiplicity.emplace_back(hitPtr->Multiplicity());
	//fHitPeakTimeTicks.emplace_back(hitPtr->PeakTime());
	//convert hit time from ticks to us
	double hitT = clockData.TPCTick2Time(hitPtr->PeakTime());
	fHitPeakTime.emplace_back(hitT);
	meanHit += hitT;
	//	if(thishit==0) firstHit=hitT;
	/*	if(static_cast<int>(iHit) == 0){
	  firstHit=hitT;
	  lastHit=hitT;
	  }*/
	//if(hitT<firstHit) firstHit=hitT;
	//if(hitT>lastHit) lastHit=hitT;
      } //hit in collection                 
    } //hit
  
    if(collpl==0){
      std::cout<<"********** no hits in the Collection Plane... ***********"<<'\n';
      //  fFirstHitTime = -1000;
      //fLastHitTime = -1000;
      fMeanHitTime = -1000;
    }else{
      // fFirstHitTime = firstHit;
      //fLastHitTime = lastHit;
      fMeanHitTime = meanHit/collpl;
      fHitDist = (fMeanHitTime * 1.6)/10.; //cm
      //std::cout << "First Hit in time is " << firstHit <<" (us), Last Hit in time is " << lastHit <<" (us)" << '\n';
      std::cout << "Number of hits: " << collpl <<" Mean hit time is " << fMeanHitTime <<" (us)" << " --> should correspond to charge travelling: "<<fHitDist<<" cm"<< '\n';    
    }
        
    // Output the total charge
    std::cout << "Total Charge from Reco Hits: " << fTotalCharge << '\n';    
    
    
   //compute the correction to the collected charge using DUNEAnaHitUtils
   std::vector<art::Ptr<recob::Hit>> HitsInColl = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsLabel),2); //Collection plane has Plane_ID = 2
   fTotalChargeCorr = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, HitsInColl);
    std::cout << "Total Charge from reco hits corrected with electron lifetime: " << fTotalChargeCorr << '\n';    
         
               
    /////////////////////////
    // Analyze the flashes //
    /////////////////////////

    // Set up some flags to fill as we loop
    // through flashes.
 
    // For every OpFlash in the vector
    fNOpDets   = geom->NOpDets();
    fNFlashes  = flashlist.size();
    std::cout << "Number of flashes " << fNFlashes << std::endl;
    
 /*   for(unsigned int i = 0; i < flashlist.size(); ++i)
    {
      // Get OpFlash and associated hits
      recob::OpFlash TheFlash = *flashlist[i];
      art::Ptr<recob::OpFlash> FlashP = flashlist[i];
      std::vector< art::Ptr<recob::OpHit> > hitFromFlash = pbt->OpFlashToOpHits_Ps(FlashP);
      std::vector< art::Ptr<recob::OpHit> > matchedHits = pbt->OpFlashToOpHits_Ps(FlashP);

      // Calculate the flash purity
      double purity = pbt->OpHitCollectionPurity(signal_trackids, matchedHits);

      // Calcuate relative detection time
      //if event is not from neutrino beam add true MC time
      double flashrealtime= TheFlash.Time() + fTrueT;
      //       double timeDiff = fDetectedT - TheFlash.Time();
      double timeDiff = fDetectedT - flashrealtime;

      if (!planeid) {
        // planeid isn't valid
        fRecoX = 0;
      }
      else {
        double ticks = clockData.Time2Tick(timeDiff);
        fRecoX = detProp.ConvertTicksToX(ticks, planeid);
      }

      // Check if this is a possible flash (w/in 1 drift window)
        //if (fBeam) { 
        if (timeDiff < -10 || timeDiff > maxT){
	  mf::LogError("ChargeAndLightAna") << "Skipping Flash: not w/in 1 drift window (timeDiff < -10 OR timeDiff > maxT)";
	  continue;
	} //} 
*/
      //*************************************************************
/*      
      // Put flash info into variables
      fFlashID     = i;
      fYCenter     = TheFlash.YCenter();
      fZCenter     = TheFlash.ZCenter();
      fYWidth      = TheFlash.YWidth();
      fZWidth      = TheFlash.ZWidth();
      fTime        = TheFlash.Time();
      fTimeWidth   = TheFlash.TimeWidth();
      fTimeDiff    = timeDiff;
      fTotalPE     = TheFlash.TotalPE();
      fPurity      = purity;
      std::cout<<"Flash: "<<fFlashID<<" PE of this Flash: "<<fTotalPE<<'\n';
      std::cout<<"Flash time: "<<fTime<<'\n';
      std::cout<<"True X: "<<fTrueX<<" RecoX from ticks "<<fRecoX<< '\n';

      // Calculate distance from MC truth vertex in the Y-Z plane
      fDistance = sqrt( pow(fTrueY-fYCenter,2) +  pow(fTrueZ-fZCenter,2) );
      std::cout<<"Flash distance from MC truth vertex in the Y-Z plane: "<<fDistance<< '\n';
      std::cout<<"TrueX+flash distance gives: "<< sqrt( pow(fTrueX,2) +  pow(fDistance,2) )<< '\n';

      // Loop through all the opdets with hits in this flash
      fPEsPerOpDetVector.clear();
          
      for(unsigned int iOD = 0; iOD < geom->NOpDets(); ++iOD){
        fPEsPerOpDetVector.emplace_back(0);
      }
      
      for(unsigned int iC=0; iC < geom->NOpChannels(); ++iC)
      {
        unsigned int iOD = geom->OpDetFromOpChannel(iC);
        fPEsPerOpDetVector[iOD] += TheFlash.PE(iC);
      }

      fNHitOpDets = 0;
      for(unsigned int iOD = 0; iOD < geom->NOpDets(); ++iOD){
        if (fPEsPerOpDetVector[iOD] > 0) ++fNHitOpDets;
      }
      fNHitOpDetVector.emplace_back(fNHitOpDets);

      // Add flash info to the tree of all possible flashes
      fFlashIDVector    .emplace_back(fFlashID);
      fYCenterVector    .emplace_back(fYCenter);
      fZCenterVector    .emplace_back(fZCenter);
      fYWidthVector     .emplace_back(fYWidth);
      fZWidthVector     .emplace_back(fZWidth);
      fTimeVector       .emplace_back(fTime);
      fTimeWidthVector  .emplace_back(fTimeWidth);
      fTimeDiffVector   .emplace_back(fTimeDiff);
      fTotalPEVector    .emplace_back(fTotalPE);
      fPurityVector     .emplace_back(fPurity);
      fDistanceVector   .emplace_back(fDistance);
      fRecoXVector      .emplace_back(fRecoX);
    } */
    
    ///////////////////////////////////////
    // Reco tracks                       //
    ///////////////////////////////////////    
    
    const std::vector<art::Ptr<recob::PFParticle>> particles =  dune_ana::DUNEAnaEventUtils::GetPFParticles(evt,fParticleModuleLabel);
    std::cout <<"---- Select Pandora tracks ---- " << std::endl; 

     int nPart = 0;
     fNTrack = 0;
     double maxL = 0;
     //double maxE = 0;
     int selTrack = 0;
     //fSelTrkEn = 0;
     fSelTrkLength = 0;   

     for (const art::Ptr<recob::PFParticle> &particle : particles){
     nPart++;

    //Loop over all Tracks in the event
     if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(particle,evt,fParticleModuleLabel,fTrackLabel)){
	const art::Ptr<recob::Track> trk = dune_ana::DUNEAnaPFParticleUtils::GetTrack(particle,evt,fParticleModuleLabel,fTrackLabel);
	fNTrack++;
	
	fTrkLengthVector.emplace_back(trk->Length());
	
	std::cout<<" RECONSTRUCTED TRACK: nPart: "<<nPart<<" track ID: "<<trk->ID()<<" particle ID: "<<trk->ParticleId() << std::endl;
	TVector3 trackvtx = trk->Vertex<TVector3>(); 
	std::cout<<" Track Vertex from recob::track : "<<trackvtx(0)<<" "<<trackvtx(1)<<" "<<trackvtx(2)<< " Length: "<<trk->Length()<<std::endl;
	
	/*for(unsigned int i=0; i<trk->Energy().size(); i++){
	  std::cout << trk->Energy().at(i)<<" "<<trk->dEdx().at(i)<<std::endl;
	  fTrkEnVector.emplace_back(trk->Energy().at(i));
	  fTrkdEdxVector.emplace_back(trk->dEdx().at(i));
	}*/
	
	//Select to longest track
        if(trk->Length() > maxL){
	  maxL = trk->Length();
	  //maxE = trk->Energy().at(2);
	  selTrack = nPart;
	}
					
      }//end if IsTrack
      else{ mf::LogWarning("ChargeAndLightAna") << "Cannot Find Recontructed track in this event. Failing" << std::endl; } 
     } //end PFParticles
     
     //Information on Selected Track (the longest in this case)
     if (fNTrack == 0 ) { std::cout<<"This event doesn't have any reco track"<<std::endl;}
     else{ //std::cout << "The Longest track info " << std::endl;
     
      int l = 0;
      int npoint_trk = 0;
       
      for (const art::Ptr<recob::PFParticle> &part : particles){
       l++;
        if (l == selTrack) {
        const art::Ptr<recob::Track> seltrk = dune_ana::DUNEAnaPFParticleUtils::GetTrack(part,evt,fParticleModuleLabel,fTrackLabel);
        
          //fSelTrkEn=seltrk->Energy().at(2);	
	  fSelTrkLength=seltrk->Length(); 
	  
	  //save track points
	  npoint_trk=seltrk->CountValidPoints();
	  std::cout << "Number of points in the selected track " << npoint_trk << std::endl;
	   for(int il=0; il<npoint_trk;il++){
            fSelTrkPointX.emplace_back(seltrk->LocationAtPoint(il).X());
            fSelTrkPointY.emplace_back(seltrk->LocationAtPoint(il).Y());
            fSelTrkPointZ.emplace_back(seltrk->LocationAtPoint(il).Z());
            } 
        
        //get selcted track Pandora Vertex
	if(dune_ana::DUNEAnaPFParticleUtils::GetVertex(part,evt,fParticleModuleLabel)->isValid()){
	  const art::Ptr<recob::Vertex> trkvtx = dune_ana::DUNEAnaPFParticleUtils::GetVertex(part,evt,fParticleModuleLabel);
	      recob::Track::Point_t trkvtxpos =trkvtx->position();
	      fPandoraVtxX=trkvtxpos.X();
	      fPandoraVtxY=trkvtxpos.Y();
	      fPandoraVtxZ=trkvtxpos.Z();
	      std::cout<<" Vertex from PandoraVertex X Y Z: "<<fPandoraVtxX<<" "<<fPandoraVtxY<<" "<<fPandoraVtxZ<<std::endl;
	    }else{ std::cout<<"Vertex fron Pandora GetVertex status is not valid"<<std::endl; }
        
        } //end if l seltrack
       }
     } //end else selcted longest track
     
     ///////////////////////////////////////////////
     // Use of BackTracker Service		   //
    ///////////////////////////////////////////////
   //following a Protodune example
              
     for(size_t iH = 0; iH < hitListHandle->size(); ++iH){
      art::Ptr<recob::Hit> hitPtr(hitListHandle, iH);
      //std::cout << "Backtracker per hit " << iH << std::endl;
     // if (isMC) {
      std::vector<const sim::IDE*> ides;
          try{          
            ides = bt_serv -> HitToSimIDEs_Ps(clockData, hitPtr);
          }
          catch(...){}
                 
      if (ides.size()>0) { 
         std::vector<double> SingleHitToXYZVector = bt_serv->HitToXYZ(clockData, hitPtr);
         fHitToXVector.emplace_back(SingleHitToXYZVector[0]);
         fHitToYVector.emplace_back(SingleHitToXYZVector[1]);
         fHitToZVector.emplace_back(SingleHitToXYZVector[2]); 
         
         //SingleHitToXYZVector.clear();     
	}
      //} //end if mc
    }//for hitlist
    
    
    /////////////////////////
    // Space Points	   //
    /////////////////////////
    std::cout<<std::endl;
    std::cout<<"****Space points**** "<< std::endl;    
        
/*
//Prova1
//Get all the spacepoints for the hits of the collection plane
   for (unsigned int i = 0; i < HitsInColl.size(); ++i){
     const std::vector<art::Ptr<recob::SpacePoint> > spacePoints(dune_ana::DUNEAnaHitUtils::GetSpacePoints(HitsInColl[i],evt,fHitsLabel,fHitToSpacePointLabel));
     for (unsigned int iSpacePoint = 0; iSpacePoint < spacePoints.size(); ++iSpacePoint){
       const art::Ptr<recob::SpacePoint> spacePoint(spacePoints[iSpacePoint]);
       double thispointx=spacePoint->XYZ()[0];
       double thispointy=spacePoint->XYZ()[1];
       double thispointz=spacePoint->XYZ()[2];
       std::cout<<"The space point is in x,y,z: "<<thispointx<<" "<<thispointy<<" "<<thispointz<<std::endl;
     }//loop on spacepoints collection
   }//loop on hits
  
*/
/*
//prova2
double sppointx, sppointy, sppointz;

//const std::vector<art::Ptr<recob::PFParticle>> particles = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt,fParticleModuleLabel);

for (const art::Ptr<recob::PFParticle> &particle : particles){

  std::vector<art::Ptr<recob::SpacePoint>> spacePoints= dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(particle,evt,fParticleModuleLabel);
	std::cout<<"Get Spacepoints from Pandora"<<std::endl;
   	if(spacePoints.size()>0){
	  std::cout<<"Found "<<spacePoints.size()<<" space points"<<std::endl;
	  for (unsigned int iSpacePoint = 0; iSpacePoint < spacePoints.size(); ++iSpacePoint){
	    //spacePoint = spacePoints[iSpacePoint];
	     const art::Ptr<recob::SpacePoint> spacePoint(spacePoints[iSpacePoint]);
	    sppointx=spacePoint->XYZ()[0];
	    sppointy=spacePoint->XYZ()[1];
	    sppointz=spacePoint->XYZ()[2];
	    std::cout<<"The space point is in x,y,z: "<<sppointx<<" "<<sppointy<<" "<<sppointz<<std::endl;
	  }//loop on spacepoints collection
	}//if spacepoints
	else{
	  std::cout<<"space point vector empty"<<std::endl;
	}
	
	}//end pfparticles        
 */   
 
 //prova 3 mt 
 int spcollp = 0;
 Float_t SPtotcharge = 0;
 Float_t spx, spy, spz, spq;
 Float_t fvis_point = 0;
 Float_t fvis_tmp = 0;
 //double binval=0;
 //int ibin = 0;
 fFvis = 0;
   
    auto SPHandle = evt.getValidHandle< std::vector<recob::SpacePoint> >(fHitToSpacePointLabel);
    art::FindManyP<recob::Hit> hitsFromSP(SPHandle, evt, fHitToSpacePointLabel);
    std::cout<<"Total number of Space Points in this event: "<< SPHandle->size()<<std::endl;
   
        for (size_t i = 0; i < SPHandle->size(); ++i){
      auto collhits = hitsFromSP.at(i);
      for(auto & h : collhits){
      //std::cout << "SP in collection " << std::endl;
	if (h->SignalType() == geo::kCollection){
	  spcollp++;
	  //std::cout << "SP numero" << spcollp<<  std::endl;
	  art::Ptr<recob::SpacePoint> hitSP(SPHandle, i);
	  spx = hitSP->XYZ()[0];
	  spy = hitSP->XYZ()[1];
	  spz = hitSP->XYZ()[2];
	  spq = h->Integral();
	  
	  fSpacePointX.emplace_back(spx);
	  fSpacePointY.emplace_back(spy);
	  fSpacePointZ.emplace_back(spz); 
	  fHitSPCharge.emplace_back(spq);
	  SPtotcharge += spq;
	  
	  //std::cout << "spx " << spx << ", spy " << spy << ", spz " << spz <<", spq " << spq << std::endl;
	  
	  
	  //std::cout << "Visibility " << std::endl;
	  //ibin = h3LightMap->FindBin(spx, spz, spy); 
          //binval = h3LightMap->GetBinContent(ibin); 
    
          //std::cout<<" ---> f is: "<<binval<< std::endl;
          //fvis_point = binval*spq;


	  fvis_point = GetFvisFromHisto(spx, spy, spz, spq);
	  //std::cout << "fvis_point = " << fvis_point << " point " << spcollp << std::endl;
	  fvis_tmp += fvis_point;
	  //std::cout << "fvis_tmp = " << fvis_tmp << " point " << spcollp << std::endl;
 
	}
      }
    }
    std::cout<<"Found: "<< spcollp <<" spacepoints associated to collection hits"<<std::endl;
    std::cout << "SP tot charge " << SPtotcharge << std::endl;
    
    std::cout<<"****Calculate Fvis from visibility map**** "<< std::endl;
    if (spcollp > 0 ) {
      fFvis = fvis_tmp/SPtotcharge; //F_vis weighted for the total Space Point Charge
      std::cout << "fvis " <<  fFvis << std::endl;
      }     
    else{ fFvis = 0;
       std::cout << "fvis " <<  fFvis << std::endl;
	mf::LogWarning("ChargeAndLightAna") << "No Space Point ---> No Fvis calculation. Failing";
	//return;
      }  


       
    ///////////////////////////////////////////////
    // Write out the ChargeLightTree and clean up //
    ///////////////////////////////////////////////

    fChargeLightTree->Fill();
    /*fFlashIDVector              .clear();
    fYCenterVector              .clear();
    fZCenterVector              .clear();
    fYWidthVector               .clear();
    fZWidthVector               .clear();
    fTimeVector                 .clear();
    fTimeWidthVector            .clear();
    fTimeDiffVector             .clear();
    fTotalPEVector              .clear();
    fNHitOpDetVector            .clear();
    fPurityVector               .clear();
    fDistanceVector             .clear();
    fRecoXVector                .clear();
    */
    fTruePxallpart              .clear();
    fTruePyallpart              .clear();
    fTruePzallpart              .clear();
    fTrueEallpart               .clear();
    fTrueAllPDG                 .clear();
    fEnergyDepositionVector     .clear();
    fHitMultiplicity            .clear();
    fHitCharge                  .clear();
    fHitPeakTime                .clear();
    fPointX                     .clear();
    fPointY                     .clear();
    fPointZ                     .clear();
    fGammaScint                 .clear();
    fStepLength		 .clear();
    fStepLCumVector		 .clear();
    fStepEdepCumVector          .clear();
    fOpHitPeakTime  		 .clear();
    fTrkLengthVector   	 .clear();
    //fTrkEnVector                .clear();
    //fTrkdEdxVector              .clear();
    fSelTrkPointX		 .clear();
    fSelTrkPointY		 .clear();
    fSelTrkPointZ		 .clear();
    fHitToXVector		 .clear();
    fHitToYVector		 .clear(); 
    fHitToZVector		 .clear();
    fSpacePointX                .clear(); 
    fSpacePointY                .clear(); 
    fSpacePointZ                .clear(); 
    fHitSPCharge                .clear();    
  }

  //-----------------------------------------------------------------------

 void ChargeAndLightAna::endJob(){
 
 
 }
 
  float ChargeAndLightAna::GetFvisFromHisto(float x, float y, float z, float q){
  
  if (q == 0) return 0;
  
    float binval=0;
    float vis=0;
    Int_t ibin = 0;
    //std::cout << "Funzione F vis" << std::endl;
    //ibin = h3LightMap->FindBin(x, z, y);
    //std::cout << "ibin" << ibin << std::endl;
    
    if (fIsVD) {
      ibin = h3LightMap->FindBin(y, z, x);  
      }
    else {
      ibin = h3LightMap->FindBin(x, z, y);
      }
    
    binval = h3LightMap->GetBinContent(ibin); 
    
    //std::cout<<" ---> f is: "<<binval<< std::endl;
    vis = binval*q;

    //std::cout<<"vis = "<<vis<<std::endl;
    return vis;

 }
 
 

} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(ChargeAndLightAna)
}
