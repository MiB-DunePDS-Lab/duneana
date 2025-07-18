#include "duneana/SolarNuAna/SolarNuAnaProcessors/SolarMCTruthProcessor.hh"

namespace solar {
  void SolarMCTruthProcessor::process() {

  }
  
  void SolarMCTruthProcessor::LinkTree( TTree* output_tree ) 
  {
    output_tree->Branch("Event", &data_event.Event, "Event/I");                                        // Event number
    output_tree->Branch("Flag", &data_event.Flag, "Flag/I");                                           // Flag used to match truth with reco tree entries
    output_tree->Branch("TruthPart", &data_truth.TPart);                                               // Number particles per generator
    output_tree->Branch("Interaction", &data_truth.TNuInteraction);                                    // True signal interaction process
    output_tree->Branch("SignalParticleE", &data_truth.SignalParticleE, "SignalParticleE/F");          // True signal energy [MeV]
    output_tree->Branch("SignalParticleP", &data_truth.SignalParticleP, "SignalParticleP/F");          // True signal momentum [MeV]
    output_tree->Branch("SignalParticleK", &data_truth.SignalParticleK, "SignalParticleK/F");          // True signal K.E. [MeV]
    output_tree->Branch("SignalParticleX", &data_truth.SignalParticleX, "SignalParticleX/F");          // True signal X [cm]
    output_tree->Branch("SignalParticleY", &data_truth.SignalParticleY, "SignalParticleY/F");          // True signal Y [cm]
    output_tree->Branch("SignalParticleZ", &data_truth.SignalParticleZ, "SignalParticleZ/F");          // True signal Z [cm]
    output_tree->Branch("SignalParticlePDG", &data_truth.SignalParticlePDG, "SignalParticlePDG/I");    // True signal PDG
    output_tree->Branch("SignalParticleTime", &data_truth.SignalParticleTime, "SignalParticleTime/F"); // True signal time [tick]
    output_tree->Branch("OpHitNum", &data_ophit.OpHitNum, "OpHitNum/I");                               // Number of OpHits
    output_tree->Branch("OpFlashNum", &data_opflash.OpFlashNum, "OpFlashNum/I");                       // Number of OpFlashes
    output_tree->Branch("HitNum", &data_cluster.HitNum);                                                 // Number of hits in each TPC plane
    output_tree->Branch("ClusterNum", &data_cluster.ClusterNum);                                         // Number of clusters in each TPC plane
    output_tree->Branch("TrackNum", &data_cluster.TrackNum, "TrackNum/I");                               // Number of PMTracks
    if (fSaveSignalDaughters)
    { // Save Signal Daughters. (Only makes sense for marley)
      output_tree->Branch("TSignalPDG", &data_truth.SignalPDGList);         // PDG of Signal marticles
      output_tree->Branch("TSignalE", &data_truth.SignalEList);             // Energy of Signal particles [MeV]
      output_tree->Branch("TSignalP", &data_truth.SignalPList);             // Energy of Signal momentum [MeV]
      output_tree->Branch("TSignalK", &data_truth.SignalKList);             // Kinetik Energy of Signal particles [MeV]
      output_tree->Branch("TSignalT", &data_truth.SignalTimeList);          // Time of Signal particles [ticks]
      output_tree->Branch("TSignalEndX", &data_truth.SignalEndXList);       // X of Signal particles [cm]
      output_tree->Branch("TSignalEndY", &data_truth.SignalEndYList);       // Y of Signal particles [cm]
      output_tree->Branch("TSignalEndZ", &data_truth.SignalEndZList);       // Z of Signal particles [cm]
      output_tree->Branch("TSignalMaxEDep", &data_truth.SignalMaxEDepList); // Energy of Signal particles [MeV]
      output_tree->Branch("TSignalX", &data_truth.SignalMaxEDepXList);      // X of Signal particles [cm]
      output_tree->Branch("TSignalY", &data_truth.SignalMaxEDepYList);      // Y of Signal particles [cm]
      output_tree->Branch("TSignalZ", &data_truth.SignalMaxEDepZList);      // Z of Signal particles [cm]
      output_tree->Branch("TSignalID", &data_truth.SignalIDList);           // TrackID of Signal particles
      output_tree->Branch("TSignalMother", &data_truth.SignalMotherList);   // TrackID of Signal mother
    }
    if (fSaveSignalEDep)
    {
      output_tree->Branch("TSignalPDGDepList", &data_truth.SignalPDGDepList);           // PDG for Energy deposited of Signal particles
      output_tree->Branch("TSignalEDepList", &data_truth.SignalEDepList);               // Energy deposited of Signal particles [MeV]
      output_tree->Branch("TSignalXDepList", &data_truth.SignalXDepList);               // X deposited of Signal particles [cm]
      output_tree->Branch("TSignalYDepList", &data_truth.SignalYDepList);               // Y deposited of Signal particles [cm]
      output_tree->Branch("TSignalZDepList", &data_truth.SignalZDepList);               // Z deposited of Signal particles [cm]
      output_tree->Branch("TSignalIDDepList", &data_truth.SignalIDDepList);             // ParentID of Signal particles
      output_tree->Branch("TSignalElectronDepList", &data_truth.SignalElectronDepList); // Number of electrons in the Signal particles
    }
    if (fSaveSignalOpHits)
    { // Save OpHits. (Can be very heavy for background productions)
      output_tree->Branch("OpHitPur", &data_ophit.SOpHitPur);         // OpHit Purity
      output_tree->Branch("OpHitPlane", &data_ophit.SOpHitPlane);     // OpHit Plane
      output_tree->Branch("OpHitPE", &data_ophit.SOpHitPE);           // OpHit PE
      output_tree->Branch("OpHitX", &data_ophit.SOpHitX);             // OpHit X
      output_tree->Branch("OpHitY", &data_ophit.SOpHitY);             // OpHit Y
      output_tree->Branch("OpHitZ", &data_ophit.SOpHitZ);             // OpHit Z
      output_tree->Branch("OpHitTime", &data_ophit.SOpHitTime);       // OpHit Time
      output_tree->Branch("OpHitChannel", &data_ophit.SOpHitChannel); // OpHit Channel
      output_tree->Branch("OpHitFlashID", &data_ophit.SOpHitFlashID); // OpHit Area
    }
    if (fSaveOpFlashInfo)
    {
      output_tree->Branch("OpFlashPur", &data_opflash.OpFlashPur);     // OpFlash Purity
      output_tree->Branch("OpFlashID", &data_opflash.OpFlashID);       // OpFlash ID
      output_tree->Branch("OpFlashPE", &data_opflash.OpFlashPE);       // OpFlash PE
      output_tree->Branch("OpFlashX", &data_opflash.OpFlashX);         // OpFlash X
      output_tree->Branch("OpFlashY", &data_opflash.OpFlashY);         // OpFlash Y
      output_tree->Branch("OpFlashZ", &data_opflash.OpFlashZ);         // OpFlash Z
      output_tree->Branch("OpFlashTime", &data_opflash.OpFlashTime);   // OpFlash Time
      output_tree->Branch("OpFlashSTD", &data_opflash.OpFlashSTD);     // OpFlash STD
      output_tree->Branch("OpFlashNHits", &data_opflash.OpFlashNHits); // OpFlash NHit
      output_tree->Branch("OpFlashPlane", &data_opflash.OpFlashPlane); // OpFlash Plane
      output_tree->Branch("OpFlashMaxPE", &data_opflash.OpFlashMaxPE); // OpFlash Max PE
    }

    return;
  }
}
