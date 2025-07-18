/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : SolarNuAnaData
 * @created     : Thursday Jul 17, 2025 09:29:24 CDT
 */

#ifndef SOLARNUANADATA_HH

#define SOLARNUANADATA_HH

#include <vector>
#include <string>

namespace solar {

struct SolarConfigData {
  std::string fHitLabel = {};
  std::string fTrackLabel = {};
  std::string fOpHitLabel = {};
  std::string fOpFlashLabel = {};
  std::string fGEANTLabel = {};
  std::string fSignalLabel = {};
  std::string fGeometry = {};
  int fDetectorSizeY = {};
  int fDetectorSizeZ = {};
  int fClusterAlgoAdjChannel = {};
  int fClusterInd0MatchTime = {};
  int fClusterInd1MatchTime = {};
  int fClusterPreselectionNHits = {};
  int fAdjOpFlashMinNHitCut = {};
  float fClusterMatchTime = {};
  float fAdjClusterRad = {};
  float fMinClusterCharge = {};
  float fClusterMatchCharge = {};
  float fAdjOpFlashX = {};
  float fAdjOpFlashY = {};
  float fAdjOpFlashZ = {};
  float fAdjOpFlashMaxPERatioCut = {};
  float fAdjOpFlashMinPECut = {};
  float fClusterMatchNHit = {};
  float fClusterAlgoTime = {};
  float fOpFlashAlgoMinTime = {};
  float fOpFlashAlgoMaxTime = {};
  float fOpFlashAlgoRad = {};
  float fOpFlashAlgoPE = {};
  float fOpFlashAlgoTriggerPE = {};
  float fOpFlashAlgoHotVertexThld = {};
  bool fClusterPreselectionSignal = {};
  bool fClusterPreselectionPrimary = {};
  bool fClusterPreselectionTrack = {};
  bool fClusterPreselectionFlashMatch = {};
  bool fGenerateAdjOpFlash = {};
  bool fFlashMatchByResidual = {};
  bool fSaveSignalDaughters = {};
  bool fSaveSignalEDep = {};
  bool fSaveSignalOpHits = {};
  bool fSaveOpFlashInfo = {};
  bool fSaveTrackInfo = {};
};

struct SolarEventData {
  int Event = {}; //!< Event number
  int Flag = {};  //!< Flag used to match truth with reco tree entries
};

struct SolarMCTruthData {
  std::string TNuInteraction = {};              //!< True signal interaction process 
  int SignalParticlePDG = {};                   //!< True signal PDG
  float SignalParticleE = {};                   //!< True signal energy [MeV]
  float SignalParticleP = {};                   //!< True signal momentum [MeV]
  float SignalParticleK = {};                   //!< True signal K.E. [MeV]
  float SignalParticleX = {};                   //!< True signal X [cm]
  float SignalParticleY = {};                   //!< True signal Y [cm]
  float SignalParticleZ = {};                   //!< True signal Z [cm]
  float SignalParticleTime = {};                //!< True signal time [tick]
  std::vector<int> TrackNum = {};               //!< Number 
  std::vector<int> TPart = {};                  //!< Number particles per generator
  std::vector<int> SignalPDGList = {};          //!< PDG of Signal particles
  std::vector<int> SignalPDGDepList = {};       //!< PDG for Energy deposited of Signal particles
  std::vector<int> SignalIDList = {};           //!< TrackID of Signal particles
  std::vector<int> SignalMotherList = {};       //!< TrackID of Signal mother
  std::vector<int> SignalIDDepList = {};        //!< ParentID of Signal particles
  std::vector<int> SignalElectronDepList = {};  //!< Number of electrons in the Signal particles
  std::vector<float> SignalEDepList = {};       //!< Energy deposited of Signal particles [MeV]
  std::vector<float> SignalXDepList = {};       //!< X deposited of Signal particles [cm]
  std::vector<float> SignalYDepList = {};       //!< Y deposited of Signal particles [cm]
  std::vector<float> SignalZDepList = {};       //!< Z deposited of Signal particles [cm]
  std::vector<float> SignalEList = {};          //!< Energy of Signal particles [MeV]
  std::vector<float> SignalPList = {};          //!< Energy of Signal momentum [MeV]
  std::vector<float> SignalKList = {};          //!< Kinetik Energy of Signal particles [MeV]
  std::vector<float> SignalTimeList = {};       //!< Time of Signal particles [ticks]
  std::vector<float> SignalEndXList = {};       //!< X of Signal particles [cm]
  std::vector<float> SignalEndYList = {};       //!< Y of Signal particles [cm]
  std::vector<float> SignalEndZList = {};       //!< Z of Signal particles [cm]
  std::vector<float> SignalMaxEDepList = {};    //!< Energy of Signal particles [MeV]
  std::vector<float> SignalMaxEDepXList = {};   //!< X of Signal particles [cm]
  std::vector<float> SignalMaxEDepYList = {};   //!< Y of Signal particles [cm]
  std::vector<float> SignalMaxEDepZList = {};   //!< Z of Signal particles [cm]

  void reset(); 
};

struct SolarOpHitData {
  int OpHitNum = {};                        //!< Number of OpHits
  std::vector<int> SOpHitPlane = {};        //!< OpHit Plane
  std::vector<float> SOpHitPur = {};        //!< OpHit Purity
  std::vector<float> SOpHitPE = {};         //!< OpHit PE
  std::vector<float> SOpHitX = {};          //!< OpHit X
  std::vector<float> SOpHitY = {};          //!< OpHit Y
  std::vector<float> SOpHitZ = {};          //!< OpHit Z
  std::vector<float> SOpHitTime = {};       //!< OpHit Time
  std::vector<float> SOpHitChannel = {};    //!< OpHit Channel
  std::vector<float> SOpHitFlashID = {};    //!< OpHit FlashID

  void reset();
};

struct SolarOpFlashData {
  int OpFlashNum = {};                      //!< Number of OpFlash
  std::vector<int> OpFlashID = {};          //!< OpFlash ID
  std::vector<int> OpFlashNHits = {};       //!< OpFlash NHit
  std::vector<int> OpFlashPlane = {};       //!< OpFlash Plane
  std::vector<float> OpFlashPur = {};       //!< OpFlash Purity
  std::vector<float> OpFlashPE = {};        //!< OpFlash PE
  std::vector<float> OpFlashMaxPE = {};     //!< OpFlash Max PE
  std::vector<float> OpFlashX = {};         //!< OpFlash X
  std::vector<float> OpFlashY = {};         //!< OpFlash Y
  std::vector<float> OpFlashZ = {};         //!< OpFlash Z
  std::vector<float> OpFlashTime = {};      //!< OpFlash Time 
  std::vector<float> OpFlashSTD = {};       //!< OpFlash STD
  std::vector<float> OpFlashDeltaT = {};
  std::vector<float> OpFlashFast = {};

  void reset();
};

struct SolarClusterData {
  bool MPrimary = {};                           //!< Cluster hasn't any adjcl with AdjClCharge > MCharge (bool)
  int TrackNum = {};                            //!< Number of PMTracks
  int MGen = {};                                //!< Main cluster generator idx
  int MTPC = {};                                //!< Main cluster TPC
  int MInd0TPC = {};                            //!< Main cluster ind0 TPC    
  int MInd1TPC = {};                            //!< Main cluster ind1 TPC
  int MInd0NHits = {};                          //!< Main cluster ind0 Hits
  int MInd1NHits = {};                          //!< Main cluster ind1 Hits
  int MMainID = {};                             //!< Main cluster main track ID
  int MMainPDG = {};                            //!< Main cluster main pdg
  int MMainParentPDG = {};                      //!< Main cluster main pdg
  float MInd0dTime = {};                        //!< Main cluster ind0 dT [Ticks]
  float MInd1dTime = {};                        //!< Main cluster ind1 dT [Ticks]
  float MInd0RecoY = {};                        //!< Main cluster ind0 reco Y [cm]
  float MInd1RecoY = {};                        //!< Main cluster ind1 reco Y [cm]
  float MRecX = {};                             //!< Main cluster reco X [cm] (from matched Flash)
  float MRecY = {};                             //!< Main cluster reco Y [cm]
  float MRecZ = {};                             //!< Main cluster reco Z [cm]
  float MPur = {};                              //!< Main cluster reco signal purity
  float MInd0Pur = {};                          //!< Main cluster ind0 reco signal purity
  float MInd1Pur = {};                          //!< Main cluster ind1 reco signal purity
  float MTime = {};                             //!< Main cluster time [ticks]
  float MCharge = {};                           //!< Main cluster charge [ADC*ticks]
  float MMaxCharge = {};                        //!< Main cluster's max TPCHit-charge [ADC*ticks]
  float MInd0Charge = {};                       //!< Main cluster ind0 MaxHit
  float MInd1Charge = {};                       //!< Main cluster ind1 MaxHit
  float MInd0MaxCharge = {};                    //!< Main cluster ind0 MaxHit
  float MInd1MaxCharge = {};                    //!< Main cluster ind1 MaxHit
  float MGenPur = {};                           //!< Main cluster reco generator purity
  float MMainE = {};                            //!< Main cluster main energy [MeV]
  float MMainP = {};                            //!< Main cluster main momentum [MeV]
  float MMainK = {};                            //!< Main cluster main kinetic energy [MeV]
  float MMainTime = {};                         //!< Main cluster main Time [ticks]
  float MMainParentE = {};                      //!< Main cluster main parent energy [MeV]
  float MMainParentP = {};                      //!< Main cluster main parent momentum [MeV]
  float MMainParentK = {};                      //!< Main cluster main parent kinetic energy [MeV]
  float MMainParentTime = {};                   //!< Main cluster main parent Time [ticks]
  std::vector<int> HitNum = {};                 //!< Number of hits in each TPC plane
  std::vector<int> ClusterNum = {};             //!< Number of clusters in each TPC plane
  std::vector<float> MSignalFrac = {};          //!< Main cluster particle contribution (electron, gamma, neutron) 
  std::vector<float> MGenFrac = {};             //!< Main cluster reco purity complete
  std::vector<double> MMainVertex = {};         //!< Main cluster main particle vertex [cm]
  std::vector<double> MEndVertex = {};          //!< Main cluster end particle vertex [cm]
  std::vector<double> MMainParentVertex = {};   //!< Main cluster parent particle vertex [cm]

  //---------------------------------------------------------------- Track info
  int MTrackNPoints = {};                       //!< Track #points     
  std::vector<double> MTrackStart = {};         //!< Track start point
  std::vector<double> MTrackEnd = {};           //!< Track end point
  float MTrackChi2 = {};                        //!< Track chi2
                                                
  //---------------------------------------------------- Adjacent Clusters info
  std::vector<int> MAdjClGen = {};              //!< Adj. clusters' generator idx
  std::vector<int> MAdjClMainID = {};           //!< Adj. clusters' main track ID
  std::vector<int> MAdjClMainPDG = {};          //!< Adj. clusters' main PDG
  std::vector<float> MAdjClMainE = {};          //!< Adj. clusters' main energy [MeV]
  std::vector<float> MAdjClMainP = {};          //!< Adj. clusters' main momentum [MeV]
  std::vector<float> MAdjClMainK = {};          //!< Adj. clusters' main K.E. [MeV]
  std::vector<float> MAdjClMainX = {};          //!< Adj. clusters' main X [cm]
  std::vector<float> MAdjClMainY = {};          //!< Adj. clusters' main Y [cm]
  std::vector<float> MAdjClMainZ = {};          //!< Adj. clusters' main Z [cm]
  std::vector<float> MAdjClTime = {};           //!< Adj. clusters' time [ticks]
  std::vector<float> MAdjClCharge = {};         //!< Adj. clusters' charge [ADC*ticks]
  std::vector<float> MAdjClInd0Charge = {};     //!< Adj. clusters' charge [ADC*ticks]
  std::vector<float> MAdjClInd1Charge = {};     //!< Adj. clusters' charge [ADC*ticks]
  std::vector<float> MAdjClMaxCharge = {};      //!< Adj. clusters' charge [ADC*ticks]
  std::vector<float> MAdjClInd0MaxCharge = {};  //!< Adj. clusters' charge [ADC*ticks]
  std::vector<float> MAdjClInd1MaxCharge = {};  //!< Adj. clusters' charge [ADC*ticks]
  std::vector<float> MAdjClNHits = {};          //!< Adj. clusters' #hits
  std::vector<float> MAdjClInd0NHits = {};      //!< Adj. clusters' #hits
  std::vector<float> MAdjClInd1NHits = {};      //!< Adj. clusters' #hits
  std::vector<float> MAdjClRecoY = {};          //!< Adj. clusters' reco Y [cm]
  std::vector<float> MAdjClRecoZ = {};          //!< Adj. clusters' reco Z [cm]
  std::vector<float> MAdjClR = {};              //!< Adj. clusters' distance [cm]
  std::vector<float> MAdjClPur = {};            //!< Adj. clusters' purity
  std::vector<float> MAdjClGenPur = {};
  std::vector<float> MAdjClEndX = {};           //!< Adj. clusters' end X [cm]
  std::vector<float> MAdjClEndY = {};           //!< Adj. clusters' end Y [cm]
  std::vector<float> MAdjClEndZ = {};           //!< Adj. clusters' end Z [cm]


  void reset();
};

struct SolarAdjClFlashData {
  std::vector<int> MAdjFlashPlane = {};         //!< Adj. flash' Plane
  std::vector<int> MAdjFlashNHits = {};         //!< Adj. flash' #hits
  std::vector<float> MAdjFlashTime = {};        //!< adj. flash' time [ticks]
  std::vector<float> MAdjFlashPE = {};          //!< Adj. flash' tot #PE [ADC*ticks]
  std::vector<float> MAdjFlashMaxPE = {};       //!< Adj. flash' max #PE [ADC*ticks]
  std::vector<float> MAdjFlashRecoX = {};       //!< Adj. flash' reco X [cm]
  std::vector<float> MAdjFlashRecoY = {};       //!< Adj. flash' reco Y [cm]
  std::vector<float> MAdjFlashRecoZ = {};       //!< Adj. flash' reco Z [cm]
  std::vector<float> MAdjFlashR = {};           //!< Adj. flash' reco distance [cm]
  std::vector<float> MAdjFlashPur = {};         //!< Adj. flash' purity
  std::vector<float> MAdjFlashSTD = {};         //!< Adj. flash' STD
  std::vector<float> MAdjFlashFast = {};        //!< Adj. flash' Fast Component
  std::vector<float> MAdjFlashResidual = {};    //!< Adj. flash' residual wrt. cluster
                                                
  void reset();
};

struct SolarMatchFlashData {
  int MFlashNHits = {};                         //!< Matched flash' #hits
  int MFlashPlane = {};                         //!< Matched flash' Plane
  float MFlashR = {};                           //!< Matched flash' reco distance [cm]
  float MFlashPE = {};                          //!< Matched flash' tot #PE [ADC*ticks]
  float MFlashMaxPE = {};                       //!< Matched flash' max #PE [ADC*ticks]
  float MFlashPur = {};                         //!< Matched flash' purity
  float MFlashFast = {};                        //!< Matched flash' Fast Component
  float MFlashTime = {};                        //!< Matched flash' time [ticks]
  float MFlashSTD = {};                         //!< Matched flash' STD
  float MFlashRecoX = {};                       //!< Matched flash' reco X [cm]    
  float MFlashRecoY = {};                       //!< Matched flash' reco Y [cm]
  float MFlashRecoZ = {};                       //!< Matched flash' reco Z [cm]
  float MFlashResidual = {};                    //!< Matched flash' residual wrt. cluster
  bool MFlashCorrect = {};                      //!< Matched flash' correctnes (bool)

  void reset();
};

}


#endif /* end of include guard SOLARNUANADATA_HH */

