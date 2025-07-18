/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : SolarNuAnaData.cc
 * @created     : Thursday Jul 17, 2025 10:38:26 CDT
 */

#include "duneana/SolarNuAna/SolarNuAnaProcessors/SolarNuAnaData.hh"

namespace solar {

  void SolarMCTruthData::reset() {
    TNuInteraction = {};
    SignalParticlePDG = {};
    SignalParticleE = {};
    SignalParticleP = {};
    SignalParticleK = {};
    SignalParticleX = {};
    SignalParticleY = {};
    SignalParticleZ = {};
    SignalParticleTime = {};
    TPart.clear();
    SignalPDGList.clear();
    SignalPDGDepList.clear();
    SignalIDList.clear();
    SignalMotherList.clear();
    SignalIDDepList.clear();
    SignalElectronDepList.clear();
    SignalEDepList.clear();
    SignalXDepList.clear();
    SignalYDepList.clear();
    SignalZDepList.clear();
    SignalEList.clear();
    SignalPList.clear();
    SignalKList.clear();
    SignalTimeList.clear();
    SignalEndXList.clear();
    SignalEndYList.clear();
    SignalEndZList.clear();
    SignalMaxEDepList.clear();
    SignalMaxEDepXList.clear();
    SignalMaxEDepYList.clear();
    SignalMaxEDepZList.clear();

    return;
  }

  void SolarOpHitData::reset() {
    OpHitNum = {};
    SOpHitPlane.clear();
    SOpHitPur.clear();
    SOpHitPE.clear();
    SOpHitX.clear();
    SOpHitY.clear();
    SOpHitZ.clear();
    SOpHitTime.clear();
    SOpHitChannel.clear();
    SOpHitFlashID.clear();
    return;
  }

  void SolarOpFlashData::reset() {
    OpFlashNum = {};
    OpFlashID.clear();
    OpFlashNHits.clear();
    OpFlashPlane.clear();
    OpFlashPur.clear();
    OpFlashPE.clear();
    OpFlashMaxPE.clear();
    OpFlashX.clear();
    OpFlashY.clear();
    OpFlashZ.clear();
    OpFlashTime.clear();
    OpFlashSTD.clear();
    OpFlashDeltaT.clear();
    OpFlashFast.clear();
    return;
  }

  void SolarClusterData::reset() {
    MPrimary = false;
    TrackNum = 0;
    MGen = 0;
    MTPC = 0;
    MInd0TPC = 0;
    MInd1TPC = 0;
    MInd0NHits = 0;
    MInd1NHits = 0;
    MMainID = 0;
    MMainPDG = 0;
    MMainParentPDG = 0;
    MInd0dTime = 0;
    MInd1dTime = 0;
    MInd0RecoY = 0;
    MInd1RecoY = 0;
    MRecX = 0;
    MRecY = 0;
    MRecZ = 0;
    MPur = 0;
    MInd0Pur = 0;
    MInd1Pur = 0;
    MTime = 0;
    MCharge = 0;
    MMaxCharge = 0;
    MInd0Charge = 0;
    MInd1Charge = 0;
    MInd0MaxCharge = 0;
    MInd1MaxCharge = 0;
    MGenPur = 0;
    MMainE = 0;
    MMainP = 0;
    MMainK = 0;
    MMainTime = 0;
    MMainParentE = 0;
    MMainParentP = 0;
    MMainParentK = 0;
    MMainParentTime = 0;
    HitNum.clear();
    ClusterNum.clear();
    MSignalFrac.clear();
    MGenFrac.clear();
    MMainVertex.clear();
    MEndVertex.clear();
    MMainParentVertex.clear();

    //---------------------------------------------------------------- Track info
    MTrackNPoints = 0;
    MTrackStart.clear();
    MTrackEnd.clear();
    MTrackChi2 = 0.0;

    //---------------------------------------------------- Adjacent Clusters info
    MAdjClGen.clear();
    MAdjClMainID.clear();
    MAdjClMainPDG.clear();
    MAdjClMainE.clear();
    MAdjClMainP.clear();
    MAdjClMainK.clear();
    MAdjClMainX.clear();
    MAdjClMainY.clear();
    MAdjClMainZ.clear();
    MAdjClTime.clear();
    MAdjClCharge.clear();
    MAdjClInd0Charge.clear();
    MAdjClInd1Charge.clear();
    MAdjClMaxCharge.clear();
    MAdjClInd0MaxCharge.clear();
    MAdjClInd1MaxCharge.clear();
    MAdjClNHits.clear();
    MAdjClInd0NHits.clear();
    MAdjClInd1NHits.clear();
    MAdjClRecoY.clear();
    MAdjClRecoZ.clear();
    MAdjClR.clear();
    MAdjClPur.clear();
    MAdjClGenPur.clear();
    MAdjClEndX.clear();
    MAdjClEndY.clear();
    MAdjClEndZ.clear();

    return;
  }

  void SolarAdjClFlashData::reset() {
    MAdjFlashPlane.clear();
    MAdjFlashNHits.clear();
    MAdjFlashTime.clear();
    MAdjFlashPE.clear();
    MAdjFlashMaxPE.clear();
    MAdjFlashRecoX.clear();
    MAdjFlashRecoY.clear();
    MAdjFlashRecoZ.clear();
    MAdjFlashR.clear();
    MAdjFlashPur.clear();
    MAdjFlashSTD.clear();
    MAdjFlashFast.clear();
    MAdjFlashResidual.clear();

    return;
  }

  void SolarMatchFlashData::reset() {
    MFlashNHits = 0;
    MFlashPlane = 0;
    MFlashR = 0;
    MFlashPE = 0;
    MFlashMaxPE = 0;
    MFlashPur = 0;
    MFlashFast = 0;
    MFlashTime = 0;
    MFlashSTD = 0;
    MFlashRecoX = 0;
    MFlashRecoY = 0;
    MFlashRecoZ = 0;
    MFlashResidual = 0;
    MFlashCorrect = false;

    return;
  }

}

