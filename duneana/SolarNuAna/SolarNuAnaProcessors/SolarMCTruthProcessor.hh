/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : SolarMCTruthProcessor.hh
 * @created     : Thursday Jul 17, 2025 07:57:07 CDT
 */

#ifndef SOLARMCTRUTHPROCESSOR_HH

#define SOLARMCTRUTHPROCESSOR_HH

#include "duneana/SolarNuAna/SolarNuAnaProcessors/SolarNuAnaData.hh"

#include "fhiclcpp/ParameterSet.h"

#include "TTree.h"
#include "TH1D.h"

namespace solar{
class SolarMCTruthProcessor {
  public: 
    inline SolarMCTruthProcessor( 
        SolarEventData& _data_event,
        SolarMCTruthData& _data_truth, 
        SolarOpHitData& _data_ophit, 
        SolarOpFlashData& _data_opflash, 
        SolarClusterData& _data_cluster) :
      data_event(_data_event),
      data_truth(_data_truth), 
      data_ophit(_data_ophit),
      data_opflash(_data_opflash),
      data_cluster(_data_cluster) {}
    void Config( const fhicl::ParameterSet& pset ) {}
    void LinkTree( TTree* output_tree ); 

    inline void process();
  
  private: 
    bool fSaveSignalDaughters = false;
    bool fSaveSignalEDep = false;
    bool fSaveSignalOpHits = false;
    bool fSaveOpFlashInfo = false;
    bool fSaveTrackInfo = false; 

  public: 
    SolarEventData& data_event;
    SolarMCTruthData& data_truth; 
    SolarOpHitData& data_ophit;
    SolarOpFlashData& data_opflash;
    SolarClusterData& data_cluster;
};
}


#endif /* end of include guard SOLARMCTRUTHPROCESSOR_HH */

