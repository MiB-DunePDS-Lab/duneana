/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : PhotoVisDump_module.cc
 * @created     : Thursday Feb 10, 2022 05:19:33 CST
 */

#ifndef PHOTOVISDUMP_MODULE_CC

#define PHOTOVISDUMP_MODULE_CC

// ROOT includes
#include "TH1D.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TTree.h"
#include "TVectorF.h"
#include "TRandom3.h"
#include "TParameter.h"

// C++ includes
#include <cstdio>
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"
#include "larsimdnn/PhotonPropagation/TFLoaderTools/TFLoader.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larcorealg/CoreUtils/counter.h"

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
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/FindManyP.h"

namespace opdet {

  class PhotoVisDump : public art::EDAnalyzer 
  {
    public:
      enum EVisModel {kSemiAnalytical = 0, kCompGraph = 1};

      PhotoVisDump(const fhicl::ParameterSet&);
      virtual ~PhotoVisDump() {};

      void beginJob();
      void endJob() {}
      void dumpMap();
      void analyze (const art::Event&);

    private: 
      std::unique_ptr<phot::SemiAnalyticalModel> fVisibilityModel;
      std::unique_ptr<phot::TFLoader> fTFGenerator; 
      size_t nOpDets; 
      std::vector<geo::Point_t> fOpDetCenter;

      EVisModel kVisModel; 

      bool fDoReflectedLight;
      bool fIncludeAnodeReflections;
      bool fIncludeBuffer; 
      bool fUseXeAbsorption;

      double fVoxelSizeX;
      double fVoxelSizeY;
      double fVoxelSizeZ;

      fhicl::ParameterSet fVUVHitsParams;
      fhicl::ParameterSet fVISHitsParams;
      fhicl::ParameterSet fTFLoaderPars;

      bool fIsDone = false;
  };

} // close opdet namespace

namespace opdet {
  DEFINE_ART_MODULE(PhotoVisDump)
}


#endif /* end of include guard VISMAPDUMP_MODULE_CC */

namespace opdet {
  PhotoVisDump::PhotoVisDump(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset)
  {
    //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    // Read inputs from fcl file
    fVoxelSizeX = pset.get<double>("voxel_dx", 10.0);
    fVoxelSizeY = pset.get<double>("voxel_dy", 10.0);
    fVoxelSizeZ = pset.get<double>("voxel_dz", 10.0);

    fDoReflectedLight = pset.get<bool>("do_refl", false); 
    fIncludeAnodeReflections = pset.get<bool>("do_include_anode_refl", false); 
    fIncludeBuffer = pset.get<bool>("do_include_buffer", false); 
    fUseXeAbsorption = pset.get<bool>("do_include_xe_absorption", false);

    fVUVHitsParams = pset.get<fhicl::ParameterSet>("vuvhitspars"); 
    fVISHitsParams = pset.get<fhicl::ParameterSet>("vishitspars"); 

    fTFLoaderPars = pset.get<fhicl::ParameterSet>("tfloaderpars");

    TString vis_model_str = pset.get<std::string>("vis_model"); 

    if (vis_model_str.Contains("compgraph")) {
      kVisModel = kCompGraph;
    }
    else if (vis_model_str.Contains("semianalytical")) {
      kVisModel = kSemiAnalytical;
    }
  }

  void PhotoVisDump::beginJob() {
    // create the photo-detector visibility model
    if (kVisModel == kSemiAnalytical) {
      printf("Creating Semi-analytical visibility model\n");
      fVisibilityModel = std::make_unique<phot::SemiAnalyticalModel>(
          fVUVHitsParams, fVISHitsParams, 
          fDoReflectedLight, fIncludeAnodeReflections, fUseXeAbsorption
          ); 
    }
    else if (kVisModel == kCompGraph) {
      fTFGenerator = art::make_tool<phot::TFLoader>(fTFLoaderPars);
      fTFGenerator->Initialization();
    }

    return;
  }

  void PhotoVisDump::analyze(const art::Event&) {

    //if (fIsDone) return

    dumpMap(); 
  }

  void PhotoVisDump::dumpMap() {
    // retrieve geometry
    //geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());
    const auto geom = art::ServiceHandle<geo::Geometry>();

    const auto photonVisService = art::ServiceHandle<phot::PhotonVisibilityService>(); 

    // open file
    art::ServiceHandle< art::TFileService > tfs;

    // setup variables
    struct  MeshPoint_t{
      double x, y, z;
      MeshPoint_t() {};
      MeshPoint_t(const double& x_, const double& y_, const double& z_) {
        x = x_; y = y_; z = z_;
      }
      MeshPoint_t(const TVector3& vv) {
        x = vv.x();
        y = vv.y();
        z = vv.z();
      }
      MeshPoint_t(const geo::Point_t& pt) {
        x = pt.x(); 
        y = pt.y(); 
        z = pt.z(); 
      }
    };
    MeshPoint_t point_;
    double visDirect = 0;
    double visReflct = 0;
    MeshPoint_t pointBuff_;
    double visDirectBuff = 0;
    double visReflctBuff = 0;

    float tpcH = 0;  
    float tpcW = 0;
    float tpcL = 0;
    float tpcPos[3]; 

    float opDetH = 0; 
    float opDetW = 0;
    float opDetL = 0;
    float opDetPos[3]; 
    std::vector<size_t> opChannel; 

    const double voxelDim[3] = {fVoxelSizeX, fVoxelSizeY, fVoxelSizeZ}; 

    // open tree
    TTree* tMap = tfs->make<TTree>("photoVisMap", "photoVisMap");
    tMap->SetNameTitle("photoVisMap", "photoVisMap"); 
    tMap->Branch("meshpoint", &point_, "x/D:y:z");
    tMap->Branch("visDirect", &visDirect, "visDirect/D");
    tMap->Branch("visReflct", &visReflct, "visReflct/D");

    TTree* tMapBuffer = tfs->make<TTree>("photoVisMapBuffer", "photoVisMapBuffer"); 
    tMapBuffer->SetNameTitle("photoVisMapBuffer", "photoVisMapBuffer"); 
    tMapBuffer->Branch("meshpoint", &pointBuff_, "x/D:y:z");
    tMapBuffer->Branch("visDirect", &visDirectBuff, "visDirect/D");
    tMapBuffer->Branch("visReflct", &visReflctBuff, "visReflct/D");

    TTree* tTPC = tfs->make<TTree>("tpcMap", "tpcMap");
    tTPC->Branch("tpcH", &tpcH); 
    tTPC->Branch("tpcL", &tpcL); 
    tTPC->Branch("tpcW", &tpcW);
    tTPC->Branch("tpcPos", &tpcPos, "tpcPos[3]/F");
    
    TTree* tOpDet = tfs->make<TTree>("opDetMap", "opDetMap");
    tOpDet->Branch("opDetH", &opDetH); 
    tOpDet->Branch("opDetL", &opDetL); 
    tOpDet->Branch("opDetW", &opDetW);
    tOpDet->Branch("opDetPos", &opDetPos, "opDetPos[3]/F");
    tOpDet->Branch("opDetCh", &opChannel); 

    geo::Point_t point_center;
    TVector3 center;
    TVector3 tpcpoint;
    TVector3 vpoint;
    std::vector<double> opdetvis_dir;
    std::vector<double> opdetvis_rfl;
    
    // store info from Geometry service
    nOpDets = geom->NOpDets();
    for (size_t i : util::counter(nOpDets)) {
      opChannel.clear(); 
      geo::OpDetGeo const& opDet = geom->OpDetGeoFromOpDet(i);
      auto center = opDet.GetCenter();
      center.GetCoordinates( opDetPos );
      opDetH = opDet.Height();
      opDetW = opDet.Width();
      opDetL = opDet.Length();

      size_t n_ch = geom->NOpHardwareChannels(i); 
      opChannel.resize(n_ch, 0); 
      for (size_t ich = 0; ich < n_ch; ich++) {
        opChannel.at(ich) = geom->OpChannel(i, ich);  
      }

      tOpDet->Fill();
    }

    

    // get cryostat dimensions 
    float cryostatMax[3] = {0}; 
    float cryostatMin[3] = {0};
    for (geo::CryostatGeo const& cryo : geom->Iterate<geo::CryostatGeo>()) {
      printf("cryostat: %u [%g, %g, %g] - [%g, %g, %g]\n", cryo.ID().getIndex(),
          cryo.MinX(), cryo.MinY(), cryo.MinZ(), cryo.MaxX(), cryo.MaxY(), cryo.MaxZ()); 
      cryostatMin[0] = cryo.MinX();   cryostatMax[0] = cryo.MaxX(); 
      cryostatMin[1] = cryo.MinY();   cryostatMax[1] = cryo.MaxY(); 
      cryostatMin[2] = cryo.MinZ();   cryostatMax[2] = cryo.MaxZ(); 
    }
    // define mesh points with an histogram helper to align tpcs and buffer
    // sampling coordinates
    TH1D* hgrid[3]; 
    for (int i=0; i<3; i++) {
      double   xc = 0.5*(cryostatMax[i] + cryostatMin[i]); 
      int    nbin = ceil( (cryostatMax[i]-cryostatMin[i]) / voxelDim[i] );
      double xmin = xc - nbin*0.5*voxelDim[i]; 
      double xmax = xc + nbin*0.5*voxelDim[i]; 
      hgrid[i] = tfs->make<TH1D>(Form("hgrid%i", i), Form("mesh points axis %i", i), 
          nbin, xmin, xmax); 
    }


    float tpcMin[3] = {0};
    float tpcMax[3] = {0};
    // loop over all TPCs to get the active volume dimensions
    for (geo::TPCGeo const& tpc : geom->Iterate<geo::TPCGeo>()) {
      point_center = tpc.GetCenter();
      //point_center.GetCoordinates( tpcPos ); 
      center = TVector3( point_center.x(), point_center.y(), point_center.z() ); 
      double hlfW = tpc.ActiveHalfWidth (); tpcW = 2*hlfW;
      double hlfH = tpc.ActiveHalfHeight(); tpcH = 2*hlfH;
      double hlfL = tpc.ActiveHalfLength(); tpcL = 2*hlfL;

      const double hlfDim[3] = {hlfW, hlfH, hlfL}; 

      double x_min[3]; double x_max[3]; 
      point_center.GetCoordinates( x_min );
      point_center.GetCoordinates( x_max ); 
      for (int idim =0; idim<3; idim++) {
        x_min[idim] -= hlfDim[idim]; 
        x_max[idim] += hlfDim[idim]; 

        if (tpcMin[idim] > x_min[idim]) tpcMin[idim] = x_min[idim]; 
        if (tpcMax[idim] < x_max[idim]) tpcMax[idim] = x_max[idim];
      }

      //for (auto const& plane : tpc.IteratePlanes() ) {
      //printf("\tplane ID: %u - %i\n", plane.ID().getIndex(), plane.View()); 
      //}

      tTPC->Fill();

      printf("TPC dimensions: hlfW %.2f - hlfH %.2f hlfL %.2f, ", 
          hlfW, hlfH, hlfL); 
      printf("TPC center:  %.2f -  %.2f - %.2f\n", 
          center.x(), center.y(), center.z()); 
    }

    bool tpc_range_x = false, tpc_range_y = false, tpc_range_z = false; 

    printf("Cryostat min-max: [%g, %g, %g] - [%g, %g, %g]\n", 
        cryostatMin[0], cryostatMin[1], cryostatMin[2], 
        cryostatMax[0], cryostatMax[1], cryostatMax[2]);
    printf("TPC min-max: [%g, %g, %g] - [%g, %g, %g]\n", 
        tpcMin[0], tpcMin[1], tpcMin[2], 
        tpcMax[0], tpcMax[1], tpcMax[2]);
    tfs->make<TParameter<Double_t>>("tpc_min_0", tpcMin[0]); 
    tfs->make<TParameter<Double_t>>("tpc_max_0", tpcMax[0]); 
    tfs->make<TParameter<Double_t>>("tpc_min_1", tpcMin[1]); 
    tfs->make<TParameter<Double_t>>("tpc_max_1", tpcMax[1]);  
    tfs->make<TParameter<Double_t>>("tpc_min_2", tpcMin[2]); 
    tfs->make<TParameter<Double_t>>("tpc_max_2", tpcMax[2]); 
    tfs->make<TParameter<Double_t>>("cryo_min_0", cryostatMin[0]); 
    tfs->make<TParameter<Double_t>>("cryo_max_0", cryostatMax[0]); 
    tfs->make<TParameter<Double_t>>("cryo_min_1", cryostatMin[1]); 
    tfs->make<TParameter<Double_t>>("cryo_max_1", cryostatMax[1]);  
    tfs->make<TParameter<Double_t>>("cryo_min_2", cryostatMin[2]); 
    tfs->make<TParameter<Double_t>>("cryo_max_2", cryostatMax[2]); 


    // here we can loop over the points of the pre-defined grid
    double x_= 0, y_= 0, z_= 0;
    for (int ix=1; ix<=hgrid[0]->GetNbinsX(); ix++) {
        x_ = hgrid[0]->GetBinCenter(ix); 
        if (x_> tpcMin[0] && x_< tpcMax[0]) tpc_range_x = true;
        else tpc_range_x = false; 

        for (int iy=1; iy<=hgrid[1]->GetNbinsX(); iy++) {
          y_= hgrid[1]->GetBinCenter(iy);            
          if (y_> tpcMin[1] && y_< tpcMax[1]) tpc_range_y = true; 
          else tpc_range_y = false; 

          for (int iz=1; iz<=hgrid[2]->GetNbinsX(); iz++) {
            z_= hgrid[2]->GetBinCenter(iz);            
            if (z_> tpcMin[2] && z_< tpcMax[2]) tpc_range_z = true; 
            else tpc_range_z = false; 

            if ( (tpc_range_x && tpc_range_y && tpc_range_z) == false) {
              if (fIncludeBuffer) {
                visDirectBuff = 0; 
                visReflctBuff = 0;

                geo::Point_t point( x_, y_, z_); 
                auto mapped_vis = photonVisService->GetAllVisibilities(point); 
                for (const auto &vis : mapped_vis) visDirectBuff += vis; 

                pointBuff_ = MeshPoint_t( point ); 
                tMapBuffer->Fill(); 
              }
            } 
            else {
              vpoint = TVector3(x_, y_, z_); 
              //printf("vpoint: %.2f - %.2f - %.2f\n", vpoint.x(), vpoint.y(), vpoint.z()); 
              visDirect = 0.;  visReflct = 0.;

              const double n_samplings = 5.0; 

              for (int i = 0; i < n_samplings; i++) {
              auto vpoint_ = vpoint; 
              vpoint_ += TVector3(gRandom->Uniform(-0.5*voxelDim[0], 0.5*voxelDim[0]), 
              gRandom->Uniform(-0.5*voxelDim[1], 0.5*voxelDim[1]), 
              gRandom->Uniform(-0.5*voxelDim[2], 0.5*voxelDim[2]) );  

              geo::Point_t xspot(vpoint_);

              if (kVisModel == kSemiAnalytical ) {
                fVisibilityModel->detectedDirectVisibilities   (opdetvis_dir, xspot);
                //fVisibilityModel->detectedReflectedVisibilities(opdetvis_rfl, xspot, fIncludeAnodeReflections);
              }
              else if (kVisModel == kCompGraph ) {
                std::vector<Double_t> pos_tmp {xspot.x(), xspot.y(), xspot.z()}; 
                fTFGenerator->Predict( pos_tmp ); 
                opdetvis_dir = fTFGenerator->GetPrediction(); 
              }
              for (const auto &vis : opdetvis_dir)  visDirect += vis;
              //for (const auto &vis : opdetvis_rfl)  visReflct += vis;
              }

              visDirect = visDirect / n_samplings; 
              //visReflct = visReflct / n_samplings; 

              // Fill tree
              point_ = MeshPoint_t(vpoint);
              tMap->Fill();
            }
          }
        }
      }

    fIsDone = true;
    return;
  }
}
