//author pf <pbutti@slac.stanford.edu>

//not sure of the meaning of this here. 
#pragma once

#include "Acts/Utilities/Units.hpp"

//Eigen definitions
#include "Acts/Utilities/Definitions.hpp"


//ACTS
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Material/MaterialProperties.hpp"

#include <memory>
#include <vector>
#include <map>

//ROOT
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "TGeoMatrix.h"
#include "TKey.h"

namespace Acts {
    class TrackingGeometry;
}

namespace TGeoUnits {
    double mm = 1000.;
    double cm = 100.;
    double dm = 10.;
    double m  = 1.;
}

class LdmxTracker {
    
public:
    //Instance of the ACTS tracking geometry
    //Made upon construction
    std::shared_ptr<const Acts::TrackingGeometry> _trackingGeo;
    
    //Geometry context
    std::shared_ptr<const Acts::GeometryContext*> _geoctx;
    
    
    LdmxTracker(const std::string & gdmlFile);
    ~LdmxTracker(){ if (_geo) delete _geo;}
    


    //Methods
    bool getPlanes(TGeoNode* tracker_node);
    TGeoNode* getNode(const std::string& nodeName);
    bool getTrackerVolumeInfo(TGeoNode* tracker_node);

    bool FillGlobalToLocalTrans(const TGeoMatrix* geoMatrix, 
                                Acts::RotationMatrix3D& rot_matrix,
                                Acts::Vector3D& trans_vec);
    
    void SetDebugMode(bool val) {_debug = val;};
    
    //TODO The override functions should go once geo is fixed!
    //Override functions
    
    //Change planes location
    void SetPlaneGlobalZPosition(int iplane, double z);
    void SetPlaneGlobalPosition(int iplane, const Acts::Vector3D& pos);
    
    //TODO FIX FIX FIX
    void SetVolumeDimension(double dx, double dy, double dz) {
        _tracker_dimension(0)=dx;_tracker_dimension(1)=dy;_tracker_dimension(2)=dz; }; 

    //Override the surface thickness
    void SetPlaneThickness(int iplane, double x){};

    //Build the surface configs
    bool BuildSurfaceConfigurations ();

    //Build the layers configs
    bool BuildLayerConfigurations ();

    //Build the volume configs
    bool BuildSubVolumesConfigurations();

    //Build Top Volume config
    bool BuildTopVolumeGeometry();
    
    //Build the tracking geometry
    bool BuildTrackingGeometry();

    //Print out the geometry that will be built
    void PrintOutGeometrySimple(); 
    
    //Material Properties for typical materials - TODO fill for other materials
    std::shared_ptr<Acts::MaterialProperties> getMaterial(const std::string& material,double thickness);
    
    //Get tracker surfaces
    std::map<uint, std::shared_ptr<const Acts::Surface>> surfaces();


private:

    struct TGeoManager_shptr : public TGeoManager {};
    
    // Global (z,x,y) coordinates of the centers of the tracking planes
    std::vector<Acts::Vector3D> planePos;
    std::vector<Acts::Vector3D> planeBoundaries; 
    std::vector<std::string> planeNames;
    std::vector<double> planeThickness;
    /// Unit vector pointing in direction of increasing local coordinate u of each
    /// tracking plane in terms of global coordinates (z,x,y).
    /// 
    std::vector<Acts::Vector3D> planeDirU;
    std::vector<Acts::Vector3D> planeDirV;
    std::vector<Acts::Vector3D> planeDirW; //needed?
    
    //TODO better storage?
    std::vector<Acts::Vector3D> global_to_local_t;
    std::vector<Acts::RotationMatrix3D > global_to_local_rot;
    
    //Shared didn't work
    //std::shared_ptr<TGeoManager> _geoManager;
    
    bool _debug{false};
    TGeoManager* _geo{nullptr};
    
    //Global Axes orientation
    Acts::Vector3D _globXhat {1,0,0};
    Acts::Vector3D _globYhat {0,1,0};
    Acts::Vector3D _globZhat {0,0,1};
    double _thickness {-999};

    //Global Tracker volume info
    Acts::Vector3D _tracker_position{-999,-999,-999};
    Acts::Vector3D _tracker_dimension{0.,0.,0.};
    std::string _tracker_name{"tagger"};

    //There will be bunch of duplicated information in this class
    std::vector<Acts::CuboidVolumeBuilder::SurfaceConfig> _SurfaceConfigs;
    std::vector<Acts::CuboidVolumeBuilder::LayerConfig>   _LayerConfigs;
    std::vector<Acts::CuboidVolumeBuilder::VolumeConfig>  _VolumeConfigs;
    
    Acts::CuboidVolumeBuilder::Config _trackerConfig;
    
    Acts::CuboidVolumeBuilder _cuboidVolumeBuilder;

    bool _surfacesBuilt{false};
    bool _layersBuilt{false};
    bool _subVolumesBuilt{false};

    Acts::GeometryContext* dummy_context;
    /// Instance of ACTS tracking geometry (used when employing ACTS machinery).
    /// 
    std::shared_ptr<const Acts::TrackingGeometry> tracker_detector;
};
