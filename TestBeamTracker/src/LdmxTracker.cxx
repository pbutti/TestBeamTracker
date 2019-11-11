//author pf <pbutti@slac.stanford.edu>

#include "TestBeamTracker/LdmxTracker.h"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"

#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"


LdmxTracker::LdmxTracker(const std::string& gdmlFile) { 
    
    //TODO remove hardcoded plane thickness!!
    _thickness = 0.320;
    
    //TODO what is this?
    dummy_context = new Acts::GeometryContext(); 
    
    _geo = new TGeoManager();
    _geo->Import(gdmlFile.c_str());
    
    //_geoManager = std::make_shared<TGeoManager>();
    //_geoManager->Import(gdmlFile.c_str());
    TGeoNodeMatrix* tracker_node = static_cast<TGeoNodeMatrix*>(getNode("tagger_1"));
    bool success = getTrackerVolumeInfo(tracker_node);
    success = getPlanes(tracker_node);
}

//Make it recurrent to find
TGeoNode* LdmxTracker::getNode(const std::string& nodeName) {
    
    TGeoNode* node = (TGeoNode*)_geo->GetListOfNodes()->FindObject(nodeName.c_str());
    return node;
}

bool LdmxTracker::getTrackerVolumeInfo(TGeoNode* tracker_node) {
    TGeoBBox* top_box = (TGeoBBox*) tracker_node->GetVolume()->GetShape();
    //TODO: Not sure if this actually work in a more general case!!
    double lo = -999;
    double hi = 999;
    //Axes:: 1 - X, 2 - Y, 3 - Z
    for (int i = 1; i<=3; ++i) {
        top_box->GetAxisRange(i,lo,hi);
        _tracker_position(i-1) = ((hi+lo) / 2.) * TGeoUnits::mm;
        //The max of the axis should be corresponding to the DX,DY,DZ which are half of the size. So multiply by 2
        _tracker_dimension(i-1) = hi*2 * TGeoUnits::mm;
    }
    
    if (_debug) {
        std::cout<<"tracker position"<<std::endl;
        std::cout<<_tracker_position.transpose()<<std::endl;

        std::cout<<"tracker dimension"<<std::endl;
        std::cout<<_tracker_dimension.transpose()<<std::endl;
    }
    
    return true;
}

bool LdmxTracker::getPlanes(TGeoNode* tracker_node) {

    int n_planes = tracker_node->GetNdaughters();
    for (int i_plane = 0; i_plane < n_planes; i_plane++) {
        TGeoNodeMatrix* plane_node = (TGeoNodeMatrix*) tracker_node->GetDaughter(i_plane);
        planeNames.push_back(plane_node->GetName());
        //Get the boxes
        TGeoBBox* plane_box = (TGeoBBox*) plane_node->GetVolume()->GetShape();
        //Get the half-sizes of the box - units of TGeo are in meters 
        

        //ATTENTION! 
        //Factor 10 wrong in DX and DY, corrected using cm instead of mm => need to be checked.
        
        Acts::Vector3D box_half_sizes{plane_box->GetDX()*TGeoUnits::cm, plane_box->GetDY()*TGeoUnits::cm, plane_box->GetDZ()*TGeoUnits::mm};
        if (_debug){
            std::cout<<"LdmxTracker::DEBUG plane box sizes [mm] = " <<std::endl;
            std::cout<< box_half_sizes<<std::endl;
        }
        //The matrices are needed to place the planes in the right position
        TGeoMatrix* plane_matrix = (TGeoMatrix*) plane_node->GetMatrix();
        //Get the global to local transformation
        Acts::RotationMatrix3D box_rot_matrix;
        Acts::Vector3D box_translation(3);
        if (!FillGlobalToLocalTrans(plane_matrix, box_rot_matrix, box_translation)) {
            std::cerr<<"LdmxTracker::ERROR in building the global to local transformations"<<std::endl;
            return false;
        }        
        
        //TODO - Check
        planeBoundaries.push_back(box_half_sizes);
        //Setup the global positions of the layers. The global position is pos = -1*translation //TODO CHECK THIS!!!
        planePos.push_back(box_translation);
        global_to_local_rot.push_back(box_rot_matrix);

        //Compute local axes in global ref frame - alternatively get the columns of the rotation matrix
        Acts::Vector3D dirU = box_rot_matrix * _globXhat;
        Acts::Vector3D dirV = box_rot_matrix * _globYhat;

        //Not necessary as it's along _globZhat
        Acts::Vector3D dirW = dirU.cross(dirV); 
        
        if (_debug) {
            std::cout<<dirU<<std::endl;
            std::cout<<dirV<<std::endl;
            std::cout<<dirW<<std::endl;
        }
        
        planeDirU.push_back(dirU);
        planeDirV.push_back(dirV);
        planeDirW.push_back(dirW);
    }//loop on planes
    
    return true;
 }

bool LdmxTracker::FillGlobalToLocalTrans(const TGeoMatrix* geoMatrix, 
                                         Acts::RotationMatrix3D& rot_matrix,
                                         Acts::Vector3D& trans_vec) {
    
    if (!geoMatrix->GetRotationMatrix() || !geoMatrix->GetTranslation())
        return false;
    
    
    //LDMX geometry rotation matrix
    //The rot matrix from geoMatrix has 9 elements: 0-8
    // 0 1 2        |0,0  0,1  0,2|
    // 3 4 5   =>   |1,0  1,1  1,2|
    // 6 7 8        |2,0  2,1  2,2|
    for (int i=0; i<9; i++) {
        rot_matrix(i / 3, i % 3) = geoMatrix->GetRotationMatrix()[i];
    }
    
    //TODO FIX THIS!
    //for (int i=0; i<3; i++) { 
    //I should use the conversion to mm s. TODO CHECK THIS
    //trans_vec(i) = geoMatrix->GetTranslation()[i] * TGeoUnits::mm;
    //}

    //Rotate the axes to (zyx)
    //TODO FIX!! - This is due to an hardcode in ACTS
    
    Acts::RotationMatrix3D AxesRotation; 
    double rotationAngle = M_PI * 0.5;
    Acts::Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Acts::Vector3D yPos(0., 1., 0.);
    Acts::Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));

    
    AxesRotation.col(0) = xPos;
    AxesRotation.col(1) = yPos;
    AxesRotation.col(2) = zPos;
    
    Acts::Vector3D originalPos = {geoMatrix->GetTranslation()[0],geoMatrix->GetTranslation()[1],geoMatrix->GetTranslation()[2]};
    trans_vec   = AxesRotation * originalPos;
    rot_matrix = AxesRotation.transpose()*(rot_matrix*AxesRotation);
        
    if (_debug) {
        std::cout<<"Transform Global to Local"<<std::endl;
        std::cout<<rot_matrix<<std::endl;
        std::cout<<trans_vec<<std::endl;
    }
    
    return true;   
}


void LdmxTracker::SetPlaneGlobalPosition(int iplane, int icoord, double val) {
    if (iplane > planePos.size()-1) {
        std::cerr<<"LdmxTracker::SetPlaneGlobalPosition::ERROR plane="<<iplane<<" not available. NPlanes "<<planePos.size()<<std::endl; 
        return;
    }
    
    if (icoord <0 || icoord > 2) {
        std::cerr<<"LdmxTracker::Coordinate must be between 0 and 2"<<std::endl;
    }
    
    planePos[iplane](icoord) = val;
    
}

void LdmxTracker::SetPlaneGlobalPosition(int iplane, const Acts::Vector3D& pos) {
    if (iplane > planePos.size()-1) {
        std::cerr<<"LdmxTracker::SetPlaneGlobalPosition::ERROR plane="<<iplane<<" not available. NPlanes "<<planePos.size()<<std::endl; 
        return;
    }
    planePos[iplane] = pos;
}

//A surface configuration needs:
//1) center position
//2) rotation
//3) bounds
//4) Material
//5) Thickness
//6) Optional detector elements (transform, rectangle bounds and thickness) => check acts-core
//TODO Remove hardcoded thickness
bool LdmxTracker::BuildSurfaceConfigurations () {
    
    //No planes were loaded from gdml
    if (planePos.size()==0)
        return false;
    
    for (int i_plane = 0 ; i_plane<planePos.size(); i_plane++) {
        Acts::CuboidVolumeBuilder::SurfaceConfig surf_cfg;
        //center position
        surf_cfg.position = planePos[i_plane];
        //surface rotation - Check if that is correct! TODO FIX FIX FIX 

        //================== Pass the identity. No  bending of the planes. All parallel! ==========///
        surf_cfg.rotation = global_to_local_rot[i_plane];
        
        
        //bounds - check is this correct?  - TODO REMOVE HARDCODING!
        surf_cfg.rBounds  = std::make_shared<const Acts::RectangleBounds>(planeBoundaries[i_plane](0),planeBoundaries[i_plane](1));
        
        
        //TODO fix the thickness from geometry!!
        std::shared_ptr<Acts::MaterialProperties> matProp = getMaterial("Si",_thickness);
        surf_cfg.surMat = std::shared_ptr<Acts::ISurfaceMaterial>(new Acts::HomogeneousSurfaceMaterial(*matProp));
        
        //TODO FIX THIS!
        surf_cfg.thickness = _thickness;
        //surf_cfg.thickness = 0.;

        //Create a detector element stub => update?
        //Trans/bounds/thickness are then taken from the surface configuration =)
        surf_cfg.detElementConstructor =
            [](std::shared_ptr<const Acts::Transform3D> trans,
               std::shared_ptr<const Acts::RectangleBounds> bounds, double thickness) {
            return new Acts::Test::DetectorElementStub(trans, bounds, thickness);
        };
        

        _SurfaceConfigs.push_back(surf_cfg);
    }    
    _surfacesBuilt = true;
    
    return _surfacesBuilt;
}
//thickness in mm
std::shared_ptr<Acts::MaterialProperties> LdmxTracker::getMaterial(const std::string & material, double thickness) {
    
    //http://pdg.lbl.gov/2011/AtomicNuclearProperties/HTML_PAGES/014.html
    if (material == "Silicon" || material == "Si" || material == "silicon") {
        return std::make_shared<Acts::MaterialProperties>(93.70,591.4,28.086,14,2.329e-3,thickness);
    }
    return std::make_shared<Acts::MaterialProperties>();
}


bool LdmxTracker::BuildLayerConfigurations() {
    
    //Without surfaces do not build the layers
    if (!_surfacesBuilt)
        return false;
    
    for (auto& surf_cfg : _SurfaceConfigs) {
        Acts::CuboidVolumeBuilder::LayerConfig lay_cfg;
        lay_cfg.surfaceCfg = surf_cfg;
        lay_cfg.active = true;
        _LayerConfigs.push_back(lay_cfg);
    }
    
    _layersBuilt = true;
    return _layersBuilt;
}


//TODO::You'll need 2 volumes? One for the tagger and one for the recoil?? For the moment I'll consider 
//a single configuration. Remove hardcodings!!
bool LdmxTracker::BuildSubVolumesConfigurations() {
    
    if (!_layersBuilt)
        return false;
    
    Acts::CuboidVolumeBuilder::VolumeConfig volumeConfig;
    volumeConfig.position = _tracker_position;
    volumeConfig.length   = _tracker_dimension;
    volumeConfig.layerCfg = _LayerConfigs;
    volumeConfig.name     = _tracker_name;
    
    _VolumeConfigs.push_back(volumeConfig);
    _subVolumesBuilt = true;
    return _subVolumesBuilt;
}


bool LdmxTracker::BuildTopVolumeGeometry() {
    
    _trackerConfig.position  = _tracker_position;
    _trackerConfig.length    = _tracker_dimension;
    _trackerConfig.volumeCfg = _VolumeConfigs;
    
    _cuboidVolumeBuilder.setConfig(_trackerConfig);
    return true;
}

bool LdmxTracker::BuildTrackingGeometry() { 
    Acts::TrackingGeometryBuilder::Config tgbCfg;
        
    tgbCfg.trackingVolumeBuilders.push_back(
        [=](const auto& dummy_context, const auto& placeholder_1, const auto& placeholder_2) 
        { return  _cuboidVolumeBuilder.trackingVolume(dummy_context,placeholder_1,placeholder_2);}
        );

    Acts::TrackingGeometryBuilder tgb(tgbCfg);
    
    std::cout<<"Create the Tracking Geometry"<<std::endl;
    
    auto tGeometry = tgb.trackingGeometry(dummy_context);
    
    tracker_detector = std::move(tGeometry);
    return true;
}

void LdmxTracker::PrintOutGeometrySimple() { 
    
    std::cout<<"-------------------------------------------------------------"<<std::endl;
    std::cout<<"Print out of the tracking geometry being built. Units in [mm]"<<std::endl;
    std::cout<<"-------------------------------------------------------------"<<std::endl;
    
    std::cout<<"Top Volume position : "<<_trackerConfig.position.transpose()<<std::endl;
    std::cout<<"Top Volume dimension: "<<_trackerConfig.length.transpose()<<std::endl;
    
    std::cout<<"-------------------------------------------------------------"<<std::endl;
    
    std::cout<<"SubVolumes:"<<std::endl;
    for (auto cfgvol : _VolumeConfigs) {
        std::cout<<"name: "<<cfgvol.name<<std::endl;
        std::cout<<"position: "<<cfgvol.position.transpose()<<std::endl;
        std::cout<<"length: "<<cfgvol.length.transpose()<<std::endl;

        for (auto cfglay : _LayerConfigs) {
            std::cout<<"---------------------------------"<<std::endl;
            std::cout<<"Layer is Active: "<<cfglay.active<<std::endl;
            std::cout<<"Surface:: "<<std::endl;
            std::cout<<"position::"<<cfglay.surfaceCfg.position.transpose()<<std::endl;
            std::cout<<"rotation::"<<std::endl;
            std::cout<<cfglay.surfaceCfg.rotation<<std::endl;
            std::cout<<"Half length x::"<<cfglay.surfaceCfg.rBounds->halflengthX()<<std::endl;
            std::cout<<"Half length y::"<<cfglay.surfaceCfg.rBounds->halflengthY()<<std::endl;
            std::cout<<"Material::"<<std::endl;
            
            std::dynamic_pointer_cast<const Acts::HomogeneousSurfaceMaterial>(cfglay.surfaceCfg.surMat)->toStream(std::cout);
            std::cout<<std::endl;
            
           
            
        }   
    }   
}


