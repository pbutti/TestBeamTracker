#include "TestBeamTracker/LdmxTracker.h"
#include "TestBeamTracker/common_includes.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <cstdlib>



//------------------------------------------------------------------------------

void usage(const string & executable) {
  cout << "Usage: " << executable << " [options]" << "\n"
       << "Options:\n"
       << "  -h/--help                      print help\n"
       << "  -i/--input FILE_NAME           path to input ROOT file\n"
       << "  -n/--nEntries N_ENTRIES        number of entries from input to process\n"
       << "\n";
}

//------------------------------------------------------------------------------

/// 
/// Function specifying what to do at each surface of detector.
/// Note: The detector consists of:
///         * Tracking volume
///         * Layers
///         * Surfaces
/// 
void surfaceVisitor(const Acts::Surface * surf) {
    /*
  if (!surf) {
    cout << "        Surface is null!\n";
    return;
  }

  Acts::Vector3D center = surf->center();
  Acts::Vector3D normal = surf->normal();
  cout << "        Surface(center=(" << center.x() << ','
                                     << center.y() << ','
                                     << center.z() << "),\n"
       << "                normal=(" << normal.x() << ','
                                     << normal.y() << ','
                                     << normal.z() << "))\n";
    */
}

//------------------------------------------------------------------------------

int main(int argc, char ** argv) {
    
    

    bool debug = false;
    float scaling_factor = 1.;
    if (argc > 1) 
        debug = (bool)std::stoi(argv[1]);
    if (argc > 2)
        scaling_factor = std::stof(argv[2]);
    
    
    LdmxTracker tracker("/nfs/slac/g/hps2/pbutti/LDMX/sw/TestBeamTracker/run/taggerLdmx_test.gdml");
    
    //Overwrite wrong parsing effects. Layers every 100 mm, separated by 3mm
    double z_layers[14] = {100,103,200,203,300,303,400,403,500,503,600,603,700,703};
    double x_layers[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double y_layers[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    //for (int i=0; i<14; i++) {
    //  tracker.SetPlaneGlobalZPosition(i,scaling_factor*z_layers[i]);
    //}

    for (int i=0; i<14; i++) {
        //Acts::Vector3D pos{z_layers[i],x_layers[i],y_layers[i]};
        //tracker.SetPlaneGlobalPosition(i,pos);
        tracker.SetPlaneGlobalPosition(i,0,z_layers[i]);
    }
    
    tracker.SetVolumeDimension(100,100,scaling_factor*1000);
    
    tracker.BuildSurfaceConfigurations();
    tracker.BuildLayerConfigurations();
    tracker.BuildSubVolumesConfigurations();
    tracker.BuildTopVolumeGeometry();
    tracker.BuildTrackingGeometry();
    if (debug)
        tracker.PrintOutGeometrySimple();
  
}
