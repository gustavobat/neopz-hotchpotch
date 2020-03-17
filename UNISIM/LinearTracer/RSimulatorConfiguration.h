//
//  RSimulatorConfiguration.hpp
//  MonophasicTest
//  Created by Jose on 27/8/19.
//

#ifndef RSimulatorConfiguration_hpp
#define RSimulatorConfiguration_hpp

#include <cstdio>

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzbiharmonic.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"
#include "TPZHybridizeHDiv.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeopyramid.h"
#include "TPZGeoLinear.h"

#include <TPZRefPattern.h>

#include "TPZMaterial.h"
#include "pzelasmat.h"
#include "pzlog.h"

#include "TPZGenGrid2D.h"

#include <ctime>
#include <cstdio>

#include <cmath>

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>

#include <set>
#include <map>
#include <vector>
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"
#include "TPZExtendGridDimension.h"

#include "../TRSLinearInterpolator.h"
#include "TPZMatLaplacian.h"
#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "pzbndcond.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzsolve.h"
#include "TPZPersistenceManager.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include  "RSimulationCase.h"
#include "TPZMixedDarcyFlow.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzl2projection.h"
#include "../TPZTracerFlow.h"

class RSimulatorConfiguration{
    
private:
    TPZGeoMesh *fGmesh;
    TPZMultiphysicsCompMesh *multmesh;
    SimulationCase fsim_case;
    
public:
    RSimulatorConfiguration();
    explicit RSimulatorConfiguration(SimulationCase sim_case);
    
    explicit RSimulatorConfiguration(TPZGeoMesh *gmesh);

    void CreateGeomesh(int nx, int ny, double l, double h, MMeshType type);
    TPZGeoMesh * GetGeomesh();
    
    TPZCompMesh * CreateFluxCmesh(TPZGeoMesh * gmesh, int order);
    TPZCompMesh * CreatePressureCmesh(TPZGeoMesh * gmesh, int order);
    TPZCompMesh * CreateTransportMesh(TPZMultiphysicsCompMesh *cmesh, int ref);
    void CreateTransportElement(int p_order, TPZCompMesh *cmesh, TPZGeoEl *gel, bool is_BC);
    TPZMultiphysicsCompMesh *CreateMultiPhysicsCompMesh(TPZGeoMesh * gmesh);
       TPZMultiphysicsCompMesh * MPTransportMesh(TPZMultiphysicsCompMesh * mixed, TPZManVector<TPZCompMesh *> meshvec);
    
    TPZAnalysis * CreateAnalysis(TPZMultiphysicsCompMesh * cmesh_c,  bool must_opt_band_width_Q, int number_threads, bool UsePardiso_Q);
    
    void PrintCmesh(int mesh_index, std::ofstream &file_name);
    void PrintGmesh( std::ofstream &file_name);
    void PosProcess();
    void Run();
    void InsertTransportInterfaceElements(TPZMultiphysicsCompMesh *cmesh);
    void InsertInterfacesBetweenElements(int transport_matid, TPZCompMesh * cmesh, std::vector<int> & cel_indexes);
    TPZAnalysis * CreateTransportAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data);
    TPZFMatrix<STATE> TimeForward(TPZAnalysis * tracer_analysis, int & n_steps, REAL & dt, TPZFMatrix<STATE> & M_diag);
};

#endif /* RSimulatorConfiguration_hpp */
