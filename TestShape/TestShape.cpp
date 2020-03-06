/**
 * @file TestShape.cpp
 * @brief Defines a Unit Test to be used with Boost to check the generation of shape functions
 * of all kind of geometries available in the library.
 *
 */

#include <iostream>

#include <Mesh/pzgmesh.h>
#include <Mesh/pzcmesh.h>
#include <Mesh/pzintel.h>

#include <Geom/pzgeotriangle.h>
#include <Geom/pzgeoquad.h>
#include <Geom/pzgeotetrahedra.h>
#include <Geom/pzgeoprism.h>
#include <Geom/TPZGeoCube.h>
#include <Geom/pzgeopyramid.h>

#include <Post/TPZVTKGeoMesh.h>
#include <Material/TPZNullMaterial.h>
#include <Material/TPZVecL2.h>
#include <Geom/pzgeotetrahedra.h>

template<class TGeo>
void AddSampleElement(TPZGeoMesh& gmesh);

template<class TGeo>
void CheckDivergenceOnInternalConnect();

int main() {

    CheckDivergenceOnInternalConnect<pzgeom::TPZGeoTriangle>();
    CheckDivergenceOnInternalConnect<pzgeom::TPZGeoQuad>();
    CheckDivergenceOnInternalConnect<pzgeom::TPZGeoTetrahedra>();
    CheckDivergenceOnInternalConnect<pzgeom::TPZGeoPrism>();
    CheckDivergenceOnInternalConnect<pzgeom::TPZGeoCube>();
    // Pyramid fails since its space is not complete
    //CheckDivergenceOnInternalConnect<pzgeom::TPZGeoPyramid>();

    std::cout << "Blob\n";
}

template<class TGeo>
void AddSampleElement(TPZGeoMesh& gmesh) {
    std::string elName = TGeo::TypeName();

    std::cout << "Creating " << elName << " element.\n";

    const int matId = 1;
    TPZManVector<REAL, 3> lowerCorner(3, 0);
    TPZManVector<REAL, 3> size(3, 1);
    TGeo::InsertExampleElement(gmesh, matId, lowerCorner, size);
    gmesh.BuildConnectivity();
    gmesh.SetDimension(gmesh.Element(0)->Dimension());
    std::ofstream fileName(elName + ".vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&gmesh, fileName);
}

template<class TGeo>
void CheckDivergenceOnInternalConnect() {

    // Creates element mesh
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    AddSampleElement<TGeo>(*gmesh);
    int dim = gmesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(2);
    cmesh->SetDimModel(dim);

    TPZNullMaterial *mat = new TPZNullMaterial(1);
    mat->SetDimension(dim);
    cmesh->InsertMaterialObject(mat);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();
    cmesh->InitializeBlock();

    const auto fluxEl = dynamic_cast<TPZInterpolatedElement *>(cmesh->Element(0));
    if (!fluxEl) DebugStop();
    if (fluxEl->Material()->Id() != 1) DebugStop();

    // Gets flux geometric element
    const auto gel = fluxEl->Reference();
    if (!gel) DebugStop();

    // Initialize material requirements
    TPZMaterialData elData;
    fluxEl->InitMaterialData(elData);

    // Gets last connect, which is the one that contains internal shape functions
    TPZConnect &con = fluxEl->Connect(fluxEl->NConnects() - 1);

    // Creates integration rule on edge
    const int pOrderIntRule = fluxEl->EffectiveSideOrder(gel->NSides() - 1) * 2;
    TPZIntPoints *intRule = gel->CreateSideIntegrationRule(gel->NSides() - 1, pOrderIntRule);

    TPZManVector<REAL, 3> xi(dim, 0);
    REAL w;

    // Stores results of the integration
    const int nInternalPhi = con.NShape();
    const int firstInternalPhi = fluxEl->NShapeF() - nInternalPhi;
    
    //TPZFMatrix<REAL> integrationResult(nInternalPhi, 1, 0);
    TPZFMatrix<REAL> integrationResult(fluxEl->NShapeF(), 1, 0);

    const int npts = intRule->NPoints();
    for (auto ipt = 0; ipt < npts; ipt++) {
        intRule->Point(ipt, xi, w);

        fluxEl->ComputeRequiredData(elData, xi);
        elData.ComputeFunctionDivergence();

       // if (ipt == 0) {
       //     elData.Print(std::cout);
       // }
        
        //for (int iPhi = 0; iPhi < nInternalPhi; iPhi++) {
        for (int iPhi = 0; iPhi < fluxEl->NShapeF(); iPhi++) {
            //integrationResult(iPhi, 0) += w * elData.divphi(iPhi + firstInternalPhi, 0);
            integrationResult(iPhi, 0) += w * elData.divphi(iPhi, 0);
        }
    }
    
    bool isZero;
    for (int i = 0; i < integrationResult.Rows(); i++) {
   //     if (IsZero(integrationResult(i, 0))) {
            std::cout << integrationResult(i, 0) << '\n';
 //       }
    }
    std::cout << std::endl;
}
