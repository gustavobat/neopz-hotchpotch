/**
 * @file TestShape.cpp
 * @brief Defines a Unit Test to be used with Boost to check the generation of shape functions
 * of all kind of geometries available in the library.
 *
 */
 
#include <iostream>
#include <Mesh/pzgmesh.h>
#include <Post/TPZVTKGeoMesh.h>

template <class TGeo>
TPZGeoMesh * CreateSampleElementMesh();

template <class TGeo>
void CheckDivergenceOnInternalConnect();

int main() {
    std::cout << "Blob\n";
}

template <class TGeo>
TPZGeoMesh * CreateSampleElementMesh() {
    std::string elName = TGeo::TypeName();
    
    std::cout<<"Starting mesh refinement of element type: "<<elName<<std::endl;
    
    const int matId = 1;
    TPZManVector<REAL,3> lowerCorner(3,0);
    TPZManVector<REAL,3> size(3,1);
    TPZGeoMesh gmesh;
    TGeo::InsertExampleElement(gmesh, matId, lowerCorner, size);
    gmesh.BuildConnectivity();
    
    std::ofstream fileName(elName+".vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&gmesh, fileName);
}

template <class TGeo>
void CheckDivergenceOnInternalConnect() {

//    // Gets H(div) mesh
//    TPZCompMesh* fluxMesh = multiphysics->MeshVector()[0];
//
//    int64_t nel = fluxMesh->NElements();
//    for (int64_t i = 0; i < nel; i++) {
//        TPZCompEl* cel = fluxMesh->Element(i);
//
//        if (cel->Material()->Id() != 1) continue;
//        if (cel->Dimension() != 2) continue;
//
//        const auto fluxEl = dynamic_cast<TPZInterpolatedElement*> (cel);
//        if (!fluxEl) DebugStop();
//
//        // Gets flux element
//        const auto gel = fluxEl->Reference();
//
//        // Initialize material requirements
//        TPZMaterialData elData;
//        fluxEl->InitMaterialData(elData);
//
//        // Gets last connect, which contains internal shape functions
//        TPZConnect& con = cel->Connect(fluxEl->NConnects() - 1);
//        const int nInternalPhi = con.NShape();
//
//        // Iterates element edges (d-1 sides)
//        const int nNodes = gel->NCornerNodes();
//
//        // Creates integration rule on edge
//        const int pOrderIntRule = 3;
//        TPZIntPoints* intRule = gel->CreateSideIntegrationRule(gel->NSides() - 1, pOrderIntRule);
//
//        TPZManVector<REAL, 3> xi(2, 0);
//        REAL w;
//
//        // Stores results of the integration
//        int nShapeF = fluxEl->NShapeF();
//        TPZFMatrix<REAL> numericalIntegration(nShapeF, 1, 0);
//
//        TPZManVector<REAL, 2> gradPhi7(2, 0);
//        const int npts = intRule->NPoints();
//        for (auto ipt = 0; ipt < npts; ipt++) {
//            intRule->Point(ipt, xi, w);
//
//            fluxEl->ComputeRequiredData(elData, xi);
//            elData.ComputeFunctionDivergence();
//
//            if (ipt == 0) {
//                elData.Print(std::cout);
//            }
//
//            TPZManVector<REAL, 3> phi6(3, 0);
//            TPZManVector<REAL, 3> phi7(3, 0);
//            REAL div6 = elData.divphi(6, 0);
//            REAL div7 = elData.divphi(7, 0);
//
//            int vec6Id = elData.fVecShapeIndex[6].first;
//            int vec7Id = elData.fVecShapeIndex[7].first;
//
//            int phi6Id = elData.fVecShapeIndex[6].second;
//            int phi7Id = elData.fVecShapeIndex[7].second;
//
//            for (int i = 0; i < 3; i++) {
//                phi6[i] = elData.fMasterDirections(i, vec6Id) * elData.phi(phi6Id, 0);
//                phi7[i] = elData.fMasterDirections(i, vec7Id) * elData.phi(phi7Id, 0);
//            }
//
//            std::cout << "xi: " << xi << " phi6 = " << phi6 << " phi7 = " << phi7 << " div6 = " << div6
//                      << " div7 = " << div7 << std::endl;
//
//            gradPhi7[0] += w * elData.dphi(0, phi7Id);
//            gradPhi7[1] += w * elData.dphi(1, phi7Id);
//
//            std::cout << "gradPhi7: " << elData.dphi(0, phi7Id) << " " << elData.dphi(1, phi7Id) << '\n';
//            std::cout << "gradPhi6: " << elData.dphi(0, phi6Id) << " " << elData.dphi(1, phi6Id) << '\n';
//
//            for (int iPhi = 0; iPhi < nShapeF; iPhi++) {
//                numericalIntegration(iPhi, 0) += w * elData.divphi(iPhi, 0);
//            }
//
//        }
//
//        for (int i = 0; i < numericalIntegration.Rows(); i++) {
//            std::cout << numericalIntegration(i, 0) << '\n';
//        }
//
//        std::cout << "gradPhiX = " << gradPhi7[0] << " gradPhiY = " << gradPhi7[1] << '\n';
//
//        std::cout << '\n';
    
    }
}
