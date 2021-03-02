/**
 * @file MultiThreadUnitTest.cpp
 * @brief Define a Unit Test using Boost for validation of multi-thread computations
 *
 */
#include <iostream>
#include <chrono>

#include "pzlog.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "pzbdstrmatrix.h"
#include "pzbstrmatrix.h"
#include "TPZSpStructMatrix.h"

#include "TPZMatLaplacian.h"
#include "TPZGeoMeshTools.h"
#include "pzbndcond.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.testmultithread"));
#endif

namespace threadTest {
  constexpr int dim{2};//aux variable
  //aux function for creating 2d gmesh on unit square
  TPZGeoMesh *CreateGMesh(const int nDiv, int& matIdVol, int& matIdBC);
  //aux function for creating 2d cmesh with laplacian mat
  TPZCompMesh* CreateCMesh(TPZGeoMesh *gmesh, const int pOrder, const int matIdVol, const int matIdBC);


  //test the stiffness matrices in serial and parallel computations
  template <class TSTMAT>
  void CompareStiffnessMatrices(const int nThreads);

  template <class TSTMAT>
  void ComparePostProcError(const int nThreads);

  static void ForcingFunction(const TPZVec<REAL>& pt, TPZVec<STATE> &result);
  static void ExactSolution(const TPZVec<REAL> &pt, TPZVec<STATE> &sol, TPZFMatrix<STATE> &solDx);
}

int main()
{
  threadTest::ComparePostProcError<TPZSkylineStructMatrix>(4);
}



TPZGeoMesh *threadTest::CreateGMesh(const int nDiv, int& matIdVol, int& matIdBC)
{
  constexpr MMeshType meshType = MMeshType::ETriangular;

  static TPZManVector<REAL,3> minX(3,0);
  static TPZManVector<REAL,3> maxX(3,1);
  maxX[2] = 0.;
  TPZVec<int> nDivs(dim,nDiv);
  TPZManVector<int,5> matIdVec(5, -1);
  matIdVec[0] = 1;
  matIdVol = matIdVec[0];
  matIdBC = matIdVec[1];
  return TPZGeoMeshTools::CreateGeoMeshOnGrid(threadTest::dim,minX,maxX,matIdVec,nDivs,meshType,true);
}

TPZCompMesh *threadTest::CreateCMesh(TPZGeoMesh *gmesh, const int pOrder, const int matIdVol, const int matIdBC)
{
  auto *cmesh = new TPZCompMesh(gmesh);
  auto *laplacianMat = new TPZMatLaplacian(matIdVol, threadTest::dim);

  TPZDummyFunction<STATE>* forcingFunction = new TPZDummyFunction<STATE>(&threadTest::ForcingFunction, pOrder);
  TPZDummyFunction<STATE>* exactSol = new TPZDummyFunction<STATE>(&threadTest::ExactSolution, pOrder);
  laplacianMat->SetForcingFunctionExact(forcingFunction);
  laplacianMat->SetForcingFunction(exactSol);

  TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
  int bctype = 0;
  val2.Zero();
  TPZBndCond *bc = laplacianMat->CreateBC(laplacianMat, matIdBC, bctype, val1, val2);
  bc->TPZMaterial::SetForcingFunction(exactSol);

  cmesh->InsertMaterialObject(laplacianMat);
  cmesh->InsertMaterialObject(bc);

  cmesh->SetDefaultOrder(pOrder);
  cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
  cmesh->AutoBuild();

  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();

  return cmesh;
}

template <class TSTMAT>
void threadTest::ComparePostProcError(const int nThreads)
{
    constexpr int nDiv{4};
    constexpr int pOrder{3};

    int matIdVol;
    int matIdBC;
    auto *gMesh = CreateGMesh(nDiv, matIdVol, matIdBC);
    auto *cMesh = CreateCMesh(gMesh, pOrder, matIdVol, matIdBC);

    constexpr bool optimizeBandwidth{false};

    auto GetErrorVec = [cMesh, optimizeBandwidth](const int nThreads){
        TPZAnalysis an(cMesh, optimizeBandwidth);
        TSTMAT matskl(cMesh);
        matskl.SetNumThreads(nThreads);
        an.SetStructuralMatrix(matskl);

        TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
        direct->SetDirect(ELDLt);
        an.SetSolver(*direct);
        delete direct;
        direct = 0;

        TPZFunction<STATE>* exactSol = new TPZDummyFunction<STATE>(threadTest::ExactSolution, pOrder);
        an.SetExact(threadTest::ExactSolution);
        an.Assemble();
        an.Solve();

        an.SetThreadsForError(nThreads);

        TPZManVector<REAL> errorVec(3, 0.);
        int64_t nelem = cMesh->NElements();
        cMesh->LoadSolution(cMesh->Solution());
        cMesh->ExpandSolution();
        cMesh->ElementSolution().Redim(nelem, 10);

        an.PostProcessError(errorVec);
        return errorVec;
    };

    auto errorVecSer = GetErrorVec(0);
    std::cout << "Serial errors: " << errorVecSer << std::endl;

    auto errorVecPar = GetErrorVec(nThreads);
    std::cout << "Parallel errors: " << errorVecPar << std::endl;

    bool pass = true;
    for (auto i = 0; i < errorVecSer.size(); i++) {
        if (!IsZero(errorVecSer[i] - errorVecPar[i])) {
            pass = false;
        }
    }

    delete gMesh;
}

static void threadTest::ForcingFunction(const TPZVec<REAL>& pt, TPZVec<STATE> &result) {
    const auto x = pt[0], y = pt[1];
    result[0] = 8 * M_PI * M_PI * std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y);
}

static void threadTest::ExactSolution(const TPZVec<REAL> &pt, TPZVec<STATE> &sol, TPZFMatrix<STATE> &solDx) {
    const auto x = pt[0], y = pt[1];
    sol[0] = std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y);
    solDx(0, 0) = 2 * M_PI * std::cos(2 * M_PI * x) * std::sin(2 * M_PI * y);
    solDx(1, 0) = 2 * M_PI * std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y);
}
