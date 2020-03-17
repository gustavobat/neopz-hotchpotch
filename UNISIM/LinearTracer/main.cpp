//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $

#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzvec.h"

#include "pzlog.h"

#include <cmath>
#include <iostream>

#include <set>
#include <map>
#include <vector>
#include "TPZVTKGeoMesh.h"
#include "TPZExtendGridDimension.h"
#include "TPZMultiphysicsCompMesh.h"
#include "RSimulatorConfiguration.h"

#include "TMRSApproxSpaceGenerator.h"
#include "TMRSDataTransfer.h"
#include "TMRSSFIAnalysis.h"

#include <libInterpolate/Interpolate.hpp>

TMRSDataTransfer SettingHDivUNISIM();

void ComputeCoarseIndices(TPZGeoMesh* gmesh, TPZVec<int64_t>& coarseindices);


TPZGeoMesh* CreateGeoMeshWithTopeAndBase(std::string& geometry_file2D, int nLayers, bool print3DMesh, bool Is3DQ);

void ModifyTopeAndBase(TPZGeoMesh* gmesh, std::string& filename);

void ReadData(std::string name, bool print_table_Q, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z);


//UNISIM
void UNISIMHDiv();

//
int main() {

    InitializePZLOG();
    UNISIMHDiv();

    return 0;
}

void UNISIMHDiv() {

    int nLayers = 3;
    bool print3DMesh = true;
    constexpr bool IS3DQ = false;
    std::string geometry_file2D = "InputUNISIM/Contorno.msh";
    TPZGeoMesh* gmesh = CreateGeoMeshWithTopeAndBase(geometry_file2D, nLayers, print3DMesh, IS3DQ);

    for (auto gel:gmesh->ElementVec()) {
        if (gel->MaterialId() == 2) {
            gel->SetMaterialId(1);
        }
    }

    TMRSApproxSpaceGenerator aspace;
    aspace.SetGeometry(gmesh);
    std::string name = "NewMesh";
    std::cout << gmesh->NElements();
    aspace.PrintGeometry(name);

    TMRSDataTransfer sim_data = SettingHDivUNISIM();
    aspace.SetDataTransfer(sim_data);
    int order = 1;
    //Revisar match (   )

    bool must_opt_band_width_Q = true;
    int n_threads = 4;
    bool UsePardiso_Q = true;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh* mixed_operator = aspace.GetMixedOperator();


    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh* transport_operator = aspace.GetTransportOperator();
    std::ofstream file("mixed.vtk");

    TPZVTKGeoMesh::PrintCMeshVTK(mixed_operator, file);
    aspace.LinkMemory(mixed_operator, transport_operator);

    TMRSSFIAnalysis* sfi_analysis = new TMRSSFIAnalysis(mixed_operator, transport_operator, must_opt_band_width_Q);
    sfi_analysis->Configure(n_threads, UsePardiso_Q);
    sfi_analysis->SetDataTransfer(&sim_data);

    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;

    // Fill time steps vector
    TPZStack<REAL, 100> reporting_times;
    reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;

    REAL sim_time = 0.0;
    int pos = 0;
    REAL current_report_time = reporting_times[pos];

    for (int it = 1; it <= n_steps; it++) {

        sfi_analysis->RunTimeStep();
        sim_time = it * dt;
        sfi_analysis->m_transport_module->SetCurrentTime(sim_time);
        if (sim_time >= current_report_time) {
            std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
            if (it == 1) {
                sfi_analysis->PostProcessTimeStep(1);
                sfi_analysis->PostProcessTimeStep(2);
            } else {
                sfi_analysis->PostProcessTimeStep(2);
            }

            pos++;
            current_report_time = reporting_times[pos];

        }
    }


}

TPZGeoMesh* CreateGeoMeshWithTopeAndBase(std::string& geometry_file2D, int nLayers, bool print3DMesh, bool Is3DQ) {

    TPZGmshReader Geometry;
    TPZGeoMesh* gmesh2d;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.1");
    gmesh2d = Geometry.GeometricGmshMesh(geometry_file2D);
    Geometry.PrintPartitionSummary(std::cout);

    TPZGeoMesh * returnedMesh = nullptr;
    if (!Is3DQ){
        std::string filename1 = "InputUNISIM/tope_unisim2.txt";
        ModifyTopeAndBase(gmesh2d,filename1 );
        returnedMesh =gmesh2d;
    }
    return returnedMesh;
}

void ModifyTopeAndBase(TPZGeoMesh* gmesh, std::string& filename) {


    std::vector<double> x, y, z;
    ReadData(filename, true, x, y, z);


    _2D::ThinPlateSplineInterpolator<double> interp;
    interp.setData(x, y, z);

    int nCoordinates = gmesh->NodeVec().NElements();
    double sum = 0.0;
    for (auto val:z) {
        sum += val;
    }
    double val_storage = sum / z.size();

    for (int icoord = 0; icoord < nCoordinates; icoord++) {
        TPZGeoNode node = gmesh->NodeVec()[icoord];
        TPZVec<REAL> co(3);
        node.GetCoordinates(co);
        double val_interp = interp(co[0], co[1]);

        if (val_interp == 0.0) {
            co[2] = val_storage;
        }
        if (val_interp > 1.0) {
            co[2] = val_interp;
        }

        gmesh->NodeVec()[icoord].SetCoord(co);
    }


}

TMRSDataTransfer SettingHDivUNISIM() {

    TMRSDataTransfer sim_data;

    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["RockMatrix"] = 1;
//    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["RockMatrix2"] = 2;


    //zero flux boundaries
    int D_Type = 0;
    int N_Type = 1;
    REAL zero_flux = 0.0;
    REAL pressure_in = 30.0;
    REAL pressure_out = 20.0;
    int bcInlet = 3;
    int bcOutlet = 4;
    int bcZeroFlux1 = 5;
    int bcZeroFlux2 = 6;
    int bcZeroFlux3 = 7;
    int bcZeroFlux4 = 8;

    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);

    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(bcInlet, D_Type, pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(bcOutlet, D_Type, pressure_out);

    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(bcZeroFlux1, N_Type, zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(bcZeroFlux2, N_Type, zero_flux);

    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 1.0;

    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(4);

    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(bcInlet, bc_inlet, sat_in);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(bcOutlet, bc_outlet, 0.0);

    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(bcZeroFlux1, bc_inlet, 0.0);


    //Relative permermeabilities
    TRSLinearInterpolator krw, kro;
    std::string name_krw("InputUNISIM/krw_linear.txt");
    std::string name_kro("InputUNISIM/krow_linear.txt");

    krw.ReadData(name_krw, true);
    kro.ReadData(name_kro, true);

    kro.SetLeftExtension(TRSLinearInterpolator::ELinear);
    krw.SetLeftExtension(TRSLinearInterpolator::ELinear);
    kro.SetRightExtension(TRSLinearInterpolator::ELinear);
    krw.SetRightExtension(TRSLinearInterpolator::ELinear);

//        kro.SetInterpolationType(TRSLinearInterpolator::InterpType::THermite);
//        krw.SetInterpolationType(TRSLinearInterpolator::InterpType::THermite);

    sim_data.mTPetroPhysics.mLayer_Krw_RelPerModel.resize(1);
    sim_data.mTPetroPhysics.mLayer_Kro_RelPerModel.resize(1);
    sim_data.mTPetroPhysics.mLayer_Krw_RelPerModel[0] = krw;
    sim_data.mTPetroPhysics.mLayer_Kro_RelPerModel[0] = kro;

    // Fractional flows composition
    auto fw = [](TRSLinearInterpolator& krw, TRSLinearInterpolator& kro, double sw, double p) {

        double mu_w = 1.0;
        double mu_o = 1.0;
        double Bw = 1.0;
        double Bo = 1.0;

        std::tuple<double, double> krw_t = krw.ValDeriv(sw);
        std::tuple<double, double> kro_t = kro.ValDeriv(sw);
        double krwv = std::get<0>(krw_t);
        double krov = std::get<0>(kro_t);
        double dkrw_dswv = std::get<1>(krw_t);
        double dkro_dswv = std::get<1>(kro_t);
        double lwv = krwv / (mu_w * Bw);
        double dlw_dswv = dkrw_dswv / (mu_w * Bw);
        double lov = krov / (mu_o * Bo);
        double dlo_dswv = dkro_dswv / (mu_o * Bo);
        double lv = lwv + lov;
        double dl_dswv = dlw_dswv + dlo_dswv;
        double fwv = lwv / lv;
        double dfw_dswv = (dlw_dswv / lv) - lwv * (dl_dswv / (lv * lv));
        std::tuple<double, double, double> fw_t(sw, 1, 0.0);
        return fw_t;
    };

    std::function<std::tuple<double, double, double>(TRSLinearInterpolator&, TRSLinearInterpolator&, double sw,
                                                     double p)> fo = [](TRSLinearInterpolator& krw,
                                                                        TRSLinearInterpolator& kro, double sw,
                                                                        double p) {

        double mu_w = 1.0;
        double mu_o = 1.0;
        double Bw = 1.0;
        double Bo = 1.0;

        std::tuple<double, double> krw_t = krw.ValDeriv(sw);
        std::tuple<double, double> kro_t = kro.ValDeriv(sw);
        double krwv = std::get<0>(krw_t);
        double krov = std::get<0>(kro_t);
        double dkrw_dswv = std::get<1>(krw_t);
        double dkro_dswv = std::get<1>(kro_t);
        double lwv = krwv / (mu_w * Bw);
        double dlw_dswv = dkrw_dswv / (mu_w * Bw);
        double lov = krov / (mu_o * Bo);
        double dlo_dswv = dkro_dswv / (mu_o * Bo);
        double lv = lwv + lov;
        double dl_dswv = dlw_dswv + dlo_dswv;
        double fov = lov / lv;
        double dfo_dswv = (dlo_dswv / lv) - lov * (dl_dswv / (lv * lv));
        std::tuple<double, double, double> fo_t(1 - sw, -1, 0.0);
        return fo_t;
    };

    std::function<std::tuple<double, double, double>(TRSLinearInterpolator&, TRSLinearInterpolator&, double sw,
                                                     double p)> lambda = [](TRSLinearInterpolator& krw,
                                                                            TRSLinearInterpolator& kro, double sw,
                                                                            double p) {


        double mu_w = 1.0;
        double mu_o = 1.0;
        double Bw = 1.0;
        double Bo = 1.0;

        std::tuple<double, double> krw_t = krw.ValDeriv(sw);
        std::tuple<double, double> kro_t = kro.ValDeriv(sw);
        double krwv = std::get<0>(krw_t);
        double krov = std::get<0>(kro_t);
        double dkrw_dswv = std::get<1>(krw_t);
        double dkro_dswv = std::get<1>(kro_t);
        double lwv = krwv / (mu_w * Bw);
        double dlw_dswv = dkrw_dswv / (mu_w * Bw);
        double lov = krov / (mu_o * Bo);
        double dlo_dswv = dkro_dswv / (mu_o * Bo);
        double lv = lwv + lov;
        double dl_dswv = dlw_dswv + dlo_dswv;
        std::tuple<double, double, double> l_t(lv, dl_dswv, 0.0);
        return l_t;
    };

    sim_data.mTMultiphaseFunctions.mLayer_fw[0] = fw;
    sim_data.mTMultiphaseFunctions.mLayer_fo[0] = fo;
    sim_data.mTMultiphaseFunctions.mLayer_lambda[0] = lambda;

    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 20;
    sim_data.mTNumerics.m_max_iter_transport = 20;
    sim_data.mTNumerics.m_max_iter_sfi = 20;
    sim_data.mTNumerics.m_res_tol_mixed = 1.0e-4;
    sim_data.mTNumerics.m_corr_tol_mixed = 1.0e-4;
    sim_data.mTNumerics.m_res_tol_transport = 1.0e-4;
    sim_data.mTNumerics.m_corr_tol_transport = 1.0e-4;
    sim_data.mTNumerics.m_n_steps = 100;
    sim_data.mTNumerics.m_dt = 100000.0;
    sim_data.mTNumerics.m_four_approx_spaces_Q = false;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;


    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator_3d.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator_3d.vtk";
    TPZStack<std::string, 10> scalnames, vecnames;

    vecnames.push_back("Flux");
    scalnames.push_back("Pressure");
//    scalnames.push_back("kappa");
//    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
//        scalnames.Push("g_average");
//        scalnames.Push("u_average");
//    }
    sim_data.mTPostProcess.m_file_time_step = 20000.0;
    sim_data.mTPostProcess.m_vecnames = vecnames;
    sim_data.mTPostProcess.m_scalnames = scalnames;

    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;
    TPZStack<REAL, 100> reporting_times;
    REAL time = sim_data.mTPostProcess.m_file_time_step;
    int n_reporting_times = (n_steps) / (time / dt) + 1;
    REAL r_time = 0.0;
    for (int i = 1; i <= n_reporting_times; i++) {
        r_time += dt * (time / dt);
        reporting_times.push_back(r_time);
    }
    sim_data.mTPostProcess.m_vec_reporting_times = reporting_times;
    return sim_data;
}

void
ReadData(std::string name, bool print_table_Q, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z) {

    std::ifstream file;
    file.open(name);

    int i = 1;


    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        char l = line[0];
        if (l != '/') {
            i = i + 1;
            int val = i % 15;
            if (val == 0) {
                double a, b, c;
                if (iss >> a >> b >> c);
                x.push_back(a);
                y.push_back(b);
                z.push_back(c);
            };
        };
    };

    if (x.empty()) {
        std::cout << "No data read." << std::endl;

        DebugStop();
    }
    if (print_table_Q) {
        std::cout << "*************************" << std::endl;
        std::cout << "Reading file... ok!" << std::endl;
        std::cout << "*************************" << std::endl;
        std::cout << x.size() << std::endl;
        std::cout << y.size() << std::endl;
        std::cout << z.size() << std::endl;
    }

}
