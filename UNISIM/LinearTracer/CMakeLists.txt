set(HEADERS RSimulationCase.h RSimulatorConfiguration.h)
set(SOURCE RSimulatorConfiguration.cpp)

add_executable(LinearTracer main.cpp ${SOURCE} ${HEADERS} ${common_files})
target_include_directories(LinearTracer PUBLIC ${common_dir})

find_package(libInterpolate REQUIRED)
target_link_libraries(LinearTracer pz libInterpolate::Interpolate)

configure_file(../gmsh/Contorno.msh InputUNISIM/Contorno.msh COPYONLY)
configure_file(../Reservoir/base_unisim2.txt InputUNISIM/base_unisim2.txt COPYONLY)
configure_file(../Reservoir/tope_unisim2.txt InputUNISIM/tope_unisim2.txt COPYONLY)
configure_file(../PetroPhysics/krw_linear.txt InputUNISIM/krw_linear.txt COPYONLY)
configure_file(../PetroPhysics/krow_linear.txt InputUNISIM/krow_linear.txt COPYONLY)
