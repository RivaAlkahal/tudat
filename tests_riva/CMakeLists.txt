set(CMAKE_BUILD_TYPE Debug)


TUDAT_ADD_EXECUTABLE(computeCovariances_test_skipping
        "computeCovariances_test_skipping.cpp"
        ${Tudat_ESTIMATION_LIBRARIES}
)
add_definitions(-DESIGN_USE_THREADS)
find_package(Threads REQUIRED)
target_link_libraries(computeCovariances_test_skipping PRIVATE Threads::Threads)
find_package(OpenMP REQUIRED)
target_link_libraries(computeCovariances_test_skipping PRIVATE OpenMP::OpenMP_CXX)


TUDAT_ADD_EXECUTABLE(computeCovariance_TestMultiSat
        "computeCovariance_TestMultiSat.cpp"
        ${Tudat_ESTIMATION_LIBRARIES}
)
add_definitions(-DESIGN_USE_THREADS)
find_package(Threads REQUIRED)
target_link_libraries(computeCovariance_TestMultiSat PRIVATE Threads::Threads)
find_package(OpenMP REQUIRED)
target_link_libraries(computeCovariance_TestMultiSat PRIVATE OpenMP::OpenMP_CXX)


TUDAT_ADD_EXECUTABLE(computeCovariances_ForEmps
        "computeCovariances_ForEmps.cpp"
        ${Tudat_ESTIMATION_LIBRARIES}
)
add_definitions(-DESIGN_USE_THREADS)
find_package(Threads REQUIRED)
target_link_libraries(computeCovariances_ForEmps PRIVATE Threads::Threads)
find_package(OpenMP REQUIRED)
target_link_libraries(computeCovariances_ForEmps PRIVATE OpenMP::OpenMP_CXX)

TUDAT_ADD_EXECUTABLE(parallelArcProcess
        "parallelArcProcess.cpp"
        ${Tudat_ESTIMATION_LIBRARIES}
)
target_include_directories(parallelArcProcess PRIVATE /home/ralkahal/tudat-bundle/external/json)

TUDAT_ADD_EXECUTABLE(parallelArcProcessVikings
        "parallelArcProcessVikings.cpp"
        ${Tudat_ESTIMATION_LIBRARIES}
)
target_include_directories(parallelArcProcessVikings PRIVATE /home/ralkahal/tudat-bundle/external/json)

TUDAT_ADD_EXECUTABLE(parallelArcProcessMultiSat
        "parallelArcProcessMultiSat.cpp"
        ${Tudat_ESTIMATION_LIBRARIES}
)
target_include_directories(parallelArcProcessMultiSat PRIVATE /home/ralkahal/tudat-bundle/external/json)

TUDAT_ADD_EXECUTABLE(parallelArcProcessEmpiricals
        "parallelArcProcessEmpiricals.cpp"
        ${Tudat_ESTIMATION_LIBRARIES}
)
target_include_directories(parallelArcProcessEmpiricals PRIVATE /home/ralkahal/tudat-bundle/external/json)


TUDAT_ADD_EXECUTABLE(generate_configs
        "generate_configs.cpp"
        ${Tudat_ESTIMATION_LIBRARIES}
)
target_include_directories(generate_configs PRIVATE /home/ralkahal/tudat-bundle/external/json)
