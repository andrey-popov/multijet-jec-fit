cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
project (jecfit CXX)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-Wall -Wextra -pedantic)


# External dependencies
find_package(Boost 1.63 COMPONENTS program_options REQUIRED)
find_package(ROOT 6 COMPONENTS Minuit2 REQUIRED)


# Main library
add_library(jecfit SHARED
    src/JetCorrDefinitions.cpp
    src/FitBase.cpp
    src/Nuisances.cpp
    src/PhotonJetBinnedSum.cpp
    src/PhotonJetRun1.cpp
    src/ZJetRun1.cpp
    src/MultijetBinnedSum.cpp
    src/MultijetCrawlingBins.cpp
    src/JetCorrConstraint.cpp
    src/Morphing.cpp
    src/Rebin.cpp
)
target_include_directories(jecfit PUBLIC include)
target_link_libraries(jecfit
    PUBLIC
        ROOT::Hist ROOT::MathCore ROOT::Matrix ROOT::RIO
)


# Auxiliary library with Python wrVappings
add_library(jecfit_pythonwrapping SHARED src/PythonWrapping.cpp)
target_link_libraries(jecfit_pythonwrapping
    PUBLIC jecfit
)


# Main application
add_executable(fit prog/fit.cpp)
target_link_libraries(fit
    PRIVATE
        jecfit
        ROOT::Minuit2
        Boost::program_options
)


# Some unit tests
add_subdirectory(tests)
