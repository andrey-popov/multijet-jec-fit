add_library(jecfit SHARED JetCorrDefinitions.cpp FitBase.cpp Nuisances.cpp
    PhotonJetBinnedSum.cpp PhotonJetRun1.cpp ZJetRun1.cpp MultijetBinnedSum.cpp
    MultijetCrawlingBins.cpp JetCorrConstraint.cpp Morphing.cpp Rebin.cpp)
target_link_libraries(jecfit ${ROOT_LIBRARIES})

add_library(jecfit_pythonwrapping SHARED PythonWrapping.cpp)
target_link_libraries(jecfit_pythonwrapping jecfit ${ROOT_LIBRARIES})
