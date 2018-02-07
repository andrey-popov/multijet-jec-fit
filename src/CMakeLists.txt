add_library(jecfit SHARED JetCorrDefinitions.cpp FitBase.cpp Nuisances.cpp
    PhotonJetBinnedSum.cpp PhotonJetRun1.cpp ZJetRun1.cpp MultijetBinnedSum.cpp Rebin.cpp)
target_link_libraries(jecfit ${ROOT_LIBRARIES})
