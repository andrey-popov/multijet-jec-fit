add_executable(test_rebin test_rebin.cpp)
target_link_libraries(test_rebin PRIVATE jecfit)

add_executable(test_lossFunc test_lossFunc.cpp)
target_link_libraries(test_lossFunc PRIVATE jecfit)

