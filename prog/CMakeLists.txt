add_executable(fit fit.cpp)
target_link_libraries(fit jecfit ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
