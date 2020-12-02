project(newmoon_test)
add_executable(test main.cpp lalgebra.cpp jd_clock.cpp)
target_link_libraries(test ${GTEST_LIB} pthread -lquadmath -lgtest)
