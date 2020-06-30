
add_executable(test main.cpp lalgebra.cpp)
target_link_libraries(test ${GTEST_LIB} pthread -lquadmath -lgtest)
