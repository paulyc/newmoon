#set(EXTERN_INSTALL_DIR ${PROJECT_SOURCE_DIR}/extern)
#set(LIBRARY_PATHS ${EXTERN_INSTALL_DIR}/lib /usr/lib /usr/local/lib)

#set(LIBRARY_PATHS /usr/lib /usr/local/lib)

#include_directories(${EXTERN_INSTALL_DIR}/include)

#find_library(GTEST_LIBRARY gtest PATHS ${LIBRARY_PATHS})
#find_library(GTEST_MAIN_LIBRARY gtest_main PATHS ${LIBRARY_PATHS})
#find_library(GMOCK_LIBRARY gmock PATHS ${LIBRARY_PATHS})
#find_library(GMOCK_MAIN_LIBRARY gmock_main PATHS ${LIBRARY_PATHS})


add_executable(test main.cpp lalgebra.cpp)
target_include_directories(test PUBLIC ../include)
target_link_libraries(test ${GTEST_LIBRARY} pthread -lquadmath -lgtest)
