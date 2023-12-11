cmake_minimum_required(VERSION 3.5)

project(ranfing LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

# path to gmsh SDK
set(GMSH_SDK_DIR "C:/gmsh-4.11.1-Windows64-sdk")

# set paths to headers and libs
include_directories(${GMSH_SDK_DIR}/include)
# include_directories("C:/msys64/mingw64/include/eigen3/Eigen/")
# include_directories("C:/msys64/mingw64/include/eigen3/")
link_directories(${GMSH_SDK_DIR}/lib)

set(HEADERS fiber.h
    body.h
    model.h
)

set(SOURCES main.cpp
    ${HEADERS}
    fiber.cpp
    body.cpp
    model.cpp
)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(ranfing ${SOURCES})

# linking gmsh library
target_link_libraries(${PROJECT_NAME} gmsh
                                        )

include(GNUInstallDirs)
install(TARGETS ranfing
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
