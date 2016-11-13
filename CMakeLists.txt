project(numericalmethods)
#set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if(APPLE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ ")
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/findQt.cmake")
#add_compile_options(" -std=c++11 -stdlib=libstdc++")

set(PATH_TO_QWT "$ENV{HOME}/source/qwt/qwt/src")
include_directories(${PATH_TO_QWT})

find_library(QWT_LIBRARY qwt HINTS "$ENV{HOME}/source/qwt/build/src")


set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/source/numericalmethods.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/source/numericalmatrix.cpp
        ${PATH_TO_QWT}
        )

set(HDS ${CMAKE_CURRENT_SOURCE_DIR}/include/numericalmethods.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/numericalmatrix.h)

set(target test)

add_library(numerical ${SRC} ${HDS} )
#set_target_properties(numerical PROPERTIES LINKER_LANGUAGE CXX)
add_executable(${target} source/${target}.cpp)
target_link_libraries(${target} numerical ${QWT_LIBRARY})
qt5_use_modules(${target}  Core Gui Widgets Network PrintSupport OpenGL OpenGLExtensions Svg)


