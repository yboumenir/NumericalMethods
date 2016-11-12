project(numericalmethods)
#set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ ")

#add_compile_options(" -std=c++11 -stdlib=libstdc++")

set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/source/numericalmethods.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/source/numericalmatrix.cpp)

set(HDS ${CMAKE_CURRENT_SOURCE_DIR}/include/numericalmethods.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/numericalmatrix.h)

    set(target test)

add_library(numerical ${SRC} ${HDS})

add_executable(${target} source/${target}.cpp)
target_link_libraries(${target} numerical)


