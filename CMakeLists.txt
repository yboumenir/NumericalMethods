project(numericalmethods)

add_compile_options(-std=c++1y )

set(SRC
        source/numericalmethods.cpp
        source/numericalmatrix.cpp)

set(HDS

        include/numericalmethods.h
        include/numericalmatrix.h)

    set(target test)

add_library(numerical ${SRC} ${HDS})

add_executable(${target} source/${target}.cpp)
target_link_libraries(${target} numerical)


