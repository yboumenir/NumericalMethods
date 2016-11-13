set(CMAKE_AUTOMOC ON)
set(QTDIR "/Users/yasserboumenir/Qt/5.7/clang_64")

add_compile_options(  -std=c++1y -pthread -MMD -fPIC -pipe -g3  ${WARNING_TOGGLES} )
if(EXISTS ${QTDIR})
  set(ENV{QTDIR}  ${QTDIR})
endif()


if(EXISTS $ENV{QTDIR})
  list(APPEND CMAKE_MODULE_PATH $ENV{QTDIR}/lib/cmake )
  list(APPEND CMAKE_PREFIX_PATH $ENV{QTDIR}/lib/cmake )
  message(STATUS "${Yellow} CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH} ${RESET}")
  find_package(Qt5 REQUIRED
    Core
    Widgets
    Gui
    Network
    OpenGL
    Test
    HINTS $ENV{QTDIR}/lib/cmake
          $ENV{QTDIR}/lib/cmake/Qt5 )
else()
  message(WARNING  "QTDIR environment variable not set... define QTDIR cmake var and re-run cmake.")
  set( QTDIR "/opt/pathToQt5/lib/cmake"  CACHE PATH "where are qt5 cmake files." )
endif()


include_directories(${QTDIR}/include)

message("QT Dir: " ${QTDIR})

