################################################################################
## Data analysis project based on J-PET Framework
## Created by J-PET Framework developers 2016-2018
##
## Description:
##   Extension of the LargeBarrelAnalysis example
##   with Physical Hit Class - a simplified Hit Class
################################################################################

cmake_minimum_required(VERSION 3.1...3.14)

if(${CMAKE_VERSION} VERSION_LESS 3.14)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
    cmake_policy(VERSION 3.14)
endif()

################################################################################
## Project name
set(projectName HitTransforming)

################################################################################
# Adjust the list below to contain your custom class names                     #
################################################################################
set(USER_DEFINED_CLASSES JPetPhysHit)

## Auxiliary files
set(AUXILIARY_FILES README.md)

################################################################################
## Binary, header and source files definitions
set(projectBinary ${projectName}.x)
project(${projectName} CXX)

## Using source files of modules from LargeBarrelAnalysis
set(use_modules_from ../LargeBarrelAnalysis)

set(HEADERS
  ${use_modules_from}/EventCategorizerTools.h
  ${use_modules_from}/HitFinder.h
  ${use_modules_from}/HitFinderTools.h
  ${use_modules_from}/SignalFinder.h
  ${use_modules_from}/SignalFinderTools.h
  ${use_modules_from}/SignalTransformer.h
  ${use_modules_from}/TimeWindowCreator.h
  ${use_modules_from}/TimeWindowCreatorTools.h
  ${use_modules_from}/UniversalFileLoader.h
  ${use_modules_from}/ToTConverter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/JPetPhysHit.h
  ${CMAKE_CURRENT_SOURCE_DIR}/HitTransformer.h
)

set(SOURCES
  ${use_modules_from}/EventCategorizerTools.cpp
  ${use_modules_from}/HitFinder.cpp
  ${use_modules_from}/HitFinderTools.cpp
  ${use_modules_from}/SignalFinder.cpp
  ${use_modules_from}/SignalFinderTools.cpp
  ${use_modules_from}/SignalTransformer.cpp
  ${use_modules_from}/TimeWindowCreator.cpp
  ${use_modules_from}/TimeWindowCreatorTools.cpp
  ${use_modules_from}/UniversalFileLoader.cpp
  ${use_modules_from}/ToTConverter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/JPetPhysHit.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/HitTransformer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)

set(DICTIONARY_LIB UserClasses1)
set(DICTIONARY_NAME ${DICTIONARY_LIB}Dict)

set(LINKDEF_NAME LinkDef.h)

# extract includes of the JPetFramework library and dependencies
# needed for dictionary generation
get_property(DICTIONARY_INCLUDE_DIRS TARGET
  JPetFramework::JPetFramework PROPERTY INTERFACE_INCLUDE_DIRECTORIES
)
foreach(dir ${DICTIONARY_INCLUDE_DIRS})
  list(APPEND DICTIONARY_INCLUDES -I${dir})
endforeach()

get_property(DICTIONARY_INCLUDE_DIRS TARGET
Boost::log PROPERTY INTERFACE_INCLUDE_DIRECTORIES
)
foreach(dir ${DICTIONARY_INCLUDE_DIRS})
  list(APPEND DICTIONARY_INCLUDES -I${dir})
endforeach()

foreach(file ${USER_DEFINED_CLASSES})
  list(APPEND USER_DEFINED_CLASS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/${file}.h)
  list(APPEND USER_DEFINED_CLASS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${file}.cpp)
endforeach()

include_directories(ROOT_BUG) # Fix bug in NewMacros Root file
ROOT_GENERATE_DICTIONARY(
  ${DICTIONARY_NAME}
  ${USER_DEFINED_CLASS_HEADERS}
  MODULE ${DICTIONARY_LIB}
  LINKDEF ${LINKDEF_NAME}
  OPTIONS ${DICTIONARY_INCLUDES}
)

# ROOT 6 requires that dictionaries' symbols be stored in a shared lib
# rather than directly in the executable
add_library(${DICTIONARY_LIB} SHARED ${DICTIONARY_NAME}.cxx
  ${USER_DEFINED_CLASS_HEADERS} ${USER_DEFINED_CLASS_SOURCES}
)
target_link_libraries(
  ${DICTIONARY_LIB}
  JPetFramework::JPetFramework
  ROOT::Core
  ROOT::Hist
  ROOT::RIO
  ROOT::Flags_CXX
)

target_link_libraries(
  ${DICTIONARY_LIB}
  JPetFramework::JPetFramework ROOT::Core ROOT::Hist ROOT::RIO ROOT::Flags_CXX
)

add_executable(${projectBinary} ${SOURCES} ${HEADERS})
target_link_libraries(${projectBinary}
  ${DICTIONARY_LIB}
  JPetFramework::JPetFramework
  Boost::program_options)

add_custom_target(clean_data_${projectName}
  COMMAND rm -f *.tslot.*.root *.phys.*.root *.sig.root)

################################################################################
## Copy the example auxiliary files
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${AUXILIARY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
