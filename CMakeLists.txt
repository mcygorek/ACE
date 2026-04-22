cmake_minimum_required(VERSION 3.18)
project(ACE VERSION 1.0 LANGUAGES CXX)

# ── C++ standard ──────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# All targets need -fPIC (shared lib + Python module)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ── Default build type ────────────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Debug Release RelWithDebInfo MinSizeRel)
endif()

# ── Dependencies ──────────────────────────────────────────────────────────────
find_package(Threads REQUIRED)

# Eigen3: try cmake config/module first, then fall back to EIGEN_HOME
find_package(Eigen3 3.4 QUIET)
if (Eigen3_FOUND)
  message(STATUS "Found Eigen3: ${Eigen3_DIR}")
endif()
if(NOT Eigen3_FOUND)
  set(_eigen_home "$ENV{EIGEN_HOME}")
  if("${_eigen_home}" STREQUAL "")
    set(_eigen_home "/usr/include/eigen3")
  endif()
  if(EXISTS "${_eigen_home}/Eigen/Eigen")
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${_eigen_home}")
    message(STATUS "Using Eigen3 from: ${_eigen_home}")
  else()
    # No local Eigen, download it
    message(STATUS "Eigen3 not found locally – fetching via FetchContent")
    include(FetchContent)
    FetchContent_Declare(eigen3
      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
      GIT_TAG        3.4.0
      GIT_SHALLOW    TRUE)
    set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_DOC     OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(eigen3)
  endif()
endif()

# ── Library ───────────────────────────────────────────────────────────────────
# Compile sources once into an OBJECT library so the same .o files can be
# used by both libACE.so (C++ consumers) and the self-contained Python module.
file(GLOB LIB_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_library(ACE_objects OBJECT ${LIB_SOURCES})

target_include_directories(ACE_objects
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(ACE_objects PUBLIC Eigen3::Eigen Threads::Threads)

target_compile_options(ACE_objects PRIVATE
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-fno-math-errno>
  $<$<NOT:$<OR:$<CXX_COMPILER_ID:Intel>,$<CXX_COMPILER_ID:MSVC>>>:-Wno-ignored-attributes>
  $<$<CXX_COMPILER_ID:MSVC>:/EHsc /utf-8>)

target_compile_definitions(ACE_objects PRIVATE
  # M_PI and friends; needed on MSVC and MinGW (MSYS2) where they are
  # not exposed without this define.
  $<$<BOOL:${WIN32}>:_USE_MATH_DEFINES>
  # Prevent windows.h from defining min/max macros that clash with Eigen.
  $<$<BOOL:${WIN32}>:NOMINMAX>)

target_precompile_headers(ACE_objects PRIVATE include/PCH.hpp)

# Shared library – thin wrapper for C++ consumers.
add_library(ACE SHARED $<TARGET_OBJECTS:ACE_objects>)
set_target_properties(ACE PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
target_include_directories(ACE
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(ACE PUBLIC Eigen3::Eigen Threads::Threads)

# ── Main executable ───────────────────────────────────────────────────────────
add_executable(ACE_bin src_exec/ACE.cpp)
set_target_properties(ACE_bin PROPERTIES OUTPUT_NAME ACE)
target_link_libraries(ACE_bin PRIVATE ACE)

# ── Tool executables ──────────────────────────────────────────────────────────
option(BUILD_TOOLS "Build tool executables (src_exec/tools/)" OFF)

if(BUILD_TOOLS)
  file(GLOB TOOL_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src_exec/tools/*.cpp")

  foreach(TOOL_SRC IN LISTS TOOL_SOURCES)
    get_filename_component(TOOL_NAME "${TOOL_SRC}" NAME_WE)
    add_executable("${TOOL_NAME}" "${TOOL_SRC}")
    target_link_libraries("${TOOL_NAME}" PRIVATE ACE)
  endforeach()
endif()

# ── Python bindings (nanobind) ────────────────────────────────────────────────
option(BUILD_PYTHON_BINDINGS "Build Python bindings with nanobind" ON)

if(BUILD_PYTHON_BINDINGS)
  if(SKBUILD)
    # ── scikit-build-core path ───────────────────────────────────────────────
    # scikit-build-core sets CMAKE_PREFIX_PATH so FindPython locates the right
    # interpreter, but nanobind still requires find_package(Python) to have run.
    find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)
    find_package(nanobind CONFIG REQUIRED)
  else()
    # ── Manual / developer path ──────────────────────────────────────────────
    # Enforce virtualenv – never silently use system Python.
    if(NOT DEFINED Python_ROOT_DIR)
      if(NOT "$ENV{VIRTUAL_ENV}" STREQUAL "")
        set(Python_ROOT_DIR "$ENV{VIRTUAL_ENV}")
        message(STATUS "Using active virtualenv: $ENV{VIRTUAL_ENV}")
      elseif(NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
        set(Python_ROOT_DIR "$ENV{CONDA_PREFIX}")
        message(STATUS "Using active conda env: $ENV{CONDA_PREFIX}")
      else()
        message(FATAL_ERROR
          "BUILD_PYTHON_BINDINGS is ON but no Python virtualenv is active.\n"
          "Either activate one first:\n"
          "  source /path/to/.venv/bin/activate && cmake -B build\n"
          "or point CMake at an existing env explicitly:\n"
          "  cmake -DPython_ROOT_DIR=/path/to/.venv -B build")
      endif()
    endif()
    set(Python_FIND_VIRTUALENV ONLY)
    find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)
    message(STATUS "Python interpreter: ${Python_EXECUTABLE}")

    # nanobind: already importable -> use it; else pip install (safe in venv).
    find_package(nanobind CONFIG QUIET)
    if(NOT nanobind_FOUND)
      execute_process(
        COMMAND "${Python_EXECUTABLE}" -c
                "import nanobind; print(nanobind.cmake_dir())"
        OUTPUT_VARIABLE _nb_cmake_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _nb_import_result
        ERROR_QUIET)
      if(NOT _nb_import_result EQUAL 0)
        message(STATUS "nanobind not in env – running: pip install nanobind")
        execute_process(
          COMMAND "${Python_EXECUTABLE}" -m pip install nanobind
          RESULT_VARIABLE _pip_result
          ERROR_VARIABLE  _pip_error)
        if(NOT _pip_result EQUAL 0)
          message(FATAL_ERROR
            "pip install nanobind failed:\n${_pip_error}\n"
            "Install it manually:  pip install nanobind")
        endif()
        execute_process(
          COMMAND "${Python_EXECUTABLE}" -c
                  "import nanobind; print(nanobind.cmake_dir())"
          OUTPUT_VARIABLE _nb_cmake_dir
          OUTPUT_STRIP_TRAILING_WHITESPACE
          RESULT_VARIABLE _nb_import_result)
        if(NOT _nb_import_result EQUAL 0)
          message(FATAL_ERROR "nanobind installed but cmake_dir() not found.")
        endif()
      endif()
      message(STATUS "Using nanobind from: ${_nb_cmake_dir}")
      list(APPEND CMAKE_PREFIX_PATH "${_nb_cmake_dir}")
      find_package(nanobind CONFIG REQUIRED)
    endif()
  endif()

  # Compiled extension is named _ACE so the ACE/ Python package can wrap it.
  nanobind_add_module(ACE_nb NB_STATIC
    pybind/pybind_nb.cpp
    $<TARGET_OBJECTS:ACE_objects>)
  set_target_properties(ACE_nb PROPERTIES OUTPUT_NAME _ACE)
  target_include_directories(ACE_nb PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
  target_link_libraries(ACE_nb PRIVATE Eigen3::Eigen Threads::Threads)
  target_compile_definitions(ACE_nb PRIVATE
    $<$<BOOL:${WIN32}>:_USE_MATH_DEFINES>
    $<$<BOOL:${WIN32}>:NOMINMAX>)

  if(SKBUILD)
    # Install the extension into the ACE/ sub-directory of the wheel.
    install(TARGETS ACE_nb LIBRARY DESTINATION ACE)
    # Install the pure-Python package files alongside it.
    install(DIRECTORY pybind/ACE/ DESTINATION ACE)
  endif()
endif()

# ── Install ───────────────────────────────────────────────────────────────────
# include(GNUInstallDirs)

# if(NOT SKBUILD)
#   install(TARGETS ACE
#     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

#   install(TARGETS ACE_bin
#     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

#   install(DIRECTORY include/
#     DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ACE
#     FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")
# endif()
