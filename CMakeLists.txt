cmake_minimum_required(VERSION 3.14)
project(RootForge)

set(CMAKE_CXX_STANDARD 14)

# 1. Automatically fetch and build GLFW from source
# 1. Automatically fetch and build GLFW from source (Bypassing Git)
include(FetchContent)
FetchContent_Declare(
        glfw
        URL https://github.com/glfw/glfw/archive/refs/tags/3.3.8.zip
)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(glfw)

# 2. Find OpenGL (Built into your OS)
find_package(OpenGL REQUIRED)

# 3. Define ImGui Source Files (Requires the 'imgui' folder in your project root)
set(IMGUI_SOURCES
        imgui/imgui.cpp
        imgui/imgui_draw.cpp
        imgui/imgui_tables.cpp
        imgui/imgui_widgets.cpp
        imgui/backends/imgui_impl_glfw.cpp
        imgui/backends/imgui_impl_opengl3.cpp
)

# 4. Your RootForge Sources
set(ROOTFORGE_SOURCES
        src/RootForge.cpp
        include/RootForgeCore.h
        include/NumericalAnalysis.h
        include/NumberTheory.h
        include/MathematicalFunctions.h
        include/RootForge.h
        include/Polynomial.h
        src/NumericalAnalysis.cpp
        src/MathematicalFunctions.cpp
        src/Polynomial.cpp
)

# 5. Compile the Executable
add_executable(RootForge ${ROOTFORGE_SOURCES} ${IMGUI_SOURCES})

# 6. Link everything together
# FIXED: Added 'include' so the compiler can resolve your #include "NumericalAnalysis.h" statements
target_include_directories(RootForge PRIVATE include imgui imgui/backends)
target_link_libraries(RootForge PRIVATE glfw OpenGL::GL)