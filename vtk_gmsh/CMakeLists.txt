cmake_minimum_required(VERSION 3.10)
project(VTKGMSHProject)

# Поиск библиотеки VTK 9.1
find_package(VTK 9.1 REQUIRED)

# Подключение GMSH
add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/mnt/c/Users/tryto/OneDrive/Рабочий стол/семинары инфа/плюсы/2 семестр/lab1/vtk_gmsh/gmsh-4.13.1-Linux64-sdk/gmsh-4.13.1-Linux64-sdk/lib/libgmsh.so.4.13.1"
    INTERFACE_INCLUDE_DIRECTORIES "/mnt/c/Users/tryto/OneDrive/Рабочий стол/семинары инфа/плюсы/2 семестр/lab1/vtk_gmsh/gmsh-4.13.1-Linux64-sdk/gmsh-4.13.1-Linux64-sdk/include"
)

# Добавление исполняемых файлов
add_executable(simple simple.cpp)
add_executable(challenging challenging.cpp)

# Подключение библиотек к исполняемым файлам
target_link_libraries(simple ${VTK_LIBRARIES})
target_link_libraries(challenging ${VTK_LIBRARIES} gmsh)

# Установка директории для выходных файлов
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)