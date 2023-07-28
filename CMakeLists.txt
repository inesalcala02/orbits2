cmake_minimum_required(VERSION 3.0)  # Versión mínima requerida de CMake

project(Orbits)  # Nombre del proyecto
#find_package(FFTW3 REQUIRED)

# Configurar el ejecutable
add_executable(orbits main.cpp)  # Nombre del ejecutable y lista de archivos fuente
add_executable(fftw fftw3.cpp)

# Opcional: Configurar las opciones de compilación
set(CMAKE_CXX_STANDARD 14)  # Establecer el estándar de C++ (por ejemplo, C++11)

# Opcional: Agregar directorios de inclusión
# include_directories(path/to/include/directory)

# Opcional: Agregar bibliotecas

target_link_libraries(fftw fftw3)

# Opcional: Configurar opciones de compilación adicionales
# target_compile_options(MyExecutable PRIVATE option1 option2)

# Opcional: Configurar opciones de enlace adicionales
# target_link_options(MyExecutable PRIVATE option1 option2)

# Target de instalación
install(TARGETS orbits
        DESTINATION bin)
        
install(TARGETS fftw
        DESTINATION bin)
        
