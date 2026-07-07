TEMPLATE = app
TARGET = AnalyticConvectionSolver
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    convection.cpp

HEADERS +=

QMAKE_CXXFLAGS_RELEASE += -O3

QMAKE_LFLAGS += -O3

# Requires freeglut (e.g. via vcpkg: vcpkg install freeglut:x64-windows).
# See README.md for the exact cl.exe build command actually used to build this project.
LIBS += -lopengl32 -lglu32 -lfreeglut

