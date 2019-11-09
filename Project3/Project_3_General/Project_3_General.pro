TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        CoordinateSystem.cpp \
        GaussianQuadrature.cpp \
        MonteCarlo.cpp \
        lib.cpp \
        main.cpp

HEADERS += \
    CoordinateSystem.h \
    GaussianQuadrature.h \
    MonteCarlo.h \
    lib.h
