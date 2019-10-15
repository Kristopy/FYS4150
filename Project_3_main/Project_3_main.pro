TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    Integration_loops.cpp \
    Laguerre.cpp \
    Legendre.cpp \
    Repulsion_function.cpp \
    main.cpp

HEADERS += \
    Integration_loops.h \
    Laguerre.h \
    Legendre.h \
    Repulsion_function.h
