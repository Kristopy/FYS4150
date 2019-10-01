TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        Jacobi.cpp \
        functions.cpp \
        initializematrix.cpp \
        main.cpp \
        test.cpp

INCLUDEPATH += /usr/local/include

LIBS += -L/usr/local/lib
LIBS += -larmadillo -llapack -lblas


HEADERS += \
    Jacobi.h \
    catch.hpp \
    functions.h \
    initializematrix.h

DISTFILES +=
QMAKE_CXXFLAGS += -O3
