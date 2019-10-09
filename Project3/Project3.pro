TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        Gauss_laguerre.cpp \
        Various_functions.cpp \
        main.cpp

HEADERS += \
    Gauss_laguerre.h \
    Various_functions.h

INCLUDEPATH += /usr/local/include

LIBS += -L/usr/local/lib
LIBS += -larmadillo -llapack -lblas

DISTFILES +=

