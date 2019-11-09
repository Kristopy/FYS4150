TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        inputoutput.cpp \
        lib.cpp \
        main.cpp \
        metropolis.cpp

HEADERS += \
    inputoutput.h \
    lib.h \
    metropolis.h

INCLUDEPATH += /usr/local/include

LIBS += -L/usr/local/lib
LIBS += -larmadillo -llapack -lblas
