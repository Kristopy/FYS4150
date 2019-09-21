TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        Function.cpp \
        main.cpp

INCLUDEPATH += /usr/local/include

LIBS += -L/usr/local/lib
LIBS += -larmadillo -llapack -lblas

DISTFILES +=

HEADERS += \
    Function.h
