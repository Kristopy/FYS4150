TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        backward_euler.cpp \
        main.cpp \
        sym_tridiag_solver.cpp

HEADERS += \
    backward_euler.h \
    sym_tridiag_solver.h

INCLUDEPATH += /usr/local/include

LIBS += -L/usr/local/lib
LIBS += -larmadillo -llapack -lblas
