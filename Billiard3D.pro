TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -O2

SOURCES += \
        Billiard.cpp \
        NumericalIntegrator.cpp \
        main.cpp

DISTFILES += \
    Billiard3D.pro.user \
    TODO.txt

HEADERS += \
    Billiard.h \
    Matrix.h \
    NumericalIntegrator.h \
    Vector.h
