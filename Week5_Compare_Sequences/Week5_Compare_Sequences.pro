TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += C++11

SOURCES += main.cpp \
    Graph.cpp

HEADERS += \
    Graph.h \
    BLOSUM62.h \
    PAM250.h

