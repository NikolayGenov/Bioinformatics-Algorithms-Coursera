TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += C++11
SOURCES += main.cpp


QMAKE_LFLAGS += -lgmp  -lcrypto -L/usr/lib/x86_64-linux-gnu -lgmp
