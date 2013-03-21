TEMPLATE = app

CONFIG = qt release 

TARGET = crystal_to_box_gui

SOURCES += \
    main.cpp \
    mainwindow.cpp \
    crystal_to_box.c

HEADERS += \
    mainwindow.h \
    crystal_to_box.h

LIBS += -llapack -lblas





