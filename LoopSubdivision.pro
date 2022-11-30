QT+=opengl
TEMPLATE = app
TARGET = LoopSubdivision
INCLUDEPATH += . 
CONFIG += c++17

# The following define makes your compiler warn you if you use any
# feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

linux-g++ | linux-g++-64 | linux-g++-32 {
    LIBPATH+= lib
    INCLUDEPATH+= include
}

# Input
HEADERS += src/ArcBall.h \
           src/ArcBallWidget.h \
           src/Cartesian3.h \
           src/Homogeneous4.h \
           src/DirectedEdgeSurface.h \
           src/Matrix4.h \
           src/Quaternion.h \
           src/RenderController.h \
           src/RenderParameters.h \
           src/RenderWindow.h \
           src/RenderWidget.h \
           src/RGBAImage.h \
           src/RGBAValue.h \
           src/SphereVertices.h

SOURCES += src/ArcBall.cpp \
           src/ArcBallWidget.cpp \
           src/Cartesian3.cpp \
           src/Homogeneous4.cpp \
           src/DirectedEdgeSurface.cpp \
           src/main.cpp \
           src/Matrix4.cpp \
           src/Quaternion.cpp \
           src/RenderController.cpp \
           src/RenderWindow.cpp \
           src/RenderWidget.cpp \
           src/RGBAImage.cpp \
           src/RGBAValue.cpp \
           src/SphereVertices.cpp
