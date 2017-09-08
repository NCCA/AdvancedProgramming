TARGET=PathTracer

#QMAKE_CXXFLAGS += -O3
QMAKE_CXXFLAGS += -O2
#DEFINES+= USEFRAMEBUFFER
SOURCES+=$$PWD/src/main.cpp \
         $$PWD/src/Util.cpp \
         $$PWD/src/Image.cpp \
         $$PWD/src/Framebuffer.cpp

HEADERS+=$$PWD/include/Vec.h \
         $$PWD/include/Ray.h \
         $$PWD/include/Sphere.h \
         $$PWD/include/Material.h \
         $$PWD/include/Util.h \
         $$PWD/include/Image.h \
         $$PWD/include/Framebuffer.h
OBJECTS_DIR=$$PWD/obj
INCLUDEPATH +=$$PWD/include
CONFIG+=c++17

macx:CONFIG -=app_bundle



# core Qt Libs to use add more here if needed.
QT+=gui opengl core
# as I want to support 4.8 and 5 this will set a flag for some of the mac stuff
# mainly in the types.h file for the setMacVisual which is native in Qt5
isEqual(QT_MAJOR_VERSION, 5) {
  cache()
  DEFINES +=QT5BUILD
}
LIBS+=-L/usr/local/lib -lglfw3
macx:INCLUDEPATH+=/usr/local/include
LIBS += -L/usr/local/lib
macx:LIBS+= -framework OpenGL -framework IOKit -framework Cocoa -framework CoreVideo
linux:LIBS+= -lGL -lX11 -lXxf86vm -lXrandr -lXi -ldl -lXinerama -lXcursor -lglfw3 -lGLEW -lpthread
macx:QMAKE_MAC_SDK = macosx10.12

# where our exe is going to live (root of project)
DESTDIR=./
