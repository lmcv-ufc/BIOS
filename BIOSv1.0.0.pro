TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TARGET = bios

DESTDIR =
OBJECTS_DIR = obj/#.obj
MOC_DIR = obj/#.moc
RCC_DIR = obj/#.rcc
UI_DIR = obj/#.ui

# OpenMP flags =================================

DEFINES += _OMP_
QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp

#===============================================

SOURCES += \
#    prob/compriser.cpp \
#    prob/riserseg.cpp \
    prob/problem.cpp \
#    prob/ccr.cpp \
#    prob/reduce.cpp \
#    prob/meshgen2d.cpp \
#    prob/riser.cpp \
    prob/metaprob.cpp \
    prob/material.cpp \
#    prob/lamriser.cpp \
    prob/metaopt.cpp \
    prob/lam.cpp \
    prob/benchmark.cpp \
    prob/lamplt.cpp \
    ctrl/lamga.cpp \
    ctrl/modpso.cpp \
    ctrl/algparam.cpp \
    ctrl/lamalg.cpp \
    ctrl/sel.cpp \
    ctrl/stdga.cpp \
    ctrl/stdabc.cpp \
    ctrl/stdais.cpp \
    ctrl/stdpso.cpp \
    ctrl/optalg.cpp \
    main/utl.cpp \
    main/input.cpp \
    main/main.cpp \
    pen/penalty.cpp \
    #sur/rbfm.cpp \
    sur/rbf.cpp \
    #sur/svm.cpp \
    sol/optsolution.cpp \
    sol/group.cpp \
    sol/food.cpp \
    sol/particle.cpp \
    sol/individual.cpp \
    sol/sampsao.cpp \
    mlib/rk4.cpp \
    mlib/vec.cpp \
    mlib/mat.cpp \
    mlib/matvec.cpp \
    mlib/sysmat.cpp \
    ctrl/modlamnsgaII.cpp \
    ctrl/modnsgaII.cpp \
    sur/surr.cpp \
    sur/krg.cpp \
    sur/samp.cpp \
    ctrl/saokrg.cpp \
    prob/fgm.cpp \
    prob/fgmplt.cpp \
    ctrl/rs.cpp \
    #ctrl/kitayamasao.cpp \
    ctrl/saorbf.cpp \
    prob/probsurr.cpp \
    ctrl/sao.cpp \
    sur/problike.cpp \
    #ctrl/grego.cpp \
    ctrl/stdde.cpp

HEADERS += \
    #prob/ccr.h \
    #prob/compriser.h \
    #prob/meshgen2d.h \
    #prob/reduce.h \
    #prob/lamriser.h \
    prob/metaprob.h \
    prob/lam.h \
    prob/material.h \
    #prob/riserseg.h \
    prob/metaopt.h \
    #prob/riser.h \
    #prob/scr.h \
    prob/benchmark.h \
    prob/lamplt.h \
    prob/problem.h \
    ctrl/modpso.h \
    ctrl/lamalg.h \
    ctrl/lamga.h \
    ctrl/algparam.h \
    ctrl/stdabc.h \
    ctrl/stdais.h \
    ctrl/sel.h \
    ctrl/stdga.h \
    ctrl/stdpso.h \
    ctrl/optalg.h \
    main/gblvar.h \
    main/utl.h \
    main/gbldef.h \
    main/input.h \
    pen/penalty.h \
    #sur/rbfm.h \
    #sur/svm.h \
    sur/rbf.h \
    sol/group.h \
    sol/food.h \
    sol/individual.h \
    sol/particle.h \
    sol/optsolution.h \
    sol/sampsao.h \
    mlib/rk4.h \
    mlib/sysmat.h \
    mlib/mat.h \
    mlib/vec.h \
    mlib/matvec.h \
    ctrl/modlamnsgaII.h \
    ctrl/modnsgaII.h \
    sur/surr.h \
    sur/krg.h \
    sur/samp.h \
    ctrl/saokrg.h \
    prob/fgm.h \
    prob/fgmplt.h \
    ctrl/rs.h \
    #ctrl/kitayamasao.h \
    ctrl/saorbf.h \
    prob/probsurr.h \
    ctrl/sao.h \
    sur/problike.h \
    #ctrl/grego.h \
    ctrl/stdde.h

INCLUDEPATH += \
   prob/ \
   ctrl/ \
   sol/ \
   pen/ \
   mlib/ \
   sur/ \
   main/
