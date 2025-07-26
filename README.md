# LADlib
hcana classes for the Hall C LAD experiments 

To build

    % cmake -B build -S . -DCMAKE_INSTALL_PREFIX=<install location>
    % cmake --build build -j4

Check that the library and classes are visible

    % hcana -l
    hcana [0] gSystem->Load("build/src/libLAD")
    (int) 0
    hcana [1] .class THcLADHodoPlane
    ===========================================================================
    class THcLADHodoPlane
    SIZE: 1736 FILE: THcLADHodoPlane.h LINE: 18
    Base classes: --------------------------------------------------------
    0x0        public THaSubDetector
       0x0        public THaDetectorBase
          0x0        public THaAnalysisObject
             0x0        public TNamed
                0x0        public TObject
    List of member variables --------------------------------------------------
    ...
Install the library

    % cmake --install build

Add the library to `LD_Library_Path`
eg.

    setenv LD_LIBRARY_PATH ${HOME}/hallc/software/LADlib/lib64:$LD_LIBRARY_PATH
