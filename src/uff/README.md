## Moving faults while unfaulting 3D seismic images

This repository contains computer programs written and used by 
[Xinming Wu](http://www.jsg.utexas.edu/wu/) 
for seismic unfaulting that is discussed in our Geophysics paper 
[Moving faults while unfaulting 3D seismic images]
(http://www.jsg.utexas.edu/wu/files/wu2016MovingFaultsWhileUnfaulting3dSeismicImages.pdf).

If you find this work helpful in your research, please cite:

    @article{doi:wu2016moving,
        author = {Xinming Wu and Simon Luo and Dave Hale},
        title = {Moving faults while unfaulting 3D seismic images},
        journal = {GEOPHYSICS},
        volume = {81},
        number = {2},
        pages = {IM25-IM33},
        year = {2016},
        doi = {10.1190/GEO2015-0381.1},
        URL = {https://library.seg.org/doi/abs/10.1190/geo2015-0381.1},
    }

This software depends on that in the [Mines Java Toolkit
(JTK)](https://github.com/dhale/jtk/). If you want to do more than browse the
source code, you must first download and build the Mines JTK using
[Gradle](http://www.gradle.org). The build process for software in
this repository is the same.

### Summary

Here are brief descriptions of key components:

#### FaultSlipConstraints
Constructs unfaulting equations from pre-estimated fault slip vectors

#### UnfaultS
Estimates 3 components of unfaulting shifts by solving the unfaulting 
equations with smoothing regularization.

---
Copyright (c) 2018, Xinming Wu. All rights reserved.
This software and accompanying materials are made available under the terms of
the [Common Public License - v1.0](http://www.eclipse.org/legal/cpl-v10.html),
which accompanies this distribution.
