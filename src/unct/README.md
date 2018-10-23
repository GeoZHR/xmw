## 3D seismic image processing for unconformities

This repository contains computer programs written and used by 
[Xinming Wu](http://www.jsg.utexas.edu/wu/) 
for salt boundary interpretation that is discussed in our Geophysics paper 
[Fault salt bounary interpretation with optimal path picking]
(http://www.jsg.utexas.edu/wu/files/wu2018FastSaltBoundaryInterpretationWithOptimalPathPicking.pdf).

If you find this work helpful in your research, please cite:

    @article{doi:wu2015unconformities,
        author = {Xinming Wu and Dave Hale},
        title = {3D seismic image processing for unconformities},
        journal = {GEOPHYSICS},
        volume = {80},
        number = {2},
        pages = {IM35-IM44},
        year = {2015},
        doi = {https://doi.org/10.1190/geo2014-0323.1},
        URL = {https://library.seg.org/doi/10.1190/geo2014-0323.1},
    }

This software depends on that in the [Mines Java Toolkit
(JTK)](https://github.com/dhale/jtk/). If you want to do more than browse the
source code, you must first download and build the Mines JTK using
[Gradle](http://www.gradle.org). The build process for software in
this repository is the same.

Like the Mines JTK, this is a toolkit for computer programmers. It is not a
complete system for seismic interpretation of geologic horizons. Others
(including commercial software companies) have built such systems using
earlier versions of one or more of the software tools provided in this
repository.

### Summary

Here are brief descriptions of key components:

#### LocalOrientFilterUM
Estimate seismic normal vectors with vertical caucal and anti-caucal filters

#### UncSurfer
Calculates unconformity likelihoods from seismic normal vectors that are 
estimated by using caucal and anti-caucal filters

Extracts unconformity surfaces from the unconformity likelihoods


---
Copyright (c) 2018, Xinming Wu. All rights reserved.
This software and accompanying materials are made available under the terms of
the [Common Public License - v1.0](http://www.eclipse.org/legal/cpl-v10.html),
which accompanies this distribution.
