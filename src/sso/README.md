## CSM structure tensors

This repository contains computer programs written and used by 
[Xinming Wu](http://www.jsg.utexas.edu/wu/) 
for multiple types of structure tensors that are discussed in our paper 
[Directional structure tensors in estimating seismic structural and 
stratigraphic orientations](http://www.jsg.utexas.edu/wu/files/wu2017Dst.pdf).

If you find this work helpful in your research, please cite:

    @article{doi:wu2017directional,
        author = {Xinming Wu and Xavier Janson},
        title = {Directional structure tensors in estimating seismic structural 
        and stratigraphic orientations},
        journal = {Geophysical Journal International},
        volume = {210},
        number = {1},
        pages = {534-548},
        year = {2017},
        doi = {doi: 10.1093/gji/ggx194},
        URL = {https://academic.oup.com/gji/article/210/1/534/3805465},
    }

This software depends on that in the [Mines Java Toolkit
(JTK)](https://github.com/dhale/jtk/). If you want to do more than browse the
source code, you must first download and build the Mines JTK using
[Gradle](http://www.gradle.org). The build process for software in
this repository is the same.

### Summary

Here are brief descriptions of key components:

#### LocalStratigraphicFilter
Constructs directional structure tensors for estimating seismic structural 
and stratigraphic features

#### StructureTensorAttribute
Estimates seismic linearity (2D) and planarity (3D) attributes by 
using eigenvalues of structure tensors

#### ShapeCoherence
Calculates seismic coherence with structure-guied shaping regularization

---
Copyright (c) 2018, Xinming Wu. All rights reserved.
This software and accompanying materials are made available under the terms of
the [Common Public License - v1.0](http://www.eclipse.org/legal/cpl-v10.html),
which accompanies this distribution.
