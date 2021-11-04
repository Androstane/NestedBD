# NestedBD

NestedBD is a BEAST 2 package for inference of evolutionary trees (topologies and branch lengths) from single-cell copy number profiles. 

### Installation 
We are working on making NestedBD availalbe through BEAST2 Package Manager. Currently, if you would like to use NestedBD it can be installed mannualy by building from the source file contained in this repository. 

#### Prerequisite 
Apache Ant nees to be installed.

#### Build from source
First, download and build the latest release of Beast2 Project: 

    git clone https://github.com/CompEvol/beast2.git
    cd beast2
    ant

Download and build NestedBD

    cd ../
    git clone https://github.com/Androstane/NestedBD.git
    cd NestedBD
    ant addon

Install NestedBD by moving the package to a place where BEAUti can find it

    cp -r release/add-on ~/.beast/2.4/NestedBD

