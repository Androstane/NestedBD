# NestedBD

NestedBD is a BEAST 2 package for inference of evolutionary trees (topologies and branch lengths) from single-cell copy number profiles. 

### Installation 
We are working on making NestedBD availalbe through BEAST2 Package Manager. Currently, if you would like to use NestedBD it can be installed mannualy by building from the source file contained in this repository. 

#### Prerequisite 
Apache Ant needs to be installed.

#### Build from source
First, download and build the latest release of Beast2 Project: 

    git clone https://github.com/CompEvol/beast2.git
    cd beast2
    ant

Download and build NestedBD. Make sure the folder beast2 and BD are under the same directory when running ant addon. 

    cd ..
    git clone https://github.com/Androstane/NestedBD.git
    cp -r NestedBD/BD/ BD/
    cd BD/
    ant addon

If build is successful, you should be able to find the zip file named BD.addon.zip contain the package under /beast2/build/dist/. You can use the following command to check whether such file exist:

    cd ..
    test -f beast2/build/dist/BD.addon.zip && echo "$FILE exists."

### Usage

NestedBD is implemented as a BEAST2 package, and requires beast2 installation to run.

#### Installing BEAST2
Please follow the instruction on http://www.beast2.org/ and download the BEAST2 version for your operating system

Navigate to where BEAST2 is installed, and use the following command to test whether Beast2 and NestedBD package have been properly installed:
    ~/.beast/bin/beast sample.xml > sample.out
