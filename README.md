# NestedBD

NestedBD is a BEAST 2 package for inference of evolutionary trees (topologies and branch lengths) from single-cell copy number profiles. 

### Build the Package 
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

### Installation

NestedBD is implemented as a BEAST2 package, and requires beast2 installation to run.

#### Installing BEAST2
Please follow the instruction on http://www.beast2.org/ and download the BEAST2 version for your operating system. Note that this is different to download and build source code of beast2 from GitHub repository. 

#### Move package to where BEAUti can find it
Then, to allow BEAST utilize the NestedBD package, we need to locate where BEAUti look for packages and move BD.addon.zip to that folder.

By Default, BEAUti would look for package at ''$HOME/.beast2/''. Note that if you have not yet run BEAUti, you may not yet have a $HOME/.beast2 directory.
or,
If you have BEAUti UI installed (look for executable file of BEAUti at where you install BEAST). Then you can check the path by:
- open BEAUti
- Open the Package Manager through File > Manage Packages. 
- Click the question mark on right corner.
 You should see the path at which BEAUti looks for packages as in this screenshot:  
![image](https://user-images.githubusercontent.com/31413803/160720018-4c5707c8-1ae2-479f-bb61-ee036f4420a2.png)


### Usage

Navigate to where BEAST2 is installed, and use the following command to test whether Beast2 and NestedBD package have been properly installed:
    ~/.beast/bin/beast sample.xml > sample.out
    


