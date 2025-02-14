# NestedBD & NestedBD-Long

## NestedBD  
![Figure1](https://github.com/Androstane/NestedBD/assets/31413803/4baea497-4bfc-4743-8d14-dfd0d3e6205e)  

NestedBD is a BEAST 2 package for inference of evolutionary trees (topologies and branch lengths) from single-cell copy number profiles.  

## NestedBD-Long  
![Figure2](https://github.com/user-attachments/assets/bb54fa25-2f44-459d-b6ac-c303ab8bf8a5)

**NestedBD-Long** is an upgraded version of NestedBD that extends its functionality to handle longitudinal data.

We are working on making NestedBD and NestedBD-Long available through the BEAST2 Package Manager. Currently, if you would like to use NestedBD, it can be installed manually using `BDaddon.zip` as described in the **Installation** section. If you prefer to build the package from source files (which are also available in this repository), see the **Build the Package** section.

## Build the Package 
### Prerequisite 
Apache Ant needs to be installed.

### Build from source
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

NestedBD is implemented as a BEAST2 package, and requires beast2 installation to run.NestedBD is currently not compatible with beast 2.7. 

#### Installing BEAST2
Please follow the instruction on http://www.beast2.org/ and download the BEAST2 version for your operating system. Note that this is different to download and build source code of beast2 from GitHub repository. 

#### Installing the Pakcage 
Please follow the instructions on http://www.beast2.org/managing-packages/#Install_by_hand and use BD.addon.zip. 
NestedBD uses CoupledMCMC to speed up execution, the CoupledMCMC package for beast2 can be installed according to instruction at: https://github.com/nicfel/CoupledMCMC.
NestedBD-Long uses Sampled Ancestor to specify the tree prior, the Sampled Ancestor package for beast2 can be installed according to instructions at: https://github.com/CompEvol/sampled-ancestors/tree/v2.6

#### Test if the package is successfully installed
Please follow the instructions to run beast2 through the command line: https://www.beast2.org/2019/09/26/command-line-tricks.html, using scripts/BD.xml. For example, for linux, run:

        ~/beast/bin/beast -seed 1 test.xml
        
If all required packages are successfully installed, you should see ''1.trees'' and ''1.log'' (along with other files) at where you execute the command. 

### Usage
#### Prepare Input Data for NestedBD
NestedBD takes an integer copy number matrix, with each row representing a bin and each column representing a cell. An example of an accepted data matrix is available at ''scripts/test.cnp''.
After the data is properly formatted, you can then generate the XML file to be used as input by: 

        python xml_generator.py -data ${path to the data matrix file you prepared} -out ${path to generated xml file}

Please ensure ''text.xml'' is in the same folder as the ''xml_generator.py'' to avoid error. 

#### Prepare Input Data for NestedBD-Long

#### Running NestedBD or NestedBD-Long
The MCMC chain length can be manually set to the desired value, such as N, by setting ''chainLength="{N}" in the generated XML file. You may also want to change the frequency of the logger by setting logEvery="{M}" in the generated XML file accordingly. 

When the XML file is ready, again follow the instructions to run beast2 through the command line: https://www.beast2.org/2019/09/26/command-line-tricks.html. For example, with Linux and default installation:

        ~/beast/bin/beast {your_input}.xml


### Parameter and Prior specification

There's no command line parameter required to be specified except the input file location when running NestedBD. The prior distribution and parameter values used for the study are specificed within the xml file. 

    


