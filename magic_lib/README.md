# MAGIC_LIB
#### A Python Library for Region-Based Image Segmentation

This project is an extension to the MAp Guided Ice Classification (MAGIC) system developed by the Vision and Image Processing lab at the University of Waterloo. MAGIC was developed primarily for the purpose of segmentation and classification of sea ice in synthetic aperture radar imagery. 

This repository implements a user-friendly python interface to the core of the MAGIC system, namely:
- creation of region adjacency graphs (RAGs), which allow efficient region-based segmentation of large images
- access to the RAG attributes such as connectivity and region statistics through the python interface
- a wrapper for the multivariate IRGS algorithm
- various useful routines such as watershed segmentation, vector gradient computation, and GLCM texture feature extraction

Many of the workflows implemented in the MAGIC GUI are not implemented here, for example "glocal" segmentation and gluing. This is because I made this library cross-platform at the expense of stripping out all the Windows-specific .NET stuff in the GU. These routines should in principle be pretty easy to add to this library - if any VIP students feel like volunteering for this task I encourage you to do so :))))

## Installation

To build and install magic_lib you will need a python installation, a c++ compiler, and cmake. Instructions for building on Windows and Linux are given below. Support for MacOS is not currently included, but could probably be added pretty easily. 

### Windows Installation:

Make sure you have a microsoft visual C++ compiler installed on your system. Any recent version should work, but I have only tested the build with MSVC 2022. It is probably also possible to compile under MinGW but I haven't tried it. 

Download the cmake installer from [here](https://cmake.org/download/) and run it. When you are prompted, make sure to select 'Add cmake to PATH for the current user'. You will need to reboot your system for the PATH environment variable to be updated.

Clone this repository into a convenient location: 
`git clone --recurse-submodules https://github.com/Max-Manning/magic_lib`

Build the library and install it:
```
cd magic_lib
python setup.py install
```

### Linux Installation (tested on Ubuntu 20.04 LTS):

Install the dependencies:
`sudo apt install build-essential cmake`

Clone this repository into a convenient location: 
`git clone --recurse-submodules https://github.com/Max-Manning/magic_lib`

Build the library and install it:
```
cd magic_lib
python setup.py install
```

## Quick Start

`from magic_lib import magic_py`

## I've Found a Bug!!
Feel free to sumbit an issue or a pull request to this repository, or else harass me by email at `miamanni@uwaterloo.ca`

## Acknowledgements

The MAGIC system contains the contributions of many students and researchers who have worked with the VIP lab over the course of the past two decades. Sadly, records of the contributions of many of these individuals have been lost to the sands of time (and apparently a server hack at one point). An incomplete list of contributors is given below.
```
Max Manning
Mohsen Ghanbari
Major Jiang
Michael Stone
Fan Li
Steven Leigh
Zhijie Wang
Peter Yu
A. K. Quin
Quiyao Yu
David Clausi

```

### Note on Licensing
It would be nice to release this library as open source at some point, but there are a couple of licensing issues which prevent that in the current embodiment. In particular, `magic_core/Matrix/src/Eigen.cpp` and `magic_core/Src/Texture/matrix_gabor.cpp` contain code from numerical recipes (whose license forbids distribution of the source code). The fft routine in `magic_core/Src/Texture/SignalProcessing.cpp` is also clearly derived from numerical recipes in C++.
There may be other uses of copyrighted work scattered throughout the codebase, but if so they are unlabeled (and checking all of them would be a significant undertaking). 



