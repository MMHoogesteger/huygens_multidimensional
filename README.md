# Multivariable Huygens' synchronization tools

[![DOI](https://zenodo.org/badge/525821288.svg)](https://zenodo.org/badge/latestdoi/525821288)

Author: M.M. Hoogesteger

Eindhoven University of Technology, Mechanical Engineering, Dynamics and control

email address: m.m.hoogesteger@tue.nl

First version: November 2018

Current version: July 2022

## Description

This repo contains all software used in the master thesis of the author:
"Huygens' Synchronization with Multivariable Coupling" 

An overview of what can be found in this repo is also given in Appendix C of this thesis.

## Folder structure
`data`

This folder is used to store simulation data among others. 
It initially only contains the camera parameters used in the image_app.

`scripts`

This folder contains all scripts that need to be ran to simulate, analyze, 
plot and so on.

`lib`

This folder contains many functions and classes that can be used by the user
and are necessary for many of the scripts in the scripts folder.
Can be added to the path by means of the prepEnv script.

`models`

Serves as the main source for models for the simulation. Generated with the
derive_models script.

## Not included
`lib/matlab2tikz`

https://www.mathworks.com/matlabcentral/fileexchange/22022-matlab2tikz-matlab2tikz

`lib/mtimesx`

https://www.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support?s_tid=srchtitle

`Cameras`

Software to obtain images from the two PS Eye cameras: 
https://codelaboratories.com/downloads/