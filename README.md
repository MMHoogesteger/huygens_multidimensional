# Multivariable Huygens'synchronization tools
Author: M.M. Hoogesteger
Eindhoven University of Technology, Mechanical Engineering, Dynamics and control
email address: m.m.hoogesteger@tue.nl
First version: November 2018
Current version: July 2022


This repo contains all software used in 
"Huygens' Synchronization with Multivariable Coupling"
DC 2018.088

An overview of what can be found in this folder is also given in Appendix C.

%% data
This folder is used to store simulation data among others. 
It initially only contains the camera parameters used in the image_app.

%% scripts
This folder contains all scripts that need to be ran to simulate, analyze, 
plot and so on.

%% lib
This folder contains many functions and classes that can be used by the user
and are necessary for many of the scripts in the scripts folder.
Can be added to the path by means of the prepEnv script.

%% models
Serves as the main source for models for the simulation. Generated with the
derive_models script.

%% Cameras
Contains the C++ program used to capture images.

