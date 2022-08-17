Author: M.M. Hoogesteger
Eindhoven University of Technology, Mechanical Engineering, Dynamics and control
email address: m.m.hoogesteger@tue.nl
November 2018

Many of the scripts are specific cases and only serve as inspiration for 
further use, contrary to the lib folder that contains files that have been 
kept as general as possible.

Note that the models, jobs and temp folder can all be empty and are used to
store generated models, simulation settings and temporary function files.


A few important files are:

prepEnv, which loads all lib folders and thereby prepares the environment for
the other scripts.

derive_models, which derives all models and stores them in the local models
folder. Manually copy them to the global models folder afterwards.

launch_image_app, opens the video processing application.

jobrunner, automatically runs simulations of all available simSettings in
the jobs folder. See the available examples and the code for further details.

Other files include:

describing_***, contains code that evaluates the describing functions for
various settings.

gatherFourier_***, loads all simulations of a certain set and gathers 
certain results.

genSimSettings_***, contains code that generates simulation settings with
one or more varying parameters.

load***, loads certain simulations.

plot***, plots certain simulations.

testEscOr, generates the data as found in Appendix H.



