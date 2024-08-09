AARAE - Audio and Acoustical Response Analysis Environment for MATLAB

AARAE is not a stand-alone program, but runs within MATLAB. It requires the
following toolboxes:
* Audio Toolbox (previously Audio System Toolbox)
* DSP System Toolbox
* Signal Processing Toolbox
* Statistics and Machine Learning Toolbox
Other toolboxes are also used by some parts of AARAE, including:
* Curve Fitting Toolbox


(Type ver into MATLAB's command line to find out what toolboxes you have
installed.)

In order to start AARAE, please make the AARAE folder your 'Current folder'
in the MATLAB path. You should find all the AARAE subfolders:

    - Analysers
    - Audio (should have REQUIRED_AUDIO within it)
    - Calculators
    - Documents
    - Framework
    - Generators
    - Log (AARAE will automatically create this folder when it runs)
    - Processors
    - Projects
    - Templates
    - Toolboxes & General Utilities
    - Utilities
    - Workflows

These folders do not need to be added to the MATLAB path in order for the
GUI interface to launch. We advise you not to add them to the MATLAB path 
to reduce the risk of function name conflicts. Instead, simply make the AARAE 
folder MATLAB’s current directory.

You should also find the aarae.m and aarae.fig files the AARAE
directory, along with the Licence for AARAE and this README.txt file. 
When you run AARAE, the file Settings.mat will be created.

Once you've made the AARAE folder your current folder in the MATLAB path,
in order to launch AARAE, please type in MATLAB's Command Window:

>> aarae

This command will launch the user interface for AARAE.

We hope you find this project useful and we would highly appreciate your
valuable feedback.

Reference:
Cabrera, D., Jimenez, D., & Martens, W. L. (2014, November).
Audio and Acoustical Response Analysis Environment (AARAE):
a tool to support education and research in acoustics.
In Proceedings of Internoise. Melbourne, Australia.
