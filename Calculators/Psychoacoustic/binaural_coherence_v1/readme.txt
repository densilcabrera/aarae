-------------------------------------------------------------------------------
Matlab source code for binaural coherence model
-------------------------------------------------------------------------------
Related paper:
M. Jeub, M. Dörbecker, P. Vary: "A Semi-Analytical Model for the Binaural
Coherence of Noise Fields", IEEE Signal Processing Letters, Vol.18, No.3,
March 2011 (DOI 10.1109/LSP.2011.2108284)
-------------------------------------------------------------------------------
Copyright (c) 2011, Marco Jeub
Institute of Communication Systems and Data Processing
RWTH Aachen University, Germany
Contact information: jeub@ind.rwth-aachen.de

All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

    * Redistributions of source code must retain the above copyright 
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright 
      notice, this list of conditions and the following disclaimer in 
      the documentation and/or other materials provided with the distribution
    * Neither the name of the RWTH Aachen University nor the names 
      of its contributors may be used to endorse or promote products derived 
      from this software without specific prior written permission.
      
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.

-------------------------------------------------------------------------------
Version 1.0
-------------------------------------------------------------------------------
The folder contains the following files:

binaural_coherence.m
	Main function to generate 3D and 2D noise fields taking into account
	head shadowing

binaural_coherence_example.m
	Example script for coherence calculation and the comparison to 
	free-field coherence models
	
Pre-calculated look-up tables with different settings are available on request	
-------------------------------------------------------------------------------
Tested under Windows 7 (64bit) and Matlab r2010b (64bit)
-------------------------------------------------------------------------------