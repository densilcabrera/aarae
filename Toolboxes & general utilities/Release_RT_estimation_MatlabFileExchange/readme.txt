-------------------------------------------------------------------------------
Matlab source code for blind reverberation time estimation from 
reverberant (denoised) speech signals.
-------------------------------------------------------------------------------

Reference:
Heinrich W. Löllmann, Emre Yilmaz, Marco Jeub and Peter Vary:
"An Improved Algorithm for Blind Reverberation Time Estimation"
International Workshop on Acoustic Echo and Noise Control (IWAENC), 
Tel Aviv, Israel, August 2010.

-------------------------------------------------------------------------------
Copyright (c) 2012, Heinrich Loellmann and Marco Jeub
Institute of Communication Systems and Data Processing
RWTH Aachen University, Germany
Contact information: loellmann@ind.rwth-aachen.de

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

The folder contains the following files and folders:

ML_RT_estimation_example.m
  Script to demonstrate the use of the function

function 
 Folter with functions for RT estimation

utilities
 Folder with utility functions for determing the RT by the Schroeder method and
 the reading of speech files

AIR
 Folder with binaural room impulse responses for the example taken from the
 AIR database
 Download link: http://www.ind.rwth-aachen.de/~air

speech_file
 Folder with speech file for the simulation example taken from the
	"TSP speech database"
	Reference: P. Kabal, "TSP speech database", Tech. Rep., Department of
	Electrical & Computer Engineering, McGill University, Montreal,
	Quebec, Canada, 2002.
	Download link: 	http://www-mmsp.ece.mcgill.ca/Documents/Data/index.html 
	
-------------------------------------------------------------------------------
Tested under Windows 7 (64bit) and Matlab r2011b (64bit)
-------------------------------------------------------------------------------