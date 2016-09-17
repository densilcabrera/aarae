function OUT = rotateHOA_aarae(IN, axisxyz, rotAng)
% This function rotates a HOA matrix around an axis (specified by
% Cartesian coordinates axisxyz) by an angle specified in degrees.
%
% The function calls a function within Nicolas Epain's HOAToolbox.

if nargin ==1 
    
    param = inputdlg({'Rotation axis theta (x,y,z)';... 
        'Angle to roate (deg)'},...
        'Rotate HOA matrix',... 
        [1 30],... 
        {'0,0,1';'90'}); % default values
    
    
    if length(param) < 2, param = []; end 
    if ~isempty(param) 
        % inputs to your function's input parameters.
        axisxyz = str2num(char(param(1)));
        rotAng = str2num(char(param(2)));
    else
        OUT = [];
        return
    end
else
    param = [];
end
if isstruct(IN) 
    audio = IN.audio; % Extract the audio data
    
    
    
elseif ~isempty(param) || nargin > 1    
    audio = IN;
end

if ~isempty(audio)
   
    [~,chans,bands,dim4,dim5,dim6]=size(audio);
    
if abs(chans^0.5 - round(chans^0.5)) >1e-20
    h=warndlg('This audio does not appear to be in HOA format. Unable to analyse with rippleplotfromHOA.','AARAE info','modal');
    uiwait(h)
    OUT = [];
    return
end

max_order=round(chans.^0.5-1);
hoaFmt = GenerateHoaFmt('res2d',max_order,'res3d',max_order);
    
    for b = 1:bands
        for d4 = 1:dim4
            for d5 = 1:dim5
                for d6 = 1:dim6
                    audio(:,:,b,d4,d5,d6) = audio(:,:,b,d4,d5,d6) * Hoa2HoaRotationMatrix(hoaFmt,axisxyz,rotAng*pi/180);
                end
            end
        end
    end
    
    if isstruct(IN)
        OUT = IN; 
        OUT.audio = audio; 
        OUT.funcallback.name = 'rotateHOA_aarae.m'; 
        OUT.funcallback.inarg = {axisxyz, rotAng}; % assign all of the 
        
    else
        
        OUT = audio;
    end
    
else
    
    OUT = [];
end

%**************************************************************************
% Copyright (c) 2014, Densil Cabrera
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%  * Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  * Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%  * Neither the name of the University of Sydney nor the names of its contributors
%    may be used to endorse or promote products derived from this software
%    without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%**************************************************************************