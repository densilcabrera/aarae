function OUT = rectroommodesynthesiser(L,s,spattern,r,rpattern,alpha,c,maxf,maxftable,dur,fs,scaling)
% This function calculates the room mode frequencies, decay rates,
% amplitude and phase offsets for a set of source and receiver positions in
% a rectangular room. A kind of impulse response (IR) is generated from this,
% and a table of data is created. The pseudo-IR is created by additive
% synthesis of exponentially decaying tones, corresponding to the calculated
% frequencies, amplitudes, phase offsets and decay rates. A high cutoff
% frequency limits the calculation (this could potentially take a very long
% time if the cutoff frequency is set high, especially if the room volume
% is large).
%
% This function can produce pseudo-IRs for multiple inputs (sources) and/or 
% multiple outputs (receivers). For consistency with AARAE audio recording,
% multiple outputs are written as multiple channels (dimension 2), and 
% multiple inputs are written in dimension 5. (Dimension 5 is used because
% dimension 3 is used by AARAE for bands, and dimension 4 is used by AARAE
% for multicycle signal measurements. This is the same as the dimension
% assignment used in AARAE's audio recording GUI for asynchronous output
% channels.)
%
% INPUTS
% L is the length, width and height of the room, in metres, as a row
% vector, for example [6.34 5.1 4]. If a single value is input, then the
% room is taken as cubic.
% 
% s is the source coordinates (in Cartesian form), or the number of source
% coordinates to generate. The coordinate system has the room origin at a
% corner, with positive values only within the room. If a single source is
% specified, then input a 3-valued row vector, e.g. [1 2.3 1.2]. If you
% wish to specify more than one source, then use a matrix with one
% source's coordinates on each row. If you wish source coordinates to be
% automatically generated, then use a single number to specify the number
% of sources to generate. The generation method is controlled by the
% parameter spattern.
%
% r is re receiver coordinates (in Cartesian form), or the number of
% receiver coordinates to generate. Its format and usage is essentially the
% same as for s (described above).
%
% alpha is the sound absorption coefficient of the surfaces in each of the
% three dimensions e.g. [0.3 0.2 0.25]. If a single value is input, then it
% is used for all six surfaces.
%
% c is the speed of sound in m/s.
%
% maxf is the maximum frequency to synthesise. Note that the modal density
% increases dramatically with frequency, so it is dangerous to set this
% value high (please experiment with intermediate values to get a sense of
% the computational resources required).
%
% maxftable is the maximum frequency to include in the output table of mode
% frequency and individual mode reverberation time.
%
% dur is the duration of the audio that is generated in seconds.
%
% fs is the audio sampling rate in Hz.
%
% scaling is the method for amplitude scaling:
%   0: no scaling (sum all of the individual modes)
%   1: normalize
%   2: divide by the number of synthesised modes (i.e., average)







if nargin == 0
    L = [6.34 5.1 4];
    s = [0 0 0];
    r = [0 0 0];
    alpha = [0.2 0.2 0.2];
    c = 343;
    dur = 2;
    fs = 48000;
    maxf = 500;
    %dialog box for settings
    prompt = {'Room length, width and height (m)';...
              'Source coordiates (x,y,z) or number of sources'; ...
              'Source distribution method (if specifying number of sources): 0 random distribution if first octant';...
              'Receiver coordiates (x,y,z) or number of receivers'; ...
              'Receiver distribution method (if specifying number of receivers): 0 random distribution if first octant';...
              'Absorption coefficient (x,y,z), or a single value';...
              'Speed of sound (m/s)';...
              'Maximum frequency to synthesise (Hz)'; ...
              'Maximum frequency in table (Hz)'; ...
              'Duration of the impulse response (s)';...
              'Sampling rate (Hz)';...
              'Scaling method: Add without scaling (0), Normalize (1), Average (2)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {num2str(L);...
        num2str(s); '0';...
        num2str(r); '0';...
        num2str(alpha);...
        num2str(c);...
        num2str(maxf); '100';...
        num2str(dur);...
        num2str(fs); '1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(answer)
        L = str2num(answer{1,1});
        s = str2num(answer{2,1});
        spattern = str2num(answer{3,1});
        r = str2num(answer{4,1});
        rpattern = str2num(answer{5,1});
        alpha = str2num(answer{6,1});
        c = str2num(answer{7,1});
        maxf = str2num(answer{8,1});
        maxftable = str2num(answer{9,1});
        dur = str2num(answer{10,1});
        fs = str2num(answer{11,1});
        scaling = str2num(answer{12,1});
    else
        OUT = [];
        return
    end
end

if length(L) ~= 3
    L = [L(1) L(1) L(1)];
end

if length(alpha) ~= 3
    alpha = [alpha(1) alpha(1) alpha(1)];
end


maxnx = ceil(2*L(1)*maxf/c);
maxny = ceil(2*L(2)*maxf/c);
maxnz = ceil(2*L(3)*maxf/c);
aprime = -log(1-alpha); % absorption exponent
len = ceil(dur*fs);



if size(s,2) ~= 3
    switch spattern
        otherwise % random distribution in first octant of the room
            % used for case 0
            dim5 = s(1);
            s = rand(s(1),3);
            s(:,1) = s(:,1) * L(1)/2;
            s(:,2) = s(:,2) * L(2)/2;
            s(:,3) = s(:,3) * L(3)/2;
    end
else
    dim5 = size(s,1);
end

if size(r,2) ~= 3
    switch rpattern
        otherwise % random distribution in first octant of the room
            % used for case 0
            chans = r(1);
            r = rand(r(1),3);
            r(:,1) = r(:,1) * L(1)/2;
            r(:,2) = r(:,2) * L(2)/2;
            r(:,3) = r(:,3) * L(3)/2;
    end
else
    chans = size(r,1);
end

t = repmat(((1:len)'-1)./fs,[1,chans,1,1,dim5]);
OUT.fs = fs;
OUT.audio = zeros(len,chans,1,1,dim5);

tablelengthestimate = ceil(4*pi*L(1)*L(2)*L(3)/3 * (maxftable/c).^3);
tablelengthestimate = ceil(tablelengthestimate * 1.5); % increase the estimate
% table columns: nx ny nz f T
tabledata = zeros(tablelengthestimate,5);

fcount = 0;
tablecount = 0;
for nx = (1:maxnx+1)-1
    if nx > 0,jx = L(1)/nx; else jx = 0;end
    for ny = (1:maxny+1)-1
        if ny > 0,jy = L(2)/ny; else jy = 0;end
        for nz = (1:maxnz+1)-1
            if nx+ny+nx > 0
                if nz > 0,jz = L(3)/nz; else jz = 0;end
                f =c/2*((nx/L(1))^2+(ny/L(2))^2+(nz/L(3))^2)^0.5;
                if f <= maxf
                    fcount = fcount+1;
                    % cosine of reflection incidence angle
                    % for each dimension
                    denominator = (jx^2+jy^2+jz^2)^0.5;
                    costhetax = jx/denominator;
                    costhetay = jy/denominator;
                    costhetaz = jz/denominator;
                    % caclulate damping factor
                    delta = c/4 *(costhetax/L(1)*2*aprime(1) +...
                        costhetay/L(2)*2*aprime(2) +...
                        costhetaz/L(3)*2*aprime(3));
                    T = 3*log(10)./delta;
                    %tau = 2*T / log(1e6); % decay constant = 1/delta
                    % amplitude specific to each source-receiver combination
                    A = repmat(cos(nx*pi*s(:,1)/L(1)).*...
                        cos(ny*pi*s(:,2)/L(2)).*...
                        cos(nz*pi*s(:,3)/L(3)),[1,chans]).*...
                        repmat(cos(nx*pi*r(:,1)'/L(1)).*...
                        cos(ny*pi*r(:,2)'/L(2)).*...
                        cos(nz*pi*r(:,3)'/L(3)),[dim5,1]);
                    % phase difference specific to each source-receiver combination
                    phi = 2*pi*c/f*(costhetax * (repmat(s(:,1),[1,chans]) - repmat(r(:,1)',[dim5,1]))+...
                    costhetay * (repmat(s(:,2),[1,chans]) - repmat(r(:,2)',[dim5,1]))+...
                    costhetaz * (repmat(s(:,3),[1,chans]) - repmat(r(:,3)',[dim5,1])));
                
                    A1 = repmat(permute(A,[5,2,3,4,1]),[len,1,1,1,1]);
                    phi1 = repmat(permute(phi,[5,2,3,4,1]),[len,1,1,1,1]);
                    modaldensityestimate = 4*pi*L(1)*L(2)*L(3)*f.^2/c.^3; % this should probably be done a different way
                    OUT.audio = OUT.audio + A1.*sin(2*pi*f.*t + phi1)...
                        ./ (modaldensityestimate.^0.5 * exp(t.*delta));
                    
                    if f <=maxftable
                        tablecount = tablecount+1;
                        tabledata(tablecount,:) = [nx, ny, nz, f, T];
                    end
                end
            end
        end
    end
end

% Table
tabledata = tabledata(1:tablecount,:);
[~,row] = sort(tabledata(:,4),1,'ascend');
for col = 1:size(tabledata,2) % a very inelegant way of sorting!
    tabledata(:,col) = tabledata(row,col);
end
fig1 = figure('Name','Rectangular Room Modes');
table1 = uitable('Data',tabledata,...
    'ColumnName',{'nx','ny','nz','f (Hz)','T (s)'},...
    'RowName',{num2str((1:tablecount)')});
[~,table] = disptables(fig1,table1);
OUT.tables = table;



switch scaling
    case 1 % normalize
        OUT.audio = OUT.audio ./ max(max(max(abs(OUT.audio))));
    case 2 % average
        OUT.audio = OUT.audio ./ fcount;
end


OUT.chanID = makechanID(chans,4,r);
OUT.properties.dim5ID = makechanID(dim5,4,s); % dim5ID is not established in aarae at the time of writing this function
OUT.properties.roomdimensions = L;
OUT.properties.absorptioncoef = alpha;
OUT.properties.speedofsound = c;
OUT.properties.maximumfrequency = maxf;
OUT.funcallback.name = 'rectroommodesynthesiser.m';
OUT.funcallback.inarg = {L,s,spattern,r,rpattern,alpha,c,maxf,maxftable,dur,fs,scaling};

%**************************************************************************
% Copyright (c) 2015, Densil Cabrera
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
%  * Neither the name of The University of Sydney nor the names of its contributors
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