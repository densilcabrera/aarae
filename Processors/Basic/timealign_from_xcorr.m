function [OUT,varargout] = timealign_from_xcorr(IN,dims,maxlag,lin_or_circ,cthreshold,wavenv)
% This function aligns the columns of audio in time based on the peak of
% their cross-correlation function with the average of all or a selection
% of columns.
%
% The dims input argument identifies dimensions that are not included in
% the averaging that is used to create the reference waveform. These
% dimensions are therefore time-aligned independently from other
% dimensions. Typically, dimension 3 (which is used in AARAE for bands) is
% not well-suited to cross-correlation-based time-alignment, and so is
% better included in the dims argument.
%
% If all of the non-singleton dimensions are listed in dims (in the most
% extreme case: 2,3,4,5,6), then there would be nothing to do. So in this
% circumstance the last column of audio in each dimension (apart from 3) is
% used as the refernce (i.e., no averaging is done). The last is used
% instead of the first because of the possibility that silent cycles are
% present in the audio (and these are first in a dimension 4 stack).
%
% The maxlag input is the maximum lag considered for time-alignment (from
% Matlab's xcorr function) - i.e. the possible time shift is from -maxlag
% to +maxlag.
%
% If linear shifting is selected, then the audio is zeropadded by maxlag at
% the start and end. If circular shifting is used, this zeropadding is not
% done. In either case Matlab's circshift function is used for time
% alignment.
%
% A correlation coefficient threshold can be assigned, such that no time
% shifting is done for audio that yields a maximum correlation coefficient
% below the threshold. This should be a value between 0 and 1 if positive
% correlations only are used for time-alignment. If it is set to 0, then
% there is no threshold. If it is set to 1 then there is no point running
% this function because no time-alignment will be done. If it is negative
% (between -1 and almost 0) then time algnment is done based on the maximum
% absolute value correlation coefficent (i.e., the sign of the correlation
% coefficient is ignored).
%
% If phase reversals between channels exist, then it might not make sense
% to use averaging in generating the reference, and so listing all
% non-singleton dimensions in dims is likely to be the best approach.
%
% However, an alternative approach is to use the envelope cross-correlation
% lag for time-alignment, rather than the waveform itself. If wavenv > 0,
% then this is done. The value of wavenv is the length of a zero phase
% smoothing filter applied to the envelope. Using as well as dealing with
% phase reversal issues, cross-correlation of envelopes is most likely to
% be useful if you are trying to time-align bands (dimension 3) with each
% other.
%
% Code by Densil Cabrera
% version 1.00 (13 September 2014)


if isstruct(IN)
    audio = IN.audio;
else
    audio = IN;
end


if ~exist('wavenv','var')
    wavenv = 0;
end


if ~exist('cthreshold','var')
    cthreshold = 0.5;
end


if ~exist('lin_or_circ','var')
    lin_or_circ = 0;
end


if ~exist('maxlag','var')
    maxlag = 500;
end


if ~exist('dims','var')
    dims = 3;
end



% input dialog
if nargin == 1
    prompt = {'List any dimensions to align independently: channels [2], bands [3], cycles [4], output channel sequence [5]';...
        'Maximum lag in samples';...
        'Circular [0] or linear [1] shift';...
        'Correlation coefficient threshold [between 0 and 1 for time-alignment with positive correlations, or between -1 and 0 for time-alignment with absolute correlations';...
        'Xcorr using the wave [0] or using the envelope with a smoothing filter [specify length in samples]'};
    dlg_title = 'Time Alignment Settings';
    num_lines = 1;
    def = {num2str(dims);...
        num2str(maxlag);...
        num2str(lin_or_circ);...
        num2str(cthreshold)...
        ;num2str(wavenv)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        OUT = [];
        return
    else
        dims = str2num(answer{1});
        maxlag = round(str2double(answer{2}));
        lin_or_circ = str2double(answer{3});
        cthreshold = str2double(answer{4});
        wavenv = round(str2double(answer{5}));
    end
end



if ~isempty(audio) && ~isempty(dims)...
        && ~isempty(maxlag) && ~isempty(lin_or_circ)...
        && ~isempty(cthreshold)
    
    [~,chans,bands,dim4,dim5,dim6] = size(audio);
    
    if chans+bands+dim4+dim5+dim6 == 5
        % audio is only one column, so there is nothing to do
        OUT = [];
        return
    end
    
    if wavenv >= 1
        [refaudio, env] = deal(abs(audio));
    else
        refaudio = audio;
    end
    
    for n = 2:6
        if isempty(find(dims == n, 1))
            refaudio = mean(refaudio,n);
        end
    end
    [~,refd2,refd3,refd4,refd5,refd6] = size(refaudio);
    
    
    
    if chans == refd2 && bands == refd3 && dim4 == refd4 ...
            && dim5 == refd5 && dim6 == refd6
        % use the last column instead
        refaudio = audio(:,end,:,end,end,end);
        [~,refd2,refd3,refd4,refd5,refd6] = size(refaudio);
    end
    
    
    % zero pad if linear shifting
    if lin_or_circ == 1
        audio = [zeros(maxlag,chans,bands,dim4,dim5,dim6);...
            audio;zeros(maxlag,chans,bands,dim4,dim5,dim6)];
        refaudio = [zeros(maxlag,refd2,refd3,refd4,refd5,refd6);...
            refaudio;...
            zeros(maxlag,refd2,refd3,refd4,refd5,refd6)];
    end
    
    % envelope smoothing filter (zero phase)
    if wavenv > 1
        b = hann(wavenv+2);
        b= b(2:end-1);
        b = b./sum(b);
        env = filtfilt(b,1,env);
        refaudio = filtfilt(b,1,env);
    end
    
    
    
    % GET LAGS AND SHIFT
    % an alternative to the following would be cross-spectrum (which would
    % not need for loops, but could blow out the calculation size if the
    % audio was large). Consider changing if this is too slow...
    [maxc,lags] = deal(zeros(chans,bands,dim4,dim5,dim6));
    for ch = 1:chans
        for b = 1:bands
            for d4 = 1:dim4
                for d5 = 1:dim5
                    for d6 = 1:dim6
                        if wavenv < 1
                            [c,lag] = xcorr(audio(:,ch,b,d4,d5,d6),...
                                refaudio(:,min([ch,refd2]),...
                                min([b,refd3]),min([d4,refd4]),...
                                min([d5,refd5]),min([d6,refd6])),...
                                maxlag, 'coeff');
                        else
                            [c,lag] = xcorr(env(:,ch,b,d4,d5,d6),...
                                refaudio(:,min([ch,refd2]),...
                                min([b,refd3]),min([d4,refd4]),...
                                min([d5,refd5]),min([d6,refd6])),...
                                maxlag, 'coeff');
                        end
                        if cthreshold >= 0
                            maxc(ch,b,d4,d5,d6) = max(c);
                            if maxc(ch,b,d4,d5,d6)>=cthreshold
                                lags(ch,b,d4,d5,d6)=lag(c==maxc(ch,b,d4,d5,d6));
                                audio(:,ch,b,d4,d5,d6) =...
                                    circshift(audio(:,ch,b,d4,d5,d6),...
                                    -lags(ch,b,d4,d5,d6));
                            end
                        else
                            maxc(ch,b,d4,d5,d6) = max(abs(c));
                            if maxc(ch,b,d4,d5,d6)>=cthreshold
                                lags(ch,b,d4,d5,d6)=lag(abs(c)==maxc(ch,b,d4,d5,d6));
                                audio(:,ch,b,d4,d5,d6) =...
                                    circshift(audio(:,ch,b,d4,d5,d6),...
                                    -lags(ch,b,d4,d5,d6));
                            end
                        end
                    end
                end
            end
        end
    end
    
    
else
    OUT = [];
    return
end




if isstruct(IN)
    OUT = IN;
    OUT.audio = audio;
    OUT.properties.maxc = maxc;
    OUT.properties.lags = lags;
    OUT.funcallback.name = 'timealign_from_xcorr.m';
    OUT.funcallback.inarg = {dims,maxlag,lin_or_circ,cthreshold,wavenv};
else
    OUT = audio;
    varargout{1} = maxc;
    varargout{2} = lags;
    
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