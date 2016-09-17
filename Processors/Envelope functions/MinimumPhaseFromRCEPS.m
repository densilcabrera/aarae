function out = MinimumPhaseFromRCEPS(in,method,thresholddB,timeshift,timereverse)
% This function returns a minimum phase reconstruction of the input audio
% using Matlab's real cepstrum (rceps) function, or using explicit code in
% AARAE's utility function minphasefreqdomain.

if nargin ==1 
    param = inputdlg({'Method: use Matlab''s RCEPS [0] or use more flexible method [1]';... 
                      'Dynamic range for cepstrum processing (if using more flexible method) in decibels';...
                      'Circular shift to zero crossing (if using more flexible method) [0 | 1]';...
                      'Time-reverse [0 | 1]'},...
                      'Window title',... 
                      [1 60],... 
                      {'1';'120';'1';'0'});
    param = str2num(char(param)); 

    if length(param) < 4, param = []; end 
    if ~isempty(param) 
        method = param(1);
        thresholddB = param(2);
        timeshift = param(3);
        timereverse = param(4);
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
end
    

[len, chans, bands,dim4,dim5,dim6] = size(in.audio);

out.audio = zeros(len,chans, bands,dim4,dim5,dim6);

if method == 0
    % Matlab's inbuilt minimum phase calculator
for ch = 1:chans
    for bnd = 1:bands
        for d4=1:dim4
            for d5 = 1:dim5
                for d6 = 1:dim6
        [~, out.audio(:,ch,bnd,d4,d5,d6)] = rceps(in.audio(:,ch,bnd,d4,d5,d6));
                end
            end
        end
    end
end
else
    for bnd = 1:bands
        for d4=1:dim4
            for d5 = 1:dim5
                for d6 = 1:dim6
                    out.audio(:,:,bnd,d4,d5,d6) =...
            ifft(minphasefreqdomain(abs(fft(in.audio(:,:,bnd,d4,d5,d6))),...
            thresholddB,timeshift));
                end
            end
        end
    end
end
if timereverse == 1
    out.audio = flip(out.audio,1);
end

out.funcallback.name = 'MinimumPhaseFromRCEPS.m';
out.funcallback.inarg = {method,thresholddB,timeshift,timereverse};