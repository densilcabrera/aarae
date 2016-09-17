function [OUT, varargout] = reduce_dims(IN, domain, operation, percent,dim,fs)
% This function reduces one of the dimensions of multidimensional audio 
% to a singleton dimension by applying a simple selection,
% arithmetic or statistical operation. Operations include:
% * just use the top layer of chosen dimension
% * mean of chosen dimension
% * trim-mean of chosen dimension
% * median of chosen dimension
% * maximum of the chosen dimension (probably most useful for complex data)
% * minimum of the chosen dimension (ditto)
% * percentile of the chosen dimension
%
% These operations can be done in the time domain, frequency domain,
% quefrency domain, Hilbert (analytic) time domain, or an intermediate 
% domain between time and frequency via the fractional Fourier transform
% (with a 0.5 fraction).


if isstruct(IN)
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    
    
elseif ~isempty(param) || nargin > 1
    audio = IN;
end

[~,chans,bands,d4,d5,d6] = size(audio);


if nargin ==1
    if d6  > 1
        dim = 6;
    elseif d5 > 1
        dim = 5;
    elseif d4 > 1
        dim = 4;
    elseif chans > 1
        dim = 2;
    elseif bands > 1
        dim = 3;
    else
        h = warndlg('Audio data is already 1-dimensional - this processor is for multidimensional data.','AARAE info','modal');
        uiwait(h)
        OUT = [];
        return
    end
        
    
    param = inputdlg({'Domain: time [0], frequency [1], quefrency [2], Hilbert [3], Fractional Fourier [4]';...
        'Top layer [0], Mean [1], Trim-mean [2], Median [3], Maximum [4], Minimum [5], Percentile [6], Sum [7], Mean difference [8]';...
        'Trim-mean or percentile percent value [0:100]';...
        'Dimension'},...
        'Window title',...
        [1 60],...
        {'0';'2';'50';num2str(dim)});
    
    param = str2num(char(param));
    
    if length(param) < 4, param = []; end
    if ~isempty(param)
        domain = param(1);
        operation = param(2);
        percent = param(3);
        dim = param(4);
    end
else
    param = [];
end




if ~isempty(audio) && ~isempty(fs)
    
    
    dim = abs(round(dim));
    if dim <=1
        warndlg('Dimension to operate on must be higher than 1','AARAE info','modal')
        OUT = [];
        return
    end
    
    if operation == 0
        % take only the first layer of the stack - discard the rest.
        % no need to consider domains in this case
        switch dim
            case 2
                audio = audio(:,1,:,:,:,:);
            case 3
                audio = audio(:,:,1,:,:,:);
            case 4
                audio = audio(:,:,:,1,:,:);
            case 5
                audio = audio(:,:,:,:,1,:);
            case 6
                audio = audio(:,:,:,:,:,1);
        end
    else
        
        % Transfrom to appropriate domain if necessary
        if domain == 1
            % frequency domain
            audio = fft(audio);
        elseif domain == 2 || domain ==3 || domain == 4
            % quefrency domain
            for ch = 1:chans
                for b = 1:bands
                    for n4 = 1:d4
                        for n5 = 1:d5
                            for n6 = 1:d6
                                if domain == 2
                                    audio(:,ch,b,n4,n5,n6) = cceps(audio(:,ch,b,n4,n5,n6));
                                elseif domain == 3
                                    audio(:,ch,b,n4,n5,n6) = hilbert(audio(:,ch,b,n4,n5,n6));
                                elseif domain == 4
                                    audio(:,ch,b,n4,n5,n6) = frft(audio(:,ch,b,n4,n5,n6),0.5);
                                end
                            end
                        end
                    end
                end
            end
            % otherwise time domain
        end
        
        if operation == 1
            % average the stack
            audio = mean(audio,dim);
        elseif operation == 2
            % trimmean
            audio = trimmean(audio,percent,'weighted',dim);
        elseif operation == 3
            % median
            audio = median(audio,dim);
        elseif operation == 4
            % max - may be useful in freq & Hilbert domains
            audio = max(audio,[],dim);
        elseif operation == 5
            % min - may be useful in freq & Hilbert domains
            audio = min(audio,[],dim);
        elseif operation == 6
            % percentile - may be useful in freq & Hilbert domains
            audio = prctile(audio,percent,dim);
        elseif operation == 7
            % sum
            audio = sum(audio,dim);
        elseif operation == 8
            % mean difference
            audio = mean(diff(audio,1,dim),dim);
        end
        
        % Go back to real time domain if necessary
        if domain == 1
            % from frequency domain
            audio = ifft(audio);
        elseif domain == 2 || domain == 4
            % from quefrency domain
            [~,chans,bands,d4,d5,d6] = size(audio);
            for ch = 1:chans
                for b = 1:bands
                    for n4 = 1:d4
                        for n5 = 1:d5
                            for n6 = 1:d6
                                if domain == 2
                                    audio(:,ch,b,n4,n5,n6) = real(icceps(audio(:,ch,b,n4,n5,n6)));
                                elseif domain == 4
                                    audio(:,ch,b,n4,n5,n6) = real(frft(audio(:,ch,b,n4,n5,n6),-0.5));
                                end
                            end
                        end
                    end
                end
            end
        elseif domain == 3
            % from Hilbert (analytic signal)
            audio = abs(audio).*cos(angle(audio));
        end
        % otherwise real time domain - no inverse transform needed
    end
    
    
    if isstruct(IN)
        OUT = IN;
        OUT.audio = audio;
    else
        
        OUT = audio;
    end
    varargout{1} = fs;
    
else
    
    OUT = [];
end

