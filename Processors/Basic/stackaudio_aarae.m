function out = stackaudio_aarae(in,method)
% This function stacks audio that was generated in AARAE using repeated
% cycles test sigals into dimension 4, or optionally in dimension 2 if it
% is singleton.
%
% The input must be an audio structure with the properties.startflag field.
%

if ~isstruct(in)
    out = [];
    return
end

if ~isfield(in,'properties')
    out = [];
    return
end

if ~isfield(in.properties,'startflag')
    out = [];
    return
end

startflag = in.properties.startflag;
dimX = length(startflag);
if dimX<2
    out = [];
    return
end
len2 = startflag(2)-startflag(1);

[len,chans,bands,d4,d5,d6] = size(in.audio);

if nargin == 1
    method = 1;
    
    if size(in.audio,2) == 1
        choice = questdlg('Dimension 2 is singleton, so do you wish to stack in dimension 2 instead of dimension 4?',...
            'Stack audio',...
            'Yes', 'No','Yes');
        if strcmp(choice,'Yes')
            method = 2;
        end
    end
    
    
    
end

if d4 ~= 1 &&  chans == 1
    method = 2;
end

if d4 ~= 1 &&  chans ~= 1 &&  d5 == 1
    method = 5; % stack in dim 5
end

if d4 ~= 1 &&  chans ~= 1 ...
        &&  d5 ~= 1 &&  d6 == 1
    method = 6; % stack in dim 6
end

if d4 ~= 1 &&  chans ~= 1 ...
        &&  d5 ~= 1 &&  d6 ~= 1
    out = []; % give up - only six dimension supported (dim3 is reserved for bands)
    return
end

% use relgain for chanIDs or similar
if isfield(in,'properties') 
    if isfield(in.properties,'relgain')
        relgain = in.properties.relgain;
    end
end

% zero-pad end if needed
if len < startflag(end)+len2-1
    in.audio = [in.audio; zeros(len2-1,chans,bands,d4,d5,d6)];
end

switch method
    case 2
        audiotemp = zeros(len2,dimX,bands,d4,d5,d6);
        for d=1:dimX
            audiotemp(:,d,:,:,:,:) = ...
                in.audio(startflag(d):startflag(d)+len2-1,1,:,:,:,:);
        end
        if isstruct(in)
            if exist('relgain','var')
                out.chanID = cellstr([repmat('chan ',size(audiotemp,2),1),...
                    num2str((1:size(audiotemp,2))'),...
                    repmat(', ',size(audiotemp,2),1),...
                    num2str(relgain'),...
                    repmat(' dB',size(audiotemp,2),1)]);
            else
                out.chanID = cellstr([repmat('chan ',size(audiotemp,2),1),...
                    num2str((1:size(audiotemp,2))')]);
            end
        end
        
    case 5
        audiotemp = zeros(len2,chans,bands,d4,dimx,d6);
        for d=1:dimX
            audiotemp(:,:,:,:,d,:) = ...
                in.audio(startflag(d):startflag(d)+len2-1,:,:,:,1,:);
        end
        
        if isstruct(in)
            if exist('relgain','var')
                out.dim5ID = cellstr([repmat('chan ',size(audiotemp,5),1),...
                    num2str((1:size(audiotemp,5))'),...
                    repmat(', ',size(audiotemp,5),1),...
                    num2str(relgain'),...
                    repmat(' dB',size(audiotemp,5),1)]);
            else
                out.dim5ID = cellstr([repmat('chan ',size(audiotemp,5),1),...
                    num2str((1:size(audiotemp,5))')]);
            end
        end
        
    case 6
        audiotemp = zeros(len2,chans,bands,d4,d5,dimx);
        for d=1:dimX
            audiotemp(:,:,:,:,:,d) = ...
                in.audio(startflag(d):startflag(d)+len2-1,:,:,:,:,1);
        end

        if isstruct(in)
            if exist('relgain','var')
                out.dim6ID = cellstr([repmat('chan ',size(audiotemp,6),1),...
                    num2str((1:size(audiotemp,6))'),...
                    repmat(', ',size(audiotemp,6),1),...
                    num2str(relgain'),...
                    repmat(' dB',size(audiotemp,6),1)]);
            else
                out.dim6ID = cellstr([repmat('chan ',size(audiotemp,6),1),...
                    num2str((1:size(audiotemp,6))')]);
            end
        end        
        
        
    otherwise      
        audiotemp = zeros(len2,chans,bands,dimX,d5,d6);
        for d=1:dimX
            audiotemp(:,:,:,d,:,:) = ...
                in.audio(startflag(d):startflag(d)+len2-1,:,:,1,:,:);
        end
        
        if isstruct(in)
            if exist('relgain','var')
                out.dim4ID = cellstr([repmat('chan ',size(audiotemp,4),1),...
                    num2str((1:size(audiotemp,4))'),...
                    repmat(', ',size(audiotemp,4),1),...
                    num2str(relgain'),...
                    repmat(' dB',size(audiotemp,4),1)]);
            else
                out.dim4ID = cellstr([repmat('chan ',size(audiotemp,4),1),...
                    num2str((1:size(audiotemp,4))')]);
            end
        end            
        
end
out.audio = audiotemp;
end

