function OUT = choose_from_higher_dimensions(IN,maxdim,method)
% This AARAE utility function addresses the difficulty of having up to
% 6-dimensional audio data but analysers (and perhaps processors) that are
% only capable of dealing with a smaller number of dimensions.
%
% The input argument maxdim is the maximum number of output dimensions;
% The input argument IN is either an AARAE audio structure or audio as a
% vector or matrix. The format of OUT is the same as the format of IN.
%
% In AARAE, dimension 1 is used for time, dimension 2 for channels,
% dimension 3 for bands, dimension 4 for cycles (from multicycle audio test
% signals), dimension 5 for output channels of multicycle sequential output
% (e.g., measuring a series of loudspeakers, one by one), and dimension 6
% does not have a defined role, but is supported in many parts of AARAE.
%
% However, many analysers are not suitable for analysing 6-dimensional
% data, and may be better suited to analysing a subset.
%
% This function offers three solutions to the problem ('method'):
%
% METHOD = 1
% If the number of dimensions of the intput audio is less than or equal to
% maxdim, then the input is returned to the output without any change.
% Otherwise either a dialog box is used for user selection
%
% METHOD = 0
% Only the first index in each higher dimension is chosen - except for
% dimension 4, for which the last index is used (to avoid selecting the
% cycle with -inf dB gain). There is no dialog box.
%
% METHOD = -1 
% Audio is reshaped to fit the dimension limitations specified by maxdim
% (i.e., all audio will be returned, rather than a subset of it). Where
% possible, reshaping is done to extend dimension 2 (channels). There is no
% dialog box. It is probably not a good idea to use this approach with
% multiband (dimension 3, if maxdim < 3) data - it may be better to mixdown
% multiband data prior to calling this function. Also, using this with
% maxdim = 1 could cause many problems (but might be of some use in special
% cases)! Therefore, in general it is best to write analysers capable of
% dealing with at least the first three dimensions of audio (so that maxdim
% >= 3).
%
% It is anticipated that a user setting will be added to AARAE that will
% override the third input argument (method) if the user chooses to do so,
% thereby taking control of the method away from the function call in the
% analyser. This is not implemented yet.



% default settings
if ~exist('method','var')
    method = 1;
end


if ~exist('maxdim','var')
    maxdim = 1;
end



% interpret input
if isstruct(IN)
    audio = IN.audio;
    
    if isfield(IN,'chanID')
        chanID = IN.chanID;
    end
    
    if isfield(IN,'cal')
        cal = IN.cal;
    end
    
    if isfield(IN,'bandID')
        bandID = IN.bandID;
    end
    
    if isfield(IN,'properties')
        if isfield(IN.properties,'startflag')
            startflag = IN.properties.startflag;
        end
        if isfield(IN.properties,'relgain')
            relgain = IN.properties.relgain;
        end
    end
    
else
    audio = IN;
end



% Return output to input if there is nothing to do
Sz = size(audio);
if length(Sz) <= maxdim
    OUT = IN;
    return
end
if maxdim == 1 && Sz(2) == 1
    OUT = IN;
    return
end
    
% Override 3rd input argument (not implemented yet)
% if strcmp(handles.Settings.extradims,'Reshape Data')
%     method = -1;
% end
    

% If we get here, then we do need to do something
handles = guidata(findobj('Tag','aarae')); % for activity logging
nonsingleton = Sz>1; % non singleton dimensions
nonsingleton = [nonsingleton,zeros(1,6-length(Sz))];
nonsingleton(1:maxdim) = 0; % don't worry about dims <= maxdim
numberofdimstochange = sum(nonsingleton);
if method == 0
    % take the first index of each excess dimension, except dimension 4,
    % where we take the last index (because of its use for variable gain in
    % multicycle tests - including the silent cycle which is always first)
    indexstring = '(:,';
    if nonsingleton(2)
        audio = audio(:,1,:,:,:,:);
        if exist('chanID','var')
            chanID = chanID{1};
        end
        if exist('cal','var')
            cal = cal(1);
        end
        indexstring = [indexstring,'1,'];
    else
        indexstring = [indexstring,':,'];
    end
    if nonsingleton(3)
        audio = audio(:,:,1,:,:,:);
        if exist('bandID','var')
            bandID = bandID(1);
        end
        indexstring = [indexstring,'1,'];
    else
        indexstring = [indexstring,':,'];
    end
    if nonsingleton(4)
        audio = audio(:,:,:,end,:,:);
        if exist('startflag','var')
            startflag = startflag(end);
        end
        if exist('relgain','var')
            relgain = relgain(end); % should be 0 dB
        end
        indexstring = [indexstring,'1,'];
    else
        indexstring = [indexstring,':,'];
    end
    if nonsingleton(5)
        audio = audio(:,:,:,:,1,:);
        indexstring = [indexstring,'1,'];
    else
        indexstring = [indexstring,':,'];
    end
    if nonsingleton(6)
        audio = audio(:,:,:,:,:,1);
        indexstring = [indexstring,'1,'];
    else
        indexstring = [indexstring,':)'];
    end
    
    
    
elseif method ==-1
    % reshape the audio into maxdim
    [d1,d2,d3,d4,d5,d6] = size(audio);
    switch maxdim
        case 1
            % This is usually not be a good idea!! It would be better to make
            % sure the analyser can deal with more dimensions (ideally at least 3)
            % Furthermore, with impulse response analysis, there could be
            % big problems with reshaping to one dimension.
            audio = reshape(audio, [d1*d2*d3*d4*d5*d6,1]);
            chanID = {'allchans'};
            if isfield(IN,'cal')
                if sum(abs(diff(cal,2))) == 0
                    cal = cal(1);
                else
                    IN = rmfield(IN,'cal');
                end
            end
            if isfield(IN,'bandID')
                IN = rmfield(IN,'bandID');
            end
            indexstring = '(DATA_RESHAPED_TO_1_DIMENSION)';
            
        case 2
            % This also is usually not a good idea for multiband audio!!
            % In many cases it would be better to mix down bands prior to
            % calling this function if the 3rd dimension cannot be used
            audio = reshape(audio,[d1,d2*d3*d4*d5*d6]);
            if isfield(IN,'chanID')
                chanID = repmat(chanID,[d3*d4*d5*d6,1]);
            end
            if isfield(IN,'cal')
                cal = repmat(cal(:),[d3*d4*d5*d6,1])';
            end
            if isfield(IN,'bandID')
                IN = rmfield(IN,'bandID');
            end
            indexstring = '(DATA_RESHAPED_TO_2_DIMENSIONS)';
        case 3
            % reshape into channels
            audio = reshape(audio,[d1,d2*d4*d5*d6,d3]);
            % need to rewrite chanIDs
            % need to deal with HOA chan format issues within analysers
            if isfield(IN,'chanID')
                chanID = repmat(chanID,[d4*d5*d6,1]);
            end
            if isfield(IN,'cal')
                cal = repmat(cal(:),[d4*d5*d6,1])';
            end
            indexstring = '(DATA_RESHAPED_TO_3_DIMENSIONS)';
        case 4
            % reshape into channels
            audio = reshape(audio,[d1,d2*d5*d6,d3,d4]);
            if isfield(IN,'chanID')
                chanID = repmat(chanID,[d5*d6,1]);
            end
            if isfield(IN,'cal')
                cal = repmat(cal(:),[d5*d6,1])';
            end
            indexstring = '(DATA_RESHAPED_TO_4_DIMENSIONS)';
        case 5
            % reshape into channels
            audio = reshape(audio,[d1,d2*d6,d3,d4,d5]);
            if isfield(IN,'chanID')
                chanID = repmat(chanID,[d6,1]);
            end
            if isfield(IN,'cal')
                cal = repmat(cal(:),[d6,1])';
            end
            indexstring = '(DATA_RESHAPED_TO_5_DIMENSIONS)';
        case 6
            % everything is fine - no need to do anything
            indexstring = [];
        otherwise
            disp('Unrecognized maxdim value in choose_from_higher_dimensions (AARAE utility)')
            indexstring = [];
    end
    
    
    
else % method = 1 (or anything else)
    % prompt user for choice of indices in higher dimensions
    prompt = cell(1,numberofdimstochange);
    def = repmat({'1'},[1,numberofdimstochange]);
    dlgtitle = 'Please select one index from each higher dimension';
    m=1;
    if nonsingleton(2)
        prompt{1,m} = ['Select one channel (1-',num2str(size(audio,2)),')'];
        m=m+1;
    end
    if nonsingleton(3)
        prompt{1,m} = ['Select one band (1-',num2str(size(audio,3)),')'];
        m=m+1;
    end
    if nonsingleton(4)
        prompt{1,m} = ['Select one dim4 index (1-',num2str(size(audio,4)),')'];
        m=m+1;
    end
    if nonsingleton(5)
        prompt{1,m} = ['Select one dim5 index (1-',num2str(size(audio,5)),')'];
        m=m+1;
    end
    if nonsingleton(6)
        prompt{1,m} = ['Select one dim6 index (1-',num2str(size(audio,6)),')'];
    end
    answer = inputdlg(prompt,dlgtitle,[1 90],def); 
    indexstring = '(:,';
    m=1;
    if nonsingleton(2)
        if isempty(answer{m})
            selection = 1;
            indexstring = [indexstring,'1,'];
        else
            selection = str2double(answer{m});
            indexstring = [indexstring,char(answer{m}),','];
        end
        audio = audio(:,selection,:,:,:,:);
        if exist('chanID','var')
            chanID = chanID{selection};
        end
        if exist('cal','var')
            cal = cal(selection);
        end
        m=m+1;
    else
        indexstring = [indexstring,':,'];
    end
    if nonsingleton(3)
        if isempty(answer{m})
            selection = 1;
            indexstring = [indexstring,'1,'];
        else
            selection = str2double(answer{m});
            indexstring = [indexstring,char(answer{m}),','];
        end
        audio = audio(:,:,selection,:,:,:);
        if exist('bandID','var')
            bandID = bandID(selection);
        end
        m=m+1;
    else
        indexstring = [indexstring,':,'];
    end
    if nonsingleton(4)
        if isempty(answer{m})
            selection = 1;
            indexstring = [indexstring,'1,'];
        else
            selection = str2double(answer{m});
            indexstring = [indexstring,char(answer{m}),','];
        end
        audio = audio(:,:,:,selection,:,:);
        if exist('startflag','var')
            startflag = startflag(selection);
        end
        if exist('relgain','var')
            relgain = relgain(selection);
        end
        m=m+1;
    else
        indexstring = [indexstring,':,'];
    end
    if nonsingleton(5)
        if isempty(answer{m})
            selection = 1;
            indexstring = [indexstring,'1,'];
        else
            selection = str2double(answer{m});
            indexstring = [indexstring,char(answer{m}),','];
        end
        audio = audio(:,:,:,:,selection,:);
        m=m+1;
    else
        indexstring = [indexstring,':,'];
    end
    if nonsingleton(6)
        if isempty(answer{m})
            selection = 1;
            indexstring = [indexstring,'1)'];
        else
            selection = str2double(answer{m});
            indexstring = [indexstring,char(answer{m}),')'];
        end
        audio = audio(:,:,:,:,:,selection);
    else
        indexstring = [indexstring,':)'];
    end
end



% Provide output
if isstruct(IN)
    OUT = IN;
    OUT.audio = audio;
    if exist('chanID','var')
        OUT.chanID = chanID;
    end
    if exist('cal','var')
        OUT.cal = cal;
    end
    if exist('bandID','var')
        OUT.bandID = bandID;
    end
    if exist('startflag','var')
        OUT.properties.startflag = startflag;
    end
    if exist('relgain','var')
        OUT.properties.relgain = relgain;
    end
else
    OUT = audio;
end

handles.choosefromhigherdims = indexstring;
guidata(findobj('Tag','aarae'),handles); % write to aarae handles