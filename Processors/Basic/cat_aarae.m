function out = cat_aarae(in,catdim)
% This function concatenates the input audio with audio selected from a
% dialog box
    [len,chans,bands,dim4,dim5,dim6] = size(in.audio);
    fs = in.fs;
    
    selection = choose_audio; % call AARAE's choose_audio function
    if ~isempty(selection)
        auxiliary_audio = selection.audio; % additional audio data
        fs2 = selection.fs; % sampling rate
        if ~(fs2 == fs)
            % match sampling rates if desired
            gcd_fs = gcd(fs,fs2); % greatest common denominator
            auxiliary_audio = resample(auxiliary_audio,fs/gcd_fs,fs2/gcd_fs);
            % note that you can also improve the accuracy of resampling
            % by specifying the filter characteristics (see Matlab's
            % help on resample)
        end
        [len2, chans2, bands2,dim42,dim52,dim62] = size(auxiliary_audio); % 2nd wave dimensions
    else
        out = [];
        return
    end
    
    dimstring1 = [num2str(len),' x ',num2str(chans),' x ',num2str(bands),' x ',...
        num2str(dim4),' x ',num2str(dim4),' x ',num2str(dim6)];
    
    dimstring2 = [num2str(len2),' x ',num2str(chans2),' x ',num2str(bands2),' x ',...
        num2str(dim42),' x ',num2str(dim42),' x ',num2str(dim62)];
    
    if ~exist('catdim','var')
        catdim = 2;
    end
    
    param = inputdlg({['ENTER THE DIMENSION TO CONCATENATE (1st audio dimensions: ,',...
        dimstring1,';     2nd audio dimensions: ,',dimstring2,')'];},...
        'Audio Concatenation',...
        [1 60],...
        {num2str(catdim)});
    
    param = str2num(char(param));

    if length(param) ~= 1, out = []; return; end 
    if ~isempty(param) 
        catdim = param(1);
    else
        out = []; return;
    end

    catdim = round(abs(catdim));
    if catdim < 1 || catdim > 6
        out = [];
        warndlg('Audio can be concatenated in dimensions 1 to 6 only','Whoops...!');
        return
    end
    
    % new wave dimensions
    if catdim == 1
        outlen = len+len2;
    else
        outlen = max([len,len2]);
    end
    
    if catdim == 2
        outchans = chans+chans2;
    else
        outchans = max([chans,chans2]);
    end
    
    if catdim == 3
        outbands = bands+bands2;
    else
        outbands = max([bands,bands2]);
    end
    
    if catdim == 4
        outdim4 = dim4+dim42;
    else
        outdim4 = max([dim4,dim42]);
    end
    
    if catdim == 5
        outdim5 = dim5+dim52;
    else
        outdim5 = max([dim5,dim52]);
    end
    
    if catdim == 6
        outdim6 = dim6+dim62;
    else
        outdim6 = max([dim6,dim62]);
    end
    
    out.audio = zeros(outlen,outchans,outbands,outdim4,outdim5,outdim6);
    
    
    out.audio(1:len,1:chans,1:bands,1:dim4,1:dim5,1:dim6) = in.audio;
    
    switch catdim
        case 1
            out.audio(len+1:end,1:chans2,1:bands2,1:dim42,1:dim52,1:dim62) = auxiliary_audio;
        case 2
            out.audio(1:len2,chans+1:end,1:bands2,1:dim42,1:dim52,1:dim62) = auxiliary_audio;
        case 3
            out.audio(1:len2,1:chans2,bands+1:end,1:dim42,1:dim52,1:dim62) = auxiliary_audio;
        case 4
            out.audio(1:len2,1:chans2,1:bands2,dim4+1:end,1:dim52,1:dim62) = auxiliary_audio;
        case 5
            out.audio(1:len2,1:chans2,1:bands2,1:dim42,dim5+1:end,1:dim62) = auxiliary_audio;
        case 6
            out.audio(1:len2,1:chans2,1:bands2,1:dim42,1:dim52,dim6+1:end) = auxiliary_audio;
    end
    
    
    
    
    
    % resolve or delete cal
    if catdim == 2 && isfield(in,'cal') && isfield(selection,'cal')
        if length(in.cal) == chans && length(selection.cal) == chans2
            out.cal = [in.cal(:);selection.cal(:)]';
        else
            out.cal = [];
        end
    else
        if isfield(in,'cal') && ~isfield(selection,'cal')
            if length(in.cal) == size(out.audio,2)
                out.cal = in.cal;
            else
                out.cal = [];
            end
            
        elseif ~isfield(in,'cal') && isfield(selection,'cal')
            if length(selection.cal) == size(out.audio,2)
                out.cal = selection.cal;
            else
                out.cal = [];
            end
            
        elseif isfield(in,'cal') && isfield(selection,'cal')
            if length(selection.cal) == size(out.audio,2) ...
                    && length(in.cal) == size(out.audio,2) ...
                    && sum(abs(selection.cal(:) - in.cal(:)))==0
                out.cal = in.cal;
            else
                out.cal = [];
            end
            
        else
            out.cal = [];
        end
    end
    
    
    
    
    
    % resolve or delete chanIDs
    if catdim == 2
        if isfield(in,'chanID') && isfield(selection,'chanID')
            if length(in.chanID) == chans && length(selection.chanID) == chans2
                out.chanID = [in.chanID;selection.chanID];
            else
                out.chanID = makechanID(size(out.audio,2),0);
            end
        else
            out.chanID = makechanID(size(out.audio,2),0);
        end
    else
        if isfield(in,'chanID')
            if size(in.chanID,1) == size(out.audio,2)
                out.chanID = in.chanID;
            elseif isfield(selection,'chanID')
                if size(selection.chanID,1) == size(out.audio,2)
                    out.chanID = selection.chanID;
                else
                    out.chanID = makechanID(size(out.audio,2),0);
                end
            else
                out.chanID = makechanID(size(out.audio,2),0);
            end
        elseif isfield(selection,'chanID')
            if size(selection.chanID,1) == size(out.audio,2)
                out.chanID = selection.chanID;
            else
                out.chanID = makechanID(size(out.audio,2),0);
            end
        else
            out.chanID = makechanID(size(out.audio,2),0);
        end
    end
    
    
    
    
    
    % resolve or delete bandIDs
    if catdim == 3
        if isfield(in,'bandID') && isfield(selection,'bandID')
            if length(in.bandID) == bands && length(selection.bandID) == bands2
                out.bandID = [in.bandID(:);selection.bandID(:)];
            else
                out.bandID = [];
            end
        else
            out.bandID = [];
        end
    else

        if isfield(in,'bandID') && ~isfield(selection,'bandID')
            if length(in.bandID) == size(out.audio,3)
                out.bandID = in.bandID;
            else
                out.bandID = [];
            end
            
        elseif ~isfield(in,'bandID') && isfield(selection,'bandID')
            if length(selection.bandID) == size(out.audio,3)
                out.bandID = selection.bandID;
            else
                out.bandID = [];
            end
            
        elseif isfield(in,'bandID') && isfield(selection,'bandID')
            if length(selection.bandID) == size(out.audio,3) ...
                    && length(in.bandID) == size(out.audio,3) ...
                    && sum(abs(selection.bandID(:) - in.bandID(:)))==0
                out.bandID = in.bandID;
            else
                out.bandID = [];
            end
            
        else
            out.bandID = [];
        end
    end
            
            
    
    % concatenate history
    
    histrow = cell(1,4);
    histrow{1,1} = datestr(now);
    histrow{1,2} = 'Concatenated';
    if isfield(in,'name')
        histrow{1,3} = in.name;
    end
    if isfield(selection,'name')
        histrow{1,4} = selection.name;
    end
    if isfield(in,'history') && isfield(selection,'history')
        out.history = [histrow; in.history; selection.history; cell(1,4)];
    elseif isfield(in,'history') && ~isfield(selection,'history')
        out.history = [histrow; in.history; cell(1,4)];
    elseif ~isfield(in,'history') && isfield(selection,'history')
        out.history = [histrow; selection.history; cell(1,4)];
    else
        out.history = histrow;
    end
    
    
    
    
    
    
    out.funcallback.name = 'cat_aarae.m';
    out.funcallback.inarg = {catdim};
end