function out = flipdim_aarae(in,flipvector)
% Dimensions of the audio can be flipped, using a list dialog to select
% dimensions, or by specifying a list of dimensions to be flipped (optional
% flipvector input). Other AARAE fields that are associated with particular
% dimensions are flipped accordingly (e.g. cal, chanID and bandID).
%
% Dimension 1 is time, so flip means time reversal
% Dimension 2 is channels - sometimes flipping channels can make them
% easier to visualise. ChanIDs are also flipped.
% Dimension 3 is bands. BandIDs are also flipped.
% Higher dimensions are not defined (although dimension 4 is often used for
% IR stacks). This function supports up to 6 dimensions.
%
% Code by Densil Cabrera 
% version 1.0 (2 August 2014)

if isstruct(in)
    audio = in.audio;
    out = in;
else
    audio = in;
end

[len,chans,bands,dim4,dim5,dim6]=size(audio);

if nargin == 1
    param = [1;2;3;4;5;6];
    menustring1 = ['Dim ' num2str(param(1)) ' Time (' num2str(len) ')'];
    menustring2 = ['Dim ' num2str(param(2)) ' Channels (' num2str(chans) ')'];
    menustring3 = ['Dim ' num2str(param(3)) ' Bands (' num2str(bands) ')'];
    menustring4 = ['Dim ' num2str(param(4)) ' (' num2str(dim4) ')'];
    menustring5 = ['Dim ' num2str(param(5)) ' (' num2str(dim5) ')'];
    menustring6 = ['Dim ' num2str(param(6)) ' (' num2str(dim6) ')'];
    menulist = {menustring1; menustring2; menustring3; menustring4; menustring5; menustring6};
        
            [S,ok] = listdlg('Name','Flip',...
                                     'PromptString','Dimension(s) to flip',...
                                     'ListString',menulist);
            flipvector = param(S);
            %disp(flipvector)
end

if isstruct(in)
    if ~isempty(find(flipvector==1, 1))
        audio = flipdim(audio,1);
    end
    
    if ~isempty(find(flipvector==2, 1))
        audio = flipdim(audio,2);
        
        if isfield(in,'cal')
            cal = in.cal;
            if ~isrow(cal), cal = cal'; end
        	out.cal = flipdim(cal,2);
        end
        
        if isfield(in,'chanID')
            chanID = in.chanID;
            if ~isrow(chanID), chanID = chanID'; end
            out.chanID = flipdim(chanID,2);
        end
    end
    
    
    if ~isempty(find(flipvector==3, 1))
        audio = flipdim(audio,3);
        
        if isfield(in,'bandID')
            bandID = in.bandID;
            if ~isrow(bandID), bandID = bandID'; end
            out.bandID = flipdim(bandID,2);
        end
    end
    
    
    if ~isempty(find(flipvector==4, 1))
        audio = flipdim(audio,4);
        
        if isfield(in,'dim4ID')
            dim4ID = in.dim4ID;
            if ~isrow(dim4ID), dim4ID = dim4ID'; end
            out.dim4ID = flipdim(dim4ID,2);
        end
    end
    
    
    if ~isempty(find(flipvector==5, 1))
        audio = flipdim(audio,5);
        
        if isfield(in,'dim5ID')
            dim5ID = in.dim5ID;
            if ~isrow(dim5ID), dim5ID = dim5ID'; end
            out.dim5ID = flipdim(dim5ID,2);
        end
    end
     
    
    if ~isempty(find(flipvector==6, 1))
        audio = flipdim(audio,6);
        
        if isfield(in,'dim6ID')
            dim6ID = in.dim6ID;
            if ~isrow(dim6ID), dim6ID = dim6ID'; end
            out.dim6ID = flipdim(dim6ID,2);
        end
    end
    
    out.audio = audio;
    out.funcallback.name = 'flipdim_aarae.m';
    out.funcallback.inarg = {flipvector};
    
else
    if ~isempty(find(flipvector==1, 1))
        audio = flipdim(audio,1);
    end
    
    if ~isempty(find(flipvector==2, 1))
        audio = flipdim(audio,2);
    end
    
    if ~isempty(find(flipvector==3, 1))
        audio = flipdim(audio,3);
    end
    
    if ~isempty(find(flipvector==4, 1))
        audio = flipdim(audio,4);
    end
    
    if ~isempty(find(flipvector==5, 1))
        audio = flipdim(audio,5);
    end
    
    if ~isempty(find(flipvector==6, 1))
        audio = flipdim(audio,6);
    end
    out = audio;
end


    

end