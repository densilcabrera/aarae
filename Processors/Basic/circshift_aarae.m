function out = circshift_aarae(in,shiftval)
% This function applies a circular shift in any or all of the first six
% dimensions of an AARAE audio field, and associated fields (such as cal,
% chanID and bandID).
%
% If an element in shiftval is positive, the values are shifted down
% (or to the right). If it is negative, the values are shifted up (or
% to the left).
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
    menustring1 = ['Shift Dim 1: Time samples (' num2str(len) ')'];
    menustring2 = ['Shift Dim 2: Channels (' num2str(chans) ')'];
    menustring3 = ['Shift Dim 3: Bands (' num2str(bands) ')'];
    menustring4 = ['Shift Dim 4 (' num2str(dim4) ')'];
    menustring5 = ['Shift Dim 5 (' num2str(dim5) ')'];
    menustring6 = ['Shift Dim 6 (' num2str(dim6) ')'];
        
                param = inputdlg({menustring1;... 
                      menustring2;...
                      menustring3;...
                      menustring4;...
                      menustring5;...
                      menustring6},...
                      'Circular shift in each dimension',... 
                      [1 60],... 
                      {'0';'0';'0';'0';'0';'0'}); 

    param = str2num(char(param)); 

    if length(param) < 6, param = []; end 
    if ~isempty(param) 
        shiftval = param;
    else
        out=[];
        return
    end
end

if isstruct(in)
    
    audio = circshift(audio,shiftval);
    
    
    
    if shiftval(2) ~= 0      
        if isfield(in,'cal')
        	out.cal = circshift(in.cal,shiftval(2));
        end
        
        if isfield(in,'chanID')
            out.chanID = circshift(in.chanID,shiftval(2));
        end
    end
    
    
    if shiftval(3) ~= 0        
        if isfield(in,'bandID')
            out.bandID = circshift(in.bandID,shiftval(3));
        end
    end
    
    
    if shiftval(4) ~= 0        
        if isfield(in,'dim4ID')
            out.dim4ID = circshift(in.dim4ID,shiftval(4));
        end
    end
    
    
    if shiftval(5) ~= 0
        if isfield(in,'dim5ID')
            out.dim5ID = circshift(in.dim5ID,shiftval(5));
        end
    end
     
    
    if shiftval(6) ~= 0        
        if isfield(in,'dim6ID')
            out.dim6ID = circshift(in.dim6ID,shiftval(6));
        end
    end
    
    out.audio = audio;
    out.funcallback.name = 'circshift_aarae.m';
    out.funcallback.inarg = {shiftval};
    
else
    out = audio;
end


    

end
