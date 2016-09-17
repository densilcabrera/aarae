function out = select_aarae(in)
% This function allows the audio to be selected via numeric indices in all
% of the dimensions present. This can be used to simultaneously truncate
% time, delete channels, delete bands, and limit higher dimensions if
% present, reducing the size of the audio leaf.
%

    indices = partial_selection(in);    
    %out.audio(indices{:}) = in.audio(indices{:}); % presumably an error
    out.audio = in.audio(indices{:});
       
    if isfield(in,'cal') && length(indices)>1
        if length(in.cal) == size(in.audio,2)
            out.cal = in.cal(indices{2});
        else
            out.cal = [];
        end      
    end
    
    if isfield(in,'chanID') && length(indices)>1
        if length(in.chanID) == size(in.audio,2)
            out.chanID = in.chanID(indices{2});
        else
            out.chanID = [];
        end      
    end
    
    
    if isfield(in,'bandID') && length(indices)>2
        if length(in.bandID) == size(in.audio,3)
            out.bandID = in.bandID(indices{3});
        else
            out.bandID = [];
        end      
    end
    
    if isfield(in,'properties') && length(indices)>3
        if isfield(in.properties,'relgain')
            if length(in.properties.relgain) == size(in.audio,4)
                out.properties.relgain = in.properties.relgain(indices{4});
            else
                out.properties.relgain = [];
            end
        end
    end
    
    % potentially include other fields for similar treatment

            
    out.funcallback.name = 'select_aarae.m';
    out.funcallback.inarg = {};
end