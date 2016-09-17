function signaldata = checkcal(signaldata)
% This function checks cal and units fields in AARAE audio structures,
% and writes neutral values to the structure.
% The cal field is set to 0 and reference value is set to 1.
% It also checks that the cal vector matches the number of channels, and 
% fixes this if there is a discrepancy.
% This function is called by AARAE framework functions in serveral places
% (whenenver there is the potential for new audio structures to be added).

nullstring = 'Val';

if isfield(signaldata,'audio')
    if ~isfield(signaldata,'cal')
        signaldata.cal = zeros(1,size(signaldata.audio,2));
    elseif length(signaldata.cal) > size(signaldata.audio,2)
        signaldata.cal = signaldata.cal(1:size(signaldata.audio,2));
    elseif length(signaldata.cal) < size(signaldata.audio,2)
        if length(signaldata.cal) == 1
            signaldata.cal = repmat(signaldata.cal,[1,size(signaldata.audio,2)]);
        else
            signaldata.cal = [signaldata.cal(:)',...
                zeros(1,size(signaldata.audio,2)-length(signaldata.cal))];
        end
    end
    
%     signaldata.cal(isnan(signaldata.cal)) = 0;
%     signaldata.cal(isinf(signaldata.cal)) = 0;

    if isfield(signaldata,'properties')
        if ~isfield(signaldata.properties,'units')
            signaldata.properties.units = nullstring;
        end
        if ~isfield(signaldata.properties,'units_ref')
            signaldata.properties.units_ref = 1;
        end
        if ~isfield(signaldata.properties,'units_type')
            signaldata.properties.units_type = 1;
        end
    else
        signaldata.properties.units = nullstring;
        signaldata.properties.units_ref = 1;
        signaldata.properties.units_type = 1;
    end
end