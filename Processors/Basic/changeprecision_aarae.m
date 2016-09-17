function OUT = changeprecision_aarae(IN,format)
% Normally audio in AARAE is in double precision format. However, this
% function allows it to be changed to another format. This function assumes
% that the input is in the range -1 to 1. If it is outside that range, peak
% clipping will occur (at least for the integer formats). 
%
% This function may be useful in studing the effects of quantization. 
%
% Currently the unsigned formats are not specifically supported by other 
% parts of AARAE (i.e. data may be interpreted as having a substantial DC 
% offset).
if nargin == 1
    prompt = {'Output choice: 1-bit unsigned (logical, without upsampling or delta modulation) [0]; 8-bit unsigned integer [1]; 8-bit signed integer (first 8 bits of a 16 bit format) [2]; 16-bit signed integer [3]; 32-bit signed integer [4]; 32-bit floating point (single) [5]; 64-bit floating point (double) [6]'};
    dlg_title = 'Audio wave precision';
    num_lines = 1;
    def = {'2'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    format = str2num(char(answer));
end

OUT = IN;

if ~isempty(format)
    switch format
        case 0
            % 1-bit (without upsampling or delta modulation!)
            OUT.audio = IN.audio > 0;
            if isfield(IN,'audio2')
                OUT.audio2 = IN.audio2 > 0;
            end
        case 1
            % 8-bit unsigned
            OUT.audio = uint8((IN.audio + 1)*127.5);
            if isfield(IN,'audio2')
                OUT.audio2 = uint8((IN.audio2 + 1)*127.5);
            end
        case 2
            % 8-bit signed (using first 8 bits of 16 bit format)
            OUT.audio = int16(IN.audio*128 - 0.5);
            if isfield(IN,'audio2')
                OUT.audio2 = int16(IN.audio2*128 - 0.5);
            end    
        case 3
            % 16-bit signed
            OUT.audio = int16(IN.audio*32768 - 0.5);
            if isfield(IN,'audio2')
                OUT.audio2 = int16(IN.audio2*32768 - 0.5);
            end
        case 4
            % 32-bit signed
            OUT.audio = int32(IN.audio*2^32 - 0.5);
            if isfield(IN,'audio2')
                OUT.audio2 = int32(IN.audio2*2^32 - 0.5);
            end
        case 5
            % 32-bit float
            OUT.audio = single(IN.audio);
            if isfield(IN,'audio2')
                OUT.audio2 = single(IN.audio2);
            end
        case 6
            % 64-bit float
            OUT.audio = double(IN.audio);
            if isfield(IN,'audio2')
                OUT.audio2 = double(IN.audio2);
            end
    end
    OUT.funcallback.name = 'changeprecision_aarae.m';
    OUT.funcallback.inarg = {format};
else
    OUT = [];
    return
end



    