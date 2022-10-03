function OUT = MTF(IR)
% This function calculates the modulation transfer function of an IR,
% without considering the effect of background noise. The function is
% intended to support a tutorial on ISO 3382-3 parameters, using a separate
% spreadsheet to caclulate the parameters.


fs = IR.fs;

[~,chans,bands,dim4,dim5,dim6] = size(IR.audio);
if bands > 1
    IR.audio = mean(IR.audio,3);
end

if dim4 > 1 && chans == 1
    IR.audio = permute(IR.audio,[1,4,3,2]);
    dim4=1;
    chans = size(IR.audio,2);
end

% zero-pad the IR so that MTF covers the full IR
IR.audio = [IR.audio; zeros(round(1.6*fs),chans,1,1,dim5,dim6)];
len = size(IR.audio,1);

% Nyquist frequency
Nyquist=fs/2;

% Time in seconds for each sample
time=((1:len)-1)'./fs;

% list of modulation frequencies
% mf = [0.63 0.8 1 1.25 1.6 2 2.5 3.15 4 5 6.3 8 10 12.5]; % nominal frequencies
mf = 10.^((-2:11)/10); % exact frequencies
fc = [125 250 500 1000 2000 4000 8000]; % octave band nominal frequencies

orderin = 12; % in-band filter pseudo-order
orderout = 12; % out-of-band filter pseudo-order
if exist('FilterStrength','var')
    orderin = orderin * FilterStrength;
    orderout = orderout * FilterStrength;
end

P_octave = octbandfilter_viaFFT(IR.audio,fs,...
    fc,[orderin,orderout],0,...
    1000,0,0,10);

% number of whole number cycles to use for each modulation frequency
Fm_cycles = floor(len .* mf./fs);
% number of samples to use for each modulation frequency
Fm_len = floor(fs.*Fm_cycles ./ mf);

for d4 = 1:dim4
    for d5 = 1:dim5
        for d6 = 1:dim6
            
            MTF = zeros(14,7,chans);
            for j = 1:length(mf) % derive MTF using Schroeder's formula
                MTF_num=abs(sum(P_octave(1:Fm_len(j),:,:,d4,d5,d6).^2.*exp(-2i*pi*mf(j)...
                    .*repmat(time(1:Fm_len(j)),[1,chans,length(fc)]))));
                MTF_den=sum(P_octave(1:Fm_len(j),:,:,d4,d5,d6).^2);
                MTF(j,:,:)=permute(MTF_num./MTF_den,[1,3,2]);
            end
            
            % output tables
            OUT.tables = [];
            for ch = 1:chans
                tablename = ['MTF Chan ' num2str(ch)];
                if dim4>1
                    tablename = [tablename ', Cycle ' num2str(d4)];
                end
                if dim5>1
                    tablename = [tablename ', OutChan ' num2str(d5)];
                end
                if dim6>1
                    tablename = [tablename ', D6 ' num2str(d6)];
                end
                f = figure('Name',tablename);
                dat1 = MTF(:,:,ch);
                cnames1 = {'125', '250', '500', '1k', '2k', '4k', '8k'};
                rnames1 = {'0.63', '0.8', '1', '1.25', '1.6', '2',...
                    '2.5', '3.15','4','5','6.3','8','10','12.5'};
                t1 =uitable('Data',dat1,'ColumnName',cnames1,'RowName',rnames1);
                [~,tables] = disptables(f,t1,{tablename});
                OUT.tables = [OUT.tables tables];
            end
        end
    end
end

OUT.funcallback.name = 'MTF.m';
OUT.funcallback.inarg = {};
