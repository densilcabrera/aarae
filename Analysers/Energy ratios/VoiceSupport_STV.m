function OUT = VoiceSupport_STV(IN,fs,cropthreshold,endtime)
% This function calculates voice support (STV) and room gain (GRG) 
% as described by:
%
% Pelegrín-García D. Comment on ?Increase in voice level and speaker
% comfort in lecture rooms? [J. Acoust. Soc. Am. 125, 2072?2082 (2009)]. J
% Acoust Soc Am 2011;129:1161-4 doi:10.1121/1.3543940.
%
% Pelegrín-García D, Brunskog J, Lyberg-Åhlander V, Löfqvist A. Measurement
% and prediction of voice support and room gain in school classrooms. The
% Journal of the Acoustical Society of America. 2012 Jan;131(1):194-204.
%

% Revised 2020 to include speech weighting & GRG & raw energy levels

if isstruct(IN)
    OUT = IN;
    audio = IN.audio;
    fs = IN.fs;
    if nargin < 3
        prompt = {'Start crop threshold (in dB)';...
            'End time'};
        dlg_title = 'STV parameters';
        num_lines = 1;
        def = {'-40','1'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        cropthreshold = str2double(char(answer{1,1}));
        endtime = str2double(char(answer{2,1}));
    end
else
    audio = IN;
end
% crop just before IR direct sound
audio = autocropstart_aarae(audio,cropthreshold,2);
[len,chans,bands,dim4,dim5,dim6] = size(audio);
% STV window
t = ((1:len)'-1)./fs;
t0 = 0.0045;
ind0 = round(t0.* fs)+1;
t1 = 0.0055;
ind1 = round(t1.* fs)+1;
T = 0.002;


figure('Name',['Voice Support ' IN.name]);
plot(t,audio(:,1,1,1,1,1),'color',[0.7 0.7 0.7]);
hold on


% Direct sound
audiodirect = audio(1:ind1,:,1,:,:,:);
audiodirect(ind0:ind1,:,:,:,:,:) = audiodirect(ind0:ind1,:,:,:,:,:) .* ...
repmat(0.5 + 0.5 .* cos(2*pi*(t(ind0:ind1)-t0)./T),[1,chans,1,dim4,dim5,dim6]);
plot(t(1:ind1),audiodirect(:,1,1,1,1,1),'color',[0 0.7 0]);

audiodirect = [zeros(fs/100,chans,1,dim4,dim5,dim6);...
    audiodirect;
    zeros(fs/100,chans,1,dim4,dim5,dim6)]; %zero-pad so that energy is not lost in filtering

% change 'audio' to IR excluding direct
audio(1:ind0-1,:,:,:,:,:) = 0;
audio(ind0:ind1,:,:,:,:,:) = audio(ind0:ind1,:,:,:,:,:) .* ...
repmat(0.5 - 0.5 .* cos(2*pi*(t(ind0:ind1)-t0)./T),[1,chans,1,dim4,dim5,dim6]);
plot(t,audio(:,1,1,1,1,1),'color',[0.5 0 0.5]);
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 0.02]);

if round(endtime*fs)<len
    audio = audio(1:round(endtime*fs),:,:,:,:,:);
end

% oct band filter
audiodirectoct = octbandfilter_viaFFT(audiodirect,fs,...
    [125,250,500,1000,2000,4000,8000,16000],[12,12],0,...
    1000,0,0,2);
audiooct = octbandfilter_viaFFT(audio,fs,...
    [125,250,500,1000,2000,4000,8000,16000],[12,12],0,...
    1000,0,0,2);


% Calculate STV (average channels, no weighting)
Lreflectoct = pow2db(mean(mean(mean(mean(sum(audiooct.^2),2),4),5),6));
Ldirectoct = pow2db(mean(mean(mean(mean(sum(audiodirectoct.^2),2),4),5),6));
STVoct = Lreflectoct-Ldirectoct;
% EnergyRatio = mean(mean(mean(mean(sum(audiooct.^2) ./ sum(audiodirectoct.^2),2),4),5),6);
% STVoct = pow2db(EnergyRatio);
STV = mean(STVoct(:,1,1:6,:,:,:),3);
Goct = pow2db(db2pow(STVoct)+1);
G = pow2db(db2pow(STV)+1);

% Calculate speech-weighted STV
% Typical speech levels at the eardrum
LDeardrum = [58.0 69.1 73.5 71.7 69.0 63.0]; % 125 Hz to 4 kHz bands
STVspeech = pow2db(sum(db2pow(LDeardrum+permute(STVoct(:,:,1:6),[1,3,2])))./sum(db2pow(LDeardrum)));
Gspeech = pow2db(db2pow(STVspeech)+1);

% 1/3-oct band filter
%thirdoctbandfilter_viaFFT(IN,fs,param,order,zeropad,minfftlenfactor,test,phasemode)
audiodirect3oct = thirdoctbandfilter_viaFFT(audiodirect,fs,...
    [100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,...
    4000,5000,6300,8000,10000,12500,16000],[12,12],0,...
    1000,0,0);
audio3oct = thirdoctbandfilter_viaFFT(audio,fs,...
    [100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,...
    4000,5000,6300,8000,10000,12500,16000],[12,12],0,...
    1000,0,0);
Lreflect3oct = pow2db(mean(mean(mean(mean(sum(audio3oct.^2),2),4),5),6));
Ldirect3oct = pow2db(mean(mean(mean(mean(sum(audiodirect3oct.^2),2),4),5),6));
%EnergyRatio = mean(mean(mean(mean(sum(audio3oct.^2) ./ sum(audiodirect3oct.^2),2),4),5),6);
%STV3oct = pow2db(EnergyRatio);
STV3oct = Lreflect3oct-Ldirect3oct;
G3oct = pow2db(db2pow(STV3oct)+1);

f=figure('Name',['Voice Support ' IN.name]);

cnames1 = {'125', '250', '500', '1k', '2k', '4k', '8k', '16k'};
rnames1 = {'STV oct (dB)','G oct (dB)','Ldirecect (dB)','Lreflect (dB)'};
t1 =uitable('Data',[permute(STVoct,[1,3,2]);permute(Goct,[1,3,2]);...
    permute(Ldirectoct,[1,3,2]);permute(Lreflectoct,[1,3,2])],...
    'ColumnName',cnames1,'RowName',rnames1);

cnames2 = {'100','125','160','200','250','315','400','500','630','800','1000','1200','1600','2000','2500','3150',...
    '4000','5000','6300','8000','10000','12500','16000'};
rnames2 = {'STV 1/3-oct (dB)','G 1/3-oct (dB)','Ldirecect (dB)','Lreflect (dB)'};
t2 =uitable('Data',[permute(STV3oct,[1,3,2]);permute(G3oct,[1,3,2]);...
    permute(Ldirect3oct,[1,3,2]);permute(Lreflect3oct,[1,3,2])],...
    'ColumnName',cnames2,'RowName',rnames2);

cnames3 = {'STV','GRG'};
rnames3 = {'Oct band mean (dB)','Speech-weighted (dB)'};
t3 =uitable('Data',[STV,G;STVspeech,Gspeech],'ColumnName',cnames3,'RowName',rnames3);
[~,OUT.tables] = disptables(f,[t3 t2 t1]);
OUT.funcallback.name = 'VoiceSupport_STV.m';
OUT.funcallback.inarg = {fs,cropthreshold,endtime};

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2017-2020, Densil Cabrera
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%  * Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  * Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%  * Neither the name of the University of Sydney nor the names of its contributors
%    may be used to endorse or promote products derived from this software
%    without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%