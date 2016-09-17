load roomir.mat

reps = 16;
N = 14;
DCCoupling = 0;
offset=0;
useimp=1;
method = 'mls';
tsplincirc = 'circ';

if(strcmp(method,'mls'))
    seq = GenerateMLSSequence(reps,N,useimp);
    h = h(:,:,9);   % 500ms only
    seqfilt = fftfilt(h,seq);
    analysed = AnalyseMLSSequence(seqfilt,offset,reps,N,DCCoupling,useimp);
elseif(strcmp(method,'tsp'))
    seq = GenerateTSPSequence(reps,N,useimp,tsplincirc);
    h = h(:,:,9);   % 500ms only
    seqfilt = fftfilt(h,seq);
    analysed = AnalyseTSPSequence(seqfilt,offset,reps,N,useimp,tsplincirc);
end

sz = size(h);

subplot(3,1,1);
plot(h);
xlim([1 sz(1)]);
title('Original Responses');
subplot(3,1,2);
plot(analysed);
xlim([1 sz(1)]);
title('Analysed responses');
subplot(3,1,3);
plot(analysed(:,1));hold on;plot(h(:,1),'r');
title('Overlaid channel 1. Blue=estimated, red=original');

[npm_lin] = npm(h,analysed(1:size(h,1),:));
npm_log=10*log10(npm_lin);
disp(['Average NPM: ' num2str(npm_log) ' dB']);