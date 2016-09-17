function out = importaudio
% import a mat or wav audio recording.
% this function is adapted from AARAE's choose_audio, but without the
% option of selecting from the uitree.
[filename,pathname] = uigetfile(...
    {'*.wav;*.mat','Calibration file (*.wav,*.mat)'});
[~,~,ext] = fileparts(filename);
if filename ~= 0
    % Check type of file. First 'if' is for .mat, second is for .wav
    if strcmp(ext,'.mat') || strcmp(ext,'.MAT')
        audio = importdata(fullfile(pathname,filename));
        fprintf(handles.fid,['%% Loaded additional audio from mat file.','\n']);
        fprintf(handles.fid,['X1 = importdata(fullfile(', pathname ',' filename '));\n']);
        if ~isstruct(audio)
            out.audio = audio;
            out.fs = inputdlg('Sampling frequency',...
                ['Please specify ' filename ' sampling frequency'],[1 50]);
            out.fs = str2double(char(out.fs))';
            fprintf(handles.fid,['X1.fs = ', num2str(out.fs) ';\n']);
            if isnan(out.fs) || out.fs <= 0
                out = [];
                warndlg('Cannot import file without a valid sampling frequency!','AARAE info')
                fprintf(handles.fid,['%% Failed to load additional audio: Invalid sampling frequency!','\n']);
            end
        else
            if isfield(audio,'audio') && isfield(audio,'fs')
                out = audio;
            else
                warndlg('Invalid AARAE file format!','AARAE info')
                fprintf(handles.fid,['%% Failed to load additional audio: Invalid AARAE file format!','\n']);
            end
        end
    elseif strcmp(ext,'.wav') || strcmp(ext,'.WAV')
        [out.audio,out.fs] = audioread(fullfile(pathname,filename));
%         fprintf(handles.fid,['%% Loaded audio from wav file.','\n']);
%         fprintf(handles.fid,['X1 = audioread(fullfile(', pathname ',' filename '));\n']);
    else
        out = [];
    end
else
    out = [];
end