function refreshplots(handles,axes)

selectedNodes = handles.mytree.getSelectedNodes;
signaldata = selectedNodes(1).handle.UserData;

% % % find dimension of plot (input or output channel)
% % switch get(handles.chandisp, 'Value')
% %     case 1 
% %         chandim = 2;
% %     case 2
% %         chandim = 5;
% % end
% %         
%if isa(signaldata.audio,'memmapfile'), signaldata.audio = signaldata.audio.Data; end
if ~ismatrix(signaldata.audio)
    % consider writing better code than this? it may take a long time for large
    % data
    if get(handles.chandisp, 'Value') == 2
        linea(:,:,:,:,1) = signaldata.audio(:,:,:,:,str2double(get(handles.IN_nchannel,'String')));
    else
        linea(:,:) = signaldata.audio(:,str2double(get(handles.IN_nchannel,'String')),:);
    end
else
    linea = signaldata.audio;
end
if size(linea,2) > handles.Settings.maxlines
    linea = linea(:,1:handles.Settings.maxlines);
end
plottype = get(handles.(matlab.lang.makeValidName([axes '_popup'])),'Value');
%if plottype == 4 || plottype == 5 || plottype > 7
    set(handles.(matlab.lang.makeValidName(['To_' axes])),'Visible','on')
    To_s = str2double(get(handles.(matlab.lang.makeValidName(['To_' axes])),'String'));
    To = floor(To_s*signaldata.fs)+1;
    set(handles.(matlab.lang.makeValidName(['Tf_' axes])),'Visible','on')
    Tf_s = str2double(get(handles.(matlab.lang.makeValidName(['Tf_' axes])),'String'));
    Tf = floor(Tf_s*signaldata.fs);
    if Tf > length(linea), Tf = length(linea); end
    linea = linea(To:Tf,:);
%else
%    set([handles.(matlab.lang.makeValidName(['To_' axes])),handles.(matlab.lang.makeValidName(['Tf_' axes]))],'Visible','off')
%end
units = '';
units_ref = 1;
units_type = 1;
if isfield(signaldata,'cal') && handles.Settings.calibrationtoggle == 1
    if size(linea,2) == length(signaldata.cal)
        if isfield(signaldata,'properties')
            if isfield(signaldata.properties,'units')
                units = signaldata.properties.units;
            else
                units = '';
            end
            if isfield(signaldata.properties,'units_ref')
                units_ref = signaldata.properties.units_ref;
            else
                units_ref = 1;
            end
            if isfield(signaldata.properties,'units_type')
                units_type = signaldata.properties.units_type;
            else
                units_type = 1;
            end
        else
            units = '';
            units_ref = 1;
            units_type = 1;
        end
        signaldata.cal(isnan(signaldata.cal)) = 0;
        if units_type == 1
            if numel(units_ref) == 1
                linea = linea * units_ref;
                signaldata.cal = signaldata.cal ./ 10.^(units_ref./20);
            elseif length(units_ref) == length(signaldata.cal)
                linea = linea .* repmat(units_ref,size(linea,1),1);
                signaldata.cal = signaldata.cal ./ 10.^(units_ref./20);
            end
        else
            if numel(units_ref) == 1
                linea = linea * units_ref.^0.5;
                signaldata.cal = signaldata.cal ./ 10.^(units_ref/10);
            elseif length(units_ref) == length(signaldata.cal)
                linea = linea .* repmat(units_ref.^0.5,size(linea,1),1);
                signaldata.cal = signaldata.cal ./ 10.^(units_ref./10);
            end
        end
        linea = linea.*repmat(10.^(signaldata.cal./20),length(linea),1);
    elseif ~ismatrix(signaldata.audio) && size(signaldata.audio,2) == length(signaldata.cal)
        if isfield(signaldata,'properties')
            if isfield(signaldata.properties,'units')
                units = signaldata.properties.units;
            else
                units = '';
            end
            if isfield(signaldata.properties,'units_ref')
                units_ref = signaldata.properties.units_ref;
            else
                units_ref = 1;
            end
            if isfield(signaldata.properties,'units_type')
                units_type = signaldata.properties.units_type;
            else
                units_type = 1;
            end
        else
            units = '';
            units_ref = 1;
            units_type = 1;
        end
        signaldata.cal(isnan(signaldata.cal)) = 0;
        if units_type == 1
            linea = linea * units_ref;
            signaldata.cal = signaldata.cal ./ 10.^(units_ref/20);
        else
            linea = linea * units_ref.^0.5;
            signaldata.cal = signaldata.cal ./ 10.^(units_ref/10);
        end
        if get(handles.chandisp, 'Value') == 2
            cal = signaldata.cal;
        else
            cal = repmat(signaldata.cal(str2double(get(handles.IN_nchannel,'String'))),1,size(linea,2));
        end
        linea = linea.*repmat(10.^(cal./20),length(linea),1);
    end
else
    units = '';
    units_ref = 1;
    units_type = 1;
end
fftlength = length(linea);
set(handles.(matlab.lang.makeValidName(['smooth' axes '_popup'])),'Visible','off');
switch handles.Settings.specmagscale;
    case {'Divided by length'}
        spectscale = 1./length(linea);
    case {'x sqrt2/length'}
        spectscale = 2.^0.5./length(linea);
    case {'x 2/length'}
        spectscale = 2./length(linea);
    otherwise
        spectscale = 1;
end
if plottype == 1, linea = real(linea); end
if plottype == 2, linea = linea.^2; end
if plottype == 3
    if units_type == 1
        if numel(units_ref) == 1
            linea = 10.*log10(linea.^2 ./ units_ref.^2);
        else
            linea = 10.*log10(linea.^2 ./...
                repmat(units_ref.^2,size(linea,1),1));
        end
    else
        if numel(units_ref) == 1
            linea = 10.*log10(linea.^2 ./ units_ref);
        else
            linea = 10.*log10(linea.^2 ./...
                repmat(units_ref,size(linea,1),1));
        end
    end
end
if plottype == 4, linea = abs(hilbert(real(linea))); end
if plottype == 5, linea = medfilt1(diff([angle(hilbert(real(linea))); zeros(1,size(linea,2))])*signaldata.fs/2/pi, 5); end
if plottype == 6, linea = abs(linea); end
if plottype == 7, linea = imag(linea); end
if plottype == 8
    if units_type == 1
        if numel(units_ref) == 1
            linea = 10*log10(abs(fft(linea).*spectscale  ./ units_ref).^2);
        else
            linea = 10*log10(abs(fft(linea).*spectscale  ./ ...
                repmat(units_ref,size(linea,1),1)).^2);
        end
    else
        if numel(units_ref) == 1
            linea = 10*log10(abs(fft(linea).*spectscale  ./ units_ref.^0.5).^2);
        else
            linea = 10*log10(abs(fft(linea).*spectscale  ./ ...
                repmat(units_ref.^0.5,size(linea,1),1)).^2);
        end
    end
end %freq
if plottype == 9, linea = (abs(fft(linea)).*spectscale).^2; end
if plottype == 10, linea = abs(fft(linea)).*spectscale; end
if plottype == 11, linea = real(fft(linea)).*spectscale; end
if plottype == 12, linea = imag(fft(linea)).*spectscale; end
if plottype == 13, linea = angle(fft(linea)); end
if plottype == 14, linea = unwrap(angle(fft(linea))); end
if plottype == 15, linea = angle(fft(linea)) .* 180/pi; end
if plottype == 16, linea = unwrap(angle(fft(linea))) ./(2*pi); end
if plottype == 17, spec = fft(linea,fftlength); linea = -diff(unwrap(angle(spec))).*length(spec)/(signaldata.fs*2*pi).*1000; end
if strcmp(get(handles.(matlab.lang.makeValidName(['smooth' axes '_popup'])),'Visible'),'on')
    smoothfactor = get(handles.(matlab.lang.makeValidName(['smooth' axes '_popup'])),'Value');
    if smoothfactor == 2, octsmooth = 1; end
    if smoothfactor == 3, octsmooth = 3; end
    if smoothfactor == 4, octsmooth = 6; end
    if smoothfactor == 5, octsmooth = 12; end
    if smoothfactor == 6, octsmooth = 24; end
    if smoothfactor ~= 1, linea = octavesmoothing(linea, octsmooth, signaldata.fs); end
end
t = (linspace(To_s,Tf_s,length(linea))).';
f = (signaldata.fs .* ((1:fftlength)-1) ./ fftlength).';
if plottype <= 7
    if ~isreal(signaldata.audio)
        set(handles.(matlab.lang.makeValidName(['complex' axes])),'Visible','on');
    else
        set(handles.(matlab.lang.makeValidName(['complex' axes])),'Visible','off');
    end
    set(handles.(matlab.lang.makeValidName(['log' axes '_chk'])),'Visible','off');
    pixels = get_axes_width(handles.(matlab.lang.makeValidName(['axes' axes])));
    [t, linea] = reduce_to_width(t, real(linea), pixels, [-inf inf]);
    t = t(:,1);
    plot(handles.(matlab.lang.makeValidName(['axes' axes])),t,real(linea)); % Plot signal in time domain
    xlabel(handles.(matlab.lang.makeValidName(['axes' axes])),'Time [s]');
    %xlim(handles.(matlab.lang.makeValidName(['axes' axes])),[0 length(signaldata.audio)/signaldata.fs])
    xlim(handles.(matlab.lang.makeValidName(['axes' axes])),[To_s Tf_s])
    set(handles.(matlab.lang.makeValidName(['axes' axes])),'XScale','linear','XTickLabelMode','auto')
    set(handles.(matlab.lang.makeValidName(['axes' axes])),'XTickLabel',num2str(get(handles.(matlab.lang.makeValidName(['axes' axes])),'XTick').'))
end
if plottype >= 8
    set(handles.(matlab.lang.makeValidName(['complex' axes])),'Visible','off')
    pixels = get_axes_width(handles.(matlab.lang.makeValidName(['axes' axes])));
    log_check = get(handles.(matlab.lang.makeValidName(['log' axes '_chk'])),'Value');
    if log_check == 0
        [f1, linea1] = reduce_to_width(f(1:length(linea)), linea, pixels, [-inf inf]);
        f1 = f1(:,1);
    else
        [f1, linea1] = reduce_to_width(log10(1:length(linea))', linea, pixels, [-inf inf]);
        f1 = f1(:,1);
        f1 = (10.^f1)./max(10.^f1).*signaldata.fs;
    end
    plot(handles.(matlab.lang.makeValidName(['axes' axes])),f1,linea1);% Plot signal in frequency domain
    xlabel(handles.(matlab.lang.makeValidName(['axes' axes])),'Frequency [Hz]');
    if ischar(handles.Settings.frequencylimits)
        xlim(handles.(matlab.lang.makeValidName(['axes' axes])),[f1(2) signaldata.fs/2])
    else
        xlim(handles.(matlab.lang.makeValidName(['axes' axes])),handles.Settings.frequencylimits)
    end
    set(handles.(matlab.lang.makeValidName(['log' axes '_chk'])),'Visible','on');
    if log_check == 1
        set(handles.(matlab.lang.makeValidName(['axes' axes])),'XScale','log')
        set(handles.(matlab.lang.makeValidName(['axes' axes])),'XTickLabel',num2str(get(handles.(matlab.lang.makeValidName(['axes' axes])),'XTick').'))
    else
        set(handles.(matlab.lang.makeValidName(['axes' axes])),'XScale','linear','XTickLabelMode','auto')
        set(handles.(matlab.lang.makeValidName(['axes' axes])),'XTickLabel',num2str(get(handles.(matlab.lang.makeValidName(['axes' axes])),'XTick').'))
    end
end

if ~exist('units','var'),units = ''; end
switch plottype
    case {1,4,6,7,10,11,12}
        if ~isempty(units) && isfield(signaldata,'cal') && handles.Settings.calibrationtoggle == 1
            if units_type == 1
                units_text = units;
            else
                units_text = ['(' units ')^0.5'];
            end
            set(handles.(matlab.lang.makeValidName([axes '_units_display'])),'String',units_text);
            set(handles.(matlab.lang.makeValidName([axes '_units_display'])),'Visible','on');
        else
            set(handles.(matlab.lang.makeValidName([axes '_units_display'])),'Visible','off');
        end
    case {2,9}
        if ~isempty(units) && isfield(signaldata,'cal') && handles.Settings.calibrationtoggle == 1
            if units_type == 2
                units_text = units;
            else
                units_text = ['(' units ')^2'];
            end
            set(handles.(matlab.lang.makeValidName([axes '_units_display'])),'String',units_text);
            set(handles.(matlab.lang.makeValidName([axes '_units_display'])),'Visible','on');
        else
            set(handles.(matlab.lang.makeValidName([axes '_units_display'])),'Visible','off');
        end
    otherwise
        set(handles.(matlab.lang.makeValidName([axes '_units_display'])),'Visible','off');
end

