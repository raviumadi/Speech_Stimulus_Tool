function speechStimulusTool
% SPEECHSTIMULUSTOOL  Interactive GUI for preparing speech segments.
%
% This tool lets you:
%   • Record or load a mono WAV file.
%   • Visually inspect the full waveform, spectrum, and spectrogram.
%   • Select a time segment for further processing and A/B compare
%     the original vs processed version.
%
% MAIN FEATURES
% -------------
% 1) Recording and I/O
%    - Select input device and sampling rate.
%    - Timed recording with a countdown shown on the Record button.
%    - Load / save full signals as WAV with peak-normalisation
%      (scaled to |x| ≤ 0.99 to avoid clipping).
%
% 2) Segment selection and playback
%    - Drag a rectangle on the waveform to define a segment.
%    - Store both the original segment (segOrig) and the current
%      processed version (segProc).
%    - Playback controls for full signal, original segment, and
%      processed segment.
%
% 3) Fade processing
%    - Apply onset/offset fades to the selected segment:
%         * Linear
%         * Raised cosine (Hann)
%         * Gaussian
%    - Fade duration specified in milliseconds.
%    - After fading, the segment is RMS-matched back to the
%      pre-fade level so loudness is preserved.
%
% 4) Band-pass filtering
%    - 4th-order Butterworth band-pass (filtfilt for zero-phase).
%    - User-defined low and high cut frequencies (constrained
%      to [0, Fs/2]).
%    - Post-filter RMS is matched to the pre-filter RMS to avoid
%      overall level changes.
%
% 5) Frequency-domain normalisation (spectral envelope shift)
%    - Goal: shift the dominant spectral envelope peak to a chosen
%      reference frequency (Freq-norm ref) while:
%          * preserving the overall spectral shape,
%          * preserving phase,
%          * preserving segment RMS.
%    - Implementation (applyFreqNormEnvelope):
%          a) Compute FFT of the current segment.
%          b) Estimate a smooth spectral envelope (Gaussian-smoothed
%             magnitude up to Nyquist).
%          c) Find the envelope peak frequency and shift this envelope
%             so that its peak aligns with the reference frequency.
%          d) Form a magnitude-shaping ratio (shifted / original
%             envelope), clipped to a reasonable range to avoid
%             extreme boosts/cuts.
%          e) Apply this ratio to the positive-frequency magnitudes
%             only, keep original phase, and reconstruct a
%             conjugate-symmetric spectrum.
%          f) Inverse FFT, trim to original length, and RMS-match
%             the result to the input segment.
%    - Effect: moves the “loudest” part of the spectrum towards the
%      reference frequency but keeps the relative structure of the
%      band and the time-domain character largely intact.
%
% 6) Visualisation
%    - Waveform: full signal or segment (original vs processed).
%    - Spectrum comparison: original vs processed segment in dB.
%    - Spectrogram: processed segment only (fixed window/overlap,
%      20*log10 magnitude, dB colormap).
%
% 7) Export and reproducibility
%    - Save processed segment as WAV (peak-limited to 0.99).
%    - Export waveform, spectrum, and spectrogram as PNG.
%    - Export current processing settings (fade type/duration,
%      band-pass limits, freq-norm ref) to a simple YAML file.
%    - Free-text notes are saved to a .txt file alongside the
%      audio and figures for later reference.

fig = uifigure('Name','Speech Stimulus Tool',...
    'Position',[100 100 1100 950],...
    'Resize','off');

% ================= LEFT: PLOTS (0–0.5 width) ========================
% Waveform (top)
axWave = uiaxes(fig,'Units','normalized',...
    'Position',[0.05 0.65 0.45 0.30]);
title(axWave,'Time Signal');
xlabel(axWave,'Time (s)');
ylabel(axWave,'Amplitude');
grid(axWave,'on');

% Spectrum (middle)
axSpec = uiaxes(fig,'Units','normalized',...
    'Position',[0.05 0.38 0.45 0.25]);
title(axSpec,'Spectrum');
xlabel(axSpec,'Frequency (kHz)');
ylabel(axSpec,'Magnitude (dB)');
grid(axSpec,'on');

% Spectrogram (bottom)
axSpecGram = uiaxes(fig,'Units','normalized',...
    'Position',[0.05 0.07 0.45 0.3]);
title(axSpecGram,'Spectrogram');
xlabel(axSpecGram,'Time (s)');
ylabel(axSpecGram,'Frequency (kHz)');
grid(axSpecGram,'on');

% ================= RIGHT: CONTROL PANELS (0.5–1 width) ==============
% Input panel (top right)
pnlInput = uipanel(fig,'Title','INPUT',...
    'Units','normalized',...
    'Position',[0.52 0.72 0.43 0.2]);

uilabel(pnlInput,'Text','Device:',...
    'Position',[20 130 50 20],...
    'HorizontalAlignment','left');

ddDevice = uidropdown(pnlInput,...
    'Position',[80 130 260 22],...
    'Items',{'Default device'},...
    'ItemsData',-1,...
    'Value',-1);

uilabel(pnlInput,'Text','Fs (Hz):',...
    'Position',[20 90 50 20],...
    'HorizontalAlignment','left');

fsOptions = [8000 16000 22050 44100 48000 96000];
fsStr = cellstr(string(fsOptions));
ddFs = uidropdown(pnlInput,...
    'Position',[80 90 80 22],...
    'Items',fsStr,...
    'Value','48000');

uilabel(pnlInput,'Text','Rec dur (s)',...
    'Position',[180 90 80 20],...
    'HorizontalAlignment','left');
editRecDur = uieditfield(pnlInput,'numeric',...
    'Position',[260 90 80 22],...
    'Value',5,'Limits',[0.1 60],...
    'RoundFractionalValues','off');

btnLoad = uibutton(pnlInput,'push','Text','Load WAV',...
    'Position',[20 50 150 30],...
    'ButtonPushedFcn',@onLoad);

btnRec  = uibutton(pnlInput,'push','Text','Record',...
    'Position',[350 25 100 100],...
    'ButtonPushedFcn',@onRecord);

btnPlayFull = uibutton(pnlInput,'push','Text','Play full',...
    'Position',[20 15 150 30],...
    'ButtonPushedFcn',@onPlayFull);

btnSaveFull = uibutton(pnlInput,'push','Text','Save full',...
    'Position',[180 15 150 30],...
    'ButtonPushedFcn',@onSaveFull);

btnResetAllSeg = uibutton(pnlInput,'push','Text','Reset segments',...
    'Position', [180 50 150 30],...
    'ButtonPushedFcn',@onResetAllSegments);

% Segment panel (under Input)
pnlSeg = uipanel(fig,'Title','SEGMENT',...
    'Units','normalized',...
    'Position',[0.52 0.55 0.43 0.15]);

btnSelectSeg = uibutton(pnlSeg,'push','Text','Select segment',...
    'Position',[20 65 180 30],...
    'ButtonPushedFcn',@onSelectSegment);

btnResetSeg = uibutton(pnlSeg,'push','Text','Reset to original',...
    'Position',[20 25 180 30],...
    'ButtonPushedFcn',@onResetSegment);

btnPlayOrigSeg = uibutton(pnlSeg,'push','Text','Play original seg',...
    'Position',[270 65 180 30],...
    'ButtonPushedFcn',@onPlayOrigSeg);

btnPlayProcSeg = uibutton(pnlSeg,'push','Text','Play processed seg',...
    'Position',[270 25 180 30],...
    'ButtonPushedFcn',@onPlayProcSeg);

% Processing panel
pnlProc = uipanel(fig,'Title','PROCESS',...
    'Units','normalized',...
    'Position',[0.52 0.28 0.43 0.25]);

% Fades
uilabel(pnlProc,'Text','Fade Shape',...
    'Position',[10 180 80 20],'HorizontalAlignment','left');

ddFadeType = uidropdown(pnlProc,...
    'Items',{'Linear','Raised cosine','Gaussian'},...
    'Value','Raised cosine',...
    'Position',[95 180 100 22]);

uilabel(pnlProc,'Text','Dur. (ms)',...
    'Position',[210 180 60 20],'HorizontalAlignment','left');

editFadeDur = uieditfield(pnlProc,'numeric',...
    'Position',[270 180 60 22],...
    'Value',25,'Limits',[1 500],...
    'RoundFractionalValues','on');

btnApplyFade = uibutton(pnlProc,'push','Text','Apply fades',...
    'Position',[350 180 100 22],...
    'ButtonPushedFcn',@onApplyFade);

% Freq-domain norm
uilabel(pnlProc,'Text','Freq-norm ref (Hz)',...
    'Position',[10 120 150 20],...
    'HorizontalAlignment','left');

editRefFreq = uieditfield(pnlProc,'numeric',...
    'Position',[140 120 80 22],...
    'Value',1000,...
    'Limits',[1 Inf],...
    'RoundFractionalValues','off');

btnFreqNorm = uibutton(pnlProc,'push',...
    'Text','Freq-domain normalise',...
    'Position',[250 120 200 22],...
    'ButtonPushedFcn',@onFreqNorm);

% --- Band-pass sub-section ---
uilabel(pnlProc,'Text','Band-pass (Hz):',...
    'Position',[10 80 120 20],'HorizontalAlignment','left');

% Low cut
uilabel(pnlProc,'Text','Low',...
    'Position',[130 80 30 18],'HorizontalAlignment','left');

editLowBP = uieditfield(pnlProc,'numeric',...
    'Position',[170 80 50 22],...
    'Value',100,...
    'Limits',[0 Inf],...       % we enforce Nyquist later
    'RoundFractionalValues','on');

% High cut
uilabel(pnlProc,'Text','High',...
    'Position',[250 80 35 18],'HorizontalAlignment','left');

editHighBP = uieditfield(pnlProc,'numeric',...
    'Position',[300 80 50 22],...
    'Value',8000,...
    'Limits',[0 Inf],...
    'RoundFractionalValues','off');

btnApplyBP = uibutton(pnlProc,'push','Text','Apply band-pass',...
    'Position',[30 15 120 25],...
    'ButtonPushedFcn',@onApplyBandpass);

btnExportSettings = uibutton(pnlProc,'push',...
    'Text','Export settings',...
    'Position',[175 15 120 25],...
    'ButtonPushedFcn',@onExportSettings);

btnLoadSettings = uibutton(pnlProc,'push',...
    'Text','Load settings',...
    'Position',[315 15 120 25],...
    'ButtonPushedFcn',@onLoadSettings);

% Notes panel (bottom right)
pnlNotes = uipanel(fig,'Title','Notes (saved with segment)',...
    'Units','normalized',...
    'Position',[0.52 0.11 0.43 0.15]);

notesArea = uitextarea(pnlNotes,...
    'Position',[10 40 450 70]);

btnSaveAll = uibutton(pnlNotes,'push','Text',...
    'Save ALL',...
    'Position',[180 10 100 22],...
    'ButtonPushedFcn',@onSaveAll);

% ================= STATE & DEVICE DISCOVERY ==========================
app.raw      = [];
app.Fs       = 48000;
app.fullTime = [];

app.segOrig  = [];
app.segProc  = [];
app.segTime  = [];

app.axWave     = axWave;
app.axSpec     = axSpec;
app.axSpecGram = axSpecGram;

% new fields for recording status
app.isRecording = false;
app.recObj      = [];
app.recTimer    = [];
app.recDur      = 0;

fig.UserData = app;

% Populate device dropdown
devNames = {'Default device'};
devIDs   = -1;
try
    info = audiodevinfo;
    if isstruct(info) && isfield(info,'input') && ~isempty(info.input)
        inputDevs = info.input;
        nIn = numel(inputDevs);
        devNames = cell(1,nIn);
        devIDs   = zeros(1,nIn);
        for k = 1:nIn
            devNames{k} = inputDevs(k).Name;
            devIDs(k)   = inputDevs(k).ID;
        end
    end
catch
end
ddDevice.Items     = devNames;
ddDevice.ItemsData = devIDs;
ddDevice.Value     = devIDs(1);

% ---------- helper: update sliders wrt Fs ---------------------------
    function updateBandpassSliderLimits(Fs)
        Nyq = Fs/2;

        % Clamp the numeric band-pass fields to [0, Nyquist]
        editLowBP.Limits  = [0 Nyq];
        editHighBP.Limits = [0 Nyq];

        % Sensible defaults
        defaultLow  = 100;
        defaultHigh = min(8000, 0.4*Nyq);

        if defaultLow >= Nyq
            defaultLow = max(0, Nyq-1);
        end
        if defaultHigh <= defaultLow
            defaultHigh = min(Nyq, defaultLow+1);
        end

        editLowBP.Value  = defaultLow;
        editHighBP.Value = defaultHigh;
    end

% ================== callbacks & helpers =============================

    function onLoad(~,~)
        [file,path] = uigetfile({'*.wav','WAV files (*.wav)'},'Select WAV file');
        if isequal(file,0), return; end
        [y,Fs] = audioread(fullfile(path,file));
        if size(y,2) > 1, y = mean(y,2); end
        y = y(:);

        app = fig.UserData;
        app.raw      = y;
        app.Fs       = Fs;
        app.fullTime = (0:length(y)-1)'/Fs;
        app.segOrig  = [];
        app.segProc  = [];
        app.segTime  = [];
        fig.UserData = app;

        updateBandpassSliderLimits(Fs);
        plotWaveformFull();
        clearSpectrum();
    end

    function onRecord(~,~)
        app = fig.UserData;
        if app.isRecording
            uialert(fig,'Recording is already in progress.','Info');
            return;
        end

        devID = ddDevice.Value;
        Fs = str2double(ddFs.Value);
        if isnan(Fs) || Fs <= 0
            uialert(fig,'Invalid sampling frequency.','Error'); return;
        end
        recDur = editRecDur.Value;

        try
            if devID < 0
                recObj = audiorecorder(Fs,16,1);
            else
                recObj = audiorecorder(Fs,16,1,devID);
            end
        catch ME
            uialert(fig,['Error creating recorder: ' ME.message],'Error');
            return;
        end

        % Start non-blocking recording
        record(recObj);

        % Update app state
        app.recObj      = recObj;
        app.isRecording = true;
        app.recDur      = recDur;
        app.Fs          = Fs;
        fig.UserData    = app;

        % Visual feedback on button: text + colour
        btnRec.Enable          = 'off'; % prevent re-clicks
        btnRec.BackgroundColor = [0.85 0.2 0.2];
        btnRec.Text            = sprintf('Recording\n%.1f s', recDur);

        % Start a timer to update countdown and stop recording when done
        t = timer( ...
            'ExecutionMode','fixedSpacing', ...
            'Period',0.2, ...
            'TimerFcn',@updateRecCountdown);

        app = fig.UserData;
        app.recTimer = t;
        fig.UserData = app;

        start(t);
    end

    function updateRecCountdown(~,~)
        % Timer callback: update Record button text and finish when done
        app = fig.UserData;
        if ~app.isRecording || isempty(app.recObj)
            return;
        end

        Fs   = app.Fs;
        recObj = app.recObj;
        % how many samples recorded so far
        currSamp = get(recObj,'CurrentSample');
        elapsed  = currSamp / Fs;
        remaining = max(0, app.recDur - elapsed);

        % Update button label with remaining time
        btnRec.Text = sprintf('Recording\n%.1f s', remaining);

        % If time is up, stop and collect data
        if remaining <= 0
            stop(recObj);
            y = getaudiodata(recObj);
            y = y(:);

            app.raw      = y;
            app.fullTime = (0:length(y)-1)'/Fs;
            app.segOrig  = [];
            app.segProc  = [];
            app.segTime  = [];
            app.isRecording = false;

            % Clean up timer
            if ~isempty(app.recTimer) && isvalid(app.recTimer)
                stop(app.recTimer);
                delete(app.recTimer);
            end
            app.recTimer = [];
            fig.UserData = app;

            % Reset button
            btnRec.Enable          = 'on';
            btnRec.BackgroundColor = [0.96 0.96 0.96];
            btnRec.Text            = 'Record';

            updateBandpassSliderLimits(Fs);
            plotWaveformFull();
            clearSpectrum();

            uialert(fig,'Recording finished.','Done');
        end
    end

    function onPlayFull(~,~)
        app = fig.UserData;
        if isempty(app.raw)
            uialert(fig,'No signal loaded or recorded yet.','Info'); return;
        end
        sound(app.raw,app.Fs);
    end

    function onSaveFull(~,~)
        app = fig.UserData;
        if isempty(app.raw)
            uialert(fig,'No signal loaded or recorded yet.','Info'); return;
        end
        [file,path] = uiputfile('*.wav','Save full signal as');
        if isequal(file,0), return; end
        sig  = app.raw(:);
        peak = max(abs(sig)) + eps;
        sig  = sig / peak * 0.99;
        audiowrite(fullfile(path,file),sig,app.Fs);
        uialert(fig,'Full signal saved.','Done');
    end

    function onResetAllSegments(~,~)
        app = fig.UserData;
        if isempty(app.raw)
            uialert(fig,'No signal loaded or recorded yet.','Info'); return;
        end
        app.segOrig = [];
        app.segProc = [];
        app.segTime = [];
        fig.UserData = app;
        plotWaveformFull();
        clearSpectrum();
    end

    function onSelectSegment(~,~)
        app = fig.UserData;
        if isempty(app.raw)
            uialert(fig,'Load or record a signal first.','Info'); return;
        end
        choice = uiconfirm(fig,...
            ['Click and drag a rectangle over the waveform to select the segment,' ...
            ' then double-click inside the rectangle to finish.'],...
            'Select segment',...
            'Options',{'Select','Cancel'},...
            'DefaultOption',1,'CancelOption',2);
        if ~strcmp(choice,'Select'), return; end

        hRect = drawrectangle(app.axWave,'DrawingArea','unlimited');
        wait(hRect);
        pos = hRect.Position;
        t1  = max(0,pos(1));
        t2  = min(app.fullTime(end),pos(1)+pos(3));

        idx1 = max(1,floor(t1 * app.Fs));
        idx2 = min(length(app.raw),ceil(t2 * app.Fs));
        seg = app.raw(idx1:idx2);
        segTime = (0:length(seg)-1)'/app.Fs;

        app.segOrig = seg;
        app.segProc = seg;
        app.segTime = segTime;
        fig.UserData = app;

        plotSegmentWaveforms();
        plotSpectrumComparison();
        plotSpectrogram();
    end

    function onResetSegment(~,~)
        app = fig.UserData;
        if isempty(app.segOrig)
            uialert(fig,'No segment selected yet.','Info'); return;
        end
        app.segProc = app.segOrig;
        fig.UserData = app;
        plotSegmentWaveforms();
        plotSpectrumComparison();
        plotSpectrogram();
    end

    function onPlayOrigSeg(~,~)
        app = fig.UserData;
        if isempty(app.segOrig)
            uialert(fig,'No segment selected yet.','Info'); return;
        end
        sound(app.segOrig,app.Fs);
    end

    function onPlayProcSeg(~,~)
        app = fig.UserData;
        if isempty(app.segProc)
            uialert(fig,'No processed segment available yet.','Info'); return;
        end
        sound(app.segProc,app.Fs);
    end

    function onApplyFade(~,~)
        app = fig.UserData;
        if isempty(app.segProc)
            uialert(fig,'Select a segment first.','Info'); return;
        end
        fadeType = ddFadeType.Value;
        fadeMs   = editFadeDur.Value;
        Fs       = app.Fs;

        segIn = app.segProc;
        fadeSamples = round(fadeMs*1e-3*Fs);
        fadeSamples = min(fadeSamples, floor(length(segIn)/2)-1);
        if fadeSamples <= 0
            uialert(fig,'Fade duration is too long for this short segment.','Warning');
            return;
        end

        env = ones(size(segIn));
        switch fadeType
            case 'Linear'
                fadeIn  = linspace(0,1,fadeSamples).';
                fadeOut = linspace(1,0,fadeSamples).';
            case 'Raised cosine'
                w = hann(2*fadeSamples);
                fadeIn  = w(1:fadeSamples);
                fadeOut = w(fadeSamples+1:end);
            case 'Gaussian'
                w = gausswin(2*fadeSamples,3);
                fadeIn  = w(1:fadeSamples);
                fadeOut = w(fadeSamples+1:end);
        end
        env(1:fadeSamples) = fadeIn;
        env(end-fadeSamples+1:end) = fadeOut;

        segOut = segIn .* env;
        rmsBefore = rms(segIn);
        rmsAfter  = rms(segOut) + eps;
        segOut    = segOut * (rmsBefore / rmsAfter);

        app.segProc = segOut;
        fig.UserData = app;
        plotSegmentWaveforms();
        plotSpectrumComparison();
        plotSpectrogram();
    end

    function onApplyBandpass(~,~)
        app = fig.UserData;
        if isempty(app.segProc)
            uialert(fig,'Select a segment first.','Info'); return;
        end

        Fs  = app.Fs;
        Nyq = Fs/2;

        % Read values from text inputs
        low  = editLowBP.Value;
        high = editHighBP.Value;

        % Enforce [0, Nyquist] and ordering
        if low  < 0,   low  = 0;   end
        if high > Nyq, high = Nyq; end

        if low >= high
            uialert(fig,'Low cut must be strictly below high cut.','Warning');
            return;
        end

        % Write back corrected values
        editLowBP.Value  = low;
        editHighBP.Value = high;

        segIn = app.segProc;

        [b,a] = butter(4,[low high]/Nyq,'bandpass');
        segOut = filtfilt(b,a,segIn);
        rmsBefore = rms(segIn);
        rmsAfter  = rms(segOut) + eps;
        segOut    = segOut * (rmsBefore / rmsAfter);

        app.segProc = segOut;
        fig.UserData = app;
        plotSegmentWaveforms();
        plotSpectrumComparison();
        plotSpectrogram();
    end

    function onFreqNorm(~,~)
        app = fig.UserData;
        if isempty(app.segProc)
            uialert(fig,'Select a segment first.','Info');
            return;
        end

        segIn = app.segProc;
        Fs    = app.Fs;
        fRef  = editRefFreq.Value;

        if fRef <= 0 || fRef > Fs/2
            uialert(fig, ...
                sprintf('Reference frequency must be between 0 and %.0f Hz.',Fs/2), ...
                'Warning');
            return;
        end

        % --- envelope-based peak shift + RMS match ---
        segOut = applyFreqNormEnvelope(segIn, Fs, fRef);

        app.segProc = segOut;
        fig.UserData = app;

        plotSegmentWaveforms();
        plotSpectrumComparison();
        plotSpectrogram();
    end

    function onExportSettings(~,~)
        cfg = getCurrentSettings();
        [file,path] = uiputfile('*.yml','Export settings to');
        if isequal(file,0), return; end
        writeSettingsYAML(fullfile(path,file),cfg);
        uialert(fig,'Settings exported.','Done');
    end

    function onLoadSettings(~,~)
        [file,path] = uigetfile('*.yml','Load settings from');
        if isequal(file,0), return; end
        cfg = readSettingsYAML(fullfile(path,file));
        if isfield(cfg,'fade_shape') && any(strcmp(ddFadeType.Items,cfg.fade_shape))
            ddFadeType.Value = cfg.fade_shape;
        end
        if isfield(cfg,'fade_duration_ms') && ~isnan(cfg.fade_duration_ms)
            editFadeDur.Value = max(min(cfg.fade_duration_ms,200),1);
        end
        if isfield(cfg,'band_low_Hz') && ~isnan(cfg.band_low_Hz)
            Fs   = fig.UserData.Fs;
            Nyq  = Fs/2;
            val  = max(min(cfg.band_low_Hz,Nyq),0);
            editLowBP.Value = val;
        end
        if isfield(cfg,'band_high_Hz') && ~isnan(cfg.band_high_Hz)
            Fs   = fig.UserData.Fs;
            Nyq  = Fs/2;
            val  = max(min(cfg.band_high_Hz,Nyq),0);
            editHighBP.Value = val;
        end
        if isfield(cfg,'freqnorm_ref_Hz') && ~isnan(cfg.freqnorm_ref_Hz)
            editRefFreq.Value = cfg.freqnorm_ref_Hz;
        end
        uialert(fig,'Settings loaded into controls. Apply fades/filters as needed.','Done');
    end

    function onSaveAll(~,~)
        app = fig.UserData;
        if isempty(app.segProc)
            uialert(fig,'No processed segment to save.','Info'); return;
        end
        [file,path] = uiputfile('*.yml',...
            'Choose base name (folder will be created to store everything)');
        if isequal(file,0), return; end
        [~,base,~] = fileparts(file);
        outFolder  = fullfile(path,base);
        if ~exist(outFolder,'dir'), mkdir(outFolder); end
        baseFull = fullfile(outFolder,base);

        seg  = app.segProc(:);
        peak = max(abs(seg)) + eps;
        seg  = seg / peak * 0.99;
        audiowrite([baseFull '_segment.wav'],seg,app.Fs);

        exportgraphics(app.axWave,    [baseFull '_waveform.png']);
        exportgraphics(app.axSpec,    [baseFull '_spectrum.png']);
        exportgraphics(app.axSpecGram,[baseFull '_spectrogram.png']);

        cfg = getCurrentSettings();
        writeSettingsYAML([baseFull '_settings.yml'],cfg);

        noteLines = notesArea.Value;
        fid = fopen([baseFull '_notes.txt'],'w');
        if fid ~= -1
            for k = 1:numel(noteLines), fprintf(fid,'%s\n',noteLines{k}); end
            fclose(fid);
        end
        uialert(fig,sprintf('Saved segment, plots, settings, and notes in:\n%s',outFolder),...
            'Done');
    end

% plotting & YAML helpers (unchanged)
    function plotWaveformFull
        app = fig.UserData;
        cla(app.axWave);
        plot(app.axWave,app.fullTime,app.raw);
        xlabel(app.axWave,'Time (s)');
        ylabel(app.axWave,'Amplitude');
        title(app.axWave,'Full signal');
        grid(app.axWave,'on');
    end

    function plotSegmentWaveforms
        app = fig.UserData;
        if isempty(app.segOrig), return; end
        cla(app.axWave);
        hold(app.axWave,'on');
        plot(app.axWave,app.segTime,app.segOrig,'DisplayName','Original');
        plot(app.axWave,app.segTime,app.segProc,'--','DisplayName','Processed');
        hold(app.axWave,'off');
        xlabel(app.axWave,'Time (s)');
        ylabel(app.axWave,'Amplitude');
        title(app.axWave,'Selected Segment: Original vs Processed');
        legend(app.axWave,'Location','best');
        grid(app.axWave,'on');
    end

    function plotSpectrumComparison
        app = fig.UserData;
        if isempty(app.segOrig), clearSpectrum; return; end
        Nfft = 4096;
        Fs   = app.Fs;
        f    = (0:Nfft-1)'*(Fs/Nfft);
        Yorig = fft(app.segOrig,Nfft);
        Yproc = fft(app.segProc,Nfft);
        magOrig = 20*log10(abs(Yorig)+eps);
        magProc = 20*log10(abs(Yproc)+eps);
        cla(app.axSpec);
        plot(app.axSpec,f(1:Nfft/2)/1000,magOrig(1:Nfft/2),'DisplayName','Original');
        hold(app.axSpec,'on');
        plot(app.axSpec,f(1:Nfft/2)/1000,magProc(1:Nfft/2),'DisplayName','Processed');
        hold(app.axSpec,'off');
        xlabel(app.axSpec,'Frequency (kHz)');
        ylabel(app.axSpec,'Magnitude (dB)');
        title(app.axSpec,'Spectrum Comparison');
        legend(app.axSpec,'Location','best');
        grid(app.axSpec,'on');
    end

    function plotSpectrogram
        app = fig.UserData;
        if isempty(app.segProc), cla(app.axSpecGram); return; end
        sig = app.segProc;
        Fs  = app.Fs;
        winLen  = round(0.03*Fs);
        win     = hamming(winLen);
        overlap = round(0.02*Fs);
        nfft    = 2048;
        [S,F,T] = spectrogram(sig,win,overlap,nfft,Fs,'yaxis');
        cla(app.axSpecGram);
        imagesc(app.axSpecGram,T,F/1000,20*log10(abs(S)+eps));
        axis(app.axSpecGram,'xy');
        xlabel(app.axSpecGram,'Time (s)');
        ylabel(app.axSpecGram,'Frequency (kHz)');
        title(app.axSpecGram,'Spectrogram - Processed Segment');
        colormap(app.axSpecGram,parula);
        colorbar(app.axSpecGram);
        caxis(app.axSpecGram,[-100 -20]);
    end

    function clearSpectrum
        cla(axSpec); cla(axSpecGram);
    end

    function segOut = applyFreqNormEnvelope(segIn, Fs, fRef)
        % APPLYFREQNORMENVELOPE
        % Shift the *envelope peak* to fRef while:
        %   - preserving spectral shape as much as possible
        %   - preserving original phase
        %   - matching RMS to segIn

        x   = segIn(:);
        N   = length(x);
        Nfft = 2^nextpow2(N);

        X = fft(x, Nfft);
        mag = abs(X);
        phs = angle(X);

        df   = Fs / Nfft;
        Npos = floor(Nfft/2) + 1;            % bins 1..Npos are 0..Fs/2
        magPos = mag(1:Npos);

        % ---- 1) Smooth magnitude to get spectral envelope ----
        sigmaHz   = 300;                     % envelope smoothing width in Hz
        sigmaBins = sigmaHz / df;
        winRad    = max(3, ceil(4*sigmaBins));
        g         = exp(-0.5*((-winRad:winRad)/sigmaBins).^2);
        g         = g(:) / sum(g);
        envOrig   = conv(magPos, g, 'same');

        % ---- 2) Find original envelope peak frequency ----
        [~, idxPeak] = max(envOrig);
        fPeak = (idxPeak-1) * df;

        % If peak already ~ fRef, just RMS-match and return
        if abs(fPeak - fRef) < df
            segOut = x;
            return;
        end

        % ---- 3) Shift the envelope peak to fRef (integer-bin shift) ----
        shiftBins = round((fRef - fPeak) / df);
        envShift = zeros(size(envOrig));

        if shiftBins >= 0
            srcIdx = 1:(Npos-shiftBins);
            dstIdx = (1+shiftBins):Npos;
        else
            shiftBinsAbs = -shiftBins;
            srcIdx = (1+shiftBinsAbs):Npos;
            dstIdx = 1:(Npos-shiftBinsAbs);
        end
        envShift(dstIdx) = envOrig(srcIdx);

        % ---- 4) Construct shaping ratio and limit extremes ----
        ratio = envShift ./ (envOrig + eps);
        % avoid crazy boosts/cuts
        ratio = max(min(ratio, 5), 0.2);

        % ---- 5) Apply ratio only to magnitudes, keep phase ----
        magPosNew = magPos .* ratio;
        XposNew   = magPosNew .* exp(1j*phs(1:Npos));

        % Build full conjugate-symmetric spectrum
        Xnew = zeros(Nfft,1);
        Xnew(1:Npos) = XposNew;
        % mirror 2..Npos-1
        Xnew(Npos+1:end) = conj(flipud(XposNew(2:end-1)));

        % ---- 6) IFFT and trim ----
        y = real(ifft(Xnew));
        y = y(1:N);

        % ---- 7) RMS match to original segment ----
        rmsIn  = sqrt(mean(x.^2)) + eps;
        rmsOut = sqrt(mean(y.^2)) + eps;
        y = y * (rmsIn / rmsOut);

        segOut = y;
    end

    function cfg = getCurrentSettings
        cfg.fade_shape       = ddFadeType.Value;
        cfg.fade_duration_ms = editFadeDur.Value;
        cfg.band_low_Hz      = editLowBP.Value;
        cfg.band_high_Hz     = editHighBP.Value;
        cfg.freqnorm_ref_Hz  = editRefFreq.Value;
    end

    function writeSettingsYAML(fname,cfg)
        fid = fopen(fname,'w'); if fid==-1, return; end
        fprintf(fid,"# Speech stimulus processing settings\n");
        fprintf(fid,"fade_shape: %s\n",cfg.fade_shape);
        fprintf(fid,"fade_duration_ms: %.3f\n",cfg.fade_duration_ms);
        fprintf(fid,"band_low_Hz: %.3f\n",cfg.band_low_Hz);
        fprintf(fid,"band_high_Hz: %.3f\n",cfg.band_high_Hz);
        fprintf(fid,"freqnorm_ref_Hz: %.3f\n",cfg.freqnorm_ref_Hz);
        fclose(fid);
    end

    function cfg = readSettingsYAML(fname)
        cfg = struct();
        fid = fopen(fname,'r'); if fid==-1, return; end
        while true
            line = fgetl(fid); if ~ischar(line), break; end
            line = strtrim(line);
            if isempty(line) || startsWith(line,'#'), continue; end
            parts = split(line,':'); if numel(parts)<2, continue; end
            key = strtrim(parts{1});
            valStr = strtrim(strjoin(parts(2:end),':'));
            switch key
                case 'fade_shape',        cfg.fade_shape        = valStr;
                case 'fade_duration_ms',  cfg.fade_duration_ms  = str2double(valStr);
                case 'band_low_Hz',       cfg.band_low_Hz       = str2double(valStr);
                case 'band_high_Hz',      cfg.band_high_Hz      = str2double(valStr);
                case 'freqnorm_ref_Hz',   cfg.freqnorm_ref_Hz   = str2double(valStr);
            end
        end
        fclose(fid);
    end
end