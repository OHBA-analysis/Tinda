% mvar_02_hcp_aal_beamform
% Updated by Mark Hymers and Andrew Quinn
% Based on an original script by Giles Colclough 24 Febouary 2015

% This function assumes that OSL has been properly added to the
% path and set up

function mvar_02_hcp_aal_beamform_CHedit(subID)

% Basic configuration
gridStep = 8;
reducedFs = 240;
recon_Fbp = [1 80];

% Directories - load from MVAR config; assumes we're running from the
% repository
conf = mvarconfig('mvarconfig.ini')

parcelInds = load(conf.FILE_AALINDS);

ROInets.make_directory(conf.DIR_FIGURES);
ROInets.make_directory(conf.DIR_VEMISC);
ROInets.make_directory(conf.DIR_RAWVEDATA);

mni_brain_2mm = fullfile(conf.DIR_EXTERNAL, 'osl/std_masks/MNI152_T1_2mm_brain.nii.gz');

sourcemodel_spec = ['MEG_anatomy_sourcemodel_3d' num2str(gridStep) 'mm'];

megFileList    = {};
headFileList   = {};
sourceFileList = {};
noiseFileList  = {};
reconFileList  = {};

conf.DIR_HCP = '/ohba/pi/mwoolrich/datasets/HCP_RAW/HCP_DATA/'
conf.DIR_VEMISC = '/ohba/pi/mwoolrich/datasets/HCP_CH/preprocessed'
conf.DIR_RAWVEDATA = '/ohba/pi/mwoolrich/datasets/HCP_CH/ve_output_rest/'
mkdir(conf.DIR_RAWVEDATA)

%megDir       = fullfile(conf.DIR_HCP, subID,  '/MEG/Resting/HCP_preprocessed/rmegpreproc/')
megDir       = fullfile(conf.DIR_HCP, subID,  '/MEG/Restin/rmegpreproc/')
anatomyDir   = fullfile(conf.DIR_HCP, subID,  '/MEG/anatomy/');
rNoiseDir    = fullfile(conf.DIR_HCP, subID, '/MEG/Noise/2_Pnoise/4D/');
%sessionList1  = osl_filelist(megDir, '*.mat');
% do manually:
sessionList = {};
for j=3:5
    sessionList = {sessionList{:},fullfile(conf.DIR_HCP, subID,'/MEG/Restin/rmegpreproc/',[subID,'_MEG_',int2str(j),'-Restin_rmegpreproc.mat'])};
end
sessionList = sessionList';

megFileList    = [sessionList]
headFileList   = [repmat({fullfile(anatomyDir, [subID, '_MEG_anatomy_headmodel.mat'])}, length(sessionList), 1)];
sourceFileList = [repmat({fullfile(anatomyDir, [subID, '_' sourcemodel_spec '.mat'])},  length(sessionList), 1)];
noiseFileList  = [repmat({fullfile(rNoiseDir,  'c,rfDC')},                                  length(sessionList), 1)];
reconFileList  = [strrep(sessionList, megDir, [conf.DIR_VEMISC, '/'])];
% Check files exist:
missing_files = [];
missing_files = [missing_files, megFileList(cellfun(   @exist, megFileList)    ~= 2)];
missing_files = [missing_files, headFileList(cellfun(  @exist, headFileList)   ~= 2)];
%missing_files = [missing_files, noiseFileList(cellfun( @exist, noiseFileList)  ~= 2)];
missing_files = [missing_files, sourceFileList(cellfun(@exist, sourceFileList) ~= 2)];
if ~isempty(missing_files)
    disp(missing_files)
    error('Files are missing!');
else
    disp('Good news! All files are valid!');
end%if

% Perform the beamforming
nSessions = length(megFileList);

for iSession = 1:nSessions,
    fprintf('\n\n%s: Source recon for session %d out of %d. #############################\n\n\n', ...
            mfilename, iSession, nSessions);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Load the data
    load(megFileList{iSession});    % loads data
    load(headFileList{iSession});   % loads headmodel
    load(sourceFileList{iSession}); % loads sourcemodel3d

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Prepare the data

    % convert to mm
    grad          = ft_convert_units(data.grad, 'mm');
    headmodel     = ft_convert_units(headmodel,'mm');
    sourcemodel3d = ft_convert_units(sourcemodel3d,'mm');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute the leadfields

    reduceLFRank    = 2;
    cfg             = [];
    cfg.headmodel   = headmodel;
    cfg.grid        = sourcemodel3d;
    cfg.grad        = grad;
    cfg.channel     = data.label;
    cfg.normalize   = 'yes'; % normalising lead fields allows us to correct for depth biases.
    cfg.reducerank  = reduceLFRank;

    % prepare leadfields
    gridLF = ft_prepare_leadfield(cfg);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Load noise
%     fprintf('Loading empty room noise data. \n');
%     % Load Raw Noise data
%     % Split Room Noise in pseudo trials of 2 seconds
%     clear pseudoTrl
%     noiseFile      = noiseFileList{iSession};
%     noiseHdr       = ft_read_header(noiseFile);
%     pseudoDur      = 2;                        % s pseudo Trial Duration
%     pseudoDurSamps = ceil(pseudoDur * noiseHdr.Fs);
%     pseudoTrl      = (1:pseudoDurSamps:(noiseHdr.nSamples-pseudoDurSamps))';
%     pseudoTrl      = [pseudoTrl pseudoTrl + (pseudoDurSamps-1)];
%     pseudoTrl      = [pseudoTrl zeros(size(pseudoTrl,1),1)];
% 
%     cfg           = [];
%     cfg.channel   = 'MEG';
%     cfg.dataset   = noiseFile;
%     cfg.trl       = pseudoTrl;
%     origNoiseData = ft_preprocessing(cfg);
% 
%     % Detrend each pseudo trial just in case there is any problematic channel with any
%     % significant slow increasing or decreasing trends. Should make little
%     % difference for channels with no problems
%     cfg         = [];
%     cfg.detrend = 'yes';
%     noiseData   = ft_preprocessing(cfg, origNoiseData);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Reshape data
    allData    = cat(2, data.trial{:});
    time       = (1:size(allData,2)) ./ data.fsample; % not real time, due to bad segments having been removed
    Fs         = data.fsample;
%     allNoise   = cat(2, noiseData.trial{:});
% 
%     % downsample noise to match data, including anti-aliasing filter
%     clear dsNoise
%     dsfactor   = round(noiseData.fsample ./ Fs);
% 
%     for i = ROInets.rows(allNoise):-1:1,
%         dsNoise(i,:) = decimate(allNoise(i,:), dsfactor);
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Bandpass filter before source localization
    fprintf('Bandpass filtering. \n');

    % use wide-band reconstruction
    N   = 5;
    filteredData  = ft_preproc_bandpassfilter(allData, data.fsample, recon_Fbp, N);
 %   filteredNoise = ft_preproc_bandpassfilter(dsNoise, Fs, recon_Fbp, N);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute covariance
    fprintf('Computing covariance. \n');
    %usedChans = ismember(noiseData.label, data.label);
    usedChans = ismember(data.label, data.label);
    C         = cov(filteredData.');
    %Ce        = cov(filteredNoise(usedChans,:).');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Write out SPM MEEG object to store the results in
    fprintf('Converting to SPM\n');

    sourcedata           = data;
    sourcedata.trial     = {filteredData};
    sourcedata.trialinfo = [nan 0];
    sourcedata.time      = {time}; % not real time, due to bad segments having been removed
    D = spm_eeg_ft2spm(sourcedata, reconFileList{iSession});

    D = osl_spmfun(@spm_eeg_downsample, struct('D', D, 'fsample_new', reducedFs));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute localization weights
    fprintf('Beamforming. \n');

    % reduce dimensionality
    icaResDir = fullfile(conf.DIR_HCP, subID, '/MEG/Restin/icaclass/');
    [~, MEGfileName] = fileparts(megFileList{iSession});
    icaClassFile = strrep(MEGfileName, 'rmegpreproc', 'icaclass.mat');
    if isfile(icaClassFile)
        tmp = load(fullfile(icaResDir, icaClassFile));
        nICAremoved = numel(tmp.comp_class.class.ecg_eog_ic);
    else
        nICAremoved = 0;
    end
    % Don't use spm_pca_order as broke - use Mark's method instead. This gets used automatically.
    dimensionality   = D.nchannels - nICAremoved - 5;
    [invC, ~] = pinv_plus(C, dimensionality);

    % compute weights

    % Figure out which nodes are inside the head model
    indices = find(sourcemodel3d.inside);
    nInside = length(indices);

    for i = nInside:-1:1,
        idx = indices(i);
        % leadfield for this dipole (already normalised)
        lf  = gridLF.leadfield{idx};

        % 3d source variances
        tmp = lf' * invC *lf;
        nn  = size(lf,2);

        % collapse to scalar
        [eta,~] = svds(real(pinv_plus(tmp(1:nn, 1:nn), reduceLFRank, 0)), 1);
        lf      = lf * eta;

        % LCMV weights
        W_nonorm{i} = lf'*invC/(lf' * invC * lf);
        W_mag(i)   = norm(W_nonorm{i});

        % pseudo-z weights are normalised by the power of the projected
        % sensor noise.
        %Z_mag(i) = W_nonorm{i} * Ce * W_nonorm{i}';
        Z_mag(i) = W_nonorm{i} * C * W_nonorm{i}';
        
        W{i} = W_nonorm{i} ./ Z_mag(i);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Montage the spm object
    fprintf('Adding weights in with a montage. \n');

    useWeights = {W_nonorm, W};
    montageName = {'sources without weights normalisation', ...
                   'sources with weights normalisation'};
    for iMontage = 1:2,
        % Set up SPM montage
        WW = useWeights{iMontage};
        D = D.montage('switch');
        D.save;

        montage = [];
        montage.labelorg = D.chanlabels;
        montage.labelnew = cell(numel(WW),1);

        % tra matrix
        montage.tra = zeros(numel(WW), size(WW{1}(1,:),2));
        for i = 1:numel(WW)
            montage.labelnew{i} = mat2str(sourcemodel3d.pos(i,:),3);
            montage.tra(i,:)    = WW{i};
        end%for

        montage.name = montageName{iMontage};
        montage.chantypenew = repmat({'VE'},   length(montage.labelnew), 1);
        montage.chanunitnew = repmat({'nA*m'}, length(montage.labelnew), 1);
        montage.chantypeorg = chantype(D, D.indchannel(montage.labelorg))';
        montage.chanunitorg = units(D,    D.indchannel(montage.labelorg))';

        % Online montage needs additional channel information
        for ch = 1:length(montage.labelnew)
            montage.channels(ch).label = montage.labelnew{ch};
            montage.channels(ch).type  = montage.chantypenew{ch};
            montage.channels(ch).units = montage.chanunitnew{ch};
            montage.channels(ch).bad   = 0;
        end

        S = [];
        S.D            = fullfile(D.path, D.fname);
        S.montage      = montage;
        S.keepsensors  = false;
        S.keepothers   = false;
        S.mode         = 'switch';

        D = spm_eeg_montage(S);
        D.save;
    end%for

    clear WW useWeights

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save variance maps
    fprintf('Computing source variance. \n');

    v = osl_source_variance(D);

    % declare empty volume of same size as grid
    vol = zeros(length(sourcemodel3d.cfg.grid.template.xgrid), ...
                length(sourcemodel3d.cfg.grid.template.ygrid), ...
                length(sourcemodel3d.cfg.grid.template.zgrid));

    % fill in volume with source variance at correct locations
    for vox = sourcemodel3d.cfg.grid.template.inside;
        vol(sourcemodel3d.cfg.grid.template.pos(vox,1) == sourcemodel3d.cfg.grid.template.xgrid,...
            sourcemodel3d.cfg.grid.template.pos(vox,2) == sourcemodel3d.cfg.grid.template.ygrid,...
            sourcemodel3d.cfg.grid.template.pos(vox,3) == sourcemodel3d.cfg.grid.template.zgrid) = ...
            v(vox == sourcemodel3d.cfg.grid.template.inside);
    end

    % perform interpolation
    Vq = hcp_resample_to_mni(vol,round(sourcemodel3d.cfg.grid.template.xgrid * 10), ...
                                 round(sourcemodel3d.cfg.grid.template.ygrid * 10), ...
                                 round(sourcemodel3d.cfg.grid.template.zgrid * 10), ...
                                 gridStep);

    % output nifti
%     niiFileName     = fullfile(conf.DIR_VEMISC, ...
%                           strrep(D.fname, '.mat', '_source_variance.nii.gz'));
%     niiFileName_2mm = fullfile(conf.DIR_VEMISC, ...
%                           strrep(D.fname, '.mat', '_source_variance_ds2mm.nii.gz'));
%     save_avw(Vq, niiFileName, 'f', [gridStep gridStep gridStep 1]);
%     nii.resample_flirt(niiFileName, niiFileName_2mm, 2, 'trilinear', mni_brain_2mm);

    % fslview(niiFileName_2mm);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save BF time series
    fprintf('Computing ROI time-series. \n');

    % Load parcellation
    p = parcellation(conf.FILE_PARCELFILE, conf.FILE_PARCELLABELS);

    % We avoid spm_eeg_copy because on various systems it results
    % in a broken set of files - answers on a postcard please
    bname = D.fname;
    bname = bname(1:end-4);
    bnameout = ['parcel', bname];

    % Copy the MAT file
    copyfile(fullfile(D.path, [bname, '.mat']),...
             fullfile(D.path, [bnameout, '.mat']));

    % And the DAT file
    copyfile(fullfile(D.path, [bname, '.dat']),...
             fullfile(D.path, [bnameout, '.dat']));

    Dparcel = spm_eeg_load(fullfile(D.path, [bnameout]));
    Dparcel = Dparcel.montage('switch', 2);
    Dparcel = ROInets.get_node_tcs(Dparcel, p.parcelflag, 'spatialBasis', 'Giles');
    Dparcel.save;

    % Save in HDF5 for export to python
    % Note that we transpose the data because MATLAB writes it out
    % transposed and we want it as (rois, timepts) in the HDF5 file
    data = Dparcel(:, :, :)';
    sample_rate = Dparcel.fsample;
    [~, fname, ~] = fileparts(D.fname);
    outfname = strcat(fullfile(conf.DIR_RAWVEDATA, prefix(fname, 'roinet')), '.mathdf5');
    save(outfname, 'data', 'sample_rate', '-v7.3');
    clear outfname data sample_rate;
end
