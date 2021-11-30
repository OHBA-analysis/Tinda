function camcan_preproc_CHedit( configfile, sessions_to_run )
%% function osl_preproc( configfile )
%
% configfile is a path to a yaml preprocessing config
%
%

%% Check compiled environment was good
% If this fails, something has gone wrong with osldir.m or the compilation
% itself

try osldir
    fprintf('OSL Environment Detected\n');
catch ME
    warning('Cannot find osldir in compiled environment');
    warning('\tlikely that compilation done without OSL initialised');
    error('OSL Environment Missing - exiting');
end

% Check that layouts are in place
if ~exist( fullfile(osldir,'layouts','neuromag306mag.lay'), 'file')
    warning( '%s', fullfile(osldir,'layouts','neuromag306mag.lay') );
    warning( 'Layout files not found in compiled environment')
    error('Stored data missing from compiled environment - exiting');
end

%% Check files and output actually exist

if ~exist( configfile, 'file' )
    error('config file list not found!');
end
addpath(genpath('/home/chiggins/Documents/MATLAB/yamltoolbox')) % toolbox must be downloaded from https://gitlab.umich.edu/lsa-ts-rsp/yaml-toolbox.git
config = yaml.ReadYaml( configfile );

%% If asked for a statfile, make it

if isfield(config,'statfile')
    stat_fid = fopen(config.statfile,'w');
    fprintf('statfile created: %s\n', config.statfile);
end

%% If asked for a logfile, make it

if isfield(config,'logfile')
    if ~exist(config.logfile)
        mkdir(fileparts(config.logfile));
        fid = fopen(config.logfile,'w+');
        fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf(fid,['New log being created: ',datestr(datetime('now')),'\n\n\n\n']);
        fclose(fid);
    end
    diary(config.logfile)
    diary on;
    fprintf('logfile created: %s\n', config.logfile);
end


%% Check output directory is specified and exists

if ~isfield( config, 'outbase' )
    error('Output directory not specified in config file');
end

if ~exist( config.outbase, 'dir' )
    error('output directory not found!');
end



%% Read and sanity check files

nfiles = length(config.datafiles);
all_good = true;
has_fif = true(nfiles,1);
has_fmri= true(nfiles,1);
for ii = 1:nfiles
    if isfield(config.datafiles{ii}, 'fif')
        if ~exist( config.datafiles{ii}.fif, 'file')
            warning('Input fif file - %s not found',config.datafiles{ii}.fif);
            %all_good = false;
            has_fif(ii) = false;
        end
    end
    if isfield(config.datafiles{ii}, 'spm')
        if ~exist( config.datafiles{ii}.spm, 'file')
            warning('Input spm file - %s not found',config.datafiles{ii}.spm);
            %all_good = false;
        end
    end
    if isfield(config.datafiles{ii}, 'mri')
        if ~exist( config.datafiles{ii}.mri, 'file')
            warning('Input mri file - %s not found',config.datafiles{ii}.mri);
            %all_good = false;
            has_fmri(ii) = false;
        elseif strcmp(config.datafiles{ii}.mri(end-1:end),'gz') && ~exist(config.datafiles{ii}.mri(1:end-3), 'file')
            gunzip(config.datafiles{ii}.mri);
        end
    end
end

if all_good == false
    error('Some files not found - exiting');
end

outbase = config.outbase;
fprintf('output will be saved to: %s\n\n', outbase);

%% Main loop

if nargin < 2 || isempty(sessions_to_run)
    sessions_to_run = 1:nfiles;
end
sessions_good = find(has_fif & has_fmri);
sessions_to_run = intersect(sessions_to_run,sessions_good');

nstages = length(config.preproc);
coregdir = '/ohba/pi/mwoolrich/datasets/CamCan_2021/MEGRS_preprocessed/references/coreg/';
mkdir(coregdir);

for ii = sessions_to_run

    fprintf('Processing file %d/%d\n', ii, nfiles);
    fprintf('\t%s\n',config.datafiles{ii}.name);

    % Start up logfile lines if requested
    if isfield(config,'statfile')
        header = 'run, fname, infile, outfile, status ';
        logline = strcat(num2str(ii), ', ', config.datafiles{ii}.name, ', ');
        if strcmp(config.preproc{1}.method, 'import')
            logline = strcat(logline, config.datafiles{ii}.fif, ',');
        elseif strcmp(  config.preproc{1}.method, 'load')
            logline = strcat(logline, config.datafiles{ii}.spm, ',');
        end
    end
    try
        for jj = 1:nstages

            if strcmp(config.preproc{jj}.method, 'import')
                % Import data from fif to spm in output dir
                %  - {method: import}

                % Get raw data file name and make spm output file name
                [~,fname,~] = fileparts(config.datafiles{ii}.fif);
                outfile = [outbase fname];

                D = osl_import(config.datafiles{ii}.fif, 'outfile', outfile );

            elseif strcmp(  config.preproc{jj}.method, 'load')
                % Load data from existing spm and copy to output dir
                %  - {method: load}
                D = spm_eeg_load( config.datafiles{ii}.spm );

                % Get raw data file name and make spm output file name
                %[~,fname,~] = fileparts(config.datafiles{ii}.spm);
                fname = [ config.datafiles{ii}.name '_rest_raw.mat' ];
                outfile = [outbase fname];

                S = [];
                S.D = D;
                S.outfile = outfile;
                D = spm_eeg_copy(S);

            elseif strcmp( config.preproc{jj}.method, 'coregister')
                % TODO
                %

%                 % Remove nose points
%                 pnt = D.fiducials.pnt;
%                 nas = D.fiducials.fid.pnt(2,:);
%                 nas(3) = nas(3) - 40; % drop nasion by ~4cm
%                 distances = sqrt((nas(1)-pnt(:,1)).^2 + (nas(2)-pnt(:,2)).^2 + (nas(3)-pnt(:,3)).^2 );
%                 keeps = distances>60; % remove points within 6cm
%                 fids = D.fiducials;
%                 fids.pnt = fids.pnt(keeps,:);
%                 D = fiducials(D,fids);
%                 D.save

                % get nasal points:
%                 fname = config.datafiles{ii}.fiducials; 
%                 fid = fopen(fname); 
%                 raw = fread(fid,inf); 
%                 str = char(raw'); 
%                 fclose(fid); 
%                 val = jsondecode(str);
%                 
%                 
%                 coreg_settings = struct;
%                 coreg_settings.D = D;
%                 coreg_settings.mri = config.datafiles{ii}.mri;
%                 coreg_settings.fid_mnicoords = 1;
%                 coreg_settings.useheadshape = false;
%                 coreg_settings.forward_meg = 'Single Shell';
%                 %coreg_settings.forward_meg = 'MEG Local Spheres';
%                 coreg_settings.use_rhino = false;
%                 coreg_settings.fid.label.nasion='Nasion';
%                 coreg_settings.fid.label.lpa='LPA';
%                 coreg_settings.fid.label.rpa='RPA';
%                 coreg_settings.fid.mnicoords.lpa = val.HeadCoilCoordinates.LPA;
%                 coreg_settings.fid.mnicoords.rpa = val.HeadCoilCoordinates.RPA;
%                 coreg_settings.fid.mnicoords.nasion = val.HeadCoilCoordinates.NAS;
%                 D = osl_headmodel(coreg_settings);
%                 coregfig = gcf;
%                 print([coregdir,config.datafiles{ii}.name],'-dpng');
                %D.save();
                mrifile = config.datafiles{ii}.mri;
                mrifile = mrifile(1:end-3)
                matlabbatch{1}.spm.meeg.source.headmodel.D = {D.fullfile};%{'/home/chiggins/Documents/MATLAB/osl/spmeeg_mf2pt2_sub-CC110033_ses-rest_task-rest_megtransdef.mat'};
                matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
                matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
                matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.mri = {mrifile};%{'/ohba/pi/mwoolrich/datasets/CamCan_2021/cc700/mri/pipeline/release004/BIDS_20190411/anat/sub-CC110033/anat/sub-CC110033_T1w.nii,1'};
                matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 1;
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'Nasion';
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'LPA';
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'RPA';
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
                matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
                matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
                spm_jobman('run',matlabbatch)
                D=spm_eeg_load(D.fullfile);
                spm_eeg_inv_checkforward(D);
                coregfig = gcf;
                print([coregdir,config.datafiles{ii}.name],'-dpng');
                close all;
            elseif strcmp( config.preproc{jj}.method, 'downsample')
                % Run downsampling
                %   - {method: downsample, new_fsample: 250}
                D = spm_eeg_downsample(struct('D',D,'fsample_new',config.preproc{jj}.new_fsample));

            elseif strcmp( config.preproc{jj}.method, 'bandpassfilter')
                % Run bandpass filter
                %   - {method: bandpassfilter, freqs: 1 90}
                freqs = split_list( config.preproc{jj}.freqs, 'num' );
                D = osl_filter(D,freqs,'prefix','');

            elseif strcmp( config.preproc{jj}.method, 'bandstopfilter')
                % Run bandstop filtre
                %   - {method: bandstopfilter, freqs: 48 52}
                freqs = split_list( config.preproc{jj}.freqs, 'num' );
                D = osl_filter(D,-freqs,'prefix','');

            elseif strcmp( config.preproc{jj}.method, 'detect_artefacts')
                % Detect artefacts automatically
                %  - {method: detect_artefacts, badtimes: true, badchannels: true, modalities: MEGMAG MEGPLANAR}

                modalities = get_opt( config.preproc{jj}, 'modalities', {'MEGMAG','MEGPLANAR'} );
                if isa( modalities, 'char' )
                    modalities = split_list( modalities, 'str' );
                end
                badchannels =  get_opt( config.preproc{jj}, 'badchannels', true );
                badtimes = get_opt( config.preproc{jj}, 'badtimes', true );

                D = osl_detect_artefacts(D, 'badchannels',badchannels,...
                    'badtimes',badtimes,...
                    'modalities',modalities);

            elseif strcmp( config.preproc{jj}.method, 'africa')
                % Run ICA denoising
                %   - {method: africa, usedmaxfilter: true, artefact_channels: [EOG, ECG], auto_artefact_chans_corr_thresh: .35}
                artefact_channels = get_opt( config.preproc{jj}, 'artefact_channels', {'EOG','ECG'} );
                if isa( artefact_channels, 'char' )
                    artefact_channels = split_list( modalities, 'str' );
                end

                used_maxfilter =  get_opt( config.preproc{jj}, 'used_maxfilter', true );
                auto_artefact_chans_corr_thresh =  get_opt( config.preproc{jj}, 'auto_artefact_chans_corr_thresh', .5 );

                D = osl_africa(D, 'used_maxfilter',used_maxfilter,...
                    'auto_artefact_chans_corr_thresh',auto_artefact_chans_corr_thresh,...
                    'artefact_channels',artefact_channels);

            elseif strcmp( config.preproc{jj}.method, 'sensor_norm' )
                % Normalise sensor modalities
                %  - {method: sensor_norm, normalise_method: min_eig, modalities: MEGMAG MEGPLANAR}

                % Process options
                modalities = get_opt( config.preproc{jj}, 'modalities', {'MEGMAG','MEGPLANAR'} );
                if isa( modalities, 'char' )
                    modalities = split_list( modalities, 'str' );
                end
                normalise_method =  get_opt( config.preproc{jj}, 'normalise_method', 'min_eig' );

                % Normalise sensor types
                S = [];
                S.D = D;
                S.modalities = modalities;
                S.do_plots = 0;
                S.samples2use = good_samples(D,D.indchantype(S.modalities,'GOOD'));
                S.trials = 1;
                S.pca_dim = 99;
                S.force_pca_dim = 0;
                S.normalise_method = normalise_method;
                D = normalise_sensor_data( S );

            elseif strcmp( config.preproc{jj}.method, 'lcmv_parcellation' )
                % Run LCMV beamformer
                %  - {method: lcmv, parcellation: fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm, pca_order: 120}

                %parc = get_opt(config.preproc{jj}, 'parcellation', 'fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm' );
                parc = get_opt(config.preproc{jj}, 'parcellation', 'fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm' );
                pca_order = get_opt( config.preproc{jj}, 'pca_order', 120);

                p = parcellation( parc );
                D = osl_inverse_model(D,p.template_coordinates,'pca_order',pca_order);

                D = ROInets.get_node_tcs(D,p.parcelflag,'spatialBasis','parcellation');

            else
                error( 'method %s not recognised', config.preproc{jj}.method);
            end
        end

        if isfield(config,'statfile')
            [out_line,out_header] = get_statfile_line( D, ii, {'MEGMAG','MEGPLANAR'} );
            logline = strcat(logline, D.fullfile, ', completed');

            if ii == 1
                fprintf(stat_fid,strcat(header,out_header));
            end
            fprintf(stat_fid,strcat(logline,out_line));
        end

    catch ME
        warning('Preproc FAILED on stage: %s', config.preproc{jj}.method);
        %warning(ME.identifier);
        warning(ME.message);
        out_line = strcat(', FAILED-', config.preproc{jj}.method ,'\n');

        fprintf(stat_fid,strcat(logline,out_line));

    end

    D.save % crucial
end


if isfield(config,'statfile')
    fclose(stat_fid);
    fprintf('statfile: %s\n', config.statfile);
end

fprintf('Preprocessing complete - exiting\n');
diary off;

end


function out = get_opt( stage, field, default )
%% function out = get_opt( stage, field, default )
%
% function for extracting an option from a field in a struct and returning a default
% value if the field is not specified

if isfield(stage,field)
    out = stage.(field);
else
    out = default;
end

end

function out = split_list( freqs, outtype )
%% function out = split_list( freqs, outtype )
%
% function taking a string containing two space-separated numbers and
% retuning a cell array of two numbers

freqs = strsplit(freqs,' ');
if strcmp(outtype,'num')
    out = [str2double(freqs{1}), str2double(freqs{2})];
else
    out = {freqs{1} freqs{2}};
end
end

%% Logfile functions %%
%

function [out,header] = get_statfile_line( D, jj, modalities )
%% function [out,header] = get_statfile_line( D, jj, modalities )
%

header = '';
out = '';
badchans = zeros(numel(modalities),1);
for m = 1:length(modalities)
    badchans(m) = length(intersect( D.badchannels, D.indchantype(modalities{m})));
    header = strcat( header, ', badchan-', modalities{m} );
    out = strcat( out, ', ', num2str(badchans(m)) );
end

baddurs = zeros(numel(modalities),1);
for m = 1:length(modalities)
    baddurs(m) = get_badtimes_duration( D, modalities{m} );
    header = strcat( header, ', badtime-', modalities{m} );
    out = strcat( out, ', ', num2str(baddurs(m)) );
end

% Get names of auto africa metrics
all_metric_names = fieldnames(D.ica.metrics);
inds = [];
for ii = 1:length(all_metric_names)
    if strfind(all_metric_names{ii},'corr_chan')
        inds(end+1) = ii;
    end
end
metric_names = all_metric_names(inds);
auto_reasons = cellfun(@num2str,D.ica.auto_reason,'UniformOutput',false);

for ii = 1:length(metric_names)
    header = strcat(header, ', IC-', metric_names{ii});
    inds = cellfun(@(x) strfind(x,metric_names{ii}), auto_reasons, 'UniformOutput', false);
    out = strcat(out, ', ', num2str(sum(cellfun(@isempty, inds ) == 0)) );
end

header = [header '\n'];
out = [out '\n'];

end

function out = get_badtimes_duration( D, modality )
%% function out = get_badtimes_duration( D, modality )
%

ev = D.events;
if isempty(ev)
    %fprintf('Bad times - none\n');
    out = 0;
else
    ev = ev(cellfun(@(x) ~isempty(strcmp('artefact',x)),{ev.type})); % Find only artefact events
    this_modality = strcmp({ev.value},modality);
    out = sum([ev(this_modality).duration]);
end
end
