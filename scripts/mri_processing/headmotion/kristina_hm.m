addpath('/projects/kg98/Sid/STAGES/STAGES_fmriprep/headmotion')

projectdir = '/home/scho0011/kg98/Sid/STAGES/STAGES_fmriprep/';
fmriprepworkdir = [projectdir,'data/fmriprep_work/'];
%change ses as required
fileID = fopen([projectdir,'list_fmri_ses-4.txt']);
ParticipantIDs = textscan(fileID,'%s');
ParticipantIDs = ParticipantIDs{1};

% compute numsubs
numSubs = length(ParticipantIDs);


% ------------------------------------------------------------------------------
% Containers
% ------------------------------------------------------------------------------
motionparams = cell(numSubs,1);

fdJenk = cell(numSubs,1);
fdJenk_mean = zeros(numSubs,1);

for i = 1:numSubs
    subject = ParticipantIDs{i};
    subjectnumber = subject(5:end);
    %change ses on next line
    motiondir = [fmriprepworkdir, subject, '/fmriprep_wf/single_subject_',subjectnumber,'_wf/func_preproc_ses_4_task_rest_wf/bold_hmc_wf/mcflirt/'];
    files = dir(motiondir);
    filesstruct = struct2cell(files);
    files = filesstruct(1,:);
    findname = strfind(files,'.nii.gz.par');
    idx_motionfile = find(~cellfun(@isempty,findname));
    filename = files{idx_motionfile};
    motionparams{i} = dlmread([motiondir,filename]);
    motionparams{i} = motionparams{i}(:,[4:6,1:3]);
    
    numVols = size(motionparams{i},1);
    
    % Get FD %kristna/power=50, lindein = 80? 
	fdJenk{i} = GetFDJenk(motionparams{i},50);
	% Calculate mean
	fdJenk_mean(i) = mean(fdJenk{i});
    
    % ------------------------------------------------------------------------------
    % Initial, gross movement exclusion
    % ------------------------------------------------------------------------------
	% 1) Exclude on mean rms displacement
	% Calculate whether subject has suprathreshold mean movement
	% If the mean of displacement is greater than 0.55 mm (Sattethwaite), then exclude
    % This value is halved for multiband data
	if fdJenk_mean(i) > 0.55
		exclude(i,1) = 1;
	else
		exclude(i,1) = 0;
    end	

    % ------------------------------------------------------------------------------
	% Stringent, multi criteria exclusion
	% ------------------------------------------------------------------------------
    % If the mean of displacement is greater than 0.2 mm (Ciric), then exclude
	if fdJenk_mean(i) > 0.2
		mean_exclude(i,1) = 1;
	else
		mean_exclude(i,1) = 0;
	end	
    
	% Calculate whether subject has >20% suprathreshold spikes
    fdThr = 0.25; 
	fdJenkThrPerc = round(numVols * 0.20);
	% If the number of volumes that exceed fdThr are greater than %20, then exclude
	if sum(fdJenk{i} > fdThr) > fdJenkThrPerc
		sum_exclude(i,1) = 1;
	else
		sum_exclude(i,1) = 0;
    end
    
    % 3) Exclude on large spikes (>5mm)
    % This value is halved for multiband data
	if any(fdJenk{i} > 5)
		spike_exclude(i,1) = 1;
	else
		spike_exclude(i,1) = 0;
    end
     
     % If any of the above criteria is true of subject i, mark for exclusion
    if mean_exclude(i,1) == 1 | sum_exclude(i,1) == 1 | spike_exclude(i,1) == 1
        exclude(i,2) = 1;
    else
        exclude(i,2) = 0;
    end
    
    % threshold for exclusion in minutes
	thresh = 4;
    TR = 2; 
	spikereg = GetSpikeRegressors(fdJenk{i},fdThr); % Spike regression exclusion
	numCVols = numVols - size(spikereg,2); % number of volumes - number of spike regressors (columns)
	NTime = (numCVols * TR)/60; % Compute length, in minutes, of time series data left after censoring
	if NTime < thresh;
		censoring_exclude(i,1) = 1;
	else
		censoring_exclude(i,1) = 0;
    end
    
    
end
    
T = table(ParticipantIDs, motionparams, fdJenk, fdJenk_mean, exclude, mean_exclude, sum_exclude, spike_exclude, censoring_exclude, exclude(:,1), exclude(:,2));
T.Properties.VariableNames = {'ParticipantIDs'  'motion_params' 'fdJenk'  'fdJenk_mean'  'exclude'  'mean_exclude'  'sum_exclude'  'spike_exclude'  'censoring_exclude'  'grossmvmt_exclude'  'stringent_exclude'};
    
