
datadir = '/projects/kg98/Sid/STAGES/STAGES_fmriprep/headmotion/hm_prams'
cd(datadir)
fileID = fopen('list_fmri_ses-1.txt');
ParticipantIDs = textscan(fileID,'%s');
ParticipantIDs = ParticipantIDs{1};

% compute numsubs
numSubs = length(ParticipantIDs);

for i = 1:202
	% load motion para
	mov = dlmread([ParticipantIDs{i},'_ses-1_task-rest_bold_mcf.nii.gz.par']);
	mov = mov(:,[4:6,1:3]);
	% Get FD
	%fdJenk{i} = GetFDJenk(headmotion.mov{i});
    fdJenk{i} = GetFDJenk(mov,80);
	% get threshold of each subject
	percentile{i} = prctile(fdJenk{i},75);
	interquartile{i} = iqr(fdJenk{i});
	thr{i} = percentile{i} + (interquartile{i} * 1.5);
end