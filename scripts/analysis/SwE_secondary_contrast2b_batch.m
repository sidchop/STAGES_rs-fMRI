% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'/projects/kg98/Sid/STAGES/STAGES_fmriprep/analyses/swe_validation/a3/med_v_placebo/batch_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
