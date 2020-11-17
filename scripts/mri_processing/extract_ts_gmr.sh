module load fsl
fmri_dir=/home/scho0011/kg98/Sid/STAGES/STAGES_fmriprep/data/derivatives/fmriprep/fmriprep

for ses in 1 2 4 ; do for s in `cat list_fmri_ses-${ses}.txt` ; do fslmeants -i $fmri_dir/${s}/ses-${ses}/dicer_lite/${s}_ses-${ses}_task-rest_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold+2P_detrended_hpf_GMR.nii.gz -o ts_c_${s}_${ses}.txt --label=restricted_atlases/Schaefer2018_300Parcels_7Networks_order_FSLMNI152_2mm_${s}_ses-${ses}_res.nii.gz --transpose ; done ; done

#The line below is also to extract the choi parc of the stri.
#for ses in 1 2 4 ; do for s in `cat list_fmri_ses-${ses}.txt` ; do fslmeants -i $fmri_dir/${s}/ses-${ses}/dicer_lite/${s}_ses-${ses}_task-rest_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold+2P_detrended_hpf_dbscan.nii.gz -o ts_stri_${s}_${ses}.txt --label=Choi2012_7Networks_MNI152_FreeSurferConformed1mm_TightMask_2mm.nii.gz --transpose ; done ; done




for ses in 1 2 4 ; do for s in `cat list_fmri_ses-${ses}.txt` ; do for k in 10 11 12 13 17 18 26 49 50 51 52 53 54 58 ; do fslmaths $fmri_dir/${s}/ses-${ses}/func/${s}_ses-${ses}_task-rest_space-MNI152Lin_desc-aseg_dseg.nii.gz -thr ${k} -uthr ${k} $fmri_dir/${s}/ses-${ses}/func/aseg_${k} ; done ; cd $fmri_dir/${s}/ses-${ses}/func/ ; fslmaths aseg_10 -add aseg_11 -add aseg_12 -add aseg_13 -add aseg_17 -add aseg_18 -add aseg_26 -add aseg_49 -add aseg_50 -add aseg_51 -add aseg_52 -add aseg_53 -add aseg_54 -add aseg_58 subcort_aseg ; rm aseg*.nii.gz -rf ; fslmeants -i /home/scho0011/kg98/Sid/STAGES/STAGES_fmriprep/data/derivatives/fmriprep/fmriprep/${s}/ses-${ses}/dicer_lite/${s}_ses-${ses}_task-rest_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold+2P_detrended_hpf_GMR.nii.gz -o /home/scho0011/kg98/Sid/STAGES/STAGES_fmriprep/timeseries/schaefer_300parc_7net/ts_sc_${s}_${ses}_temp.txt --label=subcort_aseg.nii.gz --transpose; cd  /home/scho0011/kg98/Sid/STAGES/STAGES_fmriprep/timeseries/schaefer_300parc_7net/ ; awk '!a[$0]++' ts_sc_${s}_${ses}_temp.txt >> ts_sc_${s}_${ses}.txt ; sed -i '1d' ts_sc_${s}_${ses}.txt ; rm ts_sc_${s}_${ses}_temp.txt -rf ; done ; done


for ses in 1 2 4 ; do for s in `cat list_fmri_ses-${ses}.txt` ; do cat ts_sc_${s}_${ses}.txt ts_c_${s}_${ses}.txt > ts_p300_n7/ts_${s}_${ses}.txt ;done ; done
