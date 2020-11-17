
#bash
source activate pysurfer
python 


#python



import os
import numpy as np
import nibabel as nib
from surfer import Brain

print(__doc__)
os.chdir('/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces') 
#max_degree = np.array([0,66, 115, 80, 58, 41, 10, 6, 5, 21, 15]) #now we are using the same max and min for both pos and neg analysis
max_degree = np.array([0,63, 63, 58, 58, 41, 41, 6, 6, 21, 21, 12, 12, 5, 5])
hemis = "rh lh"

for h in hemis.split():
	for x in range(1,15):
		subject_id = "fsaverage"
		maxdegree = max_degree[x]#for colourscale
		hemi = h
		surf = "inflated" #pick pial, inflated or white
		analysis = "analysis_"+h+"_"+str(x)
		analysis_name = "analysis_"+h+"_"+str(x)+"_"+surf +"_transparent_reds" #for output file name
			
		"""
		Bring up the visualization.
		"""
		brain = Brain(subject_id, hemi, surf,background="white", cortex = "grey")
		
		"""
		Read in the automatic parcellation of sulci and gyri.
		"""
		aparc_file = os.path.join('/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/data/fsaverage/label', hemi + '.Schaefer2018_300Parcels_7Networks_order.annot')
		labels, ctab, names = nib.freesurfer.read_annot(aparc_file)
		
		"""
		Make a random vector of scalar data corresponding to a value for each region in
		the parcellation.
		
		"""
		#rs = np.random.RandomState(4)
		#roi_data = rs.uniform(.5, .8, size=len(names))
		degree_file = os.path.join('/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/data/fsaverage/degree', analysis + '.txt')
		roi_data = np.loadtxt(degree_file)
		
		"""
		Make a vector containing the data point at each vertex.
		"""
		vtx_data = roi_data[labels]
		
		"""
		Handle vertices that are not defined in the annotation.
		"""
		vtx_data[labels == -1] = -1000
		
		"""
		Display these values on the brain. Use a sequential colormap (assuming
		these data move from low to high values), and add an alpha channel so the
		underlying anatomy is visible.
		"""
		brain.add_data(vtx_data, 0, maxdegree, thresh=0, colormap="Reds", 
			transparent = True, colorbar = False)
		
		
		#change view
		
		brain.save_imageset(analysis_name, ['med', 'lat'], 'jpg')





#brain.show_view('lateral')
#brain.show_view('m')
#brain.show_view('rostral')
#brain.show_view('caudal')
#brain.show_view('ve')
#brain.show_view('frontal')
#brain.show_view('par')
#brain.show_view('dor')
