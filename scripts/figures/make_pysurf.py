# This script uses pySurfer (https://pysurfer.github.io/)
# to map an c to a set of ROIs on brain surfaces. 
# 
# PySurfer is great, but can be a pain in the ass to install, the
# most fool-proof and painless method I've found so far is explained
# here: https://gist.github.com/danjgale/4f64ca81f5e91cc0669d0f744c7a9f82
#
# Once you have PySurfer installed, open up a terminal and enter:
# >$ source activate pysurfer
# >$ python 
#
# The rest of the script operates within python. 
# You will need a fsaverage .annot file. If you are using the Schaefer atlas, these can be found 
# here: https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage
# You will also need a file with metric you would like to visualise, in the same order as the annot file (you will need to add a "0"
# at the start for the medial wall
# 

import os
import numpy as np
import nibabel as nib
from surfer import Brain

print(__doc__)
os.chdir('/dir/where/you/want/output/surfaces') 

subject_id = "fsaverage"
hemi = h # l or h
surf = "inflated" # pick pial, inflated or white
image_name = "analysis_"+h+"_"+str(x)+"_"+surf #for output file name
	
"""
Bring up the visualization.
"""
brain = Brain(subject_id, hemi, surf,background="white", cortex = "grey")

"""
Read in the automatic parcellation of sulci and gyri.
"""

aparc_file = os.path.join('/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/data/fsaverage/label', hemi + '.Schaefer2018_300Parcels_7Networks_order.annot') #Add your own
labels, ctab, names = nib.freesurfer.read_annot(aparc_file)


#rs = np.random.RandomState(4)
#roi_data = rs.uniform(.5, .8, size=len(names))

#attribute:
attribute_file = os.path.join('/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/data/fsaverage/degree', analysis + '.txt')  #Add your own
roi_data = np.loadtxt(attribute)

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
See: https://pysurfer.github.io/generated/surfer.Brain.html#surfer.Brain.add_data
"""
brain.add_data(vtx_data, 0, thresh=0, colormap="YlOrRd", 
	alpha=1, transparent = True, colorbar = True)

#Save  images

brain.save_imageset(analysis_name, ['med', 'lat'], 'jpg')
# other views
#brain.show_view('lateral')
#brain.show_view('m')
#brain.show_view('rostral')
#brain.show_view('caudal')
#brain.show_view('ve')
#brain.show_view('frontal')
#brain.show_view('par')
#brain.show_view('dor')
