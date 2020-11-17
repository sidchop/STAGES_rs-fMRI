hemisphere = 'rh'

%%
addpath(genpath('/Users/sidchopra/Dropbox/Sid/matlab_functions/'));


%colourmapDir = '/home/asegal/kg98/Ashlea/code/cbrewer/cbrewer/cbrewer';
%addpath(colourmapDir)

% Get colour map
hexMap = {'9C755FFF', 'F28E2BFF', 'E15759FF', '76B7B2FF', '59A14FFF',  'EDC948FF', 'B07AA1FF'}
myColorMap = zeros(length(hexMap), 3); % Preallocate
for k = 1 : length(hexMap)
	thisCell = hexMap{k}
	r = hex2dec(thisCell(1:2))
	g = hex2dec(thisCell(3:4))
	b = hex2dec(thisCell(5:6))
	myColorMap(k, :) = [r, g, b]
end

myColorMap = myColorMap/255




% Main working directory 
wdir = '/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/data/';
cd(wdir)

%% 

% Load Stu's Template
load(['fsaverage_surface_data.mat'])
%,'lh_inflated_verts','lh_verts','lh_faces','lh_HCPMMP1','lh_aparc','lh_rand200','lh_sulc')

% Surfaces 

surface.vertices = rh_inflated_verts;
surface.faces = rh_faces;


% This just plots the ROI ID number for each ROI

Schaefer_filename = ['/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/data/fsaverage/label/rh.Schaefer2018_300Parcels_7Networks_order.txt'];
Schaefer_labels = load(Schaefer_filename);

Schaefer_net = load('net_labs_legend.txt');
Schaefer_net = Schaefer_net(151:300,1) ;




figure

%ax1 = axes('Position',[0.01 0 .3 1]);

plotSurfaceROIBoundary(surface,Schaefer_labels,Schaefer_net,'faces',myColorMap,1,2);

angle = 'lateral';
view([90 0])

angle = 'medial';
view([-90 0])        
 
   
    %view([90 0]) % -90 lateral, 90 medial
    %view([-90 0]) % -90 lateral, 90 medial (inside)

    axis off
    axis tight
    axis equal
    


