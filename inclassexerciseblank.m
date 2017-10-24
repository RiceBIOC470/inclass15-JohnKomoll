%% step 1: write a few lines of code or use FIJI to separately save the
% nuclear channel of the image Colony1.tif for segmentation in Ilastik

file = '48hColony1_DAPI.tif';
reader = bfGetReader(file);
iplane = reader.getIndex(0,0,0)+1;
img_nuc = bfGetPlane(reader, iplane);
imwrite(imadjust(img_nuc), 'Colony1.tif')

%% step 2: train a classifier on the nuclei
% try to get the nuclei completely but separate them where you can
% save as both simple segmentation and probabilities

% See saved files

%% step 3: use h5read to read your Ilastik simple segmentation
% and display the binary masks produced by Ilastik 

% Read in the data and alter it for MATLAB
segmask = h5read('Colony1_Simple Segmentation.h5', '/exported_data');
segmask = ~logical(reshape(segmask(1,:,:) - 1, [2301 2301]));
imshow(segmask)

% (datasetname = '/exported_data')
% Ilastik has the image transposed relative to matlab
% values are integers corresponding to segmentation classes you defined,
% figure out which value corresponds to nuclei

% The value of 1 corresponds to nuclei and the value of 2 corresponds to
% background. Correcting the mask, I set the value of 1 to correspond to
% nuclei and the value of 0 as background.

%% step 3.1: show segmentation as overlay on raw data

% Reverse the transpose
segmask = transpose(segmask);

% Show the overlay
overlay = cat(3, im2double(imadjust(img_nuc)), segmask, zeros(2301));
figure
imshow(overlay)

%% step 4: visualize the connected components using label2rgb
% probably a lot of nuclei will be connected into large objects

visual = label2rgb(segmask);
figure
imshow(visual)

%% step 5: use h5read to read your Ilastik probabilities and visualize
% it will have a channel for each segmentation class you defined

% Read in the data and alter it for MATLAB
probmasks = h5read('Colony1_Probabilities.h5', '/exported_data');
probmask1 = transpose(reshape(probmasks(1,:,:), [2301 2301]));
probmask2 = transpose(reshape(probmasks(2,:,:), [2301 2301]));
figure
imshow(probmask1)
figure
imshow(probmask2)

% NOTE: I performed the transpose operation to make the image coordinates
% appropriate

%% step 6: threshold probabilities to separate nuclei better
threshold = 0.99;
threshmask = probmask1 > threshold;
figure
imshow(threshmask)

%% step 7: watershed to fill in the original segmentation (~hysteresis threshold)
outside = ~imdilate(segmask, strel('disk', 1));
figure
imshow(outside)

hx = fspecial('sobel');
hy = hx';
Iy = imfilter(double(img_nuc), hy, 'replicate');
Ix = imfilter(double(img_nuc), hx, 'replicate');
gradmag = sqrt(Iy.^2 + Ix.^2);
figure
imshow(gradmag,[])

basin = imimposemin(gradmag, threshmask | outside);
L = watershed(basin);
newmask = L > 1;
figure
imshow(cat(3, newmask, im2double(imadjust(img_nuc)), outside))

%% step 8: perform hysteresis thresholding in Ilastik and compare the results
% explain the differences
Il_objects = h5read('Colony1_Object Predictions.h5', '/exported_data');
Il_objects = logical(reshape(Il_objects(1,1,:,:), [2301 2301]));
imshow(cat(3, Il_objects, im2double(imadjust(img_nuc)), outside))

% The Ilastik hysteresis thresholding did not work very well, because many
% of the cells were connected. After size filtering, only a fraction of the
% cells were left.

% The MATLAB hysteresis thresholding did not work very well either, because
% the cells were very closely clustered and the watershed could not segment
% the cells very well

%% step 9: clean up the results more if you have time 
% using bwmorph, imopen, imclose etc

