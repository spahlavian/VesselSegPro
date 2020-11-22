function [terr_images,terr_masks]=step_wise_segment(seed_voxels,mat_4D,ini_b_masks,thr_array,peak_t,search_radius,seedKernel)
% -------------------------------------------------------------------------
% Stepwise segmentation of vessel structures based on the regions grown
% from the seed_voxels

% written by: Oren Geri 
%   Tel Aviv Sourasky Medical Center
% Soroush H. Pahlavian
%   Mark and Mary Stevens Neuroimaging and Informatics Institute
%   University of Southern California
% © 2018-2020 
% ------------------------------------------------------------------------- 
[s1,s2,s3,s4]=size(mat_4D);

frames=s4;
nTHR = length(thr_array);
n_terr=size(seed_voxels,1);

hh = waitbar(0);
text = sprintf('Threshold: %d / %d',1,length(thr_array));
waitbar(1/length(thr_array),hh,text)
thr=thr_array(1);
pre_blocking_masks=ini_b_masks;
clear terr_masks;
terr_masks=false(s1,s2,s3,n_terr,frames,nTHR);
terr_mask_prev_thr = [];
[terr_masks(:,:,:,:,:,1)]=create_terr_masks(mat_4D,n_terr,frames,thr,pre_blocking_masks,seed_voxels,peak_t,search_radius,seedKernel,terr_mask_prev_thr);

if length(thr_array)>1
    for i=2:nTHR
        text = sprintf('Threshold: %d / %d',i,length(thr_array));
        waitbar(i/length(thr_array),hh,text)
        thr = thr_array(i);
        % Create a blocking mask based on the last thr at the last time-frame:
        n_masks = next_blocking_mask(terr_masks(:,:,:,:,:,i-1),frames);
        % Combine n_mask and pre_blocking_masks (mask cannot grow into where it is currently masked)
        pre_blocking_masks = bitor(pre_blocking_masks,n_masks);
        % -> Mask created based on the segmented regions from OTHER territories
        terr_mask_prev_thr = [];%terr_masks(:,:,:,:,:,i-1);
        [terr_masks(:,:,:,:,:,i)] = create_terr_masks(mat_4D,n_terr,frames,thr,pre_blocking_masks,seed_voxels,peak_t,search_radius,seedKernel,terr_mask_prev_thr);
    end
end
close(hh);



terr_images=zeros(s1,s2,s3,frames,'int16');
for f=1:frames
    for jj=1:n_terr
        terr_images(:,:,:,f)=terr_images(:,:,:,f)+jj*squeeze(int16(terr_masks(:,:,:,jj,f,nTHR)));
    end
end



