function [terr_masks]=create_terr_masks(mat_4D,n_terr,frames,thr,pre_blocking_masks,seed_voxels,peak_t,search_radius,seedKernel,terr_mask_prev_thr)
% -------------------------------------------------------------------------
% Function to create territory masks using progressive thresholding
%
% written by: Oren Geri 
%   Tel Aviv Sourasky Medical Center
% Soroush H. Pahlavian
%   Laboratory of Functional MRI Technology
%   Mark and Mary Stevens Neuroimaging and Informatics Institute
%   University of Southern California
% © 2018-2020 
% ------------------------------------------------------------------------- 
[s1,s2,s3,~] = size(mat_4D);
terr_masks = false(s1,s2,s3,n_terr,frames);


for i = 1:frames
    %         % Only keep voxels which had their peak value between 1 and the current frame:

    for j = 1:n_terr
        % Only keep voxels which had their ATT value between seedATT and the current frame:
        mat = squeeze(mat_4D(:,:,:,i));
        temp = peak_t(seed_voxels(j,1)-seedKernel:seed_voxels(j,1)+seedKernel,...
            seed_voxels(j,2)-seedKernel:seed_voxels(j,2)+seedKernel,...
            seed_voxels(j,3)-seedKernel:seed_voxels(j,3)+seedKernel);
        seedATT = nanmean(temp(:));
        T_mask = create_mask(peak_t,[1 i]);
        seedATT_mask = peak_t>=seedATT;
        % Keep voxels that are greater than the thr
        mask = create_mask(mat,thr);
        mask = bitand(mask,T_mask);
        mask = bitand(mask,seedATT_mask);
        
        % --------  b_mask does not include the current territory. --------
        if i == 1
            b_mask = squeeze(pre_blocking_masks(:,:,:,j));
        else
            % Add territories created at the current time frame to the
            % blocking masks:
            b_mask = b_mask_creation(terr_masks,j,i); 
            b_mask = bitor(b_mask,squeeze(pre_blocking_masks(:,:,:,j)));
        end
        
        % The mask in which we will grow a cluster from the seedPoint:
        mask_terr = bitand(mask,~b_mask);
        mask_terrNew = bwareaopen(mask_terr, 30);
        

        first_seed = seed_voxels(j,:);                          % Ver. 2
        
        if i>1
            mat_3D = bitor(mask_terrNew,squeeze(terr_masks(:,:,:,j,i-1)));
            [~, volIn_2] = regionGrowing(mat_3D, first_seed, 60, Inf, false, true, false,search_radius);
        else
            mat_3D = mask_terrNew;
            [~, volIn_2] = regionGrowing(mat_3D, first_seed, 60, Inf, false, true, false,search_radius);
        end
        if ~isempty(terr_mask_prev_thr)

            cluster = volIn_2;
        else
            cluster = volIn_2;
        end
        
        terr_masks(:,:,:,j,i)=cluster;
    end
    % find overlaped territories at the end of each time-step:
    % Exclude that region from its corresponding territories
    % Add that region to terr_mask(:,:,:,1,i)
    % Create seed_voxels(1,:) as the centroid of the cluster
    if i == 1
        terr_masks(:,:,:,1,i) = sum(squeeze(terr_masks(:,:,:,:,i)),4)>1;
    else
        terr_masks(:,:,:,1,i) = bitor(sum(squeeze(terr_masks(:,:,:,:,i)),4)>1,...
            squeeze(terr_masks(:,:,:,1,i-1)));
    end
    shared_mask = squeeze(terr_masks(:,:,:,1,i));
    if sum(shared_mask(:))>1
        for j = 2:n_terr
            temp = squeeze(terr_masks(:,:,:,j,i));
            temp(shared_mask==1) = 0;
            terr_masks(:,:,:,j,i) = temp;
        end
        
        [xx,yy,zz] = ind2sub(size(shared_mask),find(shared_mask==1));
        
        centroid = zeros(1,3);
        centroid(1,1) = mean(xx);
        centroid(1,2) = mean(yy);
        centroid(1,3) = mean(zz);
        P = [xx yy zz];
        seed_voxels(1,:) = P(dsearchn(P,centroid),:);
    end
end


end



