function [mask]=create_mask(img,range)
% -------------------------------------------------------------------------
% Create a bit mask according to range limits 
% range(1) - min limit
% optional range(2) - max limit

% written by: Oren Geri 
%   Tel Aviv Sourasky Medical Center
% Soroush H. Pahlavian
%   Mark and Mary Stevens Neuroimaging and Informatics Institute
%   University of Southern California
% © 2018-2020 
% ------------------------------------------------------------------------- 

 [s1,s2,s3]=size(img);
 mask=false(s1,s2,s3);
 
if length(range)==2 && range(1)==range(2) %mask using one value
    voxels = img==range(1);
    mask(voxels) = 1;
else %use lower limit
    voxels = img>=range(1);
    mask(voxels) = 1;
    if length(range)==2 %upper limit
        r_voxels = img>range(2);
        mask(r_voxels) = 0;
    end 
end

    
