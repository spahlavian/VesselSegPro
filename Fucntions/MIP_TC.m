function [mip_TC]=MIP_TC(mat_4D,dim,choose_slices)
% -------------------------------------------------------------------------
% Obtain the Maximum intensity projection of a 3D image
%
% written by: Oren Geri 
%   Tel Aviv Sourasky Medical Center
% © 2018-2020 
% ------------------------------------------------------------------------- 
if nargin<3
    choose_slices=0;
end
[s1,s2,s3,s4]=size(mat_4D);
for i=1:s4
    mip_TC(:,:,i)=MIP(squeeze(mat_4D(:,:,:,i)),dim,choose_slices);
end

    
