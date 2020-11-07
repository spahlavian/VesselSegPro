function [dim]=get_resolution(files,PathName)
% -------------------------------------------------------------------------
% Obtain the resolution of the loaded nifti data
%
% written by: Oren Geri 
%   Tel Aviv Sourasky Medical Center
% © 2018-2020 
% ------------------------------------------------------------------------- 
%  info=load_untouch_nii(get_nifti);
 info=load_untouch_nii(strcat(PathName,files{1}));
    %for TOF switch dimentions 2 and 3
    dim(1)=info.hdr.dime.pixdim(1+1);
    dim(2)=info.hdr.dime.pixdim(2+1);
    dim(3)=info.hdr.dime.pixdim(3+1);