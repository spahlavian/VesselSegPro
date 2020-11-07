function [ mat ] = read_nifti(file)
nii=load_untouch_nii(file);
mat=nii.img;
end

