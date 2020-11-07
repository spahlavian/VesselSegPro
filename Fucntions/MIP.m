function [mip_image,mip_ind]=MIP(image,dim,choose_slices)
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
[s1,s2,s3]=size(image);
switch dim
    case 3
        mip_image=zeros(s1,s2);
        mip_ind=zeros(s1,s2);
        for i=1:s1
            for j=1:s2
                if choose_slices
                    [mip_image(i,j),mip_ind(i,j)]=max(image(i,j,choose_slices));                    
                else
                     [mip_image(i,j),mip_ind(i,j)]=max(image(i,j,:));
                end
            end
        end
    case 2
        mip_image=zeros(s3,s1);
        for i=1:s1
            for j=1:s3
                if choose_slices
                     [mip_image(j,i),mip_ind(j,i)]=max(image(i,choose_slices,j));                    
                else
                     [mip_image(j,i),mip_ind(j,i)]=max(image(i,:,j));
                end
            end
        end
        temp=mip_image;
        
    case 1
        mip_image=zeros(s2,s3);
        for i=1:s2
            for j=1:s3
                if choose_slices
                     [mip_image(i,j),mip_ind(i,j)]=max(image(choose_slices,i,j));                    
                else
                     [mip_image(i,j),mip_ind(i,j)]=max(image(:,i,j));
                end
            end
        end
end