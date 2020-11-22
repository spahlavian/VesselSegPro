
function b_mask=b_mask_creation(terr_masks,terr_no,frame)
% -------------------------------------------------------------------------
% Function to create blocking masks
%
% written by: Oren Geri 
%   Tel Aviv Sourasky Medical Center
% Soroush H. Pahlavian
%   Mark and Mary Stevens Neuroimaging and Informatics Institute
%   University of Southern California
% © 2018-2020 
% ------------------------------------------------------------------------- 
  [s1,s2,s3,n_terr,~]=size(terr_masks);
  b_mask=false(s1,s2,s3);
for k=1:n_terr
    if k~=terr_no       % add newly created terr_masks ONLY from other terrs to the mask
        
        if k>terr_no % not yet created in this step 
            other_terr=terr_masks(:,:,:,k,frame-1);
        else    
            other_terr=terr_masks(:,:,:,k,frame-1);
        end
        
        b_mask=bitor(b_mask,other_terr);
   
    end
end
            
