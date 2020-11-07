function n_masks=next_blocking_mask(terr_masks,frame)
% -------------------------------------------------------------------------
% Create blocking masks according to the frame number (frame)
%
% written by: Oren Geri 
%   Tel Aviv Sourasky Medical Center
% © 2018-2020 
% ------------------------------------------------------------------------- 
  [s1,s2,s3,n_terr,~]=size(terr_masks);
  n_masks=false(s1,s2,s3,n_terr);
for i=1:n_terr
    b_mask=false(s1,s2,s3);
    for j=1:n_terr   
        if j~=i      
            other_terr=terr_masks(:,:,:,j,frame);   
            b_mask=bitor(b_mask,other_terr);
        end
    end
    n_masks(:,:,:,i)=b_mask;
end
            
