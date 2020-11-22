% -------------------------------------------------------------------------
% Vessel class containing fucntions for evaluating the vessel strcutures
% obtained from the segmentation

% written by: Soroush H. Pahlavian
% Mark and Mary Stevens Neuroimaging and Informatics Institute
% University of Southern California
% © 2018-2020 
% ------------------------------------------------------------------------- 
classdef vessel < handle
    properties(GetAccess = 'public', SetAccess = 'public')
        vol_input;
        skel;
        branch_points;
        end_points;
        seg_mask;
        vol_output;
        seed_points;
    end
    
    methods
        function obj = vessel()
        end
        function obj = plot_mip(obj)
            img_mip = MIP_TC(obj.vol_input,3);
            figure
            imagesc(img_mip)
        end
        function obj = get_skel(obj)
            obj.skel = bwskel(imbinarize(obj.vol_input,0.5));
            obj.branch_points = bwmorph3(obj.skel,'branchpoints');
            obj.end_points = bwmorph3(obj.skel,'endpoints');
        end
        function out_pt = generate_skel_points(~,skel_in)
            inds_skel_in  = find(skel_in==1);
            [xx,yy,zz] = ind2sub(size(skel_in),inds_skel_in);
            out_pt = [xx yy zz];
        end
        
        function plot_surface(obj)
            
            [faces, vertices] = isosurface(obj.vol_input,0.4);
            vert_smooth = SurfaceSmooth(vertices, faces,32);
            figure;
            p = patch('Faces',faces,'Vertices',vert_smooth);
            p.FaceColor = [0.07, 0.62, 1];
            p.EdgeColor = 'none';
            axis equal
            camlight
            lighting gouraud
        end
        
        
        function plot_skel(obj)
            if length(size(obj.skel)) == 3
                skel_pt = obj.generate_skel_points(obj.skel);
                branch_points_pt = obj.generate_skel_points(obj.branch_points);
                end_points_pt = obj.generate_skel_points(obj.end_points);
                figure;
                hold on;
                siz = 40;
                col = [0.3, 0.4, 0.7];
                scatter3(skel_pt(:,1),skel_pt(:,2),skel_pt(:,3),siz,col,'Marker','o',...
                    'MarkerFaceColor','flat','MarkerEdgeColor','k');
                if ~isempty(obj.seed_points)
                    col = [0.9, 0.4, 0.7];
                    scatter3(obj.seed_points(:,1),obj.seed_points(:,2),obj.seed_points(:,3),siz,col,'Marker','o',...
                        'MarkerFaceColor','flat','MarkerEdgeColor','k');
                end
                %                 col = [0.9, 0.4, 0.7];
                %                 scatter3(branch_points_pt(:,1),branch_points_pt(:,2),branch_points_pt(:,3),siz,col,'Marker','o',...
                %                     'MarkerFaceColor','flat','MarkerEdgeColor','k');
                %                 col = [0.3, 0.9, 0.2];
                %                 scatter3(end_points_pt(:,1),end_points_pt(:,2),end_points_pt(:,3),siz,col,'Marker','o',...
                %                     'MarkerFaceColor','flat','MarkerEdgeColor','k');
                axis('equal')
            else
                warndlg('You need to generate the skeleton first!');
            end
        end
        function region_grow(obj, seg_mask, seed_points, first_seed)
            % Add a mask around the first_seed to the block_mask.
            first_seed_mask = false(size(seg_mask));
            rad_first = 3;
            rad_first_l = max(first_seed - rad_first, 1);      % Add points within this radius to the block_mask
            rad_first_u = min(first_seed + rad_first, size(seg_mask));
            first_seed_mask(rad_first_l(1):rad_first_u(1),...
                rad_first_l(2):rad_first_u(2),...
                rad_first_l(3):rad_first_u(3)) = 1;
            %             block_mask = (obj.vol_input & ~seed_points) | first_seed_mask;
            block_mask = obj.vol_input  | first_seed_mask;
            seed_points(first_seed_mask) = 0;
            seg_mask_updated = (seg_mask & ~block_mask);
            search_radius = 2;   % Kernel around seedPoint for the region growing
            
            seed_voxels = obj.generate_skel_points(seed_points);
            
            % ------------- See only the newly segmented reion:
            %             obj.vol_output = false(size(seg_mask_updated));
            % ------------- See the final segmentation:
            obj.vol_output = obj.vol_input;
            
            for i = 1:size(seed_voxels,1)
                [~, vol_temp] = regionGrowing(seg_mask_updated, seed_voxels(i,:), 60, Inf, false, true, false,search_radius);
                if sum(vol_temp(:))>1
                    obj.vol_output = obj.vol_output | vol_temp;
                end
            end
        end
    end
end
