
function plotter(savePath,figName,hFig)
% -------------------------------------------------------------------------
% Export figure to a file

% written by: Soroush H. Pahlavian
% Mark and Mary Stevens Neuroimaging and Informatics Institute
% University of Southern California
% © 2018-2020 
% ------------------------------------------------------------------------- 
figure(hFig)
set(gcf, 'Color', 'w');
pathName = savePath;
path = sprintf('%s',pathName);
ExStr = strcat(path,figName);
export_fig(ExStr,'-m2','-opengl')
disp('done!')
%close(hFig)
end
