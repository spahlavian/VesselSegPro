function varargout = vesselSegPro(varargin)
% -------------------------------------------------------------------------
% Custom software for territorial and sub-territorial segmentation of 
% cerebral vessels based on 4D MR angiography.
%
% Developed by: Soroush H. Pahlavian, PhD
% Laboratory of Functional MRI Technology
% Mark and Mary Stevens Neuroimaging and Informatics Institute
% University of Southern California
% © 2018-2020 
% ------------------------------------------------------------------------- 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vesselSegPro_OpeningFcn, ...
                   'gui_OutputFcn',  @vesselSegPro_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before vesselSegPro is made visible.
function vesselSegPro_OpeningFcn(hObject, eventdata, handles, varargin)
global dim3Num dim4Num
global numDim3 numDim4 maxSCMAG
global numSeeds th_vis_flag
global hMin hMax delH IntpFactor
global ATTFlag search_radius search_radius_seg
global t0 delt numSeeds_seg smoothing_flag smooth_flag

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

code_dir = pwd; 
addpath(genpath(code_dir));
dim3Num = 1;
dim4Num = 1;
numDim3 = 10;
maxSCMAG = 0;
numSeeds = 4;
hMin = 120;
delH = 6;
hMax = 300;
IntpFactor = 1;
ATTFlag = 1;
t0 = 110;
delt = 113;
numSeeds_seg = 3;
search_radius = 1;   % Kernel around seedPoint for the region growing
search_radius_seg = 1;   % Kernel around seedPoint for the region growing
smoothing_flag = 0;
th_vis_flag = 0;
smooth_flag = 0;

set(handles.sl_01, 'Min', 1);
set(handles.sl_01, 'Max', numDim3);
set(handles.sl_01, 'Value', int32(numDim3/2));
set(handles.sl_01, 'SliderStep', [1/numDim3 , 1/numDim3 ]);

set(handles.sl_02, 'Min', 1);
set(handles.sl_02, 'Max', 1000);
set(handles.sl_02, 'Value', int32(1000/2));
set(handles.sl_02, 'SliderStep', [1/1000 , 1/1000 ]);

warning('off','all')



% --- Outputs from this function are returned to the command line.
function varargout = vesselSegPro_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in b_01.
function b_01_Callback(hObject, eventdata, handles)
global mat_4D mat_4DOrig
global s1 s2 s3 s4 numDim3 plotFlag
global forPlotInput maxSCMAG dims caseID

%% ~~~~~~~~~~~~~~ Load subtracted NIFTI files ~~~~~~~~~~~~~~
choice = questdlg('Select the MRI file type', ...
    'MRI file type', ...
    'Subtracted NIFTI','Un-substracted MAT','Subtracted NIFTI');

switch choice
    case 'Subtracted NIFTI'
        readFlag = 1;
    case 'Un-substracted MAT'
        readFlag = 2;
end


%read files
if readFlag ==1
    [files,PathName] = uigetfile(strcat('*.nii'),'Select subtrated NIFTI files','MultiSelect', 'on');
    mat_1=read_nifti(strcat(PathName,files{1}));
    [s1,s2,s3]=size(mat_1);
    [dims]=get_resolution(files,PathName); % get image dimensions
    
    mat_4DOrig = int16(zeros(s1,s2,s3,size(files,2)));
    
    for i=1:size(files,2)
        mat_4DOrig(:,:,:,i) = read_nifti(strcat(PathName,files{i}));
    end
    
    mat_4D = mat_4DOrig;
    [s1,s2,s3,s4]=size(mat_4D);
else
    % Read file (.mat)
    [file,PathName] = uigetfile(strcat('*.mat'),'Select the un-subtrated MAT file','MultiSelect', 'off');
    temp = load(strcat(PathName,file));
    img = temp.image_PICS;
    mat_4DOrig = 1e7*abs((squeeze(abs(img(:,:,:,2,:)))-squeeze(abs(img(:,:,:,1,:)))));
    mat_4D = zeros(size(mat_4DOrig,2),size(mat_4DOrig,1),size(mat_4DOrig,3),size(mat_4DOrig,4));
    dims = [0.9821    0.9821    1.5000];
    for i = 1:size(mat_4DOrig,4)
        mat_4D(:,:,:,i) = imrotate3(mat_4DOrig(:,:,:,i),-90,[0 0 1],'linear','loose','FillValues',0);
    end
    mat_4DOrig = mat_4D;
    [s1,s2,s3,s4]=size(mat_4D);

end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plotFlag = 1;
imgMIP = MIP_TC(mat_4D,3);
forPlotInput = imgMIP;
numDim3 = size(imgMIP,3);
imageDraw(handles,forPlotInput)

%maxSCMAG = max(mat_4D(:));
sliderInit (handles,forPlotInput)

prompt = {'Case ID/Name:',};
dlg_title = '';
defaultans = {''};
ttemp = inputdlg(prompt,dlg_title,1,defaultans);
caseID = ttemp{1,1};    
set(handles.edit8,'String',caseID);

% --- Executes on slider movement.
function sl_01_Callback(hObject, eventdata, handles)
global dim3Num
global numDim3
global forPlotInput

handles = guidata(hObject);
newVal = int32(get(hObject,'Value'));
set(hObject,'Value',newVal);
dim3Num = get(hObject,'Value');
text = strcat(num2str (dim3Num),'/',num2str(numDim3));
set (handles.slSText_01,'String',text);

imageDraw(handles,forPlotInput)
%save('MIP.mat','forPlotInput')


% --- Executes on slider movement.
function sl_02_Callback(hObject, eventdata, handles)
global th_vis plotFlag
global MtAIP_axial  th_vis_flag
global forPlotInput  imgMIP

th_vis_flag = 1;
handles = guidata(hObject);
newVal = int32(get(hObject,'Value'));
set(hObject,'Value',newVal);
th_vis = get(hObject,'Value');
text = strcat(num2str (th_vis));
set (handles.th_txt,'String',text);

if plotFlag == 1
    forPlotInput = double(imgMIP > th_vis);
elseif plotFlag == 2
    forPlotInput = double(MtAIP_axial > th_vis);
end
imageDraw(handles,forPlotInput)







% --- Executes on button press in b_02.
function b_02_Callback(hObject, eventdata, handles)
global mat_4D
global peak_v peak_t
global seed_voxels
global numSeeds

[peak_v,peak_t]=max(mat_4D,[],4); %
MtAIP_axial = max(peak_v,[],3);
forPlotInput = MtAIP_axial;
imageDraw(handles,forPlotInput)

h = [];
for i = 1:numSeeds
    h{i,1} = impoint(handles.axes1);
    color = rand(3,1);
    setColor(h{i,1},color)
end

seed_voxels = zeros(size(h,1),3);
pos = zeros(size(h,1),2);

for i = 1:size(h,1)
    pos(i,:) = getPosition(h{i,1});
    x = round(pos(i,1)); 
    y = round(pos(i,2));
    [~,z] = max(peak_v(y,x,:));
    seed_voxels(i,:) = [y,x,z];
end


% --- Executes on button press in b_03.
function b_03_Callback(hObject, eventdata, handles)
global mat_4D
global peak_v peak_t
global seed_voxels
global s1 s2 s3 s4
global ini_b_masks
global hMin hMax delH
global terr_masks
global noiseSD
global ATTFlag
global terr_images skull_mask visible_flag
global terrVol search_radius seedKernel seedKernMask segNum


set (handles.b_16,'Enable','on');
set (handles.b_17,'Enable','on');

sizePV = size(peak_v);
% skull_mask = false(s1,s2,s3); % for future code skull removal
 % create blocking zones
n_terr = size(seed_voxels,1);
% Blocking mask:
ini_b_masks = false(sizePV(1),sizePV(2),sizePV(3),n_terr);
% Add skull mask to the blocking mask for all terrs:
for i = 1:n_terr %all terr 
    ini_b_masks(:,:,:,i) = bitor(squeeze(ini_b_masks(:,:,:,i)),skull_mask);
end

% Terr # 1 is reserved fro shared regions
% Define the blocking mask for the first two terrs (#2 and 3): main branches:
ter1Area = false(sizePV(1),sizePV(2),sizePV(3));
ter2Area = false(sizePV(1),sizePV(2),sizePV(3));
seed1 = seed_voxels(2,:);
seed2 = seed_voxels(3,:);
ter1Area(1:seed1(1),:,:) = 1;
ter2Area(seed2(1):end,:,:) = 1;
ini_b_masks(:,:,:,1) = bitor(squeeze(ini_b_masks(:,:,:,i)),ter2Area);
ini_b_masks(:,:,:,2) = bitor(squeeze(ini_b_masks(:,:,:,i)),ter1Area);
% Add the region around each seedPoint to the blocking masks for the others
seedKernMask = 6;
for i = 1:n_terr
    for j = 1:n_terr
        if j~=i
            ini_b_masks(...
                max(seed_voxels(i,1)-seedKernMask,1):min(seed_voxels(i,1)+seedKernMask,s1),...
                max(seed_voxels(i,2)-seedKernMask,1):min(seed_voxels(i,2)+seedKernMask,s2),...
                max(seed_voxels(i,3)-seedKernMask,1):min(seed_voxels(i,3)+seedKernMask,s3),j) = 1;
        end
    end
end


thr_array = [hMax:-(hMax-hMin)/delH:hMin];

if ATTFlag==0
    [peak_v,peak_t]=max(mat_4D,[],4); %
else
    peak_t = nan(s1,s2,s3);
    
    for i = 1:s1
        for j = 1:s2
            for k = 1:s3
                %find the first time point at which the value is larger than 3*SD
                timeEvo = mat_4D(i,j,k,:);
                ind = find(timeEvo>3.5*noiseSD);
                if ~isempty(ind)
                    peak_t(i,j,k) = ind(1,1);
                end
            end
        end
    end
end

%segment images based first on thr, then based on frame number 
for i = 1:s4
    temp = mat_4D(:,:,:,i);
    temp (skull_mask==1) = 0;
    mat_4D(:,:,:,i) = temp;
end
peak_t(skull_mask==1) = 1000;
        

seedKernel = 0; % Take the average around the seedPoint with this kernel  

[~,terr_masks] = step_wise_segment(seed_voxels,mat_4D,ini_b_masks,thr_array,peak_t,search_radius,seedKernel);


frames=s4;
n_terr = size(seed_voxels,1);
terr_images = zeros(s1,s2,s3,frames,'int16');
for f=1:frames
    for jj=1:n_terr
        if visible_flag(jj)==1
            terr_images(:,:,:,f) = terr_images(:,:,:,f) + jj*squeeze(int16(terr_masks(:,:,:,jj,f,end)));
        else
            terr_images(:,:,:,f) = terr_images(:,:,:,f) + 0*squeeze(int16(terr_masks(:,:,:,jj,f,end)));
        end
    end
end
segNum = 0;
b_16_Callback(hObject, eventdata, handles)
% imtool3D(int16(MIP_TC(terr_images,3)))

terrVol = zeros(n_terr,frames);
for f=1:frames
    for jj=1:n_terr
        temp = terr_masks(:,:,:,jj,f,end);
        terrVol(jj,f) = sum(temp(:));
    end
end


% --- Executes on button press in b_04.
function b_04_Callback(hObject, eventdata, handles)
global terr_masks dims
global s1 s2 s3 s4 numSeeds segNum terr_masks_seg
global pathNameSave stlPath_surf stlPath_vox caseID visible_flag
global numSeeds_seg
xDim = dims(1)*s1;
yDim = dims(2)*s2;
zDim = dims(3)*s3;

% This can be modified to change the save directory:
saveRes = strcat('C:\Users\SHPro\Desktop\vessle_3D_results');

if ~7==exist(saveRes,'dir')
    mkdir(saveRes);
end
pathName =  saveRes;

rootFolder = pwd;
if ~7==exist(strcat(pathName,filesep,caseID,filesep),'dir')
        cd(pathName);
        mkdir(caseID)
        cd(rootFolder);
end
pathNameSave = strcat(pathName,filesep,caseID,filesep);


pathName =  pathNameSave;
stlPath_surf = strcat(pathName,filesep,'stl_surf',filesep);
stlPath_vox = strcat(pathName,filesep,'stl_vox',filesep);

rootFolder = pwd;
if ~7==exist(stlPath_surf,'dir')
        cd(pathName);
        mkdir stl_surf
        cd(rootFolder);
end
if ~7==exist(stlPath_vox,'dir')
        cd(pathName);
        mkdir stl_vox
        cd(rootFolder);
end

choice = questdlg('Export Main territories?', ...
    'Export Options', ...
    'Yes','No','No');
switch choice
    case 'Yes'
        main_export_flag = 1;
    case 'No'
        main_export_flag = 0;
end

if main_export_flag == 1
    hh = waitbar(0);
    for frameNum = s4 % <- Save only the last time point
        for terNum = 1:numSeeds
            text = sprintf('Frame: %d, Territory: %d / %d',frameNum,terNum,numSeeds);
            waitbar(terNum/numSeeds,hh,text)
            if visible_flag(terNum)==1
                fileName_surf = strcat(stlPath_surf,sprintf('geo3D_ter%d_f%d.stl',terNum,frameNum));
                fileName_vox = strcat(stlPath_vox,sprintf('geo3D_ter%d_f%d.stl',terNum,frameNum));
                temp_binary = logical(terr_masks(:,:,:,terNum,frameNum,end));
                temp_binary_filled = imfill(temp_binary,26, 'holes');
                mat3D = double(temp_binary_filled);

                gridX = dims(1):dims(1):xDim;
                gridY = dims(2):dims(2):yDim;
                gridZ = dims(3):dims(3):zDim;
                
                
                [~,~] = CONVERT_voxels_to_stl(fileName_vox,mat3D,gridX,gridY,gridZ,'binary');
%                 [~,~] = CONVERT_voxels_to_stl(fileName_vox,mat3DFine,gridXF,gridYF,gridZF,'binary');
                [faces, vertices] = isosurface(temp_binary_filled,0.5);
                vert_smooth = SurfaceSmooth(vertices, faces,8);
                myStlwrite(fileName_surf, faces, vert_smooth);
            end
        end
    end
    close(hh)
end

hh = waitbar(0);
if segNum~=0
    for sub_id = 1:numSeeds_seg
        text = sprintf('Sub-territory: %d / %d',sub_id,numSeeds_seg);
        waitbar(sub_id/numSeeds_seg,hh,text)
        fileName_surf = strcat(stlPath_surf,sprintf('geo3D_T%d_sub%d.stl',segNum,sub_id));
        fileName_vox = strcat(stlPath_vox,sprintf('geo3D_T%d_sub%d.stl',segNum,sub_id));
        temp_binary = logical(terr_masks_seg(:,:,:,sub_id,end,end));
        temp_binary_filled = imfill(temp_binary,26, 'holes');
        mat3D = double(temp_binary_filled);
        gridX = dims(1):dims(1):xDim;
        gridY = dims(2):dims(2):yDim;
        gridZ = dims(3):dims(3):zDim;
        [~,~] = CONVERT_voxels_to_stl(fileName_vox,mat3D,gridX,gridY,gridZ,'binary');
        [faces, vertices] = isosurface(temp_binary_filled,0.4);
        vert_smooth = SurfaceSmooth(vertices, faces,8);
        myStlwrite(fileName_surf, faces, vert_smooth);

    end
end
close(hh)



% --- Executes on button press in b_05.
function b_05_Callback(hObject, eventdata, handles)
global s1 s2 s3 s4 numSeeds pathNameSave dims
global savePathImg caseID stlPath visible_flag stlPath_surf
xDim = dims(1)*s1;
yDim = dims(2)*s2;
zDim = dims(3)*s3;

col = zeros(8,3);
col(1,:) = [164/255 199/255 71/255];
col(2,:) = [15/255 172/255 243/255];
col(3,:) = [247/255 23/255 53/255];
col(4,:) = [47/255 57/255 77/255];
col(5,:) = [164/255 199/255 71/255];
%col(5,:) = [255/255 210/255 63/255];
col(6,:) = [255/255 210/255 63/255];
% pathName =  uigetdir(pathNameSave,'select the folder containing STL files');
pathName =  stlPath_surf

rootFolder = pwd;
if ~7==exist(strcat(pathName,filesep,'savedImages',filesep),'dir')
        cd(pathName);
        mkdir savedImages
        cd(rootFolder);
end
savePathImg = strcat(pathName,filesep,'savedImages',filesep);
for frameNum = 1:s4
    try
    hFig = figure('Position', [100 100 1000 800]);
    hold on
    imgName = sprintf('frame_%d',frameNum);
    
    for terNum = 1:numSeeds
        if visible_flag(terNum)==1
            fileName = strcat(stlPath_surf,sprintf('geo3D_ter%d_f%d.stl',terNum,frameNum));
            [faces,verts] = stlread(fileName);
            cindex = ones(size(verts,1),3);
            for j = 1:size(verts,1)
                cindex(j,1:3) = col(terNum,:);
            end
            
            p = patch('Faces',faces,'Vertices',verts,'FaceVertexCData',cindex,'FaceColor','interp');
            axis equal
            p.FaceAlpha = 1;           % remove the transparency
            p.FaceColor = 'interp';    % set the face colors to be interpolated
            p.LineStyle = 'none';      % remove the lines
            colormap(winter)           % change the colormap
            lighting gouraud
            material shiny
            %         view(-28,39)     % change the orientation
            %         view(-34,27)
            view(-144,34)     % change the orientation
            l = light('Position',[-0.4 0.2 0.9],'Style','infinite');
            l.Position = [0.4 0.2 0.9];
            
            l = light('Position',[-0.4 0.2 0.9],'Style','infinite');
            l.Position = [0.3 0.2 0.9];
            
            
            p.FaceLighting = 'gouraud';
            p.AmbientStrength = 0.5;
            p.DiffuseStrength = 0.4;
            p.SpecularStrength = 0.5;
            p.SpecularExponent = 50;
            p.BackFaceLighting = 'unlit';
            material shiny
            axis equal on    % make the axes equal and invisible
            
        end
        xlim([0,xDim])
        ylim([0 yDim])
        zlim([0 zDim])
        plotter(savePathImg,imgName,hFig)
        if frameNum ~=s4
            close(hFig)
        end
    end
    end
end
axis equal on 

imgName = sprintf('frame_%d.png',1);
renderAnimation([],savePathImg,imgName,caseID)

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
global ATTFlag
ATTFlag = get(hObject,'Value');



function imageDraw(handles,forPlotInput)
global maxSCMAG mat_4D visible_flag th_vis_flag
global s1 s2 s3 dim3Num dim4Num dims smooth_flag
global plotFlag numSeeds col3D first3dFlag segNum numSeeds_seg
% plotFlag = 1 : MIP various time frames
% plotFlag = 2: Maximum temporal MIP (for selecting seedPoints)
% plotFlag = 3 : MIP of the segmented regions
% plotFlag = 4: 3D view of the segmented regions

axes(handles.axes1);
cla
if segNum==0
    num_seeds = numSeeds;
else
    num_seeds = numSeeds_seg;
end
if length(size(forPlotInput))==2
    forPlot = forPlotInput;
    maxSCMAG = 0;
elseif length(size(forPlotInput))==3
    forPlot = squeeze(forPlotInput(:,:,dim3Num));
    maxSCMAG = max(mat_4D(:));
elseif length(size(forPlotInput))==4
    forPlot = squeeze(forPlotInput(:,:,dim3Num,dim4Num));
    
end
if plotFlag==3 || plotFlag==4
    minSC = 0;
    maxSC = num_seeds;
    set(handles.sl_02, 'Enable', 'off');
    set(handles.th_txt, 'Enable', 'off');
elseif plotFlag==1 || plotFlag==2
    minSC = min(forPlot(:));
    if maxSCMAG~=0
        maxSC = maxSCMAG;
    else
        maxSC = max(forPlot(:));
    end
    set(handles.sl_02, 'Enable', 'on');
    set(handles.th_txt, 'Enable', 'on');
end
if plotFlag==4
%     col3D = rand(3,sum(visible_flag(:)));
%     col3D(:,1) = [0, 234, 255]/255;
%     col3D(:,2) = [255, 255, 0]/255;
%     col3D(:,3) = [255, 170, 0]/255;
%     col3D(:,4) = [247, 47, 47]/255;
    
    hold on
    [caz,cel] = view(handles.axes1);
    for i = 1:sum(visible_flag(:))
        try
            forPlot = squeeze(forPlotInput(:,:,:,dim3Num,i));
%             
            if smooth_flag==0
                blockPlot(forPlot,[0 0 0],'FaceColor',col3D(:,i),'EdgeColor',[0.15 0.15 0.15]);
            else
                [faces, vertices] = isosurface(forPlot,0.5);
                vert_smooth = SurfaceSmooth(vertices, faces,16);
                p = patch('Faces',faces,'Vertices',vert_smooth);
                p.FaceColor = col3D(:,i);
                p.EdgeColor = 'none';
                camlight
                lighting gouraud
            end
        end
    end
    if first3dFlag==1
        xDim = [0, dims(1)*s1];
        yDim = [0, dims(2)*s2];
        zDim = [0, dims(3)*s3];
    else
        xDim = xlim(handles.axes1);
        yDim = ylim(handles.axes1);
        zDim = zlim(handles.axes1);
    end
    first3dFlag = 0;
    axis equal
    axis tight
    xlim(xDim)
	ylim(yDim)
	zlim(zDim)
    view(handles.axes1,caz,cel)
    colorbar('off')
    legend
else
    imageX = imagesc(forPlot);
    % set(imageX,'Hittest','off');
    minSC  = min(forPlot(:));
    axis equal
    axis tight; 
    colorbar
    legend('off')
end

if th_vis_flag == 1
    minSC = 0;
    maxSC = 1;
end

caxis([minSC maxSC])
% set(gca,'DataAspectRatio',[s1/s2 1 1])
set(gca,'Color',[80/255 80/255 80/255]);
set(gca,'XTickLabel','');
set(gca,'YTickLabel','');
colormap('jet')



function sliderInit (handles,forPlotInput)
global numDim3 
set(handles.sl_01, 'Enable', 'on');
set(handles.sl_01, 'Min', 1);
set(handles.sl_01, 'Max', numDim3);
set(handles.sl_01, 'Value', int32(numDim3/2));
set(handles.sl_01, 'SliderStep', [1/numDim3 , 1/numDim3 ]);


set(handles.sl_02, 'Min', 1);
set(handles.sl_02, 'Max', 1000);
set(handles.sl_02, 'Value', int32(1000/2));
set(handles.sl_02, 'SliderStep', [1/1000 , 1/1000 ]);





function edit1_Callback(hObject, eventdata, handles)
global numSeeds visible_flag
numSeeds = str2num(get(handles.edit1,'String')) +1 ;

stText = cell(numSeeds+1,1);

for i =1:numSeeds+1
     stText(i,1) = {sprintf('%d',i)};
end
stText(numSeeds+1,1) = {'All'};

set(handles.menu_01, 'String', (stText));
visible_flag = ones(numSeeds,1);




function edit2_Callback(hObject, eventdata, handles)
global hMin
hMin = str2num(get(handles.edit2,'String'));


function edit3_Callback(hObject, eventdata, handles)
global delH
delH = str2num(get(handles.edit3,'String'));


function edit4_Callback(hObject, eventdata, handles)
global hMax
hMax = str2num(get(handles.edit4,'String'));




% --- Executes on button press in b_06.
function b_06_Callback(hObject, eventdata, handles)
global mat_4D noiseSD
global peak_v peak_t s1 s2 s3 s4

[peak_v,peak_t]=max(mat_4D,[],4); %
MtAIP_axial = max(peak_v,[],3);
forPlotInput = MtAIP_axial;
% imgMIP = MIP_TC(mat_4D,3);
% forPlotInput = imgMIP;
% forPlotInput = squeeze(mat_4D(:,:,15,:));
imageDraw(handles,forPlotInput)

sizePV = size(peak_v);
noiseMask = zeros(sizePV(1),sizePV(2),sizePV(3));
h = impoly(handles.axes1);
pos = getPosition(h);
temp = poly2mask(pos(:,1),pos(:,2),s1,s2);

for i = 1:s3%round(s4/4):round(3*s4/4)
    noiseMask(:,:,i) = temp;
end


masked4D = nan(s1,s2,s3,s4);
for i = 1:s4
    temp = double(squeeze(mat_4D(:,:,:,i)));
    temp(noiseMask==0) = NaN;    
    masked4D(:,:,:,i) = temp;
end

% figure
% imagesc(squeeze(masked4D(:,:,round(s4/2),8)))
tt = squeeze(masked4D(:,:,:,end));
noiseSD = nanstd(tt(:))


% --- Executes on button press in b_07.
function b_07_Callback(hObject, eventdata, handles)
global terrVol terr_masks terr_images
global CBF_3D peak_pos_R_3D caseID  visible_flag
global hMin delH hMax

n_terr = size(terrVol,1);
frames = size(terrVol,2);

% terrVolNew = zeros(5,size(terrVol,2));
% terrVolNew(1:4,:) = terrVol(1:4,:);
% terrVolNew(1,:) = terrVol(1,:)+terrVol(5,:);
% terrVolNew(5,:) = terrVol(6,:);

hFig = figure;
hold on
for i = 1:n_terr
    if visible_flag(i)==1
        plot(terrVol(i,:))
    end
end

% This can be modified to change the save directory:
saveRes = strcat('C:\Users\SHPro\Desktop\vessle_3D_results');

if ~7==exist(saveRes,'dir')
    mkdir(saveRes);
end
pathName =  saveRes;

rootFolder = pwd;
if ~7==exist(strcat(pathName,filesep,caseID,filesep),'dir')
        cd(pathName);
        mkdir(caseID)
        cd(rootFolder);
end
pathNameSave = strcat(pathName,filesep,caseID,filesep);
% legend(gca,'show');
% plotter(strcat(pathNameSave,filesep),'ters_Volum',hFig)
% pathNameSave = 'C:\Users\SHPro\Desktop';
savePath = strcat(pathNameSave,filesep,caseID,'.mat');
save(savePath,'terr_masks','terr_images','terrVol',...
    'hMin', 'delH', 'hMax')





% --- Executes on button press in b_08.
function b_08_Callback(hObject, eventdata, handles)
global mat_4D skull_mask
global peak_v peak_t s1 s2 s3 pos_mask

[peak_v,peak_t]=max(mat_4D,[],4); %
MtAIP_axial = max(peak_v,[],3);
forPlotInput = MtAIP_axial;
% imgMIP = MIP_TC(mat_4D,3);
% forPlotInput = imgMIP;
% forPlotInput = squeeze(mat_4D(:,:,15,:));
imageDraw(handles,forPlotInput)

sizePV = size(peak_v);
skull_mask = zeros(sizePV(1),sizePV(2),sizePV(3));
h = impoly(handles.axes1);
pos_mask = getPosition(h);
temp = poly2mask(pos_mask(:,1),pos_mask(:,2),s1,s2);

for i = 1:s3%round(s4/4):round(3*s4/4)
    skull_mask(:,:,i) = 1-temp;
end

% figure
% imagesc(skull_mask(:,:,5))



function edit5_Callback(hObject, eventdata, handles)
global IntpFactor
IntpFactor = str2num(get(handles.edit5,'String'));


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in b_11.
function b_11_Callback(hObject, eventdata, handles)
global mat_4D mat_4DOrig
global s1 s2 s3 s4 numDim3
global forPlotInput IntpFactor dims seed_voxels seed_voxels_seg
global pos_mask skull_mask

[s1OG,s2OG,s3OG,s4OG] = size(mat_4DOrig);

xDim = dims(1)*s1OG;
yDim = dims(2)*s2OG;
zDim = dims(3)*s3OG;
reverse_fact = 1;

if IntpFactor~=1
    [XCoarse,YCoarse,ZCoarse] = meshgrid(dims(2):dims(2):yDim,dims(1):dims(1):xDim,dims(3):dims(3):zDim);
    [XFine,YFine,ZFine] = meshgrid(dims(2)/IntpFactor:dims(2)/IntpFactor:yDim,...
        dims(1)/IntpFactor:dims(1)/IntpFactor:xDim,dims(3)/IntpFactor:dims(3)/IntpFactor:zDim);
    mat_4DFine=int16(zeros(IntpFactor*s1OG,IntpFactor*s2OG,IntpFactor*s3OG,s4OG));
    for i=1:s4OG
        mat_4DFine(:,:,:,i) = double(interp3(XCoarse,YCoarse,ZCoarse,double(squeeze(mat_4DOrig(:,:,:,i))),XFine,YFine,ZFine));
    end
    mat_4DFine(isnan(mat_4DFine)) = 0;
    mat_4D = mat_4DFine;
    [s1,s2,s3,s4]=size(mat_4D);
    seed_voxels = IntpFactor * seed_voxels;
    seed_voxels_seg = IntpFactor * seed_voxels_seg;
    reverse_fact = 1 / IntpFactor;
else
    mat_4D = mat_4DOrig;
    [s1,s2,s3,s4]=size(mat_4D);
    seed_voxels = reverse_fact * seed_voxels;
    seed_voxels_seg = reverse_fact * seed_voxels_seg;
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
imgMIP = MIP_TC(mat_4D,3);
forPlotInput = imgMIP;
numDim3 = size(imgMIP,3);
imageDraw(handles,forPlotInput)


function renderAnimation(mainFolder,workingDir,frameName,caseID)

% prompt = {'Enter FPS:',};
% dlg_title = '';
% defaultans = {'4'};
% ttemp = inputdlg(prompt,dlg_title,1,defaultans);
% fps = str2num(ttemp{1,1});       

fps = 4;

if nargin<4
    if ~isempty(mainFolder)
        [frameName,workingDir] = uigetfile({'*.tif';'*.jpg';'*.*'},'Select the first frame to make a rendered animation',mainFolder);
    else
        [frameName,workingDir] = uigetfile({'*.tif';'*.jpg';'*.*'},'Select the first frame to make a rendered animation');
        ind0 = strfind(workingDir,filesep);
        mainFolder = workingDir(1:ind0(end-1)-1);
    end

else
    ind0 = strfind(workingDir,filesep);
    mainFolder = workingDir(1:ind0(end-1)-1);
end


[~,~,ext] = fileparts(frameName);
ind1 = strfind(frameName,'_');
fName = frameName(1:ind1(end));



filesAndFolders = dir(workingDir);     % Returns all the files and folders in the directory
filesInDir = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory                    
stringToBeFound = fName;
numOfFiles = length(filesInDir);
numImages = 0;
i=1;
while(i<=numOfFiles)
      filename = filesInDir(i).name;                     % Store the name of the file
      found = strfind(filename,stringToBeFound);         % Search for the stringToBeFound in contentOfFile
      if ~isempty(found)
          numImages = numImages+1;
      end                                            
      i = i+1;
end

images    = cell(numImages,1);
NewImage    = cell(numImages,1);

maxWidth = 0;
maxHight = 0;
for i = 1:numImages
    imName = strcat(fName,sprintf('%d',i),ext);
    images{i} = imread(strcat(workingDir,imName));
    if size(images{i},2)>maxWidth
        maxWidth = size(images{i},2);
    end
    if size(images{i},1)>maxHight
        maxHight = size(images{i},1);
    end
end

for i = 1:numImages
    NewImage{i} = images{i};
    if maxWidth>size(NewImage{i},2)
        temp = NewImage{i};
        temp (:,size(NewImage{i},2)+1:maxWidth,:) = 255;
        NewImage{i} = temp;
    end
    if maxHight>size(NewImage{i},1)
        temp = NewImage{i};
        temp (size(NewImage{i},1)+1:maxHight,:,:) = 255;
        NewImage{i} = temp;
    end
end



saveName = strcat(caseID,'.avi');
% [file,path] = uiputfile(strcat(mainFolder,filesep,saveName),'Save rendered animation');

% create the video writer with 1 fps
%  writerObj = VideoWriter(strcat(path,file));
writerObj = VideoWriter(strcat(mainFolder,filesep,saveName));
 writerObj.FrameRate = fps;
%  % set the seconds per image
%  secsPerImage = [5 10 15];
 % open the video writer
 open(writerObj);
 % write the frames to the video
 for u=1:length(NewImage)
     % convert the image to a frame
     frame = im2frame(NewImage{u});
%      for v=1:secsPerImage(u) 
         writeVideo(writerObj, frame);
%      end
 end
 % close the writer object
 close(writerObj);
 warndlg('Video rendered successfully!')



function edit8_Callback(hObject, eventdata, handles)
global caseID
caseID = get(handles.edit8,'String');



% --- Executes on selection change in menu_01.
function menu_01_Callback(hObject, eventdata, handles)
global segFlag visible_flag isolated_segment
global numSeeds segNum terr_images

segFlag = get(hObject,'Value');
if segFlag == numSeeds+1
	segNum = 0; 
else
    segNum = segFlag;
    isolated_segment = terr_images==segFlag;
    axes(handles.axes1);
    imagesc( max(isolated_segment(:,:,:,end),[],3) )
    %     figure
    %     imagesc(int16(MIP_TC(isolated_segment,3)))
    if visible_flag(segFlag) == 1
        set(handles.checkbox2, 'Value', 1);
    else
        set(handles.checkbox2, 'Value', 0);
    end
end





% --- Executes on button press in b_13.
function b_13_Callback(hObject, eventdata, handles)
global peak_v isolated_segment skel_points
global seed_voxels_seg numSeeds_seg

choice = questdlg('Seed point selection mode', ...
    'Seed point selection mode', ...
    'MIP','3D Skeleton','3D Surface','MIP');
switch choice
    case 'MIP'
        b_18_Callback(hObject, eventdata, handles)
        h = [];
        for i = 1:numSeeds_seg
            h{i,1} = impoint(handles.axes1);
            color = rand(3,1);
            setColor(h{i,1},color)
        end
        seed_voxels_seg = zeros(size(h,1),3);
        pos = zeros(size(h,1),2);
        for i = 1:size(h,1)
            pos(i,:) = getPosition(h{i,1});
            x = round(pos(i,1));
            y = round(pos(i,2));
            [~,z] = max(peak_v(y,x,:));
            seed_voxels_seg(i,:) = [y,x,z];
        end
        seed_voxels_seg
    case '3D Skeleton'
        
        seed_voxels_seg = zeros(numSeeds_seg,3);
        v_segment = vessel();
        v_segment.vol_input = double(squeeze(isolated_segment(:,:,:,end)));
        v_segment.get_skel();
        v_segment.plot_skel();
        
        
        Prompt = {'Field'
            'Values'};
        Formats = struct('type',{'edit','edit'});
        Formats(1).enable = 'inactive'; % acts as a label
        Formats(2).format = 'vector';
        Formats = Formats.';
        % create initial answer
        DefAns = cell(size(Prompt,1),3);
        DefAns{1,1} = 'X-Coordinate';
        DefAns{1,2} = 'Y-Coordinate';
        DefAns{1,3} = 'Z-Coordinate';
        for datanr = 1:3
            DefAns{2,datanr} = randi(20,1,3);
        end
        
        DefAns = cell2struct(DefAns,Prompt,1);
        Prompt = repmat(Prompt,1,2);
        
        Title = 'Enter seed points coordinates [2:end]';
        Options.AlignControls = 'on';
        Options.CreateFcn = @(~,~,handles)celldisp(get(handles,'type'));
        Options.DeleteFcn = @(~,~,handles)celldisp(get(handles,'type'));
        
        prompt_ans = inputsdlg(Prompt,Title,Formats,DefAns,Options);
        
        seed_voxels_seg(1,:) = [2, 2, 2];
        for i = 2:numSeeds_seg
            seed_voxels_seg(i,1) = prompt_ans(1).Values(i-1);
            seed_voxels_seg(i,2) = prompt_ans(2).Values(i-1);
            seed_voxels_seg(i,3) = prompt_ans(3).Values(i-1);
        end
        
        case '3D Surface'
        numSeeds_seg = 1;
        seed_voxels_seg = 2*ones(numSeeds_seg,3);
        v_segment = vessel();
        v_segment.vol_input = double(squeeze(isolated_segment(:,:,:,end)));
        v_segment.get_skel();    
%         v_segment.plot_surface();

        skel_points = v_segment.generate_skel_points(v_segment.skel);
        vol_for_surf= double(squeeze(isolated_segment(:,:,:,end)));
       
        surface_data (vol_for_surf)
%         global AAT_segment 
%         save('vol_test.mat', 'vol_for_surf','AAT_segment')
end

function  surface_data (vol_for_surf)
global fig_class numSeeds_seg AAT_segment terr_images segFlag
global s1 s2 s3 s4
numSeeds_seg = 1;
[faces, vertices] = isosurface(vol_for_surf,0.4);
vert_smooth = SurfaceSmooth(vertices, faces,32);

fig_class.fh = figure;
fig_class.ax1 = subplot(1,2,1,...
            'buttondownfcn',{@ax_bdfcn,fig_class},...
            'nextplot','replacechildren');
% p = patch('Faces',faces,'Vertices',vert_smooth);
% p.FaceColor = [0.07, 0.62, 1];
% p.EdgeColor = 'none';
% axis equal
% axis tight
% camlight
% lighting gouraud

AAT_segment = nan(s1, s2, s3);
isolated_segment = terr_images==segFlag;
for i = 1:s1
    for j = 1:s2
        for k = 1:s3
            %find the first time point at which the value is larger than 3*SD
            timeEvo = isolated_segment(i,j,k,:);
            ind = find(timeEvo>0);
            if ~isempty(ind)
                AAT_segment(i,j,k) = ind(1,1);
            end
        end
    end
end

%
vol_ind = find(vol_for_surf == 1);
vol_points = zeros(size(vol_ind,1),3);
for ind = 1:size(vol_ind,1)
    [i,j,k] = ind2sub (size(vol_for_surf),vol_ind(ind));
    vol_points(ind,:) = [i,j,k]';
end
vert_smooth_reordered = zeros(size(vert_smooth));
vert_smooth_reordered(:,1) = vert_smooth(:,2);
vert_smooth_reordered(:,2) = vert_smooth(:,1);
vert_smooth_reordered(:,3) = vert_smooth(:,3);
dist_th = 10;
[near_ind,near_dist] = knnsearch(vol_points, vert_smooth_reordered, 'K', 50);

vert_AAT = zeros(size(vert_smooth,1),1);

for i = 1:size(vert_smooth,1)
    near_coords = sub2ind(size(AAT_segment),vol_points(near_ind(i,:),1),...
        vol_points(near_ind(i,:),2),...
        vol_points(near_ind(i,:),3));
    dist_mask = near_dist(i,:)'<dist_th;
    near_AAT = AAT_segment(near_coords);
    near_AAT(~dist_mask) = nan;
    vert_AAT(i) = nanmean(near_AAT);
end
p = patch('Faces',faces,'Vertices',vert_smooth,'FaceVertexCData',vert_AAT);
p.FaceColor = 'interp';
p.EdgeColor = 'none';
% colorbar
axis equal
axis tight
camlight
lighting gouraud
AdvancedColormap('temp')

% Enable data cursor mode
datacursormode on
dcm_obj = datacursormode(fig_class.fh);
% Set update function
set(dcm_obj,'UpdateFcn',@myupdatefcn)
disp('Clicked positioin is')
% Export cursor to workspace
info_struct = getCursorInfo(dcm_obj);
if isfield(info_struct, 'Position')
    disp('Clicked positioin is')
    disp(info_struct.Position)
end

function output_txt = myupdatefcn(~,event_obj)
global skel_points 
global seed_candidate
pos = get(event_obj, 'Position');
selected_point = [pos(2); pos(1); pos(3)]';
nearest_ind = dsearchn(skel_points,selected_point);
seed_candidate = skel_points(nearest_ind,:);
output_txt = {};% {[sprintf('%1.1f',pos(1))], [sprintf('%1.1f',pos(2))], [sprintf('%1.1f',pos(3))]};
skel_plot()

function [] = ax_bdfcn(varargin)
global seed_candidate seed_voxels_seg
global numSeeds_seg
amp = 0.02;
fs = 20500;  
duration=0.05;
values = 0:1/fs:duration;

% buttondownfcn for axes.
[h,S] = varargin{[1,3]};  % Extract the calling handle and structure.

seltype = get(S.fh,'selectiontype'); % Right-or-left click?

switch seltype
    case 'alt'
        if numSeeds_seg > 1
            seed_voxels_seg(numSeeds_seg,:) = [];
        end
        numSeeds_seg = numSeeds_seg - 1;
        if numSeeds_seg == 0
            numSeeds_seg = 1;
        end
        skel_plot()
        freq = 700;
        a = amp*sin(2*pi* freq*values);
        sound(a)
        seed_voxels_seg
        
    case 'normal'
        numSeeds_seg = numSeeds_seg + 1;
        seed_voxels_seg(numSeeds_seg,:) = seed_candidate;
        skel_plot()
        freq = 1600;
        a = amp*sin(2*pi* freq*values);
        sound(a)
        seed_voxels_seg
    otherwise
        % Do something else for double-clicks, etc.
end
function skel_plot()
global fig_class skel_points seed_candidate seed_voxels_seg numSeeds_seg
fig_class.ax2 = subplot(1,2,2,...
            'buttondownfcn',{@ax_bdfcn,fig_class},...
            'nextplot','replacechildren');
cla(fig_class.ax2)
hold on;
siz = 40;
col = [0.3, 0.4, 0.7];
scatter3(skel_points(:,1),skel_points(:,2),skel_points(:,3),siz,col,'Marker','o',...
    'MarkerFaceColor','flat','MarkerEdgeColor','k');
col = [0.9, 0.4, 0.7];
scatter3(seed_candidate(:,1),seed_candidate(:,2),seed_candidate(:,3),siz,col,'Marker','o',...
    'MarkerFaceColor','flat','MarkerEdgeColor','k');
if numSeeds_seg > 1
    col = [0.3, 0.9, 0.2];
    scatter3(seed_voxels_seg(2:end,1),seed_voxels_seg(2:end,2),seed_voxels_seg(2:end,3),siz,col,'Marker','o',...
        'MarkerFaceColor','flat','MarkerEdgeColor','k');
end
set(fig_class.ax2, 'Ydir', 'reverse')
axis('equal')
% hlink = linkprop([fig_class.ax1,fig_class.ax2],{'CameraPosition','CameraUpVector'}); 

% --- Executes on button press in b_14.
function b_14_Callback(hObject, eventdata, handles)
global peak_v peak_t
global seed_voxels_seg
global s1 s2 s3 s4
global terr_masks_seg terr_images_seg numSeeds_seg segFlag 
global terr_images 
global seedKernel_seg search_radius_seg seedKernMask AAT_segment
set(handles.edit12, 'String', num2str(numSeeds_seg))
sizePV = size(peak_v);
% Blocking mask:
ini_b_masks = false(sizePV(1),sizePV(2),sizePV(3),numSeeds_seg);
isolated_segment = terr_images==segFlag;
% Sub-territories can not get out of the main territory:
for i = 1:numSeeds_seg %all terr
    ini_b_masks(:,:,:,i) =~isolated_segment(:,:,:,end);
end
AAT_segment = nan(s1, s2, s3);
AAT_segment2 = nan(s1, s2, s3);
for i = 1:s1
    for j = 1:s2
        for k = 1:s3
            %find the first time point at which the value is larger than 3*SD
            timeEvo = isolated_segment(i,j,k,:);
            ind = find(timeEvo>0);
            if ~isempty(ind)
                AAT_segment(i,j,k) = ind(1,1);
                AAT_segment2(i,j,k) = peak_t(i,j,k);
            end
        end
    end
end

img_mip = min(AAT_segment,[],3);
figure
imagesc(img_mip)
axis equal
axis tight

img_mip = min(AAT_segment2,[],3);
figure
imagesc(img_mip)
axis equal
axis tight




AdvancedColormap('temp')

file_3d = squeeze(isolated_segment(:,:,:,end));
% save('test_3D','AAT_segment', 'file_3d');

s1 = sizePV(1);
s2 = sizePV(2);
s3 = sizePV(3);
s4 = size(isolated_segment,4);

% Add the region around each seedPoint to the blocking masks for the others
seedKernMask = 3;
for i = 1:numSeeds_seg
    for j = 1:numSeeds_seg
        if j~=i
            ini_b_masks(...
                max(seed_voxels_seg(i,1)-seedKernMask,1):min(seed_voxels_seg(i,1)+seedKernMask,s1),...
                max(seed_voxels_seg(i,2)-seedKernMask,1):min(seed_voxels_seg(i,2)+seedKernMask,s2),...
                max(seed_voxels_seg(i,3)-seedKernMask,1):min(seed_voxels_seg(i,3)+seedKernMask,s3),j) = 1;
        end
    end
end
thr_array = 1;
seedKernel_seg = 0; % Take the average around the seedPoint with this kernel  
[~,terr_masks_seg] = step_wise_segment(seed_voxels_seg,isolated_segment,ini_b_masks,thr_array,AAT_segment,search_radius_seg,seedKernel_seg);

frames=s4;
n_terr = size(seed_voxels_seg,1);
terr_images_seg = zeros(s1,s2,s3,frames,'int16');
for f=1:frames
    for jj=1:n_terr
        terr_images_seg(:,:,:,f) = terr_images_seg(:,:,:,f) + jj*squeeze(int16(terr_masks_seg(:,:,:,jj,f,end)));
    end
end

b_16_Callback(hObject, eventdata, handles)

%imtool3D(int16(MIP_TC(seg_mask,3)))

% --- Executes on button press in b_17.
function b_17_Callback(hObject, eventdata, handles)
global forPlotInput terr_images numDim3
global plotFlag maxSCMAG numSeeds first3dFlag
global col3D segNum visible_flag th_vis_flag
global terr_images_seg numSeeds_seg
h = helpdlg('Generating the graphics, please wait...');
first3dFlag = 1;
plotFlag = 4;

col3D = rand(3,10);
col3D(:,1) = [64, 64, 64]/255;
col3D(:,2) = [64, 64, 64]/255;
col3D(:,3) = [64, 64, 64]/255;
col3D(:,4) = [13, 209, 52]/255; 
col3D(:,5) = [240, 240, 81]/255;
col3D(:,6) = [64, 64, 64]/255;
col3D(:,7) = [255, 146, 74]/255;
col3D(:,8) = [112, 0, 0]/255;
col3D(:,9) = [64, 64, 64]/255;
col3D(:,10) = [64, 64, 64]/255;
if segNum==0
%     col3D = rand(3,numSeeds);
    forPlotInput = zeros(size(terr_images,1),size(terr_images,2),...
        size(terr_images,3),size(terr_images,4),sum(visible_flag(:)));
    i_counter = 1;
    for i = 1:numSeeds
        if visible_flag(i)==1
            forPlotInput(:,:,:,:,i_counter) = terr_images==i;
            i_counter= i_counter + 1;
        end
    end
    numDim3 = size(terr_images,4);
    maxSCMAG = 0;
    imageDraw(handles,forPlotInput)
    sliderInit (handles,forPlotInput)
else
%     col3D = rand(3,numSeeds_seg);
    forPlotInput = zeros(size(terr_images_seg,1),size(terr_images_seg,2),...
        size(terr_images_seg,3),size(terr_images_seg,4),numSeeds_seg);
    for i = 1:numSeeds_seg
        forPlotInput(:,:,:,:,i) = terr_images_seg==i;
    end
    numDim3 = size(terr_images_seg,4);
    maxSCMAG = 0;
    imageDraw(handles,forPlotInput)
    sliderInit (handles,forPlotInput)
end
close(h)
th_vis_flag = 0;



% --- Executes on button press in b_16.
function b_16_Callback(hObject, eventdata, handles)
global forPlotInput terr_images numDim3 
global terr_images_seg th_vis_flag
global plotFlag maxSCMAG segNum
h = helpdlg('Generating the graphics, please wait...');
plotFlag = 3;

if segNum==0
    imgMIP = MIP_TC(terr_images,3);
else
    imgMIP = MIP_TC(terr_images_seg,3);
end

forPlotInput = imgMIP;
numDim3 = size(imgMIP,3);
maxSCMAG = 0;
imageDraw(handles,forPlotInput)
sliderInit (handles,forPlotInput)
close(h)
th_vis_flag = 0;

% --- Executes on button press in b_15.
function b_15_Callback(hObject, eventdata, handles)
global forPlotInput mat_4D numDim3 imgMIP
global plotFlag th_vis_flag
h = helpdlg('Generating the graphics, please wait...');
plotFlag = 1;
imgMIP = MIP_TC(mat_4D,3);
forPlotInput = imgMIP;
numDim3 = size(imgMIP,3);
imageDraw(handles,forPlotInput)
%maxSCMAG = max(mat_4D(:));
sliderInit (handles,forPlotInput)
close(h)
th_vis_flag = 0;

function edit12_Callback(hObject, eventdata, handles)
global numSeeds_seg
numSeeds_seg = str2num(get(handles.edit12,'String'));


% --- Executes on button press in b_18.
function b_18_Callback(hObject, eventdata, handles)
global forPlotInput mat_4D numDim3 MtAIP_axial
global plotFlag th_vis_flag

h = helpdlg('Generating the graphics, please wait...');
plotFlag = 2;
[peak_v,~]=max(mat_4D,[],4); %
MtAIP_axial = max(peak_v,[],3);
forPlotInput = MtAIP_axial;
imageDraw(handles,forPlotInput)
set(handles.sl_01, 'Enable', 'off');
close(h)
th_vis_flag=0;



function edit15_Callback(hObject, eventdata, handles)
global search_radius
search_radius = str2num(get(handles.edit15,'String'));





function edit13_Callback(hObject, eventdata, handles)
global search_radius_seg
search_radius_seg = str2num(get(handles.edit15,'String'));


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
global visible_flag numSeeds segFlag segNum s4
global terr_images s1 s2 s3 terr_masks
vis_flag = get(hObject,'Value');

if segFlag ~= numSeeds+1
    visible_flag(segFlag) = vis_flag;
end
frames=s4;
terr_images = zeros(s1,s2,s3,frames,'int16');
for f=1:frames
    for jj=1:numSeeds
        if visible_flag(jj)==1
            terr_images(:,:,:,f) = terr_images(:,:,:,f) + jj*squeeze(int16(terr_masks(:,:,:,jj,f,end)));
        else
            terr_images(:,:,:,f) = terr_images(:,:,:,f) + 0*squeeze(int16(terr_masks(:,:,:,jj,f,end)));
        end
    end
end


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
global seed_voxels numSeeds seed_voxels_seg numSeeds_seg
global pathNameSave caseID hMin hMax delH

saveRes = strcat('C:\Users\SHPro\Desktop\vessle_3D_results');
if ~7==exist(saveRes,'dir')
    mkdir(saveRes);
end
pathName =  saveRes;
rootFolder = pwd;
if ~7==exist(strcat(pathName,filesep,caseID,filesep),'dir')
        cd(pathName);
        mkdir(caseID)
        cd(rootFolder);
end
pathNameSave = strcat(pathName,filesep,caseID,filesep);

fileD = strcat(caseID,'_seeds_t00','.mat');
pathD = pathNameSave;
[file,path] = uiputfile(strcat(pathD,fileD),'Save seed voxles (UNSCALED)');
    
save(strcat(path,file),'seed_voxels','numSeeds',...
    'seed_voxels_seg','numSeeds_seg',...
    'hMin', 'hMax', 'delH');
helpdlg('Seed voxles were saved! - take note of the Interpolation Factor')
    
% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
global seed_voxels numSeeds seed_voxels_seg numSeeds_seg
global hMin hMax delH caseID peak_t peak_v mat_4D IntpFactor
saveRes = strcat('C:\Users\SHPro\Desktop\vessle_3D_results');
pathNameload = strcat(saveRes,filesep,caseID,filesep);
if 7==exist(pathNameload)
    [file,PathName] = uigetfile(strcat('*.mat'),'Select seed MAt file',...
        strcat(pathNameload,caseID,'_seeds_t00','.mat'),...
        'MultiSelect', 'off');
else
    [file,PathName] = uigetfile(strcat('*.mat'),'Select seed MAt file','MultiSelect', 'off');
end
temp = load(strcat(PathName,file));
seed_voxels = temp.seed_voxels;
numSeeds = temp.numSeeds;
seed_voxels_seg = temp.seed_voxels_seg;
numSeeds_seg = temp.numSeeds_seg;
hMin = temp.hMin;
hMax = temp.hMax;
delH = temp.delH;

set(handles.edit1, 'String', num2str(numSeeds-1));
set(handles.edit2, 'String', num2str(hMin));
set(handles.edit3, 'String', num2str(delH));
set(handles.edit4, 'String', num2str(hMax));
set(handles.edit12, 'String', num2str(numSeeds_seg)); % The first seed is reserved for overlapped
edit1_Callback(hObject, eventdata, handles)
[peak_v,peak_t]=max(mat_4D,[],4); %

seed_voxels = IntpFactor * seed_voxels;
seed_voxels_seg = IntpFactor * seed_voxels_seg;
helpdlg('Seed voxles were loaded and scaled based on Interpolation factor!')

% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
global smooth_flag
smooth_flag = get(hObject,'Value');



function th_txt_Callback(hObject, eventdata, handles)
global th_vis plotFlag
global MtAIP_axial
global imgMIP th_vis_flag forPlotInput

th_vis_flag = 1;
th_vis = str2num(get(handles.th_txt,'String'));

if plotFlag == 1
    forPlotInput = double(imgMIP > th_vis);
elseif plotFlag == 2
    forPlotInput = double(MtAIP_axial > th_vis);
end

imageDraw(handles,forPlotInput)
