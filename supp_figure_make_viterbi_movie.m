clear
close all

addpath(genpath('./lib'))

% Load data
readPath = '.\data\burst_analysis_data\';

% load spot data
load([readPath 'spot_struct.mat'],'spot_struct')

% set path to Dropbox 
DropboxFolder = 'S:\Nick\Dropbox (Personal)\OptoDropbox\';
DataFolder = 'P:\Optogenetics\Data\PreProcessedData\';

% specify list of embryos to use
prefixList = {'2021-06-15-optoknirps_eve4_6_embryo27','2021-06-15-optoknirps_eve4_6_embryo28',...
              '2021-06-16-optoknirps_eve4_6_embryo30','2021-06-17-optoknirps_eve4_6_embryo35'};
%%
prefix_id_vec = 5:8;
gr = [213 108 85]/256;
cmap = brewermap([],'Set2');
rd  = cmap(2,:);
bl  = cmap(3,:);
mcp_ch = 2;
pt_ch = 1;
% set defaults
MaxmRNAFluo = 3e4;
MaxProteinFluo = 1.65e4;
CurrentNC = 14;
K = 2; %number of states to display 


for p = 1:length(prefixList)

    MaxRadius = 22;       
  
    %Embryo ID and Folder location
    Prefix = prefixList{p};
    FISHPath = [DataFolder Prefix filesep];
    PrefixDropboxFolder = [DropboxFolder prefixList{p}]; 

    %Load the data
    load([PrefixDropboxFolder,'\CompiledParticles.mat'])
    load([PrefixDropboxFolder,'\' Prefix,'_lin.mat'])
    load([PrefixDropboxFolder,'\Ellipses.mat'])
    load([PrefixDropboxFolder,'\FrameInfo.mat'])

    % load histone array
    histone_array = imreadStack([FISHPath Prefix '-His.tif']);

    % set write path 
    % Set write directory
    if K == 2
        subdir = 'two_state_movie';   
    elseif K == 3
        subdir = 'three_state_movie';    
    end
    
    writePath = [PrefixDropboxFolder filesep];
    mkdir(writePath);
        
    StackPath = [writePath '\' subdir '\'];
    mkdir(StackPath)        
            
    % get dim info
    xDim = FrameInfo.PixelsPerLine;
    yDim = FrameInfo.LinesPerFrame;

    [px, py] = meshgrid(1:xDim,1:yDim);
    FrameRange = nc14:length(ElapsedTime);

    %%% Make some useful ref cells    
    master_id_vec = [spot_struct.masterID];
    masterID = prefix_id_vec(p);
    set_ids = find(master_id_vec==masterID);
    setID = spot_struct(set_ids(1)).setID;
    particle_id_vec = [spot_struct.particleID];
    
%     v_particle_vec = [viterbi_fit_struct.ParticleID];
%     setID = 10;
    if K == 2
        color_cell = {rd,brighten([115 143 193]/256,-0.15),[115 143 193]/256};
%         color_cell = {brighten([233 108 85]/256,-0.25),brighten([115 143 193]/256,-0.15),[115 143 193]/256};
    elseif K == 3
        color_cell = {[115 143 193]/256, [122 169 116]/256, [213 108 85]/256};
    end
    extant_particle_cell = cell(1,length(FrameRange));
    for i = set_ids
        ParticleID = spot_struct(i).particleID;        
        if isnan(spot_struct(i).z_viterbi_interp)
            continue
        end    
        start_i = find(~isnan(spot_struct(i).fluo),1);
        p_frames = spot_struct(i).frames(start_i:end);
        for f = 1:length(p_frames)
            f_filter = FrameRange==p_frames(f);
            if sum(f_filter) > 0 
                extant_list = extant_particle_cell{f_filter};
                extant_list = [extant_list ParticleID];
                extant_particle_cell{f_filter} = extant_list;
            end
        end
    end
    
    %Iterate Through Frames
    for CurrentFrame = FrameRange    
        %Highlight Active Regions with Viterbi State    
        %Track pixel assignments
        NucleusStateMat = zeros(yDim,xDim,3);
        NucleusIDMat = zeros(yDim,xDim);
        NucleusDistMat = ones(yDim,xDim)*MaxRadius;    
        % Loop through ALL nuclei (inactive and active) and assign pixels
        % to a nucleus
        for s = 1:length(schnitzcells)
            CurrentEllipse= schnitzcells(s).cellno(...
                            schnitzcells(s).frames==...
                            CurrentFrame);
            x =  Ellipses{CurrentFrame}(CurrentEllipse,1)+1;
            y =  Ellipses{CurrentFrame}(CurrentEllipse,2)+1;        
            if isempty(x)
                continue
            end 
            distances = ((px-x).^2 + (py-y).^2).^.5;
            candidate_indices = NucleusDistMat > distances; 
            %Record Fluorescence        
            NucleusIDMat(candidate_indices) = s;
            NucleusDistMat(candidate_indices) = distances(candidate_indices);
        end
        % reference frame index cell
        extant_particles = extant_particle_cell{CurrentFrame-nc14 + 1};
        ParticlesToShow = [];
        for i = 1:length(extant_particles)
            ParticleID = extant_particles(i);
            Nucleus = spot_struct(particle_id_vec==ParticleID).nucleusID;
            frames = spot_struct(particle_id_vec==ParticleID).frames;
%             frame_key = round(linspace(1,length(spot_struct(particle_id_vec==ParticleID).fluo_interp),length(frames)));
            z_vec = spot_struct(particle_id_vec==ParticleID).z_viterbi_raw;                        
            z_state = z_vec(frames==CurrentFrame);
            try
                CurrentEllipse=...
                    schnitzcells(Nucleus).cellno(...
                    schnitzcells(Nucleus).frames==...
                    CurrentFrame);                    
                %Record Fluorescence
                filter = NucleusIDMat==Nucleus;
                color_vec = fliplr(color_cell{z_state});
                for k = 1:3
                    nc_slice = NucleusStateMat(:,:,k);
                    nc_slice(filter) = color_vec(k);
                    NucleusStateMat(:,:,k) = nc_slice;
                end
                ParticlesToShow = [ParticlesToShow Nucleus];
            catch
                warning('inconsistent indexing')
            end
        end
    
        %Now Draw Nucleus Borders for active nuclei
        NucleusBorderMat = zeros(size(NucleusIDMat));
        window = 2;
        for i = ParticlesToShow
            %get coordinates of nucleus patch
            if sum(sum(NucleusIDMat==i)) > 0
                x_vec = reshape(px(NucleusIDMat==i),[],1);
                y_vec = reshape(py(NucleusIDMat==i),[],1);
                for j = 1:length(x_vec)
                    metric = sum(sum(NucleusIDMat(max(1,y_vec(j)-window):min(y_vec(j)+window,yDim),...
                             max(1,x_vec(j)-window):min(x_vec(j) + window,xDim))));
                    if metric~= i*(2*window+1)^2
                        NucleusBorderMat(y_vec(j),x_vec(j)) = 1;
                    end
                end
            end
        end
        %Prevent overlap between fluroescence mask and borders
        for k = 1:3
            nc_slice = NucleusStateMat(:,:,k);
            nc_slice(NucleusBorderMat==0) = 0;
            NucleusStateMat(:,:,k) = nc_slice;
        end
        %Make a maximum projection of the mRNA channel
        n_char = length(num2str(CurrentFrame));
        
%         index_string = [z_string(1:3-n_char) num2str(CurrentFrame)];
        D=dir([FISHPath,'\',Prefix,'_',sprintf('%03d',CurrentFrame),'_ch' sprintf('%02d',mcp_ch) '*']);
        %Do not load the first and last frame as they are black
        ImageTemp = imreadStack([D(1).folder filesep D(1).name]);
        ImageTemp = ImageTemp(:,:,2:end-1);
       
        mRNAImage = double(max(ImageTemp,[],3));        
        mRNAImage = mRNAImage/MaxmRNAFluo;
        
        D=dir([FISHPath,'\',Prefix,'_',sprintf('%03d',CurrentFrame),'_ch' sprintf('%02d',pt_ch) '*']);
        
        %Do not load the first and last frame as they are black
        ImageTemp = imreadStack([D(1).folder filesep D(1).name]);
        ImageTemp = ImageTemp(:,:,2:end-1);
       
        proteinImage = imgaussfilt(double(max(ImageTemp,[],3)),2);        
        proteinImage = proteinImage/MaxProteinFluo;

        %Load the corresponding histone image
        HistoneImage= histone_array(:,:,CurrentFrame); 
        HistoneImage = double(HistoneImage)/max(double(HistoneImage(:)));
        %find most likely state
%         NucleusStateMat = NucleusStateMat;
        %Overlay all channels    
        greenImage =   proteinImage ;%NucleusStateMat(:,:,2);    % proteinImage ++
        MCPChannel = mat2gray(greenImage,[0,1.1]);
        HistoneImage =  HistoneImage; % NucleusStateMat(:,:,3))
        HistoneChannel=mat2gray(HistoneImage, [0,1.1]);
           
        BlueChannel = mat2gray(mRNAImage ,[0 1.1]);            %+ NucleusStateMat(:,:,1)/2
        
        ImOverlay=cat(3,HistoneChannel,MCPChannel,BlueChannel);

        % try adding protein info in borders
        for k = 1:3
            im_slice = ImOverlay(:,:,k);
            nc_slice = NucleusStateMat(:,:,3-k+1);
            im_slice(nc_slice~=0) = nc_slice(nc_slice~=0);%gr(k)*proteinImage(nc_slice==0)*2;
            ImOverlay(:,:,k) = im_slice;
        end
        
    
        OverlayFig = figure('Visible','off');
        clf
        imshow(fliplr(ImOverlay))   
        
        write_time = round(ElapsedTime(CurrentFrame)-ElapsedTime(FrameRange(1)));
        write_time_prev = round(ElapsedTime(CurrentFrame-1)-ElapsedTime(FrameRange(1)));
        
        xlim([0,xDim])    
        ylim([16,yDim-16])
        
        if ismember(write_time,[15 20 36]) && ~ismember(write_time_prev,[15 20 36])
            saveas(gcf,[writePath '\nc',num2str(CurrentNC),...
                '-',sprintf('%03d',CurrentFrame), '_' num2str(write_time) 'min.tif']);   
        end
        
        text(20,yDim-50,[sprintf('%02d',write_time),...
            ' min'],'Color','k','FontSize',16,'BackgroundColor',[1,1,1,.5])
        drawnow            
    
        saveas(gcf,[StackPath '\nc',num2str(CurrentNC),...
            '-',sprintf('%03d',CurrentFrame), '.tif']);   
        close all
               
    end
end