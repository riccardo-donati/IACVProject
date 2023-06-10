%%
% In this file we visualize and save the Kmeans results after the process of
% Clustering in an interactive manner where we can choose the clustering technique
% and the embedings obtained stored in Results/MoSeg/Archivio
%%
%% Check the subset constrained motion segmentation results
clear all;
close all;

max_NumHypoPerFrame = 500;
FrameGap = 1;

gamma_range = [1e-2];

%% Load Seq Information
temp = load('../../../Data/SeqList.mat');
SeqList = temp.SeqList;

seq_range = 1:length(SeqList);

error_avg = [];
error_med = [];
valid_gamma = [];

    
%Input the desired model
model_types = ["Subset", "SubsetOnlyHF","Fundamental","Homography","FundamentalA","FundamentalT"];
map = containers.Map([1,2,3,4,5,6],model_types);
prompt = "What is desired the model [1 = Subset, 2 = SubsetHF, 3 = Fundamental, 4 = Homography, 5 = FundamentalA, 6 = FundamentalT] ?  ";
x = input(prompt);
model_type = map(x);
fprintf("%s selected\n",model_type);

%Input the desired experiment
prompt = "What is desired experiment (from 1 to 8) ? ";
exp = input(prompt);
fprintf("%d selected\n",exp);

%Load the best alpha for the desired model
model_alphas = load("../../CheckPerf/best_model_alphas.mat").model_alphas;
Alpha = model_alphas(model_type);

if model_type == "Fundamental" || model_type == "Homography" || model_type == "FundamentalA" || model_type == "FundamentalT"
    gamma = -1;
else
    gamma = gamma_range(1);
end
Para = [];
Para.MaxItr = 8;   % max iteration
Para.gamma = gamma;

result_path = fullfile('../../../Results/MoSeg/Archivio/',int2str(exp),'/',model_type);

%load the desired segmentation file
if gamma == -1
    result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g.mat',...
    max_NumHypoPerFrame,Alpha));
else
    result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g_gamma-%g.mat',...
    max_NumHypoPerFrame,Alpha,gamma));
end

temp = load(result_filepath);

error = temp.error;
ClusterIdx = temp.ClusterIdx;

%set the position of the figure
figure('Position', [50, 50, 1024, 512]);

ax1 = subplot(2,1,1);
pos1 = get(ax1,'Position');

ax2 = subplot(2,1,2);
pos2 = get(ax2,'Position');

for seq = 1:length(SeqList)
    SeqName = SeqList{seq};
    
    %load the informations about the given sequence
    gt_filepath = fullfile('../../../Data/',[SeqName,'_Tracks.mat']);
    temp = load(gt_filepath);
    Data = temp.Data;
    
    %extract the useful informations
    num_frames = Data.nFrames;
    ySparse = Data.ySparse;
    %% Define the VideoWriter in order to save the segmented frames
    videoPath = strcat('SegmentationVideos/',SeqName,'.avi');
    v = VideoWriter(videoPath);
    v.FrameRate = 5;
    v.Quality = 100;
    open(v) %open the writer
    
    %iterate on the frames of the sequence
    for i=1:Data.nFrames
        ind = int2str(i);
        if strlength(ind)==2
            image_path = fullfile('../../../../OriginalSequence/',[SeqName,'/0000',ind,'.png']);
        else
            image_path = fullfile('../../../../OriginalSequence/',[SeqName,'/00000',ind,'.png']);
        end
        %%
        % read the single frame
        image = imread(image_path);
        
        % take the figure coordinates of the points in the given frame
        x = ySparse(1,:,i);
        y = ySparse(2,:,i);
        
        pred = ClusterIdx{seq}; % Clustering prediction
        true = Data.GtLabel; %Ground Truth
        
        %show the same image 2 times, the first for the prediction, the
        %second for the ground truth
        ax1 = subplot(2,1,1);
        
        imshow(image);
        title(strcat('Kmeans Result of',{' '},SeqName,' with',{' '},model_type,' frame:', int2str(i)),'Interpreter','none');
        hold(ax1,'on')

        ax2 = subplot(2,1,2);
     
        imshow(image);
        title(strcat('Ground Truth Result of',{' '},SeqName,' with',{' '},model_type,' frame:', int2str(i)),'Interpreter','none');
        hold(ax2,'on')
        
        %overimpose the segmentations on the 2 images
        gscatter(ax1,x,y,pred,'rgbym','...',10,'off',"","")
        gscatter(ax2,x,y,true,'rgbym','...',10,'off',"","");
        set(ax1,'Position',pos1)
        set(ax1,'FontSize',12)
        set(ax2,'Position',pos2)
        set(ax2,'FontSize',12)

        %extract the scattered image from the figure
        frame = getframe(gcf).cdata;
        
        %write the frame on the video
        writeVideo(v,frame)
        if i == 5 %save only the first of each sequence
            if ~exist("SegmentationImages", 'dir')
                mkdir("SegmentationImages")
            end
            if ~exist("SegmentationImages/"+model_type, 'dir')
                mkdir("SegmentationImages/"+model_type);
            end
            
            saveas(ax1, "SegmentationImages/"+model_type+"/"+SeqName+"KmeansResult.png");
        end
        
        hold(ax1,'off')
        hold(ax2,'off')
        hold off
    end
    %close the video writer
    close(v)
end

fprintf('Overall Miss Classification Rate = %.2f\n',100*mean(error(end,:)));





