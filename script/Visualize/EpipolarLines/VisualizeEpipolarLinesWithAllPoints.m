%%
% In this file we visualize and save the epipolar lines associated with all the points
% taken from the sequences in consecutive frames
%%
clear all;
close all;

addpath(genpath('../../Tools/'));


%% Para
FrameGap = 1; % gap between a pair of frames

%% Load Seq Information
temp = load('../../../Data/SeqList.mat');
SeqList = temp.SeqList;

model_type = lower('fundamental');   % Model name

seq_range = 1:length(SeqList);

if ~exist("EpipolarLinesVideosWithAllPoints",'dir')
    mkdir("EpipolarLinesVideosWithAllPoints");
end
for s_i = 1:length(SeqList)

    SeqName = SeqList{s_i}; % sequence name


    %%% Load Hypotheses
    save_path = fullfile('../../../Results/Hypotheses/',model_type);

    gt_filepath = fullfile('../../../Data/',[SeqName,'_Tracks']);
    temp = load(gt_filepath);
    Data = temp.Data;
    %% Define the VideoWriter in order to save the segmented frames
    videoPath = strcat('EpipolarLinesVideosWithAllPoints/',SeqName,'.avi');
    vid = VideoWriter(videoPath);
    vid.FrameRate = 5;
    vid.Quality = 100;
    open(vid) %open the writer
    %%
    for f_i = 1:Data.nFrames-FrameGap

        %% Select points visible on both frames
        visible_pts_ind = Data.visibleSparse(:,f_i) & Data.visibleSparse(:,f_i+1);

        matchedPoints1 = Data.ySparse(1:2,visible_pts_ind,f_i)';
        matchedPoints2 = Data.ySparse(1:2,visible_pts_ind,f_i+FrameGap)';

        %% estimate F with all the points
        [F,supp] = estimateFundamentalMatrix(matchedPoints1,matchedPoints2);
        %%
        ind = int2str(f_i);
        if strlength(ind)==2
            image_path = fullfile('../../../../OriginalSequence/',[SeqName,'/0000',ind,'.png']);
        else
            image_path = fullfile('../../../../OriginalSequence/',[SeqName,'/00000',ind,'.png']);
        end
        I1 = imread(image_path);
        ax1 = subplot(211);
        imshow(I1); 
        title('Inliers and Epipolar Lines in First Image'); hold on;
        plot(matchedPoints1(supp,1),matchedPoints1(supp,2),'go')

        epiLines = epipolarLine(F',matchedPoints2(supp,:));
        points = lineToBorderPoints(epiLines,size(I1));
        line(points(:,[1,3])',points(:,[2,4])');

        ind = int2str(f_i+1);
        if strlength(ind)==2
            image_path = fullfile('../../../../OriginalSequence/',[SeqName,'/0000',ind,'.png']);
        else
            image_path = fullfile('../../../../OriginalSequence/',[SeqName,'/00000',ind,'.png']);
        end
        I2 = imread(image_path);
        ax2 = subplot(212); 
        imshow(I2);
        title('Inliers and Epipolar Lines in Second Image'); hold on;
        plot(matchedPoints2(supp,1),matchedPoints2(supp,2),'go')

        epiLines = epipolarLine(F,matchedPoints1(supp,:));
        points = lineToBorderPoints(epiLines,size(I2));
        line(points(:,[1,3])',points(:,[2,4])');
        truesize;
        
        %extract the scattered image from the figure
        frame = getframe(gcf).cdata;

        %write the frame on the video
        writeVideo(vid,frame)
        hold(ax1,'off')
        hold(ax2,'off')

    end
    close(vid)
end

