%%
% In this file we visualize and save the epipolar lines associated with the minimum number of points
% taken from all the points in the sequences in consecutive frames
%%
clear all;
close all;

addpath(genpath('../../../Tools/'));


%% Para
FrameGap = 1; % gap between a pair of frames
max_NumHypoPerFrame = 500;  % Max number of hypotheses sampled from each frame pair

%% Load Seq Information
temp = load('../../../Data/SeqList.mat');
SeqList = temp.SeqList;

model_type = lower('fundamental');   % Model name
[ fitfn resfn degenfn psize numpar ] = getModelParam(model_type);

seq_range = 1:length(SeqList);

if ~exist("EpipolarLinesVideosOfBestF",'dir')
    mkdir("EpipolarLinesVideosOfBestF");
end
for s_i = 1:length(SeqList)

    SeqName = SeqList{s_i}; % sequence name

    %%% Load Hypotheses
    save_path = fullfile('../../../Results/Hypotheses/',model_type);

    hypo_filepath = fullfile(save_path,sprintf('Hypo_RandSamp_Sparse_seq-%s_nHypo-%d.mat',SeqName,max_NumHypoPerFrame));
    temp = load(hypo_filepath,'Hypos');

    Model = temp.Hypos;

    gt_filepath = fullfile('../../../Data/',[SeqName,'_Tracks']);
    temp = load(gt_filepath);
    Data = temp.Data;

    %% Define the VideoWriter in order to save the segmented frames
    videoPath = strcat('EpipolarLinesVideosOfBestF/',SeqName,'.avi');
    vid = VideoWriter(videoPath);
    vid.FrameRate = 5;
    vid.Quality = 100;
    open(vid) %open the writer
    %%
    for f_i = 1:Data.nFrames-FrameGap

        mdl_idx = find(Model.r == f_i)';    % all hypotheses in current frame pair
        Res = [];   % residual w.r.t. hypotheses
        %%
        for h_i = mdl_idx

            r = Model.r(h_i);   % first frame
            v = Model.v(h_i);   % second frame

            %% Select points visible on both frames
            visible_pts_ind = Data.visibleSparse(:,f_i) & Data.visibleSparse(:,f_i+1);

            y1 = Data.ySparse(:,visible_pts_ind,r);
            y2 = Data.ySparse(:,visible_pts_ind,v);

            %% Normalise raw correspondences.
            dat_img_1 = normalise2dpts(y1);
            dat_img_2 = normalise2dpts(y2);
            normalized_data = [ dat_img_1 ; dat_img_2 ];

            %% Calculate Residual
            Res = [Res feval(resfn,Model.H(:,h_i),normalized_data)]; %xFx' is residual (actually samson error i think)

        end
        %% Take the best model ( lower sum of residual errors )
        sumOfRes = sum(Res,1);
        [argvalue, argmin] = min(sumOfRes);
        %% prova prendendendo il min F dei primi due frame con i suoi 8 punti(supp)
        F = Model.H(:,argmin);
        F = reshape(F,3,3);
        supp = Model.supp(:,argmin);
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
        plot(Data.ySparse(1,supp,f_i)',Data.ySparse(2,supp,f_i)','go')

        epiLines = epipolarLine(F',Data.ySparse(1:2,supp,f_i+1)');
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
        plot(Data.ySparse(1,supp,f_i+1),Data.ySparse(2,supp,f_i+1),'go')

        epiLines = epipolarLine(F,Data.ySparse(1:2,supp,f_i)');
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



