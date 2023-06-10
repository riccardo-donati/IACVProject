%%
% In this file we just build the video trend of the silhouettes of (GRIC, GRIC_Tuned and RobustResiduals)
% that we found in ModelSelectionSilhouettePlotsNormalized.m
% The video is saved in e.g. for GRIC ModelSelectionSilhouettesNormalized/GRIC/{Seq Name}
%%
%% Build silhouette videos from images of silhouettes frame by frame.
clear all;
close all;

temp = load('../../../Data/SeqList.mat');
SeqList = temp.SeqList;
models = ["GRIC","GRIC Tuned", "RobustResiduals"];

%% build videos 
for model=1:length(models)
    model_type = models(model);
    for seq = 1:length(SeqList)
        seqData =  load(strcat("../../../Data/",SeqList{seq},'_Tracks.mat'));
        Data = seqData.Data;
        seqName = SeqList(seq);
        if model_type ~= "GRIC Tuned" && model_type ~="RobustResiduals" && model_type ~="GRIC"
            seqPath = 'ModelSilhouettesNormalized/'+model_type+'/'+seqName+'/';
        else
            seqPath = 'ModelSelectionSilhouettesNormalized/'+model_type+'/'+seqName+'/';
        end
        videoPath =  strcat(seqPath, seqName,'.avi');
        v = VideoWriter(videoPath);
        v.FrameRate = 5;
        v.Quality = 100;
        open(v) %open the writer
        for f_i = 1:Data.nFrames-1
            %% Select file by file and compose the video
            filename = "frames_"+f_i+"-"+(f_i+1)+".png";
            img = imread(strcat(seqPath, filename));
            imshow(img);
            title(strcat(model_type,': silhouettes of',{' '},seqName), 'Interpreter','none');
            hold on;

            %extract the image 
            frame = getframe(gcf).cdata;
            
            %write the frame on the video
            writeVideo(v,frame)
            hold off;
        end
        fprintf("%s: video %s created successfully\n", model_type, videoPath);
        close(v)
    end
end