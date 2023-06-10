%%
% In this file, starting from the models selected by ModelSelection phase (seq_best_models.mat)
% we save the silhouette plots obtained for the criteria of interest (GRIC, GRIC_Tuned and RobustResiduals)
% we also visualize and save the silhouette score for the entire sequence
% both divided by motion and of all the clip
%%
%% Obtain silhouette graphs and tables for all clips, for all models
clear all;
close all;

max_NumHypoPerFrame = 500;
motions = ["motion1", "motion2","motion3","motion4", "motion5"];
n1 = 1;
n2 = 7;

%% Load Seq Information
temp = load('../../../Data/SeqList.mat');
SeqList = temp.SeqList;
FrameGap = 1;
addpath '..\..\..\Tools\multigs\model_specific';
%% get model preferences 
models = load("seq_best_models.mat").seq_best_models;

scoresMotion = []; 
silhouetteMotion = [];
silhouetteMotionAvg = [];
silhouetteAvg = [];

model_selection = ["GRIC","GRIC Tuned", "RobustResiduals"];
mod_sel = containers.Map([1, 2], [2, 8]);
for modelIndex=1:length(mod_sel)
    i = mod_sel(modelIndex); % 2 is GRIC Tuned, 8 is RobustResiduals
    for seq = 1:length(SeqList)
        preferences = models{i,seq};
        seqData =  load(strcat("../../../Data/",SeqList{seq},'_Tracks.mat'));
        Data = seqData.Data;
        gt = seqData.Data.GtLabel;
        nMotions = max(gt);
        silhouettes = [];
        sizes = [];
        for f_i = 1:Data.nFrames-FrameGap
            %% Select points for every motion visible on both frames
            visible_pts_ind = Data.visibleSparse(:,f_i) & Data.visibleSparse(:,f_i+1);
            Models = [];
            matchedPoints1 = [];
            matchedPoints2 = [];
            % obtain a binary array containing indices of point divided per
            % motions and for every motion estimate the F
            for m = 1:nMotions
                pointsInMotion = visible_pts_ind & (gt==m);
                % matched points of motion m visible in frames f_i and f_i+1
                matchedPoints1{m} = Data.ySparse(:,pointsInMotion,f_i);
                matchedPoints2{m} = Data.ySparse(:,pointsInMotion,f_i+FrameGap);
                if size(matchedPoints1{m}, 2) > 1
                    matchedPoints1{m} = normalise2dpts(Data.ySparse(:,pointsInMotion,f_i));
                    matchedPoints2{m} = normalise2dpts(Data.ySparse(:,pointsInMotion,f_i+FrameGap));
                end

                sizes = [sizes size(matchedPoints1{m},2)];
               
                model_type = preferences(f_i,m);
                try
                    % if model is "None", NaN
                    [fitModel, resModel, ~, ~, ~] = getModelParam(model_type);
                    Models{m} = feval(fitModel, [matchedPoints1{m}; matchedPoints2{m}]);
                catch 
                    Models{m} = NaN;
                end
            end
            
            sizeMin =max(8,min(sizes));
            for m=1:nMotions
                try
                    model_type = preferences(f_i,m);
                    [fitModel, resModel, ~, ~, ~] = getModelParam(model_type);
    
                    Models{m} = feval(fitModel, [matchedPoints1{m}(:,1:sizeMin); matchedPoints2{m}(:,1:sizeMin)]);
                catch 
                    Models{m} = NaN;
                end
            end
            % now compute residuals of every point wrt the model for every motion.
            % This is needed to compute a and b for every point
            % ReswrtMotion: matrix m x n: every row corresponds to the residuals of
            % a motion wrt all the m motions. e.g. 2 motions => 2x2:
            % [res(points_motion_1,model_1) res(points_motion_1,model_2);
            % res(points_motion_2,model_1) res(points_motion_2,model_2) ]
            ReswrtMotion = [];
            for m=1:nMotions
                for n=1:nMotions

                    if isnan(Models{n})
                        ReswrtMotion{m}{n} = NaN(size(matchedPoints1{m},2),1);
                    else
                        model_type = preferences(f_i,n);
                        [~, resModel, ~, ~, ~] = getModelParam(model_type);
                        ReswrtMotion{m}{n} = feval(resModel, Models{n},[matchedPoints1{m}; matchedPoints2{m}]);
                    end
                end
            end
            
            %now compute a and b for every point and its silhouette score     
            for m=1:nMotions
                b = [];
                for n=1:nMotions
                    if m==n % vector of a's for motion m
                        a = ReswrtMotion{m}{n};                    
                    else %residuals of points of motion m wrt motion n
                        if isempty(b)
                            b = ReswrtMotion{m}{n};
                        else
                            b = min(b, ReswrtMotion{m}{n});
                        end
                    end
                end
                % we have a and b for every point of motion m. Now silhouette
                scoresMotion{seq}{f_i}{m} = (b - a) ./ max(a, b);
                silhouetteMotion{seq}{f_i}{m} = mean(scoresMotion{seq}{f_i}{m});
                silhouettes(m, f_i) = mean(scoresMotion{seq}{f_i}{m});
                frames(f_i) = "frames "+f_i+" "+(f_i+1);
            end
        end
        silhouetteMotionAvg{modelIndex}{seq} = nanmean(silhouettes,2); 
        silhouetteAvg{modelIndex}{seq} = mean(silhouetteMotionAvg{modelIndex}{seq}) ;
        
    end
    
    
    %% plot every frame silhouette of all clips for every model
    imageFolder = strcat("ModelSelectionSilhouettesNormalized/",model_selection(modelIndex));
    for seq=1:length(SeqList)
        clip = scoresMotion{seq};
        nMotions = size(clip{1}, 2);
        imgName = imageFolder+"/"+SeqList{seq};
        if ~exist(imgName, 'dir')
            mkdir(imgName);
        end
    
        for frame=1:size(clip, 2)
            f = figure('Name', "frames "+frame+" "+(frame+1), ...
                'NumberTitle','off', 'visible', 'off');
            maxLength = 0;        
            colors = [];
            for motion=1:nMotions
                colors = [colors; 0.2 0.5 0.8];
                if size(clip{frame}{motion}, 1) > maxLength
                    maxLength = size(clip{frame}{motion}, 1);
                end        
            end
            y = zeros(maxLength,nMotions);
            for motion=1:nMotions
                temp = clip{frame}{motion};
                if size(temp, 1) < maxLength 
                    temp((size(temp,1)+1):maxLength) = NaN;
                end
                %ascend if horizontal graph descend if vertical
                y(:,motion) = sort(temp,"ascend");
            end
            n =  size(y, 1);  % Number of rows
            m =  size(y, 2);    % Number of columns

            % Create the increasing matrix
            vec = 1:n;
            X = vec + (0:n:(m-1)*n).';
            % Plotting the bar plot
            barh(fliplr(X'),y,10)
            categoryNames = {'Motion 1', 'Motion 2', 'Motion 3','Motion 4','Motion 5','Motion 6'};
            legend(categoryNames)
            
            colormap(colors);
            %groups = motions(1:nMotions);
            %set(gca, 'XTick', 1:numel(y), 'XTickLabel', groups);
            title("frames "+frame+" "+(frame+1), 'Interpreter','none');
            saveas(f, imgName+"/frames_"+frame+"-"+(frame+1)+".png");
            fprintf("%s: saved frames %d-%d for %s successfully\n",model_selection(modelIndex), frame, (frame+1), SeqList{seq})
        end
    end
    
    %% plot avg silhoette 
    f = figure();
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    avgs = [];
    avgsTot = [];
    for i=1:22
        tmp = silhouetteMotionAvg{modelIndex}{i}';
        tmpAvg = mean(silhouetteMotionAvg{modelIndex}{i}');
        while size(tmp) ~= 5
            tmp = [tmp 0];
        end
        avgsTot = [avgsTot tmpAvg];
        avgs = [avgs; tmp];
    end
    bar(categorical(SeqList),avgs',1);
    title(strcat("avg silhouette score with ",model_selection(modelIndex)));
    legend("motion1","motion2","motion3","motion4","motion5");
    set(gca,'TickLabelInterpreter','none')
    saveas(f, "ModelSelectionSilhouettesNormalized/"+model_selection(modelIndex)+"AvgSilhouetteScorePerMotion.png");
    
    f = figure();
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    bar(categorical(SeqList),avgsTot,1);
    title(strcat("avg silhouette score for model ",model_selection(modelIndex)));
    set(gca,'TickLabelInterpreter','none')
    saveas(f, "ModelSelectionSilhouettesNormalized/"+model_selection(modelIndex)+"AvgSilhouetteScore.png");
end
