%%
% In this file, we save the silhouette plots obtained from choosing always 
% the same model among "fundamental", "fundamentala", "fundamentalt", "homography"
% we also visualize and save the silhouette score for the entire sequence
% both divided by motion and of all the 
% This is the version in which the fit of the model is not normalized in
% the number of the points
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
%% models silhouette
%Input the desired geometric model to analyze
models = ["fundamental", "fundamentala", "fundamentalt", "homography"];
scoresMotion = []; 
silhouetteMotion = [];
silhouetteMotionAvg = [];
silhouetteAvg = [];
for model=1:length(models)
    fprintf("%s selected\n",models(model));
    model_type = models(model);
    [fitModel, resModel, ~, ~, ~] = getModelParam(model_type);
    for seq = 1:length(SeqList)
        seqData =  load(strcat("../../../Data/",SeqList{seq},'_Tracks.mat'));
        Data = seqData.Data;
        gt = seqData.Data.GtLabel;
        nMotions = max(gt);
        silhouettes = [];
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
                % in case points of a certain motion are not visible in the
                % pair of frames or correspondences are less than 8, 
                % model is [0 0 0; 0 0 0; 0 0 1]
                if size(matchedPoints1{m}, 2) < 8
                    fprintf("Warning: %s in frames pair (%d-%d) for motion %d has %d correspondences\n", ...
                        SeqList{seq}, f_i, f_i+FrameGap, m, size(matchedPoints1{m},2));
                end
                %%%Models{m} = estimateFundamentalMatrix(matchedPoints1,matchedPoints2);
                try
                    Models{m} = feval(fitModel, [matchedPoints1{m}; matchedPoints2{m}]);
                catch 
                    Models{m} = [0 0 0; 0 0 0; 0 0 1];
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
                    ReswrtMotion{m}{n} = feval(resModel, Models{n},[matchedPoints1{m}; matchedPoints2{m}]);
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
        silhouetteMotionAvg{model}{seq} = nanmean(silhouettes,2); 
        silhouetteAvg{model}{seq} = mean(silhouetteMotionAvg{model}{seq}) ;
        
        % table with silhouette coefficient for every motion in every frame
        T = array2table(silhouettes);
        T.Properties.VariableNames = frames(1:Data.nFrames-FrameGap);
        T.Properties.RowNames = motions(1:nMotions);
    
        cell = strcat('A',int2str(n1),':W',int2str(n2));
        n1 = n1 + 7;
        n2 = n2 + 7;
        writetable(T,"ModelSilhouettes/"+model_type+"/silhouetteScores_"+model_type+".xls",'WriteRowNames',true,'Sheet',1,'Range',cell)
    end
    

    %% plot every frame silhouette of all clips for every model
   
    imageFolder = strcat("ModelSilhouettes/",models(model));
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
            for motion=1:nMotions
                if size(clip{frame}{motion}, 1) > maxLength
                    maxLength = size(clip{frame}{motion}, 1);
                end        
            end
            y = zeros(maxLength, nMotions);
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

            title("frames "+frame+" "+(frame+1), 'Interpreter','none');
            saveas(f, imgName+"/frames_"+frame+"-"+(frame+1)+".png");
            fprintf("%s: saved frames %d-%d for %s successfully\n",models(model), frame, (frame+1), SeqList{seq});
            
        end
    end
        %% plot svg silhoette 
    f = figure();
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    avgs = [];
    avgsTot = [];
    for i=1:22
    tmp = silhouetteMotionAvg{model}{i}';
    tmpAvg = mean(silhouetteMotionAvg{model}{i}');
        while size(tmp) ~= 5
            tmp = [tmp 0];
        end
        avgsTot = [avgsTot tmpAvg];
        avgs = [avgs; tmp];
    end
    bar(categorical(SeqList),avgs',1);
    title(strcat("avg silhouette score for model ",models(model)));
    legend("motion1","motion2","motion3","motion4","motion5");
    set(gca,'TickLabelInterpreter','none')
    saveas(f, "ModelSilhouettes/"+models(model)+"AvgSilhouetteScorePerMotion.png");
    
    f = figure();
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    bar(categorical(SeqList),avgsTot,1);
    title(strcat("avg silhouette score for model ",models(model)));
    set(gca,'TickLabelInterpreter','none')
    saveas(f, "ModelSilhouettes/"+models(model)+"AvgSilhouetteScore.png");
    
end

