function [] = Neuro_WaveAnalysis(tempWaves)

% FD(n) = s(n) - s(n - 1)
% SD(n) = FD(n) - FD(n - 1)

% s = spike waveform points
% n = n number of points

FSDEvals = zeros(length(tempWaves),3);
for wavi = 1:length(FSDEvals)
    tempFD = firstDer(tempWaves(:,wavi));
    tempSD = secondDer(tempFD);
    tempFDmin = min(tempFD);
    tempSDmin = min(tempSD);
    tempSDmax = max(tempSD);
    
    FSDEvals(wavi,1) = tempFDmin;
    FSDEvals(wavi,2) = tempSDmin;
    FSDEvals(wavi,3) = tempSDmax;
    
    
end

% s = tempWaves(:,1);


function FD = firstDer(s)
FD = zeros(length(s),1);
for n = 2:length(s)
    FD(n) = s(n) - s(n - 1);
end

return


function SD = secondDer(s)
SD = zeros(length(s),1);
for n2 = 2:length(FD)
    SD(n2) = FD(n2) - FD(n2 - 1);
end

return

% hold on
% plot(s,'.--')
% plot(FD,'r.--')
% plot(SD,'k.--')
% line([0 32],[0 0],'LineStyle','--','Color','k');

FDmin = min(FD);
FDmax = max(FD);
SDmin = min(SD);
SDmax = max(SD);

princomp([FDmin,SDmin,SDmax]);

% Method 4: FDmax, SDmin, SDmax
X = FSDEvals;
opts = statset('Display','final');

[idx,~,~,~] = kmeans(X,2,...
                    'Distance','cityblock',...
                    'Replicates',10);

test1 = mahal(X((idx == 1),:),X)
test2 = mahal(X((idx == 2),:),X)

malot = [test1;test2];

clust1 = X((idx ==1),1)
clust2 = X((idx ==2),1)

[h,p] = kstest2(clust1,clust2)

[s,~] = silhouette(X,idx,'sqeuclid') 

group = [repmat({'First'}, length(test1), 1); repmat({'Second'}, length(test2), 1)];
boxplot(malot, group,'orientation','horizontal')

figure()


plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
plot(ctrs(:,1),ctrs(:,2),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot(ctrs(:,1),ctrs(:,2),'ko',...
     'MarkerSize',12,'LineWidth',2)
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')





% Calculate SDmin, SDmax, and FDmax FDmin

% Run Kmeans function

% Bray-Curtis similarity index (Lian et al., 2010: Signal Process 90:684-8)

x = tempWaves(:,1);
y = tempWaves(:,2);

% for i = 1:length(x)
%     tempNum(i) = abs(x(i) - y(i))
%     tempDen(i) = abs(x(i)) + abs(y(i))
% end
% 
% out = 1 - (sum(tempNum)/sum(tempDen))

out2 = 1 - (sum(abs(x - y)) / sum(abs(x) + abs(y)));



return