
ScSize = get(0, 'ScreenSize');
CheckClusterFigHandle = figure('Name', 'test', 'NumberTitle', 'Off', 'Position', [ ScSize(3)/4 1 ScSize(3)/2 ScSize(4)*0.9]);
text_font_size = 7;

[acf,lags,bounds] = autocorr(spktms)

acorr_bin_msec = 4;
[histvals,x] = AutoCorr(spktms, acorr_bin_msec, 250); % This is peter's C AutoCorr

%acorr(floor(length(acorr)/2+1)) = 0;    % set 0 lag to 0 for scaling
% Cannot use bar-- there is some matlab bug that prevents it's use in this subplot.


plot(acf,lags);			% show acorr

figure(CheckClusterFigHandle);
xlabel(['msec (' num2str(acorr_bin_msec) 'msec binsize)'],'FontSize',text_font_size);
ylabel('rate','FontSize',text_font_size);
h = title('Autocorrelation','FontSize',text_font_size);
drawnow;
figure(CheckClusterFigHandle);
set(gca,'FontSize',text_font_size);



acorr_bin_msec = 1;
[histvals,x] = AutoCorr(spktms, acorr_bin_msec, 50); % This is peter's C AutoCorr

%acorr(floor(length(acorr)/2+1)) = 0;    % set 0 lag to 0 for scaling
% Cannot use bar-- there is some matlab bug that prevents it's use in this subplot.

bar(x,histvals);			% show acorr

figure(CheckClusterFigHandle);
xlabel(['msec (' num2str(acorr_bin_msec) 'msec binsize)'],'FontSize',text_font_size);
ylabel('rate','FontSize',text_font_size);
% h = title('Autocorrelation','FontSize',text_font_size);
drawnow;
figure(CheckClusterFigHandle);
set(gca,'FontSize',text_font_size);



% create a spike train
sf=1000; % sampling frequency
x=0:.25*sf; % sample for 1/4 sec
y=zeros(size(x));
l=length(y);
r=randperm(l);
y(r(1:fix(.25*l)))=1; % spike train
% compute interspike interval and
% instantaneous firing frq
ix=find(y);
v=[0 sf./diff(ix)];
% ... and create a simple plot
stairs(x(ix),v,'r');
hold on;
stem(x,.2*sf*y,'.'); % scale for disp
xlabel('time [ms]');
ylabel('events/inst firing frq [hz]');
