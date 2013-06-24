function [bin_times, isih] = interv_hist(t, ch, Fs, max_interv, gate)
%INTERV_HIST - Compute interspike interval histogram
% Usage:  [bin_times, isih] = interv_hist(t, ch, Fs, max_interv, gate)
%         t				event times (msec)
%	  	  ch			event channels
%         Fs            histogram sampling rate in Hz (default 20000)
%	      max_interv	maximum interval in histogram (msec) (default 20)
%	      gate			optional [onset end] of peristimulus time gate (msec)
%         bin_times     Vector of interspike intervals (msec)
%         isih          Interval histogram
%
% INTERV_HIST with no output argument plots the histogram
%

if nargin < 3, Fs = 20000; end
if nargin < 4, max_interv = 40; end
spike_chan = 1;

bw=1000/Fs;    % bin width in msec
bin_times=[0:bw:max_interv-bw/2];

if nargin >= 5,		% gate
	pst = pstimes(t, ch);
	t(find(ch~=spike_chan | pst < gate(1) | pst >= gate(2))) = NaN;
else 	t(find(ch~=spike_chan)) = NaN;
end

if nargout == 0,
	hist(diff(t), bin_times+bw/2);
	xlabel('Interspike Histogram (msec)')
	ylabel('Number of Intervals')
else
	isih = hist(diff(t), bin_times+bw/2);
end
