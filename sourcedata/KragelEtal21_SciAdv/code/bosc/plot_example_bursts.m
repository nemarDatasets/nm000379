function [hl, ho] = plot_example_bursts(theta_res)
% PLOT_EXAMPLE_BURSTS plots example high and low theta bursts
%
% INPUTS:  theta_res - structure, contains oscillation information from
%                      subject S1
%
% OUTPUTS: hl - handle for the reference line
%
%          ho - handle for the oscillation line

% Copyright (C) 2021  James Kragel

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

if ~exist('theta_res','var')
    load('theta_res.mat');
end

low_color = [0 102 204]/255;
high_color = [204 102 51]/255;

% high theta, ~8.5 Hz
fr = [2 1];
b = 1;
c = 2;
s = 1;
li = [7 10];


% load bosc info
res = get_bosc('S1', '..\..\scratch_data\');

chans = unique(theta_res(s).chans);

freqs = unique(theta_res(s).freqs( ...
    theta_res(s).chans == chans(c)));

is_osc = res(b,c).is_osc;

% load the data from this block and make an ied filter
fname = ['..\..\data\eeg\S1\S1_sl_block' num2str(b) '_bp.h5'];

pairs = h5read(fname, '/pairs');
fs = h5read(fname, '/srate');
is_hc = return_hc('S1', pairs);
hc_idx = find(is_hc);
hc_idx = hc_idx(c);

% load raw timeseries for this block
ieeg = h5read(fname, ...
    '/eeg', ...
    [hc_idx 1], ...
    [1 inf]);

% find ied timepoints to exclude
ied_times = round(res(b,c).ied_out.pos(res(b,c).ied_out.chan==hc_idx)*fs);
is_ied = false(1, size(is_osc,2));
is_ied(ied_times) = true;

tc = ones(1, fs);
is_ied = conv(is_ied, tc, 'same');
is_ied(is_ied>0)=1;

is_osc(:,is_ied==1) = false;

fig('FontSize', 16, ...
    'Font', 'Arial', ...
    'border','on', ...
    'units', 'centimeters', ...
    'width', 20, ...
    'height', 8);

set(gcf,'Position',[22.3308    9.7631    9.6573   12.2767]);

% find and plot contiguous segments
for f = 1:length(fr)
    subplot(2,1,f);
    
    f_idx = res(1,1).freqs == freqs(fr(f));
    this_osc = is_osc(f_idx,:);
    
    L = bwlabel(this_osc); %oscillation bursts
    [~, x] = find(L == li(f));
    
    [hl, ho] = plot_oscillation(ieeg, x, fs);
    hold on;

    
    if f == 1
        
        xtp = [ho.XData(1) ho.XData(1) + .5];
        ytp = [-200 -100];
        href_line = plot(xtp, [ytp(1) ytp(1)]);
        href_line.Color = 'k';
        href_line.LineWidth = 2;
        
        href_line = plot([xtp(1) xtp(1)], ytp);
        href_line.Color = 'k';
        href_line.LineWidth = 2;
        
        ht = text(mean(xtp), mean(ytp)-110, '500 ms');
        ht.HorizontalAlignment = 'center';
        ht.VerticalAlignment = 'baseline';
        
        ht = text(xtp(1)-.3, mean([ytp(1) ytp(1)+100]), '100 \muV');
        ht.HorizontalAlignment = 'left';
        ht.VerticalAlignment = 'middle';
        
        ho.Color =  high_color;

    end
    
    if f == 2
       
        ho.Color = low_color;

    end
                         
end


end

function [hl, ho] = plot_oscillation(ieeg, c, fs)

% ieeg is the raw timeseries

% c, the columns (samples) for the current oscillation

% fs, sample rate, for plotting a reference for 1s

% make an x, which contains c with 250 ms

padlen = fs/4;

x = [(min(c) - padlen):(min(c)-1) c (max(c) + 1):(max(c)+padlen)];
x(x<1) = []; x(x>size(ieeg,2))=[];
% background trace
hl = plot(x/fs, ieeg(x));
hl.Color = 'k';
hold on;
ho = plot(c/fs, ieeg(c));
ho.Color = 'r';
ho.LineWidth = 2;

set(gca, 'XLim', ...
    [min(x/fs) max(x/fs)], ...
    'Visible', 'off', ...
    'YLim', [-200 200]);

end