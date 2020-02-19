%-------------------------------------------------------------------------------
% 
% This program plots distributions of monthly mean surface air temperature
% anomalies over the Australian mainland and Tasmania. 
% Data supplied by NASA GISS.
%
% Written for Octave v4.0.3 by @ShaneWhiteEng
%
% Method:
% 1. Monthly mean temp anomaly observations are imported, organised by month
%    into five time periods labelled T1 to T5, and missing data is removed
%    (represented by 9999). 
% 2. Means are calculated for each month for each time period and displayed,
%    allowing confirmation that means over a time period are zero in each month.
% 3. The time period T2 (defined as years 1951-80) is set as the reference 
%    period.
% 4. The mean and standard deviation of the reference period is calculated.
% 5. The Z-score, relative to the mean and standard deviation of the reference 
%    period, is calculated for each month in each time period and organised 
%    into seasons.
% 6. Histograms are computed for each month in each time period and organised
%    into seasons.
% 7. Histograms are smoothed by a single pole low pass filter.
% 8. Histograms are plotted for each season displaying the reference time 
%    periods and all subsequent time periods.
%
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
% Tidy up
%-------------------------------------------------------------------------------

hold off;
close all;
clear all;


%-------------------------------------------------------------------------------
% Definitions
%-------------------------------------------------------------------------------

%define arbitrary time periods T1 to T5
global T1a = 1880; %T1 = start of 1880 to end of 1909
global T1b = 1909; 
global T2a = 1951; %T2 = start of 1951 to end of 1980
global T2b = 1980;
global T3a = 1990; %T3 = start of 1990 to end of 1999
global T3b = 1999;
global T4a = 2000; %T4 = start of 2000 to end of 2009 
global T4b = 2009;
global T5a = 2010; %T5 = start of 2010 to end of 2019
global T5b = 2019;

%define bin arrays for histograms (i.e x-axis)
%bin array interval
global bin_interval = 0.01;
%bin array extremities
max_z_score = 6;
%construct bin array
global bins = [(-max_z_score + bin_interval) : bin_interval : max_z_score]'; 

%define time constant of low pass filter
global lpf_time_constant = 0.1;

%define months
jan = 1;
feb = 2;
mar = 3;
apr = 4;
may = 5;
jun = 6;
jul = 7;
aug = 8;
sep = 9;
oct = 10;
nov = 11;
dec = 12;

%define function for importing tempoerature anomaly observation data, %organising it time periods T1 to T5, and removing any missing values
function [T1_month T2_month T3_month T4_month T5_month] = sep_period...
(current_month);
global T1a;
global T1b;
global T2a;
global T2b;
global T3a;
global T3b;
global T4a;
global T4b;
global T5a;
global T5b;
global years;
index1 = index2 = index3 = index4 = index5 = 1;
 for i = 1 : size(years, 1)
  if current_month(i, 1) != 9999
   if years(i, 1) >= T1a && years(i, 1) <= T1b
    T1_month(index1, 1) = current_month(i, 1);
    index1 = index1 + 1;
   endif
   if years(i, 1) >= T2a && years(i, 1) <= T2b
    T2_month(index2, 1) = current_month(i, 1);
    index2 = index2 + 1;
   endif
   if years(i, 1) >= T3a && years(i, 1) <= T3b
    T3_month(index3, 1) = current_month(i, 1);
    index3 = index3 + 1;
   endif
   if years(i, 1) >= T4a && years(i, 1) <= T4b
    T4_month(index4, 1) = current_month(i, 1);
    index4 = index4 + 1;
   endif
   if years(i, 1) >= T5a && years(i, 1) <= T5b
    T5_month(index5, 1) = current_month(i, 1);
    index5 = index5 + 1; 
   endif
  endif
 endfor
endfunction

%define function to calculate Z-score
function [z_score_value] = z_score(current_month, ref_mean, ref_std)
 z_score_value = (current_month .- ref_mean) ./ ref_std;
endfunction

%define function to compute histograms and scale so area under curve is unity
function [scaled_histogram] = histogram(z_score_val);
 global bins;
 unscaled_histogram = histc(z_score_val, bins);
 scaled_histogram = unscaled_histogram ./ trapz(bins, unscaled_histogram);
endfunction

%define function to determine share of observations very hot (z-score >= 3),
%hot (z-score >= 1 and < 3), mild (z-score < 1 and > -1) and 
%cold (z-score <= -1)
function [cold_share, mild_share, hot_share, vhot_share] = shares...
(season_pdf, bin_array);
 colds = zeros(size(season_pdf, 1), 1);
 milds = zeros(size(season_pdf, 1), 1);
 hots = zeros(size(season_pdf, 1), 1);
 vhots = zeros(size(season_pdf, 1), 1);    
 index1 = index2 = index3 = index4 = 1;
 for n = 1 : size(bin_array, 1)
  if bin_array(n, 1) <= -1
   colds(index1, 1) = season_pdf(n, 1);
   index1 = index1 + 1;
  endif
  if bin_array(n, 1) > -1 && bin_array(n, 1) < 1
   milds(index2, 1) = season_pdf(n, 1);
   index2 = index2 + 1;
  endif
  if bin_array(n, 1) >=1 && bin_array(n, 1) < 3
   hots(index3, 1) = season_pdf(n, 1);
   index3 = index3 + 1;
  endif
  if bin_array(n, 1) >= 3
   vhots(index4, 1) = season_pdf(n, 1);
   index4 = index4 + 1;
  endif
 endfor
 total = trapz(season_pdf);
 cold_share = trapz(colds) / total * 100;
 mild_share = trapz(milds) / total * 100;
 hot_share = trapz(hots) / total * 100;
 vhot_share = trapz(vhots) / total * 100;
endfunction

%provide switch to turn low pass filter on and off
global lpf_on = 1;

%define low pass filter function
function[output] = lpf(input);
 global bin_interval;
 global lpf_time_constant;
 global lpf_on;
 if lpf_on == 1
  output = zeros(size(input, 1), 1);
  for n = 2:size(input, 1)
   output(n, 1) = bin_interval / lpf_time_constant * (input(n, 1) - ...
   output((n - 1), 1)) + output((n - 1), 1);
  endfor
 else
  output = input;
 endif
endfunction

%define chart plotting parameters
plot_axes = [-6, 6, 0, 0.5];
year_posx = -5.5;
year_posy = 0.52;
info_text_posx = 6.2;
info_text_posy = 0.4;
xticks = -6 : 1 : 6;


%-------------------------------------------------------------------------------
% Program
%-------------------------------------------------------------------------------

% Import data
[lat1, lat2, long1, long2, year, jan_temp, feb_temp, mar_temp, apr_temp, ...
may_temp, jun_temp, jul_temp, aug_temp, sep_temp, oct_temp, nov_temp, ...
dec_temp] = textread("Australia.SBBX1880.Ts.GHCNv4.250.txt", "%f %f %f %f %f \
%f %f %f %f %f %f %f %f %f %f %f %f","headerlines", 2);
%define years global
global years = year;


% Uncomment to print lat and long extremeties of observations
%corner1_lat_long_lat_long = [min(lat1) min(long1) min(lat2) min(long2)]
%corner2_lat_long_lat_long = [max(lat1) min(long1) max(lat2) min(long2)]
%corner3_lat_long_lat_long = [max(lat1) max(long1) max(lat2) max(long2)]
%corner4_lat_long_lat_long = [min(lat1) max(long1) min(lat2) max(long2)]


% Organise imoported observations by month into time periods T1 to T5
[T1_jan T2_jan T3_jan T4_jan T5_jan] = sep_period(jan_temp);
[T1_feb T2_feb T3_feb T4_feb T5_feb] = sep_period(feb_temp);
[T1_mar T2_mar T3_mar T4_mar T5_mar] = sep_period(mar_temp);
[T1_apr T2_apr T3_apr T4_apr T5_apr] = sep_period(apr_temp);
[T1_may T2_may T3_may T4_may T5_may] = sep_period(may_temp);
[T1_jun T2_jun T3_jun T4_jun T5_jun] = sep_period(jun_temp);
[T1_jul T2_jul T3_jul T4_jul T5_jul] = sep_period(jul_temp);
[T1_aug T2_aug T3_aug T4_aug T5_aug] = sep_period(aug_temp);
[T1_sep T2_sep T3_sep T4_sep T5_sep] = sep_period(sep_temp);
[T1_oct T2_oct T3_oct T4_oct T5_oct] = sep_period(oct_temp);
[T1_nov T2_nov T3_nov T4_nov T5_nov] = sep_period(nov_temp);
[T1_dec T2_dec T3_dec T4_dec T5_dec] = sep_period(dec_temp);


% Compute and print means
means_T1 = [mean(T1_jan) mean(T1_feb) mean(T1_mar) mean(T1_apr) mean(T1_may) ...
mean(T1_jun) mean(T1_jul) mean(T1_aug) mean(T1_sep) mean(T1_oct) ...
mean(T1_nov) mean(T1_dec)]

means_T2 = [mean(T2_jan) mean(T2_feb) mean(T2_mar) mean(T2_apr) mean(T2_may) ...
mean(T2_jun) mean(T2_jul) mean(T2_aug) mean(T2_sep) mean(T2_oct) ...
mean(T2_nov) mean(T2_dec)]

means_T3 = [mean(T3_jan) mean(T3_feb) mean(T3_mar) mean(T3_apr) mean(T3_may) ...
mean(T3_jun) mean(T3_jul) mean(T3_aug) mean(T3_sep) mean(T3_oct) ...
mean(T3_nov) mean(T3_dec)]

means_T4 = [mean(T4_jan) mean(T4_feb) mean(T4_mar) mean(T4_apr) mean(T4_may) ...
mean(T4_jun) mean(T4_jul) mean(T4_aug) mean(T4_sep) mean(T4_oct) ...
mean(T4_nov) mean(T4_dec)]

means_T5 = [mean(T5_jan) mean(T5_feb) mean(T5_mar) mean(T5_apr) mean(T5_may) ...
mean(T5_jun) mean(T5_jul) mean(T5_aug) mean(T5_sep) mean(T5_oct) ...
mean(T5_nov) mean(T5_dec)]


% Set reference period
ref_period_means = means_T2;
ref_period_stds = [std(T2_jan) std(T2_feb) std(T2_mar) std(T2_apr) ...
std(T2_may) std(T2_jun) std(T2_jul) std(T2_aug) std(T2_sep) std(T2_oct) ...
std(T2_nov) std(T2_dec)];


% Uncomment to display mean and standard deviation of reference period, in
% units of degress Celsius
%T2 = vertcat(T2_jan, T2_feb, T2_mar, T2_apr, T2_may, T2_jun, T2_jul, ...
%T2_aug, T2_sep, T2_oct, T2_nov, T2_dec);
%T2_mean = mean(T2)
%T2_std = std(T2)


% Calculate the Z-score of observations for each month in each period
Z_T1_jan = z_score(T1_jan, ref_period_means(jan), ref_period_stds(jan));
Z_T1_feb = z_score(T1_feb, ref_period_means(feb), ref_period_stds(feb));
Z_T1_mar = z_score(T1_mar, ref_period_means(mar), ref_period_stds(mar));
Z_T1_apr = z_score(T1_apr, ref_period_means(apr), ref_period_stds(apr));
Z_T1_may = z_score(T1_may, ref_period_means(may), ref_period_stds(may));
Z_T1_jun = z_score(T1_jun, ref_period_means(jun), ref_period_stds(jun));
Z_T1_jul = z_score(T1_jul, ref_period_means(jul), ref_period_stds(jul));
Z_T1_aug = z_score(T1_aug, ref_period_means(aug), ref_period_stds(aug));
Z_T1_sep = z_score(T1_sep, ref_period_means(sep), ref_period_stds(sep));
Z_T1_oct = z_score(T1_oct, ref_period_means(oct), ref_period_stds(oct));
Z_T1_nov = z_score(T1_nov, ref_period_means(nov), ref_period_stds(nov));
Z_T1_dec = z_score(T1_dec, ref_period_means(dec), ref_period_stds(dec));

Z_T2_jan = z_score(T2_jan, ref_period_means(jan), ref_period_stds(jan));
Z_T2_feb = z_score(T2_feb, ref_period_means(feb), ref_period_stds(feb));
Z_T2_mar = z_score(T2_mar, ref_period_means(mar), ref_period_stds(mar));
Z_T2_apr = z_score(T2_apr, ref_period_means(apr), ref_period_stds(apr));
Z_T2_may = z_score(T2_may, ref_period_means(may), ref_period_stds(may));
Z_T2_jun = z_score(T2_jun, ref_period_means(jun), ref_period_stds(jun));
Z_T2_jul = z_score(T2_jul, ref_period_means(jul), ref_period_stds(jul));
Z_T2_aug = z_score(T2_aug, ref_period_means(aug), ref_period_stds(aug));
Z_T2_sep = z_score(T2_sep, ref_period_means(sep), ref_period_stds(sep));
Z_T2_oct = z_score(T2_oct, ref_period_means(oct), ref_period_stds(oct));
Z_T2_nov = z_score(T2_nov, ref_period_means(nov), ref_period_stds(nov));
Z_T2_dec = z_score(T2_dec, ref_period_means(dec), ref_period_stds(dec));

Z_T3_jan = z_score(T3_jan, ref_period_means(jan), ref_period_stds(jan));
Z_T3_feb = z_score(T3_feb, ref_period_means(feb), ref_period_stds(feb));
Z_T3_mar = z_score(T3_mar, ref_period_means(mar), ref_period_stds(mar));
Z_T3_apr = z_score(T3_apr, ref_period_means(apr), ref_period_stds(apr));
Z_T3_may = z_score(T3_may, ref_period_means(may), ref_period_stds(may));
Z_T3_jun = z_score(T3_jun, ref_period_means(jun), ref_period_stds(jun));
Z_T3_jul = z_score(T3_jul, ref_period_means(jul), ref_period_stds(jul));
Z_T3_aug = z_score(T3_aug, ref_period_means(aug), ref_period_stds(aug));
Z_T3_sep = z_score(T3_sep, ref_period_means(sep), ref_period_stds(sep));
Z_T3_oct = z_score(T3_oct, ref_period_means(oct), ref_period_stds(oct));
Z_T3_nov = z_score(T3_nov, ref_period_means(nov), ref_period_stds(nov));
Z_T3_dec = z_score(T3_dec, ref_period_means(dec), ref_period_stds(dec));

Z_T4_jan = z_score(T4_jan, ref_period_means(jan), ref_period_stds(jan));
Z_T4_feb = z_score(T4_feb, ref_period_means(feb), ref_period_stds(feb));
Z_T4_mar = z_score(T4_mar, ref_period_means(mar), ref_period_stds(mar));
Z_T4_apr = z_score(T4_apr, ref_period_means(apr), ref_period_stds(apr));
Z_T4_may = z_score(T4_may, ref_period_means(may), ref_period_stds(may));
Z_T4_jun = z_score(T4_jun, ref_period_means(jun), ref_period_stds(jun));
Z_T4_jul = z_score(T4_jul, ref_period_means(jul), ref_period_stds(jul));
Z_T4_aug = z_score(T4_aug, ref_period_means(aug), ref_period_stds(aug));
Z_T4_sep = z_score(T4_sep, ref_period_means(sep), ref_period_stds(sep));
Z_T4_oct = z_score(T4_oct, ref_period_means(oct), ref_period_stds(oct));
Z_T4_nov = z_score(T4_nov, ref_period_means(nov), ref_period_stds(nov));
Z_T4_dec = z_score(T4_dec, ref_period_means(dec), ref_period_stds(dec));

Z_T5_jan = z_score(T5_jan, ref_period_means(jan), ref_period_stds(jan));
Z_T5_feb = z_score(T5_feb, ref_period_means(feb), ref_period_stds(feb));
Z_T5_mar = z_score(T5_mar, ref_period_means(mar), ref_period_stds(mar));
Z_T5_apr = z_score(T5_apr, ref_period_means(apr), ref_period_stds(apr));
Z_T5_may = z_score(T5_may, ref_period_means(may), ref_period_stds(may));
Z_T5_jun = z_score(T5_jun, ref_period_means(jun), ref_period_stds(jun));
Z_T5_jul = z_score(T5_jul, ref_period_means(jul), ref_period_stds(jul));
Z_T5_aug = z_score(T5_aug, ref_period_means(aug), ref_period_stds(aug));
Z_T5_sep = z_score(T5_sep, ref_period_means(sep), ref_period_stds(sep));
Z_T5_oct = z_score(T5_oct, ref_period_means(oct), ref_period_stds(oct));
Z_T5_nov = z_score(T5_nov, ref_period_means(nov), ref_period_stds(nov));
Z_T5_dec = z_score(T5_dec, ref_period_means(dec), ref_period_stds(dec));


% Organise Z-scores into seasons
Z_T1_djf = vertcat(Z_T1_dec, Z_T1_jan, Z_T1_feb);
Z_T2_djf = vertcat(Z_T2_dec, Z_T2_jan, Z_T2_feb);
Z_T3_djf = vertcat(Z_T3_dec, Z_T3_jan, Z_T3_feb);
Z_T4_djf = vertcat(Z_T4_dec, Z_T4_jan, Z_T4_feb);
Z_T5_djf = vertcat(Z_T5_dec, Z_T5_jan, Z_T5_feb);

Z_T1_mam = vertcat(Z_T1_mar, Z_T1_apr, Z_T1_may);
Z_T2_mam = vertcat(Z_T2_mar, Z_T2_apr, Z_T2_may);
Z_T3_mam = vertcat(Z_T3_mar, Z_T3_apr, Z_T3_may);
Z_T4_mam = vertcat(Z_T4_mar, Z_T4_apr, Z_T4_may);
Z_T5_mam = vertcat(Z_T5_mar, Z_T5_apr, Z_T5_may);

Z_T1_jja = vertcat(Z_T1_jun, Z_T1_jul, Z_T1_aug);
Z_T2_jja = vertcat(Z_T2_jun, Z_T2_jul, Z_T2_aug);
Z_T3_jja = vertcat(Z_T3_jun, Z_T3_jul, Z_T3_aug);
Z_T4_jja = vertcat(Z_T4_jun, Z_T4_jul, Z_T4_aug);
Z_T5_jja = vertcat(Z_T5_jun, Z_T5_jul, Z_T5_aug);

Z_T1_son = vertcat(Z_T1_dec, Z_T1_oct, Z_T1_nov);
Z_T2_son = vertcat(Z_T2_dec, Z_T2_oct, Z_T2_nov);
Z_T3_son = vertcat(Z_T3_dec, Z_T3_oct, Z_T3_nov);
Z_T4_son = vertcat(Z_T4_dec, Z_T4_oct, Z_T4_nov);
Z_T5_son = vertcat(Z_T5_dec, Z_T5_oct, Z_T5_nov);


% Compute seasonal histograms
hist_T1_djf = histogram(Z_T1_djf);
hist_T2_djf = histogram(Z_T2_djf);
hist_T3_djf = histogram(Z_T3_djf);
hist_T4_djf = histogram(Z_T4_djf);
hist_T5_djf = histogram(Z_T5_djf);

hist_T1_mam = histogram(Z_T1_mam);
hist_T2_mam = histogram(Z_T2_mam);
hist_T3_mam = histogram(Z_T3_mam);
hist_T4_mam = histogram(Z_T4_mam);
hist_T5_mam = histogram(Z_T5_mam);

hist_T1_jja = histogram(Z_T1_jja);
hist_T2_jja = histogram(Z_T2_jja);
hist_T3_jja = histogram(Z_T3_jja);
hist_T4_jja = histogram(Z_T4_jja);
hist_T5_jja = histogram(Z_T5_jja);

hist_T1_son = histogram(Z_T1_son);
hist_T2_son = histogram(Z_T2_son);
hist_T3_son = histogram(Z_T3_son);
hist_T4_son = histogram(Z_T4_son);
hist_T5_son = histogram(Z_T5_son);


% Compute shares of Z-scores cold, mild, hot and vhot
[cold_T1_djf, mild_T1_djf, hot_T1_djf, vhot_T1_djf] = shares(hist_T1_djf, bins);
[cold_T2_djf, mild_T2_djf, hot_T2_djf, vhot_T2_djf] = shares(hist_T2_djf, bins);
[cold_T3_djf, mild_T3_djf, hot_T3_djf, vhot_T3_djf] = shares(hist_T3_djf, bins);
[cold_T4_djf, mild_T4_djf, hot_T4_djf, vhot_T4_djf] = shares(hist_T4_djf, bins);
[cold_T5_djf, mild_T5_djf, hot_T5_djf, vhot_T5_djf] = shares(hist_T5_djf, bins);

[cold_T1_mam, mild_T1_mam, hot_T1_mam, vhot_T1_mam] = shares(hist_T1_mam, bins);
[cold_T2_mam, mild_T2_mam, hot_T2_mam, vhot_T2_mam] = shares(hist_T2_mam, bins);
[cold_T3_mam, mild_T3_mam, hot_T3_mam, vhot_T3_mam] = shares(hist_T3_mam, bins);
[cold_T4_mam, mild_T4_mam, hot_T4_mam, vhot_T4_mam] = shares(hist_T4_mam, bins);
[cold_T5_mam, mild_T5_mam, hot_T5_mam, vhot_T5_mam] = shares(hist_T5_mam, bins);

[cold_T1_jja, mild_T1_jja, hot_T1_jja, vhot_T1_jja] = shares(hist_T1_jja, bins);
[cold_T2_jja, mild_T2_jja, hot_T2_jja, vhot_T2_jja] = shares(hist_T2_jja, bins);
[cold_T3_jja, mild_T3_jja, hot_T3_jja, vhot_T3_jja] = shares(hist_T3_jja, bins);
[cold_T4_jja, mild_T4_jja, hot_T4_jja, vhot_T4_jja] = shares(hist_T4_jja, bins);
[cold_T5_jja, mild_T5_jja, hot_T5_jja, vhot_T5_jja] = shares(hist_T5_jja, bins);

[cold_T1_son, mild_T1_son, hot_T1_son, vhot_T1_son] = shares(hist_T1_son, bins);
[cold_T2_son, mild_T2_son, hot_T2_son, vhot_T2_son] = shares(hist_T2_son, bins);
[cold_T3_son, mild_T3_son, hot_T3_son, vhot_T3_son] = shares(hist_T3_son, bins);
[cold_T4_son, mild_T4_son, hot_T4_son, vhot_T4_son] = shares(hist_T4_son, bins);
[cold_T5_son, mild_T5_son, hot_T5_son, vhot_T5_son] = shares(hist_T5_son, bins);


% Uncommment to print sum of integrations to check bin resolution is adequate
sum([cold_T1_djf, mild_T1_djf, hot_T1_djf, vhot_T1_djf])
sum([cold_T2_djf, mild_T2_djf, hot_T2_djf, vhot_T2_djf])
sum([cold_T3_djf, mild_T3_djf, hot_T3_djf, vhot_T3_djf])
sum([cold_T4_djf, mild_T4_djf, hot_T4_djf, vhot_T4_djf])
sum([cold_T5_djf, mild_T5_djf, hot_T5_djf, vhot_T5_djf])


% Smooth distribution plots
smoothed_T1_djf = lpf(hist_T1_djf);
smoothed_T2_djf = lpf(hist_T2_djf);
smoothed_T3_djf = lpf(hist_T3_djf);
smoothed_T4_djf = lpf(hist_T4_djf);
smoothed_T5_djf = lpf(hist_T5_djf);

smoothed_T1_mam = lpf(hist_T1_mam);
smoothed_T2_mam = lpf(hist_T2_mam);
smoothed_T3_mam = lpf(hist_T3_mam);
smoothed_T4_mam = lpf(hist_T4_mam);
smoothed_T5_mam = lpf(hist_T5_mam);

smoothed_T1_jja = lpf(hist_T1_jja);
smoothed_T2_jja = lpf(hist_T2_jja);
smoothed_T3_jja = lpf(hist_T3_jja);
smoothed_T4_jja = lpf(hist_T4_jja);
smoothed_T5_jja = lpf(hist_T5_jja);

smoothed_T1_son = lpf(hist_T1_son);
smoothed_T2_son = lpf(hist_T2_son);
smoothed_T3_son = lpf(hist_T3_son);
smoothed_T4_son = lpf(hist_T4_son);
smoothed_T5_son = lpf(hist_T5_son);


% Define ideal standard normal probability distribution function
ideal_pdf = normpdf(bins, 0, 0.95);


% Plot austral summer histograms 
%plot reference period
figure;subplot(2, 2, 1);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T2_djf, "color", "red");
box off;
xlabel("Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, ["Reference period " num2str(T2a) " - " ...
num2str(T2b)], 'fontsize', 12, 'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T2_djf, "%3.1f") "%\nMild = " num2str(mild_T2_djf, " \
%3.1f") "%\nHot = " num2str(hot_T2_djf, "%3.1f") "%\n\
VHot = " num2str(vhot_T2_djf, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%plot period T3
subplot(2, 2, 2);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T3_djf, "color", "red");
box off;
xlabel("Reference period Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, [num2str(T3a) " - " num2str(T3b)], ...
'fontsize', 12, 'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T3_djf, "%3.1f") "%\nMild = " num2str(mild_T3_djf, " \
%3.1f") "%\nHot = " num2str(hot_T3_djf, "%3.1f") "%\n\
VHot = " num2str(vhot_T3_djf, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%plot period T4
subplot(2, 2, 3);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T4_djf, "color", "red");
box off;
xlabel("Reference period Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, [num2str(T4a) " - " num2str(T4b)], ...
'fontsize', 12, 'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T4_djf, "%3.1f") "%\nMild = " num2str(mild_T4_djf, " \
%3.1f") "%\nHot = " num2str(hot_T4_djf, "%3.1f") "%\n\
VHot = " num2str(vhot_T4_djf, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%plot period T5
subplot(2, 2, 4);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T5_djf, "color", "red");
box off;
xlabel("Reference period Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, [num2str(T5a) " - " num2str(T5b)], ...
'fontsize', 12, 'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T5_djf, "%3.1f") "%\nMild = " num2str(mild_T5_djf, " \
%3.1f") "%\nHot = " num2str(hot_T5_djf, "%3.1f") "%\n\
VHot = " num2str(vhot_T5_djf, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%add page text
ha = axes("Position", [0 0 1 1], "Xlim", [0 1], "Ylim", [0 1], "Box", "off", ...
"Visible", "off", "Units", "normalized", "clipping", "off");
text(0.52, 0.98,"Monthly mean temperature distributions over the Australian \
mainland and Tasmania during austral summer", ...
"HorizontalAlignment", "center", "VerticalAlignment", "top",'fontsize', 16);
text(0.43, 0.54,"Cold months are between -1 and -3\nMild months are between -1 and 1\nHot months are between 1 and 3\nVery Hot (VHot) months are greater than or equal to 3", ...
"HorizontalAlignment", "left", "VerticalAlignment", "top",'fontsize', 12);


% Plot austral autumn histograms
%plot reference period
figure;subplot(2, 2, 1);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T2_mam, "color", "red");
box off;
xlabel("Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, ["Reference period " num2str(T2a) " - " ...
num2str(T2b)], 'fontsize', 12, 'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T2_mam, "%3.1f") "%\nMild = " num2str(mild_T2_mam, " \
%3.1f") "%\nHot = " num2str(hot_T2_mam, "%3.1f") "%\n\
VHot = " num2str(vhot_T2_mam, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%plot period T3
subplot(2, 2, 2);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T3_mam, "color", "red");
box off;
xlabel("Reference period Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, [num2str(T3a) " - " num2str(T3b)], ...
'fontsize', 12, 'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T3_mam, "%3.1f") "%\nMild = " num2str(mild_T3_mam, " \
%3.1f") "%\nHot = " num2str(hot_T3_mam, "%3.1f") "%\n\
VHot = " num2str(vhot_T3_mam, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%plot period T4
subplot(2, 2, 3);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T4_mam, "color", "red");
box off;
xlabel("Reference period Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, [num2str(T4a) " - " num2str(T4b)], ...
'fontsize', 12, 'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T4_mam, "%3.1f") "%\nMild = " num2str(mild_T4_mam, " \
%3.1f") "%\nHot = " num2str(hot_T4_mam, "%3.1f") "%\n\
VHot = " num2str(vhot_T4_mam, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%plot period T5
subplot(2, 2, 4);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T5_mam, "color", "red");
box off;
xlabel("Reference period Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, [num2str(T5a) " - " num2str(T5b)], ...
'fontsize', 12, 'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T5_mam, "%3.1f") "%\nMild = " num2str(mild_T5_mam, " \
%3.1f") "%\nHot = " num2str(hot_T5_mam, "%3.1f") "%\n\
VHot = " num2str(vhot_T5_mam, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%add page text
ha = axes("Position", [0 0 1 1], "Xlim", [0 1], "Ylim", [0 1], "Box", "off", ...
"Visible", "off", "Units", "normalized", "clipping", "off");
text(0.52, 0.98,"Monthly mean temperature distributions over the Australian \
mainland and Tasmania during austral autumn", ...
"HorizontalAlignment", "center", "VerticalAlignment", "top",'fontsize', 16);
text(0.43, 0.54,"Cold months are between -1 and -3\nMild months are between -1 and 1\nHot months are between 1 and 3\nVery Hot (VHot) months are greater than or equal to 3", ...
"HorizontalAlignment", "left", "VerticalAlignment", "top",'fontsize', 12);


% Plot austral winter histograms
%plot reference period
figure;subplot(2, 2, 1);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T2_jja, "color", "red");
box off;
xlabel("Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, ["Reference period " num2str(T2a) " - " ...
num2str(T2b)], 'fontsize', 12, 'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T2_jja, "%3.1f") "%\nMild = " num2str(mild_T2_jja, " \
%3.1f") "%\nHot = " num2str(hot_T2_jja, "%3.1f") "%\n\
VHot = " num2str(vhot_T2_jja, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%plot period T3
subplot(2, 2, 2);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T3_jja, "color", "red");
box off;
xlabel("Reference period Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, [num2str(T3a) " - " num2str(T3b)], 'fontsize', ...
12,'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T3_jja, "%3.1f") "%\nMild = " num2str(mild_T3_jja, " \
%3.1f") "%\nHot = " num2str(hot_T3_jja, "%3.1f") "%\n\
VHot = " num2str(vhot_T3_jja, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%plot period T4
subplot(2, 2, 3);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T4_jja, "color", "red");
box off;
xlabel("Reference period Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, [num2str(T4a) " - " num2str(T4b)], ...
'fontsize', 12, 'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T4_jja, "%3.1f") "%\nMild = " num2str(mild_T4_jja, " \
%3.1f") "%\nHot = " num2str(hot_T4_jja, "%3.1f") "%\n\
VHot = " num2str(vhot_T4_jja, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%plot period T5
subplot(2, 2, 4);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T5_jja, "color", "red");
box off;
xlabel("Reference period Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, [num2str(T5a) " - " num2str(T5b)], ...
'fontsize', 12,'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T5_jja, "%3.1f") "%\nMild = " num2str(mild_T5_jja, " \
%3.1f") "%\nHot = " num2str(hot_T5_jja, "%3.1f") "%\n\
VHot = " num2str(vhot_T5_jja, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%add page text
ha = axes("Position", [0 0 1 1], "Xlim", [0 1], "Ylim", [0 1], "Box", "off", ...
"Visible", "off", "Units", "normalized", "clipping", "off");
text(0.52, 0.98,"Monthly mean temperature distributions over the Australian \
mainland and Tasmania during austral winter", ...
"HorizontalAlignment", "center", "VerticalAlignment", "top",'fontsize', 16);
text(0.43, 0.54,"Cold months are between -1 and -3\nMild months are between -1 and 1\nHot months are between 1 and 3\nVery Hot (VHot) months are greater than or equal to 3", ...
"HorizontalAlignment", "left", "VerticalAlignment", "top",'fontsize', 12);


% Plot austral spring histograms
%plot reference period
figure;subplot(2, 2, 1);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T2_son, "color", "red");
box off;
xlabel("Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, ["Reference period " num2str(T2a) " - " ...
num2str(T2b)], 'fontsize', 12, 'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T2_son, "%3.1f") "%\nMild = " num2str(mild_T2_son, " \
%3.1f") "%\nHot = " num2str(hot_T2_son, "%3.1f") "%\n\
VHot = " num2str(vhot_T2_son, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%plot period T3
subplot(2, 2, 2);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T3_son, "color", "red");
box off;
xlabel("Reference period Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, [num2str(T3a) " - " num2str(T3b)], ...
'fontsize', 12,'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T3_son, "%3.1f") "%\nMild = " num2str(mild_T3_son, " \
%3.1f") "%\nHot = " num2str(hot_T3_son, "%3.1f") "%\n\
VHot = " num2str(vhot_T3_son, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%plot period T4
subplot(2, 2, 3);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T4_son, "color", "red");
box off;
xlabel("Reference period Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, [num2str(T4a) " - " num2str(T4b)], ...
'fontsize', 12, 'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T4_son, "%3.1f") "%\nMild = " num2str(mild_T4_son, " \
%3.1f") "%\nHot = " num2str(hot_T4_son, "%3.1f") "%\n\
VHot = " num2str(vhot_T4_son, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%plot period T5
subplot(2, 2, 4);
plot(bins, ideal_pdf, "color", "blue", bins, smoothed_T5_son, "color", "red");
box off;
xlabel("Reference period Z-score");
ylabel("Relative Frequency");
axis(plot_axes, "manual", "tic", "square"); 
current_axes = gca(); 
set(current_axes, "XTick", xticks);
text(year_posx, year_posy, [num2str(T5a) " - " num2str(T5b)], ...
'fontsize', 12, 'fontweight', 'bold');
text(info_text_posx, info_text_posy, ["Share of\nmonths -\n\n\
Cold = " num2str(cold_T5_son, "%3.1f") "%\nMild = " num2str(mild_T5_son, " \
%3.1f") "%\nHot = " num2str(hot_T5_son, "%3.1f") "%\n\
VHot = " num2str(vhot_T5_son, "%3.1f") "%"], 'fontname', 'arial', ...
'fontsize', 12);
grid on;

%add page text
ha = axes("Position", [0 0 1 1], "Xlim", [0 1], "Ylim", [0 1], "Box", "off", ...
"Visible", "off", "Units", "normalized", "clipping", "off");
text(0.52, 0.98,"Monthly mean temperature distributions over the Australian \
mainland and Tasmania during austral spring", ...
"HorizontalAlignment", "center", "VerticalAlignment", "top",'fontsize', 16);
text(0.43, 0.54,"Cold months are between -1 and -3\nMild months are between -1 and 1\nHot months are between 1 and 3\nVery Hot (VHot) months are greater than or equal to 3", ...
"HorizontalAlignment", "left", "VerticalAlignment", "top",'fontsize', 12);
