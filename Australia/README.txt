
 
The program GISS_temp_Aus.m plots distributions of monthly mean surface air temperature anomalies over the Australian mainland and Tasmania. 

Data supplied by NASA GISS.

Written for Octave v4.0.3 by @ShaneWhiteEng

Method:
1. Monthly mean temp anomaly observations are imported, organised by month into five time periods labelled T1 to T5, and 
missing data is removed(represented by 9999).

2. Means are calculated for each month for each time period and displayed, allowing confirmation that means over a time 
period are zero in each month.

3. The time period T2 (defined as years 1951 to 1980) is set as the reference period.

4. The mean and standard deviation of the reference period is calculated.

5. The Z-score, relative to the mean and standard deviation of the reference period, is calculated for each month in each
time period and organised into seasons.

6. Histograms are computed for each month in each time period and organised into seasons.

7. Histograms are smoothed by a single pole low pass filter.

8. Histograms are plotted for each season displaying the reference time periods and all subsequent time periods.



