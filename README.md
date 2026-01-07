
Prereqs: MATLAB with "wavelet" and "signal processing" toolboxes.

To start using this app, make sure you have MATLAB installed with the above toolboes. If you need to install the toolboxes, start matlab and click the "apps" tab at the top, select "get more apps", search for "wavelet" and "signal processing" toolboxes and select install. Download the .mapp file from this site and double click to run it in matlab.

Instructions:

1. Choose file to import, and use the dialog to select the *.xlsx file containing oscillations data with a time column labeled "Time [s]" or "Time [min]" and columns of signal for one or more regions. The file may contain multiple sheets.


Control panel steps/options

2. Select sheet within the file that you want to use from the dropdown under 'Control' panel.

3. You should see the data from the sheet populate the top table. You can click on a cell or column in the table and use the mouse or keyboard arrow keys to navigate the columns. The plot on the right will update with the currently selected column. If you hold shift, you may select multiple columns to see multiple plots which may be useful for selecting time windows next.

4. Select of unselect the "remove outliers" checkbox if your data contains spike artifacts. This will detect large amplitude, narrow spikes less than 2 time steps wide and ignore them.   

5. Update time window start and end values. 

6. Right click table to add / del time windows. You can right click and select "add row" to ad multiple time windows. For now, you can only delete the last time window. Make sure you set up your time windows before spending a lot of time modifying the analysis values since changing the time windows will reset that thresholds. 

7. Checking/unchecking 'detrend' will toggle between removing the slope as determined by the detected troughs.

8. [TODO] Selecting the sensor from the dropdown will eventually change how silent phase is calculated, appropriate to each sensor type. This is not yet implemented.

9. Peakdetect (default) uses Sdao's algorythm for detecting peaks and troughs. Findpeaks() uses a MATLAB function that may eventually be more robust, but requires more testing.

10. Update tHresholds: The threshold for peak detection can be modified individually for each region/row below, but the initial values are automatically determined for each region based on the standard deviation. The default is to mset the paakdetect threshold to 1 standard deviation. You can change the default to a different factor by entering it and clicking "update threshold". Warning: this will reset all the thresholds, so try to do this first before spending time upating indiviual rows. 


Edit analysis parameters:

11. Select row(s)/region(s) in the analysis table to plot. You can select multiple rows. Update the threshold value to detect desired peaks/troughs and inspect them on the plot. 

12. Check the 'ignore' box for any regions you don't wish to analyze (i.e. noisy/bad regions). This can also help avoid errors that are generated during analysis if few/no peaks are detected.

13. View the wavelet tab for each region if desired.

12. Click Analyze to refresh output table with calculated values. 

12. Save to output file
