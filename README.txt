*************************************************************************

This code processes 3D Kinematics data for the Apraxia Battery of the 
Cognition and Action Laboratory (Director, Laurel J Buxbaum). 
This code was
 written by Aaron L Wong (03/05/2019, last modified 10/31/2019).

This code makes certain assumptions made as to 
the nature of the tasks being analyzed and the setup of the kinematic trackers. 
In general, we assume that there are 8 kinematic markers, placed as follows:
   1. Thumb
   2. Index Finger
   3. Middle Finger
   4. Back of the hand
   5. Wrist
   6. Elbow
   7. Shoulder
   8. Reference shoulder
. 
To change these assumptions, see ABload.m. To change assumptions about the
 task structure and the file name labels, see decodetrials.m.


*************************************************************************

To analyze any subset of data from a single subject, call ABBLOCKCOMPILE. This function analyzes all data within a single folder, and sorts the data by block. Each trial's data file will be loaded and compiled based on its filename. It is expected that trial names will include the following information: block type, subject number, and item/trial number. It is also assumed that the data come from the same subject, who is one of ("Control", "Patient", or "Model") groups. If the data files do not meet these assumptions, the data may not be loaded properly. Finally, it is assumed that any data files that are empty or have less than 1 second of data are "bad" and are ignored, and that if there are multiple reptitions of the same trial, that the last redo is the one to store (that is, it just processes the data alphabetically and if there are multiple trials with the same trial/item number, it keeps replacing the old data with the newer one). Because of how this function works, the most efficient means of processing all data for a single subject is just to select the "Raw" data folder.

The output of this function is a set of data (.mat) files, one per block, each containing a structure with data sorted by trial along with some metadata.

NOTE: As part of the ABBLOCKCOMPILE analysis, the user may be prompted to identify sign discontinuities in the dataset. These could arise because the wrong hemisphere was specified in the trakSTAR system, and so when the tracker crossed the axis of the transmitter the sign inappropriately changed along all axes. Discontinuities during the trial should be automatically detected, but user input is required to confirm that these points have been properly identified (particularly if the start or the end of the trial should have been marked, e.g., if the marker started off on the wrong side of the transmitter axis). Markers will need to be verified for each tracker for which at least one discontinuity has been detected.

This initial simplistic marking program can be controlled using a combination of the keyboard and mouse. The mouse inputs are: 
	left-click to add a mark
	right-click to remove a mark
The useful keyboard inputs are:
	"n" to scroll forward in the data trace
	"b" to scroll backward in the data trace
	"0" to request a mark at the start of the data trace
	"9" to request a mark at the end of the data trace
	"x" to quit once all the desired marks have been identified.

ALSO NOTE: ABBLOCKCOMPILE (and possibly other functions below) calls a function called DECODETRIALS to identify the block/trial name. Thus this file contains in it a hard-coded pairing of block/trial number to block/trial name. 
***This pairing is current as of 10/31/2019. If the item order is ever changed, this is the file that must be modified.***

~~~~

When there are a large number of data folders to be processed, BATCHCOMPILE can handle this. However, use of the BATCHCOMPILE function dependss on the file structure of the folders that you want to process. When all desired folders exist within a single parent folder, leaving the "paths" variable empty will result in GUI that allows you to select multiple data folders (within a single parent folder) to be processed. When data span multiple parent folders (e.g., processing the "Raw" folder for multiple participants), the full set of paths will need to be hard-coded in the "path" variable. Regardless of how the folders are chosen, the script will automatically run ABBLOCKCOMPILE on all indicated folders. Data can be saved in one of two ways; either by specifying the data folder when prompted, or if no path is chosen (i.e, you hit "cancel" in the GUI), the data for each folder will be automatically put into a folder one level up the file path (e.g., instead of putting it in the "Raw" folder it will put it in the subject folder). 


*************************************************************************


~~~FOR ANALYSIS WITH MARKING OF "FULL" ACTION~~~


To mark the data for a single subject/block, run the ABFINDFULLACTIONSUBJ function. When this script is run, a GUI will prompt the selection of the data file to be analyzed. The program will automatically process all trials in the block and display each trial in turn for marking. To mark the data, the each trial is plotted individually and the user must manually identify the start and end of the "full" action -- that is, from movement start (lift-off) to movement end (return to the table). Suggested marks are provided to the user based on times when the velocity drops below a set threshold. It is the responsibility of the user to verify these marks, as they will not necessarily be appropriate in all cases. Data marking is handled by the MARKDATAGUI data-marking sub-function (called automatically by ABFINDFULLACTIONSUBJ); specific instructions for use are provided elsewhere. Be sure to exit out of the GUI after marking each trial to allow the next trial to be marked. Finally, you will be prompted to save the data with a GUI dialog box. Navigate to the folder you wish to save the data to and, if desired, change the file name (although the default datafile name is strongly recommended), then click "save". If the "cancel" button is pressed, the file will not be saved; if this is in error, the data must be saved manually by clicking the "Save Workspace" button in Matlab.


*************************************************************************
*************************************************************************

Please do not modify any functions, particularly those that have not been described above. Some of the functions (listed or unlisted) are actively called by the functions named above, while others remain unused. It is best to just ignore but not delete any unnamed functions.

