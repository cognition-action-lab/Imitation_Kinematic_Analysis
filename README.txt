*************************************************************************

This code processes 3D Kinematics data for the Apraxia Battery of the
Cognition and Action Laboratory (Director, Laurel J Buxbaum). This code was
written by Aaron L Wong (03/05/2019), with certain assumptions made as to 
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

To change these assumptions, see ABload.m. To change assumptions about the
task structure, see decodetrials.m.

*************************************************************************

When there are a large number of data folders to be processed, call BATCHCOMPILE. This function allows you to select (either via hard-coding or by using a selection GUI) a set of paths to data folders containing trial-level raw data to be analyzed. When all desired folders exist within a single parent folder, leaving the "paths" variable empty will result in GUI that allows you to select multiple data folders (within a single parent folder) to be processed. When data span multiple parent folders (e.g., "controls" and "patients"), either separate calls to BATCHCOMPILE must be called or the full set of paths can be hard-coded. Regardless of how the folders are chosen, the script will automatically run ABBLOCKCOMPILE on all indicated folders followed by ABSUBJECTCOMPILE to put all the processed data into a single data file for further analysis. For details on how these sub-functions work or to call them separately, see below.

>BATCHCOMPILE will ask you for the folders of the individual data sheets you want to combine into one subject/block file. It will then call ABSUBJECTCOMPILE, which will automatically know what files you are compiling and ask which “block” file (with multiple subjects), if any, you want to append this data to. If making a new file, hit cancel, otherwise select the file to be appended. 
Following this step, data across individuals may be analyzed using ABFINDCOREACTION or ABFINDFULLACTION (see below).

~~~~~

To analyze a single dataset, first call ABBLOCKCOMPILE. This function allows you to select a single folder of data, for which each trial will be loaded and compiled based on its filename. It is expected that trial names will include the following information: block type, subject number, and trial number (i.e., "meaningful_unnamed_sub1_trial1.xls"). It is also assumed that the data are sorted by Group (e.g., one of "Control", "Patient", or "Model") and that all the data in the selected folder come from a single block and a single subject. If the data files do not meet these assumptions, the data may not be loaded properly. The output of this function is a data (.mat) file containing a structure with data sorted by block and trial number.

NOTE: As part of the ABBLOCKCOMPILE analysis, the user may be prompted to identify sign discontinuities in the dataset. These arose because the wrong hemisphere was specified in the trakSTAR system, and so when the tracker crossed the axis of the transmitter the sign inappropriately changed along all axes. Discontinuities during the trial should be automatically detected, but user input is required to confirm that these points have been properly identified (particularly if the start or the end of the trial should have been marked, e.g., if the marker started off on the wrong side of the transmitter axis). Markers will need to be verified for each tracker for which at least one discontinuity has been detected.

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
***This pairing is current as of 02/14/2019. If the item order is ever changed, this is the file that must be modified.***

~~~~~

Once all of the data have been initially processed, they should be compiled. Run the ABSUBJECTCOMPILE function to compile all blocks of the same type across subjects. Note that all folders (model, control and patient) are to be selected using the GUI. This function takes the individual data from each subject, sorts it into the appropriate group, and compiles it UN-MODIFIED into a single data structure, which is then saved. This becomes the core data file for all future analysis.

NOTE: Replacing subjects in an existing compiled dataset will also work using the ABSUBJECTCOMPILE function. New subjects will be appended to the existing dataset without overwriting any data that has already been processed, by selecting only the folders containing block-processed data for the subjects to be added. If a selected data folder contains data by a subject already existing in the compiled data file for that block, you will be prompted if you want to overwrite completely (o), replace but preserve existing data markings (r), create a new/duplicate entry (n), or skip (s) that subject. For example, you may need to use the replace option if any preprocessing is modified at the block-level analysis. Because only the this allows you to update only the kinematic tracker data (without changing any of the data from subsequent steps, e.g., marking the "core action"). Note that under certain conditions this may be inappropriate and re-marking the data may be required, such as if you replace a subject with a dataset recorded from a completely different session rather than the same dataset but preprocessed differently; in that case you need to remember to update the "core" marks for that person as well (see ABFINDCOREACTION function below).

~~~~~

If there are multiple files that have been created with ABSUBJECTCOMPILE and you wish to merge them, use ABMERGEPROCESSEDDATA. Select the original processed data file, and the secondary data file to be merged in, and follow any prompts that appear (e.g., when subjects appear in both data files). 


*************************************************************************

~~~FOR ANALYSIS WITH MARKING OF "FULL" ACTION~~~

After all new subjects have been compiled into the existing data file, run the ABFINDFULLACTION function. When this script is run, a GUI will prompt the selection of the data file to be analyzed. A dialog box then pops up to prompt the user to select which blocks, items, controls, and patients to run. The user can input -1 to select all existing, or input a list (comma-delimited) to choose a subset. For blocks, the user can enter block indices (numbers only) or block names (text). For items, the user can enter item numbers (numbers only) or item names (text). For subjects, the user can input a set of indices into the data structure (e.g., the 1st, 2nd, and 5th subjects), or can identify them by subject ID number using an "S" designation (e.g., S001, S280, etc). If no subjects are to be analyzed, enter 0 to skip this group. If you want the program to analyze all existing subjects in the data file, leave the default input or enter -1. Note, in all cases consistency is required; it is NOT possible to enter a mixture of indexes and names. Also, set whether you want the script to automatically clear all existing markings in the data set or not by setting the "Clear Existing Marks" option to 1 (clear all marks and regenerate them) or 0 (keep existing marks to be displayed).

Once these options are select, the script will proceed to display those blocks/items/subjects for marking. To mark the data, the each trial is plotted individually and the user must manually identify the start and end of the "full" action -- that is, from movement start (lift-off) to movement end (return to the table). Suggested marks are provided to the user based on times when the velocity drops below a threshold. It is the responsibility of the user to verify these marks, as they will not necessarily be appropriate in all cases. Data marking is performed by calling the MARKDATAGUI data-marking function; specific instructions for use are provided elsewhere. Be sure to exit out of the GUI after marking each trial to allow the next trial to be marked. Finally, you will be prompted to resave the data back to the same file at the end of each call to ABFINDFULLACTION. If you desire to save the data into a different file, select not to save the data (hit "0" at the prompt) and then save the workspace manually.



~~~FOR ANALYSIS WITH MARKING OF "CORE" ACTION AND SUBMOVEMENTS~~~

After all new subjects have been compiled into the existing data file, run the ABFINDCOREACTION function. When this script is run, a GUI will prompt the selection of the data file to be analyzed. A dialog box then pops up to prompt the user to select which blocks, items, controls, and patients to run. The user can input -1 to select all existing, or input a list (comma-delimited) to choose a subset. For blocks, the user can enter block indices (numbers only) or block names (text). For items, the user can enter item numbers (numbers only) or item names (text). For subjects, the user can input a set of indices into the data structure (e.g., the 1st, 2nd, and 5th subjects), or can identify them by subject ID number using an "S" designation (e.g., S001, S280, etc). If no subjects are to be analyzed, enter 0 to skip this group. If you want the program to analyze all existing subjects in the data file, leave the default input or enter -1. Note, in all cases consistency is required; it is NOT possible to enter a mixture of indexes and names. 

Once these options are select, the script will proceed to display those blocks/items/subjects for marking. To mark the data, the each trial is plotted individually and the user must manually identify the start and end of the "core" part of each action -- that is, the part that is really task-relevant -- as well as any sub-actions. Instructions for how to identify the Core Action for each item are described elsewhere. If no marks exist for the current individual/gesture, suggested marks are provided to the user based on times when the velocity drops below a threshold. It is critical that the user verify these marks, as they will not be appropriate for all items. Data marking is performed by calling the MARKDATAGUI data-marking function; specific instructions for use are provided elsewhere. Be sure to exit out of the GUI after marking each trial to allow the next trial to be marked. Finally, you will be prompted to resave the data back to the same file at the end of each call to ABFINDCOREACTION. If you desire to save the data into a different file, select not to save the data (hit "0" at the prompt) and then save the workspace manually.


*************************************************************************
*************************************************************************

Please do not modify any functions, particularly those that have not been described above. Some of the functions (listed or unlisted) are actively called by the functions named above, while others remain unused. It is best to just ignore but not delete any unnamed functions.

