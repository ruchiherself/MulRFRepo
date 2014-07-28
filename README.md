MulRFRepo
=========

Repository of MulRF tool

Compilation steps for four projects---InputStats, MulRFScorer, MulRFSupertree, and MulRFGUI---are as follows: 

i) Compile InputStats: 1) download InputStats folder, 2) open Makefile (available inside InputStats folder) and uncomment two statements under "For Windows", "For Mac", or "For Linux" (depending your computer), 3) run command "make" from InputStats folder.

ii) Compile MulRFScorer: 1) download MulRFScorer folder, 2) open Makefile (available inside MulRFScorer folder) and uncomment two statements under "For Windows", "For Mac", or "For Linux" (depending your computer), 3) run command "make" from MulRFScorer folder.

iii) Compile MulRFSupertree: 1) download MulRFSupertree folder, 2) open Makefile (available inside MulRFSupertree folder) and uncomment two statements under "For Windows", "For Mac", or "For Linux" (depending your computer), 3) run command "make" from MulRFSupertree folder.

iv) Compile MulRFGUI: 1) download MulRFGUI, 2) download & install NetBeans, 3) import MulRFGUI, 4) run 'Build' command under 'Run', 5) the 'MulRF.jar' file  along with 'lib' directory is available in the 'dist' folder


Preparing MulRF to run after compilation:

MulRF runs by storing 'MulRF.jar', 'lib', 'outputData', 'inputData', and 'executables' in a folder which can be called MulRF1.2. The three executables that were produced by i - iii needs to be stored in the 'executables' folder, follow these steps for that: 

i) For Linux: Rename the new executable InputStats as InputStatsLin, MulRFSupertree as MulRFSupertreeLin, and MulRFScorer as MulRFScorerLin.

ii) For Windows: Rename the new executable InputStats.exe as InputStatsWin.exe, MulRFSupertree.exe as MulRFSupertreeWin.exe, and MulRFScorer.exe as MulRFScorerWin.exe.

iii) For Mac: Rename the new executable InputStats as InputStatsMac, MulRFSupertree as MulRFSupertreeMac, and MulRFScorer as MulRFScorerMac.

MulRF is now ready to be executed by clicking MulRF.jar or running the command "java -jar MulRF.jar"

