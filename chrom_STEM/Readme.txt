Steps of ChromSTEM Image analysis:

1. Run Domain_Centers.ijm on ImageJ and save domain center list in csv.

2. Run DomainAnalysis.m with the image and the center list.

3. Run DomainAnalysisResults.m in the output folder of the domains*.mat

4. Run DomainGroups.m to label domains into groups of nascent, small mature, decaying, or large mature domains

5. Run PlotandCompareCVC.m to compare the CVC trends for each subgroups.

---
Others:
MakeProjectionFromTiff.m: To create a 2D projection figure from the 3D volume