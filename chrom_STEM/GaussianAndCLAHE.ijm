resetMinAndMax;
run("16-bit");
run("Reslice [/]...", "output=1.000 start=Top avoid");

Stack.getStatistics(voxelCount,mean,min,max,stdDev)

CLAHE_blocksize = 120;
Gaussian_r = 1;

for (i=1; i<=nSlices; i++) {
	setSlice(i);
	run("Enhance Local Contrast (CLAHE)", "blocksize = CLAHE_blocksize, histogram = max-min, maximum = 3.0");
}

run("Reslice [/]...", "output=1.000 start=Top avoid");

run("Duplicate...", "duplicate");

Stack.getStatistics(voxelCount,mean,min,max,stdDev)
run("Min...", "value="+(mean-1.25*stdDev)+" stack");
run("Max...", "value="+(mean+3*stdDev)+" stack");
Stack.getStatistics(voxelCount,mean,min,max,stdDev)
run("Subtract...","value="+min+" stack");
Stack.getStatistics(voxelCount,mean,min,max,stdDev)
run("32-bit");
run("Gaussian Blur 3D...", "x=1 y=1 z=1");
run("Divide...", "value="+max+" stack");
