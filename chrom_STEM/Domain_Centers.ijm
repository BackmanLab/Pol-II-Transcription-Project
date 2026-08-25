Gaussian_radius = 5;
CLAHE_blocksize = 120;
Prominence = 1.5; //Criteria for Domain Center

run("Z Project...", "projection=[Average Intensity]");
run("Gaussian Blur...", "sigma=Gaussian_radius");
run("Enhance Local Contrast (CLAHE)", "blocksize=CLAHE_blocksize histogram=256 maximum=3 mask=*None*");

getStatistics(area, mean, min, max, std, histogram);

run("Find Maxima...", "prominence="+(Prominence*std)+" strict exclude output=[Point Selection]");
//run("Find Maxima...", "prominence=0.3 strict exclude output=[Point Selection]");
run("Measure");