These MATLAB scripts were developed to help automate data analysis from ACL mechanical testing 

Scripts:

1. CalculatingCSA.m -- Uses .stl file from 3D scanning. Lets user crop out fixtures and keep only the ligament. Using the ligament mesh only, cross sectional area is calculated with a 2D best-fit ellipse.

2. FitCurve.m -- Uses raw Acumen DAQ force & displacement data. User must input CSA and ligament length. Script computes stress & strain then fits a standard linear solid model to stress-relaxation curve over the auto-detected (or manually selected) region.

3. StrainPlot_Automatic_Simplified.m -- Uses raw Acumen displacement data and DIC data to create a strain vs time plot, specifically for the stress-relaxation test. The user must input the frequency that they want the Acumen data resampled at. Using the data from each input it detects the jump and hold from the stress-relaxation test to differentiate it from preconditioning and ramp-test and plots only that.
