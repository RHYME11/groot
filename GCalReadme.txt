## Calibration Protocol

Measure A Known Source

Start with a source that has well-known gamma-ray peaks. For example, Co-60 has two strong peaks at approximately: 1173.2 keV and 1332.5 keV

Then open the raw uncalibrated plot --> Locate the peaks and their channel position

create a text file of the format :

measured_channel known_energy
=================
eg;

429.8    1173.2
488.6    1332.5

=================

then open groot and
''''groot
MakeCalibration("chanelenergy.dat", "calname.cal", 1)



1 = linear calibration   // energy = C0 + C1 * channel
2 = quadratic calibration   // energy = C0 + C1 * channel + C2 * channel^2

then to apply calibration 
''''zsh
./bin/groot -g calname.cal ~/path/to/spectrum.txt3
























