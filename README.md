# A fast parameter sampler and optimiser based on differential evolution. 
Diver is a differential evolution sampler, designed primarily for mapping and optimising likelihood functions, but is easily applicable to any optimisation problem.

A full description of the package and the statistical framework behind it can be found in the GAMBIT ScannerBit paper:

G Martinez, J McKay, B Farmer, P Scott, E Roebber, A Putze & J Conrad 2017, EPJC 77 (2017) no.11, 761, arXiv:1705.07959
If you write a paper that uses Diver, please cite this paper.

Diver releases can be obtained as tarballs from Hepforge. The latest and greatest version, along with a full revision history, can always be found in the git repository. Compilation and usage instructions, as well as several example programs, can be found in the code release.

Authors: Elinore Roebber (e.roebber@bham.ac.uk) & Pat Scott (p.scott@imperial.ac.uk)
#### Source: https://diver.hepforge.org/

# What this code does:
I don't know Fortran. However, I wanted to modify it as much as possible to trace out the degenerate local minima/maxima of any multidimensional function. This code is one of the most naive possible ways I could think of to finish the task. If this code also helps you to accomplish your task for any research project, I would be very happy and glad that I could help.

# Steps to run it:
  1. Download and install Diver. ( Website: https://diver.hepforge.org/ )
  2. Then, replace the "src" folder of the installed Diver with the "src" folder given here.
  3. Then, replace the "example_c" folder with the "example_c" folder given here. (Don't worry too much about the "simpson_integration.h" file. It's a bonus for you, 😄😄)
  4. Go to "example_c" and open "example_c.c" with a text editor.
  5. Edit it as if you are using "Diver".
  6. After that come back to the previous Directory. Download "diver_looper.sh" from here to this directory.
  7. Open "diver_looper.sh" with a text editor. In line number 5, change the number to the desired times you want to run the diver. Then save it.
  8. Copy line number 12 without the "#". Then open the terminal from that directory and paste it.
  9. Then type "./diver_looper.sh" in the terminal and type enter. If everything is alright, it should run.
  10. The output will be saved as a text file with the name "output.txt". It is still not usable. I am adding another piece of code which will refine it.
# These files include a demo.
  If you find any error, please email me. I would be happy to fix it. If you have any feedback you want to send, you can send me as well. Any kind of feedback is welcome. If you have or know a much more efficient way to do so, please send me the link to the code or the code as well. It would be much help.
  My email: ehsanzaman.k@gmail.com
