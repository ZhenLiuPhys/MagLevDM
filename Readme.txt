This repository contains all the code used to generate the plots and figures in arXiv:xxxxx.xxxx.
If you have any questions, please contact the author: Saarik Kalia (kalias@umn.edu)

Each figure in the paper is produced by the following code.
- Fig. 1: Figures/scpfig.py and Figures/scpfig2.py (Note that scpfig.py creates the left half of the figure, while scpfig2.py creates the right half.  These must be combined manually.  Also note that scpfig2.py will not display the SCP if the back wall is present.  An image with the SCP and an image with the wall can each be created separately and then the SCP can be cropped onto the image with the wall.)
- Fig. 2: Plots/noise_sources.py
- Fig. 3: Plots/axion_sensitivity_alt.py and Plots/DPDM_sensitivity_alt.py (The data files for existing constraints are also present in the Plots folder.)
- Fig. B-1: Figures/loop.py
- Fig. B-2: Figures/shield.py
- Fig. B-3: Figures/omegas.py
- Fig. B-4: Figures/v1int.py

Additionally, axion_analytic.py computes the magnetic field signal of axion DM at the center of the quadrupole trap, as outlined in Appendix B.  Along with various physical parameters of the cavity/levitation apparatus (R, h, I, Lx, Ly, Lz, x0, y0, z0; these are all defined in Appendix B), the code also has an input parameter lnum, which determines the number of cavity modes over which the code will sum.