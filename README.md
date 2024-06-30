# OT_LS_FWI
Ground Penetrating Radar crosshole full-waveform inversion algorithm using the optimal-transport and least-squares distances. 

General Information: 
--------------------

This packages contains the source code of the OT-LS_FWI algorithm  to invert crosshole Ground Penetrating Radar data using the full waveform inversion algorithm based on optimal-transport and least-squares distances. It is published in the following paper:
J. Hunziker, G. Meles, and N. Linde (2024), Crosshole ground-penetrating radar
full-waveform inversion by combining optimal-transport and least-squares distances,
Journal of Applied Geophysics, in submitted.
Source: https://github.com/jhunziker/OT_LS_FWI

This package uses the optimal-transport algorithm of Malcolm Sambridge
(Australian National University).
Paper:
M. Sambridge, A. Jackson, and A. P. Valentine (2022), Geophysical inversion and optimal
transport, Geophysical Journal International, 231(1), 172–198.
doi: https://doi.org/10.1093/gji/ggac151
Source code:
https://github.com/msambridge/waveform-ot

This package uses also the 2D GPR crosshole finite difference time domain (FDTD)
solver of James Irving (University of Lausanne).
Paper:
J. Irving and R. Knight (2006), Numerical modeling of ground-penetrating radar
in 2-D using MATLAB, Computers & Geosciences, 32(9), 1247–1258.
doi: https://doi.org/10.1016/j.cageo.2005.11.006.

How to run the code: 
--------------------

Run the file run_OT_LS_FWI.m in Matlab. Note, it will take several hours to complete depending on your machine. 

How to plot the results: 
------------------------

To plot the results run the file plot_results_multisrc_multirec.m in Matlab.

License: 
--------

This file is part of the OT-LS-FWI package. OT-LS-FWI is free software: 
you can redistribute it and/or modify it under the terms of the GNU General 
Public License as published by the Free Software Foundation, either 
version 3 of the License, or (at your option) any later version. 
OT-LS-FWI is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
more details. You should have received a copy of the GNU General Public 
License along with OT-LS-FWI. If not, see <https://www.gnu.org/licenses/>. 
