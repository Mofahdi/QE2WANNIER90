# QE2WANNIER90
## Functionality
This script is designed to give good initial guesses for the trial localized orbitals projected onto the Bloch functions to obtain the Wannier functions in real space
as proposed by N. Marzari and D. Vanderbilt work in the Ref. below<br />
N. Marzari and D. Vanderbilt, Phys. Rev. B. 56, 12847 (1997)
***

## Motivation
creating an automated method to run WANNIER90 to obtain Maximally-localized Wannier functions (MLWF)s with less human intervention is challenging. Building a tool that can generate WANNIER90 input files has become more imperative due to the emergence of high-throughput (HT) workflows. Recently, the SCDM methodolgy had been proposed to be used to automatically generate MLWFs for HT frameworks in article entitled "Automated high-throughput Wannierisation" that can be accessed in this [link](https://www.nature.com/articles/s41524-020-0312-y). However, the SCDM methodolgy can only be implemented in Quantum Espresso (QE). This script can help generate WANNIER90 input files from VASP runs with a reasonable guess of the trial projections that can help the user perform HT calculations. This script enables the user to control the projections on each specie
***

## Features and Capabilities:
1. the code checks whether the POTCAR and POSCAR files correspond to each other to make sure your computations are consistent

2. The code provides the trial projected orbitals onto each specie and site along with the associated number of Wannier functions. 
   * if you select (projection_type="orbitals"), the code will give projected orbitals based on the species in s, p, d, and f format 
   * if you select (projection_type="quantum_numbers"), the code will give projected orbitals based on the species in l=#, ml=-l..l quantum numbers format 

3. The code can output the number of Wannier functions for spin-orbit coupling (SOC) calculations or non-SOC calculations
   * you have to set the argument "SOC" as True or False based on your calculations. The argument "SOC" is in "write_wan_projections" function

4. The code can also output the high-symmetry path for your structure if you include "bands_plot" and set it equal to "true" in "other_commands" argument in "write_wan_projections" function.
Please look at the example available in 'wan_projections_vasp.py' at the bottom of the script.

5. Pay attention that I only included the word 'true' in the script, but you can actually type '.true.' or 'True' in WANNIER90 input file. Please read WANNIER90 documentation for more info

6. if you run the code as it is, you will create a file named "wannier90_qe.win" which is uploaded here. However, you can always change the name of the file in the argument "file_name" inside the function "write_wan_projections". If you don't specify a name for your input file, the file outputted is named "wannier90.win" by default.
***

## Required Packages
* The script assumes that you had previously installed the following python packages <br />
<code> jarvis-tools==2022.9.16</code><br />
<code> ase==3.22.1</code><br />

* I believe installing the exact version is not necessary, so you can install other versions that are different from the above and still be able to run the script!

<dl>
<dt><code>jarvis-tools</code></dt>
<dd>The method I implemented to install jarvis-tools is through running the following in the command line:<br />
<code> pip install -U jarvis-tools </code>

for more information on how to install jarvis-tools with the source code, please use this [link](https://github.com/usnistgov/jarvis)</dd>
</dl>

## Usage and Development
* to copy this repository to your own computer please run the following in the command line: <br />
<code>git clone https://github.com/Mofahdi/QE2WANNIER90 </code>

* if you have any questions or would like to see more functionalities in this script, please do not hesistate to let me know!
* also please consider reading my published work in Google Scholar using this [link](https://scholar.google.com/citations?user=5tkWy4AAAAAJ&hl=en&oi=ao) thank you :)
