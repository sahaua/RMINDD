# RMINDD

RMINDD - Radiation-Matter Interaction and Damage calculation using Nuclear Data
RMINDD is a package of computer programs written in Python programming language. It can read the recent basic evaluated nuclear data libraries such as ENDF/B-VII.1, ENDF/B-VIII.0, TENDL-2023, etc. to estimate the quantities that determine the primary radiation damage in structural materials such as Fe, Ni, Cr, Si, W, etc. used in nuclear applications. The basic features of this program so far developed are distributed here for academic and scientific use. Future developments are possible and will be available here once developed and tested. This README file also be updated as when required.

Please cite the following article when you use this program:
U. Saha, "Primary radiation damage due to neutron interactions using inexplicit evaluated nuclear data: a case study in isotopes of tungsten using ENDF/B-VIII.0 and TENDL-2019", Pramana – J. Phys. (2024) 98:5, https://doi.org/10.1007/s12043-023-02682-2

The input cards / keywords / commands must be given through an input file. The keywords (available at present) are listed in the distributed file "Input_RMINDD.txt". Any of the modules can be repeatedly run for same or different ENDF/B file for same or different isotope as per requirement as long as it is numbered properly in creasing order in the sequence in which it is run. For example, if run module EngdepU first with Module_number = 1 and again run it after two modules (say, RecedU and TransmU), then before starting with the module specific keywords for EngdepU number it as Module_number = 4, and so on. 
