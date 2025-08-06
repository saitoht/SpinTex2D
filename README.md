# SpinTex2D
Plot 2D spin texture from QE data

## Notice
This program uses a linear interpolation method to get smoother 2D Fermi surface and spin texture. 

## Usage
First, you will generate 2D k-point mesh for nscf calculation by setting plt_2Dfermi=plt_arrow=plt_spol=False
```shell-session
python 2D-spintex.py
```
The 2D k-point mesh will be written in Kpoints-2D.txt.

After conducting nscf with the 2D k-point mesh, you will activate plt_* switches to see the results.
In this stage, you should also set the Fermi energy, file name, spin component, etc.

## Input parameters
- nbnds: number of bands
- nkx,nky: number of kx/ky division
- kcx,kcy: center position of kx/ky in the unit of reciprocal lattice vector
- kxmax,kymax: maximum of kx/ky from center in the unit of reciprocal lattice vector, the k-range will be (kxc-kxmax,kxc+kxmax), etc.
- kcxp,kcyp: kcx/kcy for plotting
- kxmp,kymp: kxmax/kymax for plotting
- ef: the Fermi energy in scf in eV
- ef_shift: shift of the Fermi energy in eV
- fn: output file name
- spin_direction: component of spin expectation, "x", "y", "z"
- plt_2Dfermi: This should be True when you plot the spin texture
- plt_arrow: plot arrows (sx,sy) in the case of spin_direction="z"
- plt_spol: plot spin polarization values with color bar or not
- efshmin,efshmax: min/max of ef_shift if sw_search=True
- ediv: energy division in eV
- sw_search: search appropriate energy level for plotting or not
