

HYPERSONIC AERODYNAMICS                                          \hyper\readme.txt

The files for this program are in the directory \hyper on the CD-ROM  and in the
archive file hyper.zip that may be downloaded from the PDAS web site.

  readme.txt      general description
  hyper.f90       the source code in modern Fortran
  methods.f90     a major module needed by hyper.f90
  input.txt       instructions for preparing input & interpreting output

Sample cases for this program are:
  tmx1242.inp     input for wing-body configuration of NASA TM X-1242
  tmx1242.out     output from running hyper on tmx1242.inp
  tmx1242.dbg     log file from running tmx1242.inp
  tmx1242.mak     input file for program makewgs to make the tmx1242.wgs input file
  tmx1242.wgs     the geometry file describing the vehicle in TM X-1242

  tnd6480.*       5 files, similar to above, but for NASA TN D-6480


The reference documents for this program may be accessed
from the web page http://www.pdas.com/hyperrefs.html. 

To compile this program for your computer, use the command
   gfortran  hyper.f90 -o hyper.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.
The module methods.f90 will be included automatically.


This program is a completely new program based on the procedures in the famous
program called "The Supersonic/Hypersonic Arbitrary-Body Program" written by
Arvel Gentry, Doug Smyth and Wayne Oliver of Douglas Aircraft for the Wright
Laboratories of the USAF. The program is still "in progress" and not quite ready 
for general release.

The essential idea behind this program is that hypersonic flow is characterized
by a lack of interference between components of a vehicle and that the local
surface pressure coefficient is determined solely by the incidence of the local
surface to the freestream and by the Mach number.

In addition to the files mentioned above, there is a directory called \hyper\igor
with a number of programs written by Igor Polykov. Igor worked in the USSR 
airplane design industry for many years and these programs are front-end and 
post-processing programs for the Gentry program. I have not had time to check
these programs out, but they should prove useful and informative.



