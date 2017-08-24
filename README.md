# slick
Software to invert fault slip data for a uniform stress field

The references for this package are:
Michael, A.J. (1984), Determination of stress from slip data: faults and folds, JGR, v. 89, 11,517-11526.

This paper covers the linear stress inversion technique used in the programs slick and slfast.

Michael, A.J. (1987a), The use of focal mechanisms to determine stress: a control study, JGR, v. 92, 357-368.

This paper covers the bootstrap statistics and the problems associated with using focal mechanisms instead of geological field data.

Michael, A.J. (1987b), Stress rotation during the Coalinga aftershock sequence, JGR, v. 92, 7963-7979.

This paper covers the addition of the plane flipping ability to the bootstrap statistics in order to account for uncertainties in the picking of the fault plane from the two nodal planes.

You will also need to obtain the stereonet plotting program onnet from
http://quake.wr.usgs.gov/~michael/software/onnet
I still need to update onnet and upload it to github.  You can used slick without onnet by plotting some other way from the output files.


CONTENTS
HOWTO the documentation without appendices in troff with ms macros format
howtoslick-with-figs.pdf documentation with annotated figures and file printouts in pdf format
bootboth.c runs the bootstrap statistics for both plane analysis
bootgrid.c runs the bootstrap statistics for the grid search inversion
bootslickw.c runs the bootstrap statistics for the linear inversion,
	and flips planes randomly
bothplanes.c makes a data file with both planes
controls include file that controls grid for gridfix.c and gridstrap.c
dirplg.c computes trend and plunge from a vector
dixie sample data file
dixieg sample data file in grid search format
draw_far.c subroutine for plotboots
eigen.c 3 by 3 eigenvalue and vector solver (subroutine)
gridfix.c grid search program with 1-cos(beta) error criterion
gridstrap.c fast version with truncated output of gridfix
leasq.c least squares solver (subroutine)
makeboth makefile for bothplanes
makebtbo makefile for bootboth
makebtgr makefile for bootgrid
makebtslw makefile for bootslick
makeplata makefile for plata
makeplbtg makefile for plotbootg
makeplbts makefile for plotboots
makeplbtso makefile for plotbootso
makeslfast makefile for slfast
makeslick makefile for slick
makeswitch makefile for switcher
myrand.c random number generator and seeder (subroutines)
plata.c program to plot data
plotbootg.c evaluates and plots confidence limits for grid search
plotboots.c evaluates and plots confidence limits for linear inversion,
   version that only plots some of the points to showing outer limits
   but not require lots of time to plot
plotbootso.c evaluates and plots confidence limits for linear inversion,
   version that plots all points
slfast.c fast version of slick with truncated output
slfast_sub.c fast version of slick with truncated output, subroutine version
slick.c linear inversion program
sort.c sorts a array of numbers (subroutine)
stridip.c finds strike and dip of a plane from the normal vector (subroutine)
switcher.c finds a new data set using the auxiliary planes of old data set
switchsub.c finds the auxiliary plane  and rake for one focal mechanism (subroutine)

The following files are object and executables for Mac OS X 10.6 corresponding to the *.c files above
switcher
switcher.o
slick
slick.o
slfast
slfast.o
plotbootso
plotbootso.o
plotboots
plotboots.o
draw_far.o
plotbootg
plotbootg.o
sort.o
plata
plata.o
bootslickw
bootslickw.o
dirplg.o
eigen.o
leasq.o
slfast_sub.o
switchsub.o
bootgrid
bootgrid.o
bootboth
bootboth.o
myrand.o
bothplanes
