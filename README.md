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

See the CONTENTS file for a list of what is in each file.
See the file howtoslick-with-figs.pdf for a tutorial
