[![View RegGrow on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/72944-reggrow)

This function implements a region growing algorithm for 2D, 3D and ND.  
In short this means that this function looks at all the pixels/voxels surrounding the segmented region and adds those that have a value that are within maxDiff of the mean of the already selected region. The region starts out with the seed position, and the loop continues until no more voxels are added. This function tests the entire shell at once.

This is slower than a mex implementation would be, but this should be compatible with any release and will return the same result on any release (including GNU Octave).

An HTML doc with examples is included. A tester function is also included.

Licence: CC by-nc-sa 4.0

![Example](https://github.com/thrynae/RegGrow/raw/master/html/RegGrow_show_default.png)
