# 3DSMARTER
Code to analyze data from 3D-SMARTER tracking of silver nanoparticles on live cells
Top-level functions
bayesProcessing
This function takes raw tracking data at performs recursive Bayesian regression to calculate the particle trajectory. One microsecond binned data can be found in "200304 TR003 raw tracking data 590_650s.mat." Load this file into MATLAB workspace before running the function. This function generates the data found in Fig. 4.
Function call:
[xOut,yOut,zOut]=bayesProcessing(ktPos,TAGPos,n,xStage,yStage,zStage)
Inputs: 
ktPos – xy laser coordinates
TAGPos – z laser coordinate
n - # of photons
xStage – X piezo readout
yStage – Y piezo readout
zStage – Z piezo readout
Outputs:
xOut, yOut, zOut – X, Y, Z trajectory data at 1 microsecond sampling

cyl_fit_70_120s
This script imports 3D coordinates (in this example: data from 70 to 120s of the entire trajectory) from MHz Bayesian reconstruction, downsample to 1 kHz, and fits the sparse positional data into a cylinder. To ensure correct data import, this script is executed with the assumption that the current directory of MATLAB is same as script directory. The user can also add necessary .mat files to MATLAB path and remove path dependency in the load syntax. Such path dependency also applies to get_2d_dens_on_cyl_surf, get_2d_dens_top_view, show_all_central_axis, and residence_time_analysis.

get_2d_dens_on_cyl_surf
This script loads MHz reconstructed Bayesian data and cylinder parameters obtained from the above template and converts 3D Cartesian coordinates to cylindrical coordinates. The median of r for all spots in a given segment is then used as the radius of the cylinder. All spots are then projected onto the surface. The spots are treated as 2D data scattered on an unfwrapped cylinder surface and converted to a density map with 10 nm grid. This script generates data necessary for Fig. 5 a-g and can serve as a framework for Fig. S14.

get_2d_dens_top_view
This script is similar to the previous code, which loads position data and cylinder parameters and converts 3D Cartesian coordinates to cylindrical coordinates. The spots were then projected onto the bottom (or an arbitrary cross section) of the surface and counted to create a top-down view density map, which shows the structure was largely hollow inside. This script generates data necessary for Fig. 5 h-j.

show_all_central_axis
This script imports cylinder parameters from 13 manually identified segments of 50-60 seconds and obtains the central axes of the cylinders using the fitting method in the previous script, and draws them in one figure. This script generates data necessary for Fig. S11b.

residence_time_analysis
This script imports data from manually labeled 'hot spot' data and performs a series of analyses to visualize features of such residence events and remove datasets that display increased frequency in a region but are unlikely to indicate actual structures (e.g. transient passing from multiple visits). This script generates data necessary for Fig. S12.

EvaluateBayesianPostProcessing 
Uses pseudo-random number generation to simulate Kalman tracking of a diffusive particle at the experimental bin time. First, particle positions and observed photons are defined at a finer time resolution to post-process data. Then, photon arrival information, laser positions, expected diffusion coefficient, and (for Bayesian post-processing only) expected count rates for both the particle and background are used to estimate the particle's position during each "sub-bin recursively." The actual particle positions are then used to estimate and compare the error in each post-processing algorithm to the error in stage positions observed during tracking, as seen in Supplementary Figure 7.

cylSurfProps
This function calculates the diffusion coefficients, force vectors, and surface potential for particles diffusing on a cylindrical surface. When called, the function will ask the user to find "unwrapped 2d traj 590-650s.mat." This function is used to generate the data in Fig. 6.
Function call:
[Dmap,VmapX,VmapY]=cylSurfProps
Outputs:
Dmap – gridded diffusion data
VmapX – force vector data along X
VmapY – force vector data along Y

Utilities
parsedCylFit
This function takes 3D Cartesian coordinates and returns fitting parameters. The function invokes fmincon, a built-in MATLAB function that performs nonlinear optimization with constraints. 

rcyl
This class is used to create a query object that takes two points as a central axis and an additional point as a reference point and performs conversions between the cartesian, cylindrical, and spherical coordinates.
cyl = rcyl(P,Q,R) creates a rcyl object; P, Q, and R are 3x1 arrays; P is the original point in the cylinder. PQ defines the central axis of the cylinder, and R is a reference point outside line PQ to help determine the angle.
cyl_pos = cyl.cart2cyl(cart_pos) converts cartesian coordinates to cylindrical coordinates
cart_pos = cyl.cyl2cart(cyl_pos) converts cylindrical coordinates to cartesian coordinates
cart_pos is an n by 3 matrix of Cartesian coordinates, all the conversions are vectorized and computationally efficient
cyl_pos is an n by 3 matrix of cylindrical coordinates, the 3 columns are r, theta and h
similarly, you can call rcyl.cart2sph, rcyl.sph2cart, rcyl.cyl2sph, rcyl.sph2cyl to convert from or to spherical coordinates 
how does the conversion work:
P,Q define the central axis of the cylinder from arbitrary points in cartesian space. In rcyl.cart2cyl, to obtain r, any point is projected to line PQ and the length of the projected path is calculated. To obtain h, a query point is first projected to PQ and the distance from projected point to origin P is calculated. Q only defines the normal direction. So there will be negative values for h. To determine theta (angle), a random point M will be projected to PQ and get projected point M'; The reference point R is also projected to PQ to get R'; Then the angle between MM' and RR' is calculated. Note that this only returns an angle between 0~pi. To assign a sign to the angle, the function of the plane defined by Ax+By+Cz+1=0 (or Ax+By+Cz=0, depending on whether [0,0,0] is on the plane) is calculated and the sign of Ax+By+Cz+1 (or Ax+By+Cz+0) will be assigned to the angle from arccos calculation. Angle will be in the range of -180 to 180. So it's in degrees. You can easily transfer between degree to radians. The spherical system is similar, which comes with (r, theta, phi). Here the range of theta is from 0 to 180 and phi is from 0 to 360.
As the angle increments from -180 to 180, the ring along with the normal direction defined by PQ forms a right-handed helix. All conversions are guaranteed to be reproducible. cyl.cyl2cart(cyl.cart2cyl([x, y, z])) will return the same Cartesian coordinates.

pt_proj_on_line
This function calculates the projected positions of input spots on a line. This is necessary for Cartesian to cylindrical conversions. 

ct_3
This function performs counting of spots (Cartesian coordinates) in 3D with given individual 1D grids. The function can also be used in 2D counting with decreased number of inputs.

get_2d_pdata_from_ct3
This function extracts 2D density information from 3D counting and removes faces in XZ planes or YZ planes.

ct_idx3
This function is called in ct_3 and is used to count unique pixels. It serves as an efficient way to handle sparsely filled voxels in 3D data.

draw_rec
This function draws a rectangle.

drawCyl
This function draws a cylinder.

ezgrid
This function evaluates input array and gives a grid that guarantees all spots in input array is between two consecutive spots of the returned grid.

traj3D
This function draws 3D trajectory with changing color to visualize the diffusion with time.

lf_obj
This function handles a series of X/Y data and performs linear fit using only a consecutive segment of the data.

movmean2
This function computes the median of every two consecutive elements in an array.

get_img_dims
This function analyzes vertex of plot data and calculates minimal canvas range.

plasma
This function returns plasma colormap in matplotlib.

Paired
This function returns paired colormap in matplotlib.

cab
This function closes all except user specified figure windows in MATLAB.
Karl (2022). Close all figures except those listed (https://www.mathworks.com/matlabcentral/fileexchange/24420-close-all-figures-except-those-listed), MATLAB Central File Exchange. Retrieved December 3, 2022.

arrow3
This function is used to create 3D arrow objects in MATLAB.
Tom Davis (2022). Arrow3 (https://www.mathworks.com/matlabcentral/fileexchange/14056-arrow3), MATLAB Central File Exchange. Retrieved December 3, 2022.

Plot
This function provides a wrapped interface of MATLAB figure properties.
K M Masum Habib (2022). PlotPub - Publication Quality Graphs in MATLAB (https://github.com/masumhabib/PlotPub), GitHub. Retrieved December 3, 2022.
