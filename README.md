# PlotP2 - v.0.5#
PlotP2 Plots a P2 Triangulated surface.

***Requires xfigure***
## Syntax ##
    H = PlotP2(tri, X, np, C)
	H = PlotP2(tri, X, np, C, fignum)

	tri  	triangulation, m-by-6
	X  		coordinates, n-by-3
	C       Color data, n-by-1
	nE      number of elements in each P2 element

If the total number of elements exceed 30 000, a warning message will be displayed. That warning message can be disabled by the command:

	warning('off','PLOTP2:manyelements')

## Properties ##
	H.light
	H.patch
	H.edge

	H.EdgeColor
   	H.FaceColor