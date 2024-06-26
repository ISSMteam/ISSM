def plotdoc():
    '''PLOTDOC - plot documentation
    %As of now it is more a TODO list
    %   Usage:
    %      plotdoc()
    '''
    pydata = {'quiver': ' quiver plot give data and a vector array [Vx, Vy]',
              'mesh': ' draw mesh using trisurf',
              'BC': ' this will draw all the boundary conditions (Dirichlet and Neumann).',
              'elementnumbering': ' numbering of elements (matlab indices)',
              '3D disclaimer': '3D is implemented with plot3d for now this is not optimal and may change to mayavi at some point. The impelementation is on the development side for now so expect some issue and question your plotting before you results.'}
    # TODOdata = {'basal_drag': ' plot the basal drag on the bed (in kPa) based on the velocity in md.initialization',
    #             'basal_dragx or basal_dragy': ' plot a component of the basal drag on the bed (in kPa)',
    #             'boundaries': ' this will draw all the segment boundaries to the model, including rifts.',
    #             'icefront': ' this will show segments that are used to define the icefront of the model (Neumann boundary conditions).',
    #             'deviatoricstress_tensor': ' plot the components of the deviatoric stress tensor (tauxx, tauyy, tauzz, tauxy, tauxz, tauyz, if computed',
    #             'deviatoricstress_principal': ' plot the deviatoricstress tensor principal axis and principal values',
    #             'deviatoricstress_principalaxis1': ' arrow plot the first principal axis of the deviatoricstress tensor(replace 1 by 2 or 3 if needed)',
    #             'driving_stress': ' plot the driving stress (in kPa)',
    #             'elements_type': ' model used for each element',
    #             'highlightvertices': ' to highlight vertices (use highlight option to enter the vertex list',
    #             'referential': ' stressbalance referential',
    #             'riftvel': ' velocities along rifts',
    #             'riftrelvel': ' relative velocities along rifts',
    #             'riftpenetration': ' penetration levels for a fault',
    #             'riftfraction': ' fill fractions for every node of the rifts',
    #             'rifts': ' plot mesh with an offset so that rifts are visible',
    #             'strainrate_tensor': ' plot the components of the strainrate tensor (exx, eyy, ezz, exy, exz, eyz) if computed',
    #             'strainrate_principal': ' plot the strainrate tensor principal axis and principal values)',
    #             'strainrate_principalaxis1': ' arrow plot the first principal axis of the strainrate tensor(replace 1 by 2 or 3 if needed)',
    #             'stress_tensor': ' plot the components of stress tensor (sxx, syy, szz, sxy, sxz, syz) if computed',
    #             'stress_principal': ' plot the stress tensor principal axis and principal values',
    #             'stress_principalaxis1': ' arrow plot the first principal axis of the stress tensor(replace 1 by 2 or 3 if needed)',
    #             'transient_results': ' this will printlay all the time steps of a transient run (use steps to specify the steps requested)',
    #             'transient_vel': ' vel can be by any field of the transient results (vx, vy, vz, vel, temperature, melting, pressure, bed, thickness, surface)',
    #             'transient_field': ' dynamic plot of results. specify ''steps'' option, as fell as ''field'' (defaults are all steps, for ''Vel'' field)',
    #             'transient_movie': ' this will printlay the time steps of a given field of a transient run',
    #             'transient_movie_field': ' field to be printlayed when doing  transient_movie data printlay',
    #             'transient_movie_output': ' filename if output is desired for movie',
    #             'transient_movie_time': ' time for each image (default 2 seconds)',
    #             'thermaltransient_results': ' this will printlay all the time steps of a thermal transient run',
    #             'qmuhistnorm': ' histogram normal distribution. needs option qmudata',
    #             'qmumean': ' plot of mean distribution in sampling analysis with scaled response. needs option qmudata for descriptor',
    #             'qmustddev': ' plot of stddev distribution in sampling analysis with scaled response. needs option qmudata for descriptor',
    #             'part_hist': ' partitioning node and area histogram'}

    pyoptions = {'axis': " show ('on') or hide ('off') axes",
                 'caxis': " modify  colorbar range. (array of type [a b] where b >= a)",
                 'colorlevels': " N, number of levels to use",
                 'colorbar': " add colorbar (string 'on', 'off' or 'one')",
                 'axes_pad': " spacing between axes (default is 0.25)",
                 'colorbartitle': " colorbar title (string)",
                 'colorbarticks': " set colorbar ticks manually (list)",
                 'colorbarfontsize': " specify colorbar fontsize",
                 'colormap': " change the default colormap ('viridis' is the default)",
                 'contourlevels': " N or [value1, ...] add the contours of the specified values or N contours",
                 'streamlines': " TOFIX argument does nothing",
                 'edgecolor': " color of mesh edges. RGB tuple or standard string",
                 'fontsize': " fontsize for the title",
                 'fontweight': " fontweight for the title 'normal', 'bold'",
                 'fontcolor': " TODO",
                 'highlight': " highlights certain nodes or elements when using 'vetrexnumbering' or 'elementnumbering' or 'highlightvertices ' or 'highlightelements' option",
                 'title': " subplot title (string)",
                 'xlim': " limits of X axis (all subplots) (ex:  [0, 500])",
                 'ylim': " limits of Y axis (all subplots) (ex:  [0, 500])",
                 'xlabel': " X axis title",
                 'ylabel': " Y axis title",
                 'scaling': " scaling factor used by quiver plots.",
                 'quivercol': " color of quiver arrows, 'values' give value colored arrows",
                 'text': " print string or list of strings",
                 'textposition': " [x, y] position of text, list if several texts (position betwee 0 and 1)",
                 'textsize': " text fontsize TOFIX ",
                 'textweight': " text fontweight",
                 'textcolor': " text color",
                 'textrotation': " text rotation angle",
                 'mask': " condition. Only 'true' values are plotted ",
                 'log': " cutoff value for log",
                 'backgroundcolor': " plot background color. RGB tuple or standard string",
                 'expdisp': " path (or list of paths) to the exp file to be plotted ",
                 'explinewidth': " linewidth ",
                 'explinestyle': " matplotlib linestyle string ",
                 'explinecolor': " matplotlib color string ",
                 'expfill': " (True / False) fill a closed contour ",
                 'expfillcolor': " Color for a filled contour, only used if expfill is True ",
                 'expfillalpha': " alpha transparency for filled contour ",
                 'overlay': " True / False. Overlay a georeferenced image (radar / visible) ",
                 'overlay_image': " path to overlay image ",
                 'overlayhist': " plot a histogram of overlay image, used for setting overlaylims ",
                 'overlaylims': " normalized limits to clip and stretch contrast of overlay image (in [0, 1], ex. [0.25, 0.75]) ",
                 'alpha': " set transparency of plotted data (in [0, 1]) ",
                 'vertexnumbering': ' numbering of vertices',
                 'elementnumbering': ' numbering of elements (matlab indices)',
                 'highlightelements': ' to highlight elements to highlight the element list',
                 'layer': "number of the layer to display for 3D runs"}

    # TODOoptions = {'basin': " zoom on a given basin ('pineislandglacier', 'ronneiceshelf', use isbasin to identify a basin",
    #                'figurebackgroundcolor': " figure background color. (default is 'none', ",
    #                'coord': "  'xy' (default) or 'latlon'",
    #                'colorbarpos': " [x, y, dx, dy] where x, y, dx and dy are within [0 1]",
    #                'colorbarcornerposition': " 'West', 'North', etc ...",
    #                'colorbartitlerotation': " - 90, etc ...",
    #                'colorbarwidth': " multiplier (default 1) to the default width colorbar",
    #                'colorbarheight': " multiplier (default 1) to the default height colorbar",
    #                'contourticks': " 'on' or 'off' to printlay the ticks of the contours",
    #                'contouronly': " 'on' or 'off' to printlay the contours on a white background",
    #                'contourcolor': " ticks and contour color",
    #                'density': " density of quivers (one arrow every N nodes, N integer)",
    #                'inset': " add an inset (zoom) of the current figure if 1 (use 'insetx', 'insety' and 'insetpos' to determine the inset position and content)",
    #                'insetx': " [min(x) max(x)] where min(x) and max(x) are values determining the inset content",
    #                'insety': " [min(y) max(y)] where min(y) and max(y) are values determining the inset content",
    #                'insetpos': " [x, y, dx, dy] where x, y, dx and dy are within [0 1]",
    #                'resolution': " resolution used by section value (array of type [horizontal_resolution vertical_resolution])",
    #                'showsection': " show section used by 'sectionvalue' (string 'on' or a number of labels)",
    #                'sectionvalue': " give the value of data on a profile given by an Argus file (string 'Argusfile_name.exp', ",
    #                'profile': " give the value of data along a vertical profile ([xlocation ylocation])",
    #                'smooth': " smooth element data (string 'yes' or integer)",
    #                'view': " same as standard matlab option (ex:  2, 3 or [90 180]",
    #                'zlim': " same as standard matlab option",
    #                'xticklabel': " specifiy xticklabel",
    #                'yticklabel': " specifiy yticklabel",
    #                'contrast': " (default 1) coefficient to add contrast to the radar amplitude image used in overlays",
    #                'highres': " resolution of overlayed radar amplitude image (default is 0, high resolution is 1).",
    #                'alpha': " transparency coefficient (the higher, the more transparent). Default is 1.5",
    #                'scaling': " scaling factor used by quiver plots. Default is 0.4",
    #                'autoscale': " set to 'off' to have all the quivers with the same size. Default is 'on'",
    #                'linewidth': " line width for expprint plot (use a cell of strings if more than one)",
    #                'border': " size of printlay border (in pixels). active only for overlay plots",
    #                'nan': " value assigned to NaNs (convenient when plotting BC)",
    #                'partition: " a partion vector. generates overlay plot of partition edges",
    #                'latlon': " 'on' or {latstep lonstep [resolution [color]]} where latstep, longstep and resolution are in degrees, color is a [r g b] array",
    #                'latlonnumbering': " 'on' or {latgap longap colornumber latangle lonangle} where latgap and longap are pixel gaps for the numbers",
    #                'latlonclick': " 'on' to click on latlon ticks positions colornumber is a [r g b] array and latangle and lonangle are angles to flip the numbers",
    #                'northarrow': " add an arrow pointing north, 'on' for default value or [x0 y0 length [ratio width fontsize]] where (x0, y0) are the coordinates of the base, ratio = headlength / length",
    #                'offset': " mesh offset used by 'rifts', default is 500",
    #                'scaleruler': " add a scale ruler, 'on' for default value or [x0 y0 length width numberofticks] where (x0, y0) are the coordinates of the lower left corner",
    #                'showregion': " show domain in Antarctica on an inset, use 'insetpos' properties",
    #                'visible': " 'off' to make figure unvisible, default is 'on'",
    #                'wrapping': " repeat 'n' times the colormap ('n' must be an integer)",
    #                'unit': " by default, in m, otherwise, 'km' is available",
    #                'legend_position': " by default, 'NorthEasth'",
    #                'qmudata': " ",
    #                'figposition': " position of figure:  'fullscreen', 'halfright', 'halfleft', 'portrait', 'landscape', ... (hardcoded in applyoptions.m)",
    #                'offsetaxispos': " offset of current axis position to get more space (ex:  [ -0.02 0  0.04 0])",
    #                'axispos': " axis position to get more space",
    #                'hmin': " (numeric, minimum for histogram)",
    #                'hmax': " (numeric, maximum for histogram)",
    #                'hnint': " (numeric, number of intervals for histogram)",
    #                'ymin1': " (numeric, minimum of histogram y - axis)",
    #                'ymax1': " (numeric, maximum of histogram y - axis)",
    #                'ymin2': " (numeric, minimum of cdf y - axis)",
    #                'ymax2': " (numeric, maximum of cdf y - axis)",
    #                'cdfplt': " (char, 'off' to turn off cdf line plots)",
    #                'cdfleg': " (char, 'off' to turn off cdf legends)",
    #                'segmentnumbering': " ('off' by default)",
    #                'kmlgroundoverlay': " ('off' by default)",
    #                'kmlfilename': " ('tempfile.kml' by default)",
    #                'kmlroot': " ('./' by default)",
    #                'kmlimagename': " ('tempimage' by default)",
    #                'kmlimagetype': " ('png' by default)",
    #                'kmlresolution': " (1 by default)",
    #                'kmlfolder': " ('Ground Overlay' by default)",
    #                'kmlfolderdescription': " ('' by default)",
    #                'kmlgroundoverlayname': " ('' by default)",
    #                'kmlgroundoverlaydescription': "N/A by default')"}

    print("   Plot usage:  plotmodel(model, varargin)")
    print("   plotting is done with couples of keywords values, the type and style of data to display is given by one (or several) of the followings")
    print("   Options:  ")
    print("     'data' :  and a model field or one of the following options.")
    for key in list(pydata.keys()):
        print((" - {} :  {}".format(key, pydata[key])))
    print("")
    print("   The general look of the plot is then given by the following keywords")
    for key in sorted(pyoptions):
        print((" - {} :  {}".format(key, pyoptions[key])))
    print("       any options (except 'data') can be followed by '  #i' where 'i' is the subplot number, or '  #all' if applied to all plots")
