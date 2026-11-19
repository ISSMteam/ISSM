#!/usr/bin/env python3
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy

from mesh2d import mesh2d
from expread import expread
from processmesh import processmesh
from processdata import processdata
from plot_unit import plot_unit
from SectionValues import SectionValues
from applyoptions import applyoptions

def plot_section(md,data,options,fig,axgrid,gridindex):
    '''
    PLOT_SECTION - plot a given field on a section

    Usage:
        plot_section(md,data,options,nlines,ncols,i)

    See Also: PLOTMODEL
    '''
    DEBUG=True

    #How many subplots?
    if options.exist('showsection'):
        raise Exception('Error: ''showsection'' option does not work in Python. Not implemented yet.')
    	#Compute the indexes of the 2 plots (one for the sectionvalue and one for showsection
        # upperplots=np.floor((gridindex)/ncols)
        # if upperplots==0:
        #     leftplots=gridindex
        # else:
        #     leftplots=i-ncols*upperplots-1
        # index1=4*ncols*upperplots+2*leftplots+1
        # index2=index1+1
        # ncols=2*ncols
    else:
        index1=gridindex
        ax=axgrid[gridindex]

    #process model
    x_m, y_m, z_m, elements_m, is2d, isplanet=processmesh(md,[],options)

    #Get number of curves and generate random colors
    if np.shape(data) == (md.mesh.numberofvertices,):
        data = np.reshape(data,(-1,1))
    numcurves=np.shape(data)[1]
    if is2d:
        #NOTE: It seems that this colormap is only used for 2D dataset.
        colorm=options.getfieldvalue(options,'colormap',mpl.cm.tab20)

    #replug x and y onto model so that SectionValue treats the problem correctly
    md3d=copy.deepcopy(md)
    if options.exist('layer'):
        #NOTE: Python works differently copying object. We copy md using "copy.deepcopy".
        md=copy.deepcopy(md)
        md.mesh.x=md.mesh.x2d
        md.mesh.y=md.mesh.y2d
        md.mesh.elements=md.mesh.elements2d
        md.mesh=mesh2d(md.mesh)

    #read contours: 
    profiles=expread(options.getfieldvalue('sectionvalue'))
    numprofiles=len(profiles)

    if DEBUG:
        print(f'number of profiles: {numprofiles}')

    #Loop over number of profiles: 
    for profile_i in range(numprofiles):
        profile=profiles[profile_i]

        #Loop over number of curves
        for i in range(numcurves):
            
            datai, datatype=processdata(md3d,data[:,i],options)

            #resolution
            if options.exist('resolution'):
                resolution=options.getfieldvalue('resolution');
            else: #Default resolution
                if is2d:
                    resolution=[1000, 1]
                else:
                    resolution=[1000, 10*md.mesh.numberoflayers]
                
                print('plot_section warning: no resolution specified, use default resolution: [horizontal_resolution vertical_resolution]=[' + str(resolution)  + ']')

            #Compute section value
            elements,x,y,z,s,data_s=SectionValues(md,datai,profile,resolution)

            if options.getfieldvalue('sectionmean',0)==1:
                print('Mean value of data along section: ' + str(np.mean(data_s)))
                print('Median value of data along section: ' + str(np.median(data_s)))
                print('Standard deviation of data along section: ' + str(np.std(data_s)))

            #units
            if options.exist('unit'):
                unit=options.getfieldvalue('unit')
                x=x*unit
                y=y*unit
                z=z*unit
                s=s*unit

            #caxis
            if not options.exist('caxis'):
                caxis_opt = [min(data_s), max(data_s)]
                options.changefieldvalue('caxis',caxis_opt)
            else:
                caxis_opt = options.getfieldvalue('caxis')
                data = np.maximum(data, caxis_opt[0]) # min
                data = np.minimum(data, caxis_opt[1]) # max

            #2D
            if is2d:
    # %		%plot section value
    # %		hold on;
    # %		subplot(nlines,ncols,index1)
    # %		%subplot(1,3,[2 3])
    # %		plot(s,data_s,'color',color(i,:),'LineWidth',getfieldvalue(options,'linewidth',1))
    # %		%3D
    # %	else
    # %		%plot section value
    # %		%if user requested view2: 2d plot with curvilinear coordinate
    # %		if (getfieldvalue(options,'view',3)==2 )

                #Show Section if requested by user
                if options.exist('showsection'):
                    raise Exception('Error: "showsection" does not work in python.')
                    # %compute number of labels
                    # numlabels=min(getfieldvalue(options,'showsection'),length(s));
                    # shift=fix(length(s)/numlabels);

                    # %plot labels on current graph
                    # hold on
                    # text(s(1),data_s(1),'1','backgroundcolor',[0.8 0.9 0.8])
                    # for i=2:numlabels-1
                    #     text(s(1+(i-1)*shift),data_s(1+(i-1)*shift),num2str(i),'backgroundcolor',[0.8 0.9 0.8])
                    # end
                    # text(s(end),data_s(end),'end','backgroundcolor',[0.8 0.9 0.8])

                    # %plot section only with labels
                    # subplot(nlines,ncols,index2)
                    # plot_unit(x_m,y_m,z_m,elements_m,data(:,i),is2d,isplanet,datatype,options)
                    # hold on
                    # text(x(1),y(1),'1','backgroundcolor',[0.8 0.9 0.8])
                    # for i=2:numlabels-1
                    #     text(x(1+(i-1)*shift),y(1+(i-1)*shift),num2str(i),'backgroundcolor',[0.8 0.9 0.8])
                    # end
                    # text(x(end),y(end),'end','backgroundcolor',[0.8 0.9 0.8])
                    # plot(x,y,'-r')
                    # axis([min(md.mesh.x)-1 max(md.mesh.x)+1 min(md.mesh.y)-1 max(md.mesh.y)+1])
                    # view(2)

                #plot section value
                if(i==0):
                    ax=axgrid[gridindex]
                ax.plot(s,data_s,'color',color[i,:],'LineWidth',options.getfieldvalue('linewidth',1))

                #3D
            else:
                if DEBUG:
                    print('Catch section value')
                #plot section value
                #if user requested view2: 2d plot with curvilinear coordinate
                if (options.getfieldvalue('view',2)==2):

                    #Show Section if requested by user
                    if options.exist('showsection'):

                        #compute number of labels
                        numlabels=np.minimum(options.getfieldvalue('showsection'),len(s))
                        shift=np.fix(len(s)/numlabels)

                        #plot labels on current graph
                        ax.text(s[0],z[0],'1',backgroundcolor=[0.8,0.9,0.8])
                        for i in range(1,numlabels):
                            ax.text(s[1+i*shift],z[1+i*shift],str(i+1),backgroundcolor=[0.8,0.9,0.8])

                        ax.text(s[-1],z[-1],'end','backgroundcolor',[0.8,0.9,0.8])

                        #plot section only with labels
                        # subplot(nlines,ncols,index2)
                        ax=axgrid[index2]
                        plot_unit(x_m,y_m,z_m,elements_m,data[:,i],is2d,isplanet,datatype,options,fig,axgrid,index2)

                        ax.text(x[0],y[0],'1',backgroundcolor=[0.8,0.9,0.8])
                        for i in range(1,numplots):
                            ax.text(x[1+i*shift],y[1+i*shift],str(i+1),backgroundcolor=[0.8,0.9,0.8])

                        ax.text(x[-1],y[-1],'end',backgroundcolor=[0.8, 0.9, 0.8])
                        ax.plot(x,y,color='r',linestyle='-')
                        ax.set_xlim([min(md.mesh.x)-1, max(md.mesh.x)+1])
                        ax.set_ylim([min(md.mesh.y)-1, max(md.mesh.y)+1])
                        #view(2)

                    #NOTE: Python does not support "patch" like function. Therefore, this work should be different from Matlab.
                    #A=elements[:,0]; B=elements[:,1]; C=elements[:,2];  D=elements[:,3]
                    #patch( 'Faces', [A B C D], 'Vertices', [s z zeros(length(s),1)],'FaceVertexCData',data_s,'FaceColor','interp','EdgeColor','none');
                    ax = axgrid[gridindex]
                    if 0:
                        #vertices = np.vstack([s, z])
                        polygons = np.stack([s[elements-1], z[elements-1]],axis=2)
                        print(polygons.shape)

                        if not options.exist('caxis'):
                            options.addfielddefault('caxis',[min(data_s), max(data_s)])
                        caxis=options.getfieldvalue('caxis')
                        norm=mpl.colors.Normalize(vmin=caxis[0],vmax=caxis[1])
                        #options.addfielddefault('colornorm',norm)

                        collection = mpl.collections.PolyCollection(polygons,
                                                                    array=np.mean(data_s[elements-1],axis=1),
                                                                    cmap=options.getfieldvalue('colormap'),
                                                                    norm=norm,
                                                                    rasterized=True,
                                                                    edgecolors='none')
                        ax.add_collection(collection)
                    else:
                        triangles = np.vstack([
                            elements[:,[0,1,2]]-1,
                            elements[:,[0,2,3]]-1
                            ])
                        tri = mpl.tri.Triangulation(s, z, triangles)

                        if not options.exist('colormap'):
                            options.addfielddefault('colormap',mpl.cm.viridis)

                        #NOTE: tricontourf does not support nan value in dataset.
                        data_s[np.isnan(data_s)]=0.0

                        ax.tricontourf(tri, data_s, edgecolor='none',
                                    levels=256,
                                    cmap=options.getfieldvalue('colormap'),
                                    vmin=caxis_opt[0],vmax=caxis_opt[1])
                        #ax.tripcolor(tri, data_s, edgecolor='none',
                        #             cmap=options.getfieldvalue('colormap'))
                                    #shading='interp',
                else:
                    raise Exception('Error: ''view''==3 is yet not supported in Python.')

                    #Show Section if requested by user
#                    if option.exist('showsection'):
#
#                        #compute number of labels
#                        numlabels=min(options.getfieldvalue('showsection'),len(s))
#                        shift=np.fix(len(x)/numlabels)
#
#                        #plot labels on current graph
#                        #hold on
#                        ax=axgrid[gridindex]
#                        ax.text(x(1),y(1),z(1),'1',backgroundcolor=[0.8,0.9,0.8])
#                        for i in range(1,numlabels):
#                            ax.text(x(1+(i-1)*shift),y(1+(i-1)*shift),z(1+(i-1)*shift),num2str(i),'backgroundcolor',[0.8 0.9 0.8])
#                        ax.text(x[-1],y[-1],z[-1],'end',backgroundcolor=[0.8,0.9,0.8])
#
#                        #plot section only with labels
#                        subplot(nlines,ncols,index2)
#                        plot_unit(x_m,y_m,z_m,elements_m,data,is2d,datatype,options)
#                        hold on
#                        text(x(1),y(1),'1','backgroundcolor',[0.8 0.9 0.8])
#                        for i=2:numlabels-1
#                            text(x(1+(i-1)*shift),y(1+(i-1)*shift),num2str(i),'backgroundcolor',[0.8 0.9 0.8])
#                        end
#                        text(x(end),y(end),'end','backgroundcolor',[0.8 0.9 0.8])
#                        plot(x,y,'-r')
#                        axis([min(md.mesh.x)-1 max(md.mesh.x)+1 min(md.mesh.y)-1 max(md.mesh.y)+1])
#                        view(2)
#
#                    #subplot(nlines,ncols,index1)
#                    A=elements(:,1); B=elements(:,2); C=elements(:,3);  D=elements(:,4); 
#                    patch( 'Faces', [A B C D], 'Vertices', [x y z],'FaceVertexCData',data_s,'FaceColor','interp','EdgeColor','none');
#                    #view(3)

    #apply options
    options.addfielddefault('title','Section value')
    if md.mesh.dimension() == 2:
        options.addfielddefault('colorbar',0)
    if (md.mesh.dimension()==2) | (options.getfieldvalue('view',2)==2):
        options.addfielddefault('xlabel','Curvilinear coordinate')
        options.addfielddefault('axis',True)
    if (md.mesh.dimension()==3) | (options.getfieldvalue('view',2)==2):
        options.addfielddefault('ylabel','z')
        options.addfielddefault('axis',True)
    applyoptions(md3d, data, options, fig, axgrid, gridindex)
