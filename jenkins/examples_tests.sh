#!/bin/bash
################################################################################
# This script runs the examples tests (i.e. contents of examples directory, 
# which are implementations of the tutorials found at 
# https://issm.jpl.nasa.gov/documentation/tutorials/). It is intended to be 
# called from jenkins/jenkins.sh.
#
# runme files are modified as needed to fill in statements that would otherwise 
# be added by user.
#
# NOTE:
# - Indentation of replacement string literals (e.g. 'STEP_EIGHT') is set to 
#	nest cleanly in this file, but will result in unclean nesting in the runme 
#	files (which should not be an issue)
# - Single-line string replacements in runme.m can effectively be performed 
#	using sed. When performing multi-line replacements, perl is a better 
#	option.
#
# TODO:
# - Figure out how to remove \ and \n\ from multiline string variables while 
#	preserving formatting when value is printed to file.
################################################################################

## Constants
#
RUNME_FILE='runme.m'
RUN_EXAMPLE=0
STATUS_HANDLING="\
		disp('SUCCESS');\n\
	catch me\n\
		message=getReport(me);\n\
		fprintf('%s',message);\n\
		disp('FAILURE');\n\
	end\n\
	exit\n\
"

cd $ISSM_DIR/examples

for dir in ./* ; do
	if [ -d "${dir}" ]; then
		# # Temporary short circuit to check single example
		# example="./AMR"
		# if [ "${dir}" != "${example}" ]; then
		# 	continue
		# fi

		# Some of the examples are incomplete (on purpose). As such, we will 
		# have to populate the missing steps in order to make sure that 
		# everything is working.

		cd ${dir}
		if [ "${dir}" == "./AMR" ]; then
			sed -i.bak -e '1 s|^.*$|try\n\n&|' $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./Data" ]; then
			echo 'Directory contains datasets only; no example to run.'
			RUN_EXAMPLE=0
		elif [ "${dir}" == "./EsaGRACE" ]; then
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:5\];\n\ntry\n|' $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./EsaWahr" ]; then
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:7\];\n\ntry\n|' $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./Functions" ]; then
			echo "Directory contains functions only; no example to run."
			RUN_EXAMPLE=0
		elif [ "${dir}" == "./Greenland" ]; then
			# STEP_SEVEN #{{{
			STEP_SEVEN="\
				if any(steps==7)\n\
					disp('   Step 7: Historical Relaxation run');\n\
					md = loadmodel('./Models/Greenland.Control_drag');\n\
					\n\
					load smbbox\n\
					\n\
					%convert mesh x,y into the Box projection\n\
					[md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,39,71);\n\
					[xi,yi]= ll2xy(md.mesh.lat,md.mesh.long,+1,45,70);\n\
					\n\
					%Interpolate and set surface mass balance\n\
					index = BamgTriangulate(x1(:),y1(:));\n\
					smb_mo = InterpFromMeshToMesh2d(index,x1(:),y1(:),smbmean(:),xi,yi);\n\
					smb = smb_mo*12/1000*md.materials.rho_freshwater/md.materials.rho_ice;\n\
					md.smb.mass_balance = [smb;1 ];\n\
					\n\
					%Set transient options, run for 20 years, saving every 5 timesteps\n\
					md.timestepping.time_step=0.2;\n\
					md.timestepping.final_time=200;\n\
					md.settings.output_frequency=5;\n\
					\n\
					%Additional options\n\
					md.inversion.iscontrol=0;\n\
					md.transient.requested_outputs={'IceVolume','TotalSmb', ...\n\
						'SmbMassBalance'};\n\
					md.verbose=verbose('solution',true,'module',true);\n\
					\n\
					%Go solve\n\
					md.cluster=generic('name',oshostname,'np',2);\n\
					md=solve(md,'Transient');\n\
					\n\
					save ./Models/Greenland.HistoricTransient_200yr md;\n\
				end\n\
			"
			#}}}
			# STEP_EIGHT #{{{
			STEP_EIGHT="\
				if any(steps==8)\n\
					%Load historic transient model\n\
					md=loadmodel('./Models/Greenland.HistoricTransient_200yr');\n\
					\n\
					%Create Line Plots of relaxation run. Create a figure.\n\
					figure;\n\
					\n\
					%Save surface mass balance, by looping through 200 years (1000 steps)\n\
					%Note, the first output will always contain output from time step 1\n\
					surfmb=[];\n\
					for i=2:201;\n\
						surfmb=[surfmb md.results.TransientSolution(i).SmbMassBalance];\n\
					end\n\
					\n\
					%Plot surface mass balance time series in first subplot\n\
					subplot(3,1,1);\n\
					plot([1:200],mean(surfmb));\n\
					\n\
					%Title this plot Mean surface mass balance\n\
					title('Mean Surface mass balance');\n\
					\n\
					%Save velocity by looping through 200 years\n\
					vel=[];\n\
					for i=2:201;\n\
						vel=[vel md.results.TransientSolution(i).Vel];\n\
					end\n\
					\n\
					%Plot velocity time series in second subplot\n\
					subplot(3,1,2);\n\
					plot([1:200],mean(vel));\n\
					\n\
					%Title this plot Mean Velocity\n\
					title('Mean Velocity');\n\
					\n\
					%Save Ice Volume by looping through 200 years\n\
					volume=[];\n\
					for i=2:201;\n\
						volume=[volume md.results.TransientSolution(i).IceVolume];\n\
					end\n\
					\n\
					%Plot volume time series in third subplot\n\
					subplot(3,1,3);\n\
					plot([1:200],volume);\n\
					\n\
					%Title this plot Mean Velocity and add an x label of years\n\
					title('Ice Volume');\n\
					xlabel('years');\n\
				end\n\
			"
			#}}}
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:8\];\n\ntry\n|' $RUNME_FILE
			perl -0755 -p -i -e "s|if any\(steps==7\).*% step 7 end|${STEP_SEVEN}|s" $RUNME_FILE
			perl -0755 -p -i -e "s|if any\(steps==8\).*% step 8 end|${STEP_EIGHT}|s" $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./IceBridge" ]; then
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:5\];\n\ntry\n|' $RUNME_FILE
			perl -0755 -p -i -e "s|\n\t%Mesh greenland without.*return;\n||s" $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./IceflowModels" ]; then
			sed -i.bak -e '1 s|^.*$|try\n\n&|' $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./Inversion" ]; then
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:4\];\n\ntry\n|' $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./ISMIP" ]; then
			# TODO:
			# - Run test again with ISMIPF configuration (will likely need to 
			#	add conditional after 'RUN_EXAMPLE -eq 1' block)
			#

			# RUNME #{{{
			RUNME="\
				try\n\
					%which steps to perform; steps are from 1 to 8\n\
					%step 7 is specific to ISMIPA\n\
					%step 8 is specific to ISMIPF\n\
					\n\
					steps=[1:7]; %ISMIPA\n\
					%steps=[1:6,8]; %ISMIPF\n\
					\n\
					% parameter file to be used, choose between IsmipA.par or IsmipF.par\n\
					ParamFile='IsmipA.par';\n\
					%ParamFile='IsmipF.par';\n\
					\n\
					%Run Steps\n\
					\n\
					%Mesh Generation #1\n\
					if any(steps==1)\n\
						%initialize md as a new model #help model\n\
						%->\n\
						md=model();\n\
						% generate a squaremesh #help squaremesh\n\
						% Side is 80 km long with 20 points\n\
						%->\n\
						if(ParamFile=='IsmipA.par'),\n\
							md=squaremesh(md,80000,80000,20,20);\n\
						elseif(ParamFile=='IsmipF.par'),\n\
							md=squaremesh(md,100000,100000,30,30);\n\
						end\n\
						% plot the given mesh #plotdoc\n\
						%->\n\
						plotmodel(md,'data','mesh')\n\
						% save the given model\n\
						%->\n\
						save ./Models/ISMIP.Mesh_generation md;\n\
					end\n\
					\n\
					%Masks #2\n\
					if any(steps==2)\n\
						% load the preceding step #help loadmodel\n\
						% path is given by the organizer with the name of the given step\n\
						%->\n\
						md = loadmodel('./Models/ISMIP.Mesh_generation');\n\
						% set the mask #help setmask\n\
						% all MISMIP nodes are grounded\n\
						%->\n\
						md=setmask(md,'','');\n\
						% plot the given mask #md.mask to locate the field\n\
						%->\n\
						plotmodel(md,'data',md.mask.ocean_levelset);\n\
						% save the given model\n\
						%->\n\
						save ./Models/ISMIP.SetMask md;\n\
					end\n\
					\n\
					%Parameterization #3\n\
					if any(steps==3)\n\
						% load the preceding step #help loadmodel\n\
						% path is given by the organizer with the name of the given step\n\
						%->\n\
						md = loadmodel('./Models/ISMIP.SetMask');\n\
						% parametrize the model # help parameterize\n\
						% you will need to fil-up the parameter file defined by the\n\
						% ParamFile variable\n\
						%->\n\
						md=parameterize(md,ParamFile);\n\
						% save the given model\n\
						%->\n\
						save ./Models/ISMIP.Parameterization md;\n\
					end\n\
					\n\
					%Extrusion #4\n\
					if any(steps==4)\n\
						\n\
						% load the preceding step #help loadmodel\n\
						% path is given by the organizer with the name of the given step\n\
						%->\n\
						md = loadmodel('./Models/ISMIP.Parameterization');\n\
						% vertically extrude the preceding mesh #help extrude\n\
						% only 5 layers exponent 1\n\
						%->\n\
						md=extrude(md,5,1);\n\
						% plot the 3D geometry #plotdoc\n\
						%->\n\
						plotmodel(md,'data',md.geometry.base)\n\
						% save the given model\n\
						%->\n\
						save ./Models/ISMIP.Extrusion md;\n\
					end\n\
					\n\
					%Set the flow computing method #5\n\
					if any(steps==5)\n\
						\n\
						% load the preceding step #help loadmodel\n\
						% path is given by the organizer with the name of the given step\n\
						%->\n\
						md = loadmodel('./Models/ISMIP.Extrusion');\n\
						% set the approximation for the flow computation #help setflowequation\n\
						% We will be using the Higher Order Model (HO)\n\
						%->\n\
						md=setflowequation(md,'HO','all');\n\
						% save the given model\n\
						%->\n\
						save ./Models/ISMIP.SetFlow md;\n\
					end\n\
					\n\
					%Set Boundary Conditions #6\n\
					if any(steps==6)\n\
						\n\
						% load the preceding step #help loadmodel\n\
						% path is given by the organizer with the name of the given step\n\
						%->\n\
						md = loadmodel('./Models/ISMIP.SetFlow');\n\
						% dirichlet boundary condition are known as SPCs\n\
						% ice frozen to the base, no velocity	#md.stressbalance\n\
						% SPCs are initialized at NaN one value per vertex\n\
						%->\n\
						md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);\n\
						%->\n\
						md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);\n\
						%->\n\
						md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);\n\
						% extract the nodenumbers at the base #md.mesh.vertexonbase\n\
						%->\n\
						basalnodes=find(md.mesh.vertexonbase);\n\
						% set the sliding to zero on the bed\n\
						%->\n\
						md.stressbalance.spcvx(basalnodes)=0.0;\n\
						%->\n\
						md.stressbalance.spcvy(basalnodes)=0.0;\n\
						% periodic boundaries have to be fixed on the sides\n\
						% Find the indices of the sides of the domain, for x and then for y\n\
						% for x\n\
						% create maxX, list of indices where x is equal to max of x (use >> help find)\n\
						%->\n\
						maxX=find(md.mesh.x==max(md.mesh.x));\n\
						% create minX, list of indices where x is equal to min of x\n\
						%->\n\
						minX=find(md.mesh.x==min(md.mesh.x));\n\
						% for y\n\
						% create maxY, list of indices where y is equal to max of y\n\
						% but not where x is equal to max or min of x\n\
						% (i.e, indices in maxX and minX should be excluded from maxY and minY)\n\
						% but not where x is equal to max or min of x\n\
						%->\n\
						maxY=find(md.mesh.y==max(md.mesh.y) & md.mesh.x~=max(md.mesh.x) & md.mesh.x~=min(md.mesh.x));\n\
						% create minY, list of indices where y is equal to max of y\n\
						%->\n\
						minY=find(md.mesh.y==min(md.mesh.y) & md.mesh.x~=max(md.mesh.x) & md.mesh.x~=min(md.mesh.x));\n\
						% set the node that should be paired together, minX with maxX and minY with maxY\n\
						% #md.stressbalance.vertex_pairing\n\
						%->\n\
						md.stressbalance.vertex_pairing=[minX,maxX;minY,maxY];\n\
						if (ParamFile=='IsmipF.par')\n\
							% if we are dealing with IsmipF the solution is in\n\
							% masstransport\n\
							md.masstransport.vertex_pairing=md.stressbalance.vertex_pairing;\n\
						end\n\
						% save the given model\n\
						%->\n\
						save ./Models/ISMIP.BoundaryCondition md;\n\
					end\n\
					\n\
					%Solving #7\n\
					if any(steps==7)\n\
						% load the preceding step #help loadmodel\n\
						% path is given by the organizer with the name of the given step\n\
						%->\n\
						md = loadmodel('./Models/ISMIP.BoundaryCondition');\n\
						% Set cluster #md.cluster\n\
						% generic parameters #help generic\n\
						% set only the name and number of process\n\
						%->\n\
						md.cluster=generic('name',oshostname(),'np',2);\n\
						% Set which control message you want to see #help verbose\n\
						%->\n\
						md.verbose=verbose('convergence',true);\n\
						% Solve #help solve\n\
						% we are solving a StressBalanc\n\
						%->\n\
						md=solve(md,'Stressbalance');\n\
						% save the given model\n\
						%->\n\
						save ./Models/ISMIP.StressBalance md;\n\
						% plot the surface velocities #plotdoc\n\
						%->\n\
						plotmodel(md,'data',md.results.StressbalanceSolution.Vel)\n\
					end\n\
					\n\
					%Solving #8\n\
					if any(steps==8)\n\
						% load the preceding step #help loadmodel\n\
						% path is given by the organizer with the name of the given step\n\
						%->\n\
						md = loadmodel('./Models/ISMIP.BoundaryCondition');\n\
						% Set cluster #md.cluster\n\
						% generic parameters #help generic\n\
						% set only the name and number of process\n\
						%->\n\
						md.cluster=generic('name',oshostname(),'np',2);\n\
						% Set which control message you want to see #help verbose\n\
						%->\n\
						md.verbose=verbose('convergence',true);\n\
						% set the transient model to ignore the thermal model\n\
						% #md.transient \n\
						%->\n\
						md.transient.isthermal=0;\n\
						% define the timestepping scheme\n\
						% everything here should be provided in years #md.timestepping\n\
						% give the length of the time_step (4 years)\n\
						%->\n\
						md.timestepping.time_step=4;\n\
						% give final_time (20*4 years time_steps)\n\
						%->\n\
						md.timestepping.final_time=4*20;\n\
						% Solve #help solve\n\
						% we are solving a TransientSolution\n\
						%->\n\
						md=solve(md,'Transient');\n\
						% save the given model\n\
						%->\n\
						save ./Models/ISMIP.Transient md;\n\
						% plot the surface velocities #plotdoc\n\
						%->\n\
						plotmodel(md,'data',md.results.TransientSolution(20).Vel)\n\
					end\n\
			"
			#}}}
			# PAR_A #{{{
			PAR_A="\
				%Parameterization for ISMIP A experiment\n\
				\n\
				%Set the Simulation generic name #md.miscellaneous\n\
				%->\n\
				\n\
				%Geometry\n\
				disp('   Constructing Geometry');\n\
				\n\
				%Define the geometry of the simulation #md.geometry\n\
				%surface is [-x*tan(0.5*pi/180)] #md.mesh\n\
				%->\n\
				md.geometry.surface=-md.mesh.x*tan(0.5*pi/180.);\n\
				%base is [surface-1000+500*sin(x*2*pi/L).*sin(y*2*pi/L)]\n\
				%L is the size of the side of the square #max(md.mesh.x)-min(md.mesh.x)\n\
				%->\n\
				L=max(md.mesh.x)-min(md.mesh.x);\n\
				md.geometry.base=md.geometry.surface-1000.0+500.0*sin(md.mesh.x*2.0*pi/L).*sin(md.mesh.y*2.0*pi/L);\n\
				%thickness is the difference between surface and base #md.geometry\n\
				%->\n\
				md.geometry.thickness=md.geometry.surface-md.geometry.base;\n\
				%plot the geometry to check it out\n\
				%->\n\
				plotmodel(md,'data',md.geometry.thickness);\n\
				\n\
				disp('   Defining friction parameters');\n\
				\n\
				%These parameters will not be used but need to be fixed #md.friction\n\
				%one friciton coefficient per node (md.mesh.numberofvertices,1)\n\
				%->\n\
				md.friction.coefficient=200.0*ones(md.mesh.numberofvertices,1);\n\
				%one friciton exponent (p,q) per element\n\
				%->\n\
				md.friction.p=ones(md.mesh.numberofelements,1);\n\
				%->\n\
				md.friction.q=ones(md.mesh.numberofelements,1);\n\
				\n\
				disp('   Construct ice rheological properties');\n\
				\n\
				%The rheology parameters sit in the material section #md.materials\n\
				%B has one value per vertex\n\
				%->\n\
				md.materials.rheology_B=6.8067e7*ones(md.mesh.numberofvertices,1);\n\
				%n has one value per element\n\
				%->\n\
				md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);\n\
				\n\
				disp('   Set boundary conditions');\n\
				\n\
				%Set the default boundary conditions for an ice-sheet \n\
				% #help SetIceSheetBC\n\
				%->\n\
				md=SetIceSheetBC(md);\n\
			"
			#}}}
			# PAR_F #{{{
			PAR_F="\
				%Parameterization for ISMIP F experiment\n\
				\n\
				%Set the Simulation generic name #md.miscellaneous\n\
				%->\n\
				\n\
				%Geometry\n\
				disp('   Constructing Geometry');\n\
				\n\
				%Define the geometry of the simulation #md.geometry\n\
				%surface is [-x*tan(3.0*pi/180)] #md.mesh\n\
				%->\n\
				md.geometry.surface=-md.mesh.x*tan(3.0*pi/180.0);\n\
				%base is [surface-1000+100*exp(-((x-L/2).^2+(y-L/2).^2)/(10000.^2))]\n\
				%L is the size of the side of the square #max(md.mesh.x)-min(md.mesh.x)\n\
				%->\n\
				L=max(md.mesh.x)-min(md.mesh.x);\n\
				%->\n\
				md.geometry.base=md.geometry.surface-1000.0+100.0*exp(-((md.mesh.x-L/2.0).^2.0+(md.mesh.y-L/2.0).^2.0)/(10000.^2.0));\n\
				%thickness is the difference between surface and base #md.geometry\n\
				%->\n\
				md.geometry.thickness=md.geometry.surface-md.geometry.base;\n\
				%plot the geometry to check it out\n\
				%->\n\
				plotmodel(md,'data',md.geometry.thickness);\n\
				\n\
				disp('   Defining friction parameters');\n\
				\n\
				%These parameters will not be used but need to be fixed #md.friction\n\
				%one friciton coefficient per node (md.mesh.numberofvertices,1)\n\
				%conversion from year to seconds with #md.constants.yts\n\
				%->\n\
				md.friction.coefficient=sqrt(md.constants.yts/(1000*2.140373*10^-7))*ones(md.mesh.numberofvertices,1);\n\
				%one friciton exponent (p,q) per element\n\
				%->\n\
				md.friction.p=ones(md.mesh.numberofelements,1);\n\
				%->\n\
				md.friction.q=zeros(md.mesh.numberofelements,1);\n\
				\n\
				disp('   Construct ice rheological properties');\n\
				\n\
				%The rheology parameters sit in the material section #md.materials\n\
				%B has one value per vertex\n\
				%->\n\
				md.materials.rheology_B=(1/(2.140373*10^-7/md.constants.yts))*ones(md.mesh.numberofvertices,1);\n\
				%n has one value per element\n\
				%->\n\
				md.materials.rheology_n=1*ones(md.mesh.numberofelements,1);\n\
				\n\
				disp('   Set boundary conditions');\n\
				\n\
				%Set the default boundary conditions for an ice-sheet \n\
				% #help SetIceSheetBC\n\
				%->\n\
				md=SetIceSheetBC(md);\n\
				\n\
				disp('   Initializing velocity and pressure');\n\
				\n\
				%initialize the velocity and pressurefields of #md.initialization\n\
				%->\n\
				md.initialization.vx=zeros(md.mesh.numberofvertices,1);\n\
				%->\n\
				md.initialization.vy=zeros(md.mesh.numberofvertices,1);\n\
				%->\n\
				md.initialization.vz=zeros(md.mesh.numberofvertices,1);\n\
				%->\n\
				md.initialization.pressure=zeros(md.mesh.numberofvertices,1);\n\
			"
			#}}}
			perl -0755 -p -i'.bak' -e "s|^.*$|${RUNME}|s" $RUNME_FILE
			perl -0755 -p -i'.bak' -e "s|^.*$|${PAR_A}|s" IsmipA.par
			perl -0755 -p -i'.bak' -e "s|^.*$|${PAR_F}|s" IsmipF.par
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./Jakobshavn" ]; then
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:4\];\n\ntry\n|' $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./LcurveAnalysis" ]; then
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:4\];\n\ntry\n|' $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./Mesh" ]; then
			# NOTE: Cannot test exptool region selection without GUI
			#

			# RUNME #{{{
			RUNME="\
				try\n\
					steps=[1:7];\n\
					\n\
					if any(steps==1) % Model\n\
						md=model;\n\
					end\n\
					\n\
					if any(steps==2) % squaremesh\n\
						md=squaremesh(md,100,200,15,25);\n\
						plotmodel(md,'data','mesh');\n\
					end\n\
					\n\
					if any(steps==3) % roundmesh\n\
						md=roundmesh(model,100,10);\n\
						plotmodel(md,'data','mesh');\n\
					end\n\
					\n\
					if any(steps==4) % triangle\n\
						md=triangle(model,'Square.exp',.2);\n\
						plotmodel(md,'data','mesh');\n\
					end\n\
					\n\
					if any(steps==5) % bamg\n\
						md=bamg(model,'domain','Square.exp','hmax',.05);\n\
						plotmodel(md,'data','mesh');\n\
					end\n\
					\n\
					if any(steps==6) % Non-Uniform mesh\n\
						hvertices=[0.2;0.2;0.005;0.2];\n\
						md=bamg(md,'domain','Square.exp','hvertices',hvertices);\n\
						plotmodel(md,'data','mesh');\n\
					end\n\
					\n\
					if any(steps==7) % Mesh adaptation\n\
						md=bamg(model,'domain','Square.exp','hmax',.05);\n\
						vel=shock(md.mesh.x,md.mesh.y);\n\
						plotmodel(md,'data',vel,'edgecolor','w');\n\
						\n\
						md=bamg(model,'domain','Square.exp','hmax',.005);\n\
						vel=shock(md.mesh.x,md.mesh.y);\n\
						md=bamg(md,'field',vel,'err',0.05,'hmin',0.005,'hmax',0.3);\n\
						vel=shock(md.mesh.x,md.mesh.y);\n\
						plotmodel(md,'data',vel,'edgecolor','w');\n\
						\n\
						md=bamg(model,'domain','Square.exp','hmax',.005);\n\
						vel=shock(md.mesh.x,md.mesh.y);\n\
						md=bamg(md,'field',vel,'err',0.03,'hmin',0.005,'hmax',0.3,'gradation',3);\n\
						vel=shock(md.mesh.x,md.mesh.y);\n\
						plotmodel(md,'data',vel,'edgecolor','w');\n\
						\n\
						md=bamg(model,'domain','Square.exp','hmax',.005);\n\
						vel=shock(md.mesh.x,md.mesh.y);\n\
						md=bamg(md,'field',vel,'err',0.03,'hmin',0.005,'hmax',0.3,'gradation',1.3,'anisomax',1);\n\
						vel=shock(md.mesh.x,md.mesh.y);\n\
						plotmodel(md,'data',vel,'edgecolor','w');\n\
					end\n\
					\n\
					if any(steps==8) % Mesh refinement in a specific region\n\
						plotmodel(md,'data',vel,'edgecolor','w');\n\
						exptool('refinement.exp');\n\
						h=NaN*ones(md.mesh.numberofvertices,1);\n\
						in=ContourToNodes(md.mesh.x,md.mesh.y,'refinement.exp',1);\n\
						h(find(in))=0.02;\n\
						plotmodel(md,'data',in,'edgecolor','w');\n\
						\n\
						vel=shock(md.mesh.x,md.mesh.y);\n\
						md=bamg(md,'field',vel,'err',0.03,'hmin',0.005,'hmax',0.3,'hVertices',h);\n\
						vel=shock(md.mesh.x,md.mesh.y);\n\
						plotmodel(md,'data',vel,'edgecolor','w');\n\
					end\n\
			"
			#}}}
			touch $RUNME_FILE
			perl -0755 -p -i'.bak' -e "s|^.*$|${RUNME}|s" $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./Pig" ]; then
			# STEP_SIX #{{{
			STEP_SIX="\
				if any(steps==6)\n\
					% Load Model\n\
					md = loadmodel('./Models/PIG_Control_drag');\n\
					% Disable inversion\n\
					md.inversion.iscontrol=0;\n\
					% Extrude Mesh\n\
					disp('   Extruding mesh');\n\
					number_of_layers=3;\n\
					md=extrude(md,number_of_layers,1);\n\
					% Set Flowequation\n\
					disp('   Using HO Ice Flow Model');\n\
					md=setflowequation(md,'HO','all');\n\
					% Solve\n\
					md=solve(md,'Stressbalance');\n\
					% Save Model\n\
					save ./Models/PIG_ModelHO md;\n\
				end\n\
			"
			#}}}
			mv ./DomainOutline.bkp ./DomainOutline.exp > /dev/null 2>&1
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:7\];\ntry\n|' $RUNME_FILE
			perl -0755 -p -i -e "s|if any\(steps==6\).*% step 6 end|${STEP_SIX}|s" $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./Pig2" ]; then
			STEP_NINE="\n disp('Needs work!'); exit"
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:9\];\n\ntry\n|' $RUNME_FILE
			perl -0755 -p -i -e "s|if any\(steps==9\).*% step 9 end|${STEP_NINE}|s" $RUNME_FILE
			RUN_EXAMPLE=0
		elif [ "${dir}" == "./PigSensitivity" ]; then
			# STEP_FOUR # {{{
			STEP_FOUR="\
				if any(steps==4)\n\
					%Load model\n\
					md = loadmodel('./Models/PIG_Transient');\n\
					\n\
					%Change external forcing basal melting rate and surface mass balance)\n\
					md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);\n\
					md.basalforcings.floatingice_melting_rate=25*ones(md.mesh.numberofvertices,1);\n\
					md.smb.mass_balance=2*md.smb.mass_balance;\n\
					\n\
					%Define time steps and time span of the simulation\n\
					md.timestepping.time_step=0.1;\n\
					md.timestepping.final_time=10;\n\
					\n\
					%Request additional outputs\n\
					md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};\n\
					\n\
					%Solve\n\
					md=solve(md,'Transient');\n\
					\n\
					%Plot\n\
					plotmodel(md, 'data', md.results.TransientSolution(1).Vel,...\n\
						'title#1', 'Velocity t=0 years (m/yr)',...\n\
						'data', md.results.TransientSolution(end).Vel,...\n\
						'title#2', 'Velocity t=10 years (m/yr)',...\n\
						'data', md.results.TransientSolution(1).MaskOceanLevelset,...\n\
						'title#3', 'Floating ice t=0 years',...\n\
						'data', md.results.TransientSolution(end).MaskOceanLevelset,...\n\
						'title#4', 'Floating ice t=10 years',...\n\
						'caxis#1',([0 4500]),'caxis#2',([0 4500]),...\n\
						'caxis#3',([-1,1]),'caxis#4',([-1,1]));\n\
					\n\
					%Save model\n\
					save ./Models/PIG_SMB md;\n\
				end\n\
			"
			#}}}
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:4\];\n\ntry\n|' $RUNME_FILE
			perl -0755 -p -i -e "s|if any\(steps==4\).*% step 4 end|${STEP_FOUR}|s" $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./shakti" ]; then
			sed -i.bak -e 's|steps=\[1:3\];|steps=\[1:3\];\n\ntry\n|' $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./SlrFarrell" ]; then
			# TODO: Convert from md.slr
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:5\];\n\ntry\n|' $RUNME_FILE
			RUN_EXAMPLE=0
		elif [ "${dir}" == "./SlrGRACE" ]; then
			# TODO: Convert from md.slr
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:7\];\n\ntry\n|' $RUNME_FILE
			RUN_EXAMPLE=0
		elif [ "${dir}" == "./SlrGRACE_NIMS" ]; then
			# TODO: Convert from md.slr
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:8\];\n\ntry\n|' $RUNME_FILE
			RUN_EXAMPLE=0
		elif [ "${dir}" == "./SquareIceShelf" ]; then
			sed -i.bak -e '1 s|^.*$|try\n\n&|' $RUNME_FILE
			RUN_EXAMPLE=1
		elif [ "${dir}" == "./UncertaintyQuantification" ]; then
			sed -i.bak -e 's|steps=\[1\];|steps=\[1:7\];\n\ntry\n|' $RUNME_FILE
			RUN_EXAMPLE=1
		else
			echo "Not implemented yet!"
			exit 1
		fi

		if [ $RUN_EXAMPLE -eq 1 ]; then
			echo "Testing example: $(basename $dir)"
			LOG_RUNME_FILE="matlab_log_$(basename $dir)_examples.log"
			echo -e ${STATUS_HANDLING} >> ${RUNME_FILE}
			$MATLAB_PATH/bin/matlab -nodisplay -nosplash -r "addpath $ISSM_DIR/src/m/dev; devpath; addpath $ISSM_DIR/nightlylog/; runme" -logfile $ISSM_DIR/nightlylog/$LOG_RUNME_FILE
			echo "starting: $(basename $dir)" >> $ISSM_DIR/nightlylog/matlab_log_examples.log
			cat $ISSM_DIR/nightlylog/$LOG_RUNME_FILE >> $ISSM_DIR/nightlylog/matlab_log_examples.log
			echo "finished: $(basename $dir)" >> $ISSM_DIR/nightlylog/matlab_log_examples.log
			mv -f ${RUNME_FILE}.bak ${RUNME_FILE}
		fi

		# Extra clean up
		if [ "${dir}" == "./ISMIP" ]; then
			mv -f IsmipA.par.bak IsmipA.par
			mv -f IsmipF.par.bak IsmipF.par
		fi

		if [ "${dir}" == "./Mesh" ]; then
			rm -f $RUNME_FILE
		fi

		if [ "${dir}" == "./Pig" ]; then
			mv -f DomainOutline.exp DomainOutline.bkp
		fi

		cd ..
	fi
done
