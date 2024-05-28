if 7 in steps:
    #  Step 7: Historical Relaxation run {{{
    print('   Step 7: Historical Relaxation run')
    md = loadmodel('./Models/Greenland.Control_drag.nc')

    with open('./smbbox.npz', "rb") as smbFile:
        SMB = np.load(smbFile)
        x1 = np.squeeze(SMB['x1'])
        y1 = np.squeeze(SMB['y1'])
        smbmean = np.squeeze(SMB['smbmean'])

    #convert mesh x, y into the Box projection
    [md.mesh.lat, md.mesh.long] = xy2ll(md.mesh.x, md.mesh.y, + 1, 39, 71)
    [xi, yi] = ll2xy(md.mesh.lat, md.mesh.long, + 1, 45, 70)

    #Interpolate and set surface mass balance
    x1 = x1.flatten()
    y1 = y1.flatten()
    smbmean = smbmean.flatten()
    index = BamgTriangulate(x1, y1)
    smb_mo = InterpFromMeshToMesh2d(index, x1, y1, smbmean, xi, yi)
    smb = smb_mo * 12 / 1000 * md.materials.rho_freshwater / md.materials.rho_ice
    md.smb.mass_balance = np.append(smb, 1)

    #Set transient options, run for 20 years, saving every timestep
    md.timestepping.time_step = 0.2
    md.timestepping.final_time = 200
    md.settings.output_frequency = 5

    #Additional options
    md.inversion.iscontrol = 0
    md.transient.requested_outputs = ['IceVolume', 'TotalSmb', 'SmbMassBalance']
    md.verbose = verbose('solution', True, 'module', True)

    #Go solve
    md.cluster = generic('name', gethostname(), 'np', 2)
    md = solve(md, 'Transient')

    export_netCDF(md, './Models/Greenland.HistoricTransient_200yr.nc')
    # }}}

if 8 in steps:
    # Step 8: Plotting exercise {{{
    print('   Step 8: Plotting exercise')
    md = loadmodel('./Models/Greenland.HistoricTransient_200yr.nc')

    #Create Line Plots of relaxation run. Create a figure.
    fig = plt.figure(tight_layout=True)

    #Save surface mass balance, by looping through 200 years (1000 steps)
    #Note, the first output will always contain output from time step 1

    Timer = np.arange(0.2, 200.2)
    surfmb = []
    for i in range(0, 200):
        surfmb.append(md.results.TransientSolution[i].SmbMassBalance)

    #Plot surface mass balance time series in first subplot
    ax = fig.add_subplot(311)
    ax.plot(Timer, np.nanmean(surfmb, axis=1))

    #Title this plot Mean surface mass balance
    ax.set_title('Mean Surface mass balance')

    #Save velocity by looping through 200 years
    vel = []
    for i in range(0, 200):
        vel.append(md.results.TransientSolution[i].Vel)

    #Plot velocity time series in second subplot
    ay = fig.add_subplot(312)
    ay.plot(Timer, np.nanmean(vel, axis=1))

    #Title this plot Mean Velocity
    ay.set_title('Mean Velocity')

    #Save Ice Volume by looping through 200 years
    volume = []
    for i in range(0, 200):
        volume.append(md.results.TransientSolution[i].IceVolume)

    #Plot volume time series in third subplot
    az = fig.add_subplot(313)
    az.plot(Timer, volume)

    #Title this plot Mean Velocity and add an x label of years
    az.set_title('Ice Volume')
    az.set_xlabel('years')
    # }}}
