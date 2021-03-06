
### 4) Flagging overlapping circles

## 4a) We sort the circles by galaxy density, and divide the circles in two samples: underdense and overdense.

# Sorting circles by density (rho) from lowest to highest
sort_rho = np.argsort(rhotheta_tot)
# Sorting overdense circles (second half of the sorted list) from highest to lowest
sort_rho_crab = [np.append( sort_rho[theta, 0:int(Ngrid_tot/2)], \
                    np.flip(sort_rho[theta, int(Ngrid_tot/2):Ngrid_tot], 0) ) for theta in range(Ntheta)]


# This array will contain the selected troughs (for all fields)
Stheta_tot = np.zeros([Ntheta, Ngrid_tot])
Ptheta_tot = np.zeros([Ntheta, Ngrid_tot])

for theta in range(Ntheta):

    ## 4b) Definition: Two circles are "overlapping" when their centers are closer than 0.5 \theta.
    thetamax = 0.5 * thetalist[theta]
    print('theta=', thetalist[theta]*60.)

    # This list will contain the selected (1) and flagged (0) troughs for this radius theta (for all fields)
    selected_theta = np.ones(Ngrid_tot)
    
    # Sort the troughs by density
    selID = gridID_tot[sort_rho[theta]]
    selcoords = gridcoords_tot[sort_rho[theta]]

    # Creating smaller grid point samples to avoid memory overload
    sampindex_tot = np.array(np.append(np.arange(0., Ngrid_tot, 1e4), Ngrid_tot), dtype=int)
    
    # Looping over the smaller samples of grid points (to avoid memory overload)
    for s in range(len(sampindex_tot)-1):
    #for s in range(1):
    
        print('Removing troughs within theta_max =' , thetamax*60., ', Grid sample: %i of %i'%(s+1, len(sampindex_tot)-1))
        
        # Defining the smaller grid point sample
        selsamp = selcoords[sampindex_tot[s]:sampindex_tot[s+1]]
        
        # 4c/d) For the underdense/overdense circles: We start from the lowest/highest density circle and flag all overlapping circles.
        
        # Find grid points within thetamax (= 0.5 theta) of each other
        sampxgrid, gridoverlap, d2d, d3d = selcoords.search_around_sky(selsamp, thetamax*u.deg)
        
        # Looping over grid points
        for g in range(len(selsamp)):
            
            # Index of the grid point in the flag list
            tot_index = sampindex_tot[s]+g
            
            # If the point is not flagged (not flagged = 1), flag all overlapping grid points (flagged = 0).
            if selected_theta[tot_index] == 1:
                removeindex = gridoverlap[(sampxgrid==g)&(gridoverlap!=tot_index)]
                selected_theta[removeindex] = 0.
            
            # If the point is already flagged, do nothing.
            else:
                pass
    
    # Sort the grid points back to their original order
    sort_ID = np.argsort(selID)

    # The list Stheta contains the selected (= not flagged = 1) grid points for the non-overlapping sample.
    Stheta_tot[theta] = selected_theta[sort_ID]

