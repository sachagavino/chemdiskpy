import numpy as np
from ..constants.constants import c

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'same') / w

def dust_temperature(temp_radmc3d, rchem, zchem, d, theta, hg):
    nbspecies, nph, nt, nr = temp_radmc3d.shape
    print('temp shape coupl ', nbspecies, nr, nt, nph)
    hhg, zz = np.meshgrid(hg, zchem, indexing='ij')
    zz = hhg*zz
    temp_naut = np.ones((nbspecies, len(rchem), len(zchem)))
    temp_naut_smooth = np.ones((nbspecies, len(rchem), len(zchem)))
    for size_id in range(nbspecies):
        for idx, r in enumerate(rchem):
            for alt in range(len(zchem)):
                d_pt = np.sqrt(r**2 + zz[idx, alt]**2)  #convert from cartesian to spherical
                theta_pt = np.arccos(zz[idx, alt]/d_pt) #convert from cartesian to spherical
                closest_d = min(enumerate(d), key=lambda x: abs(x[1]-d_pt)) #find closest grid point
                closest_t = min(enumerate(theta), key=lambda x: abs(x[1]-theta_pt)) #find closest grid point
                #temp_naut[size_id, idx, alt] = temp_radmc3d[size_id, closest_d[0], closest_t[0], 0] 
                temp_naut[size_id, idx, alt] = temp_radmc3d[size_id, 0, closest_t[0], closest_d[0]] 

    #SMOOTHING TEMPERATURE PROFILE
    for idx in range(len(rchem)):
        for size_id in range(nbspecies):
            temp_naut_smooth[size_id, idx, :] = moving_average(temp_naut[size_id, idx, :], 5) #average the values over a rolling window of 5 points.
    temp_naut_smooth[:, :, 0:2] = temp_naut[:, :, 0:2]  #clean the boundary effects
    temp_naut_smooth[:, :, -2:] = temp_naut[:, :, -2:]    

    return temp_naut_smooth  

def local_field():
    pass

def av_z(field_radmc3d, lam, rchem, zchem, d, theta, hg):

    nlam, nph, nt, nr = field_radmc3d.shape
    lamuv = np.where((lam >= 0.5) & (lam <= 0.6)) # extract the ~ visible
    freq = c/(lam*1e-6)
    fieldint = np.zeros((nt, nr))

    # Integrate over uv frequencies:
    for i in lamuv[0]:
        fieldint += field_radmc3d[i, 0, :, :]*freq[i]
    fieldint[fieldint==0.0] = 1e-30 #provide arbitrary low values in order to avoid division by zero.
    # Convert from spherical to nautilus grid:
    hhg, zz = np.meshgrid(hg, zchem, indexing='ij')
    zz = hhg*zz
    field_naut, field_naut_smooth = np.ones((len(rchem), len(zchem))), np.ones((len(rchem), len(zchem)))
    avz = np.ones((len(rchem), len(zchem)))

    for idx, r in enumerate(rchem):
        for z in range(len(zchem)):
            d_pt = np.sqrt(r**2 + zz[idx, z]**2)  #convert from cartesian to spherical
            theta_pt = np.arccos(zz[idx, z]/d_pt) #convert from cartesian to spherical
            closest_d = min(enumerate(d), key=lambda x: abs(x[1]-d_pt)) #find closest grid point
            closest_t = min(enumerate(theta), key=lambda x: abs(x[1]-theta_pt)) #find closest grid point
            field_naut[idx, z] = fieldint[closest_t[0], closest_d[0]]
    # Smoothing vertical profiles
    for idx in range(len(rchem)):
        field_naut_smooth[idx, :] = moving_average(field_naut[idx, :], 5) #average the values over a rolling window of 5 points.
    field_naut_smooth[:, 0:2] = field_naut[:, 0:2]  #clean the boundary effects
    field_naut_smooth[:, -2:] = field_naut[:, -2:]   

    for idx in range(len(rchem)):
        avz[idx, :] = abs(-1.086*np.log(field_naut[idx,:]/field_naut[idx, 0]))   #field0[idx]))
        avz[idx,1:][avz[idx,1:]==0.0] = np.trim_zeros(avz[idx])[0]
    
    # avz_df = pd.DataFrame(data=avz.transpose())
    # avz_df = avz_df.rolling(window=5, center=True, min_periods=2).mean()
    return avz