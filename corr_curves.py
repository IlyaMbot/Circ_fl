from astropy.io import fits
import numpy as np
import os, glob
import antlib

def get_data_for_freq(filenames : str, frequency : int = 0) -> [np.ndarray, np.ndarray, np.int32]:
    '''
    Get antennas' corr R+L,
    time of measurement in seconds (UT) from FITS files.

    Parameters
    ----------
    filename : str
        Path to the FITS files.
    frequency : int, optional
        Index of the frequency. The default is 0.
    '''

    amps = np.array([])
    time = np.array([])

    for filename in filenames:
        print(f"{filename} is opened")
        with fits.open( filename, memmap = False ) as f:

            f.verify('silentfix')

            amp     = f[1].data[ frequency ][ "vis_RCP" ] + f[1].data[0][ "vis_LCP" ]
            freq    = f[1].data[ frequency ][ "FREQUENCY" ]
            time20  = f[1].data[ frequency ][ "TIME" ]

        time = np.append(time, time20)
        amps = np.append(amps, np.reshape( amp, (20, 8128) ) )

    if amps.shape[0] == 0:
        print("ERROR")
        return(None, None, None)
    
    amps = np.reshape( amps, (len(filenames) * 20, 8128 ) ).swapaxes(0, 1)
    time = np.reshape( time, (len(filenames) * 20 ) )

    return(amps, time, freq)


filenames = glob.glob('./event*/srh_20210522T06*.fit')
filenames = sorted(filenames, key=os.path.basename)

ants, time, _ = get_data_for_freq(filenames)

ants = np.array(np.real(ants))
time = np.array(time)

ants = np.append( np.array([time]), ants, axis = 0)

antlib.save_fits_raw(ants, '2021_05_22_T06', "fluxes")

