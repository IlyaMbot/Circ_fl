# Universal constant
c = 3e10          #[cm/s] speed of light
kb = 1.38*10.^(-16)  # [cgs] Bolzmann constant
cgs2sfu = 10 ** 19 # SFU/[erg/s/cm^2/Hz]
cgs2jy  = 10 ** 23 # Jy/[erg/s/cm^2/Hz]
flux_factor = cgs2sfu
pixsz = index[0].CDELT1

omega_per_pix = (pixsz / 3600.) ** 2 # in radians

ff = index.OBS_D$FREQ  # from fits file 5600 GHz

fi = 2 * kb * (ff * 10 ** 9) ** 2 /c ** 2 * omega_per_pix * flux_factor # coefficient 

spec0 = (total(srh0[ind0])) * fi # result in sfu