import os
import astropy.io.fits as pf

# Writes cross-correlation function and Gaussian parameters to a .fits file.
def CrossCorWriteFits(header,ccfs,gaus,gauspath):
    os.chdir(gauspath)

    output_file_name=header['file']
    empty_primary=pf.PrimaryHDU(header=header)
    
    ccfs_hdu                  =pf.ImageHDU(ccfs,name='CCFS' )
    ccfs_hdu.header['CONTENT']='CCFs'
    
    gaus_hdu                  =pf.ImageHDU(gaus,name='FITS')
    gaus_hdu.header['CONTENT']='Fit_to_CCFS'
    
    hdu=pf.HDUList([empty_primary,ccfs_hdu,gaus_hdu])
        
    hdu.writeto(output_file_name.split('.fits')[0]+'_gaussparams.fits',overwrite=True)

# Writes normalized values and cross-correlation function to a .fits file.
def NormalizationWriteFits(header,matrices,normpath,ccfspath):
    os.chdir(normpath)

    output_file_name=header['file']
    empty_primary=pf.PrimaryHDU(header=header)
    
    flux_hdu                  =pf.ImageHDU(matrices[0],name='Normalized Flux' )
    flux_hdu.header['CONTENT']='Normalized Flux'
    
    lam_hdu                  =pf.ImageHDU(matrices[1],name='Interpolated Wavelength')
    lam_hdu.header['CONTENT']='Interpolated Wavelength'

    template_hdu                  =pf.ImageHDU(matrices[2],name='Interpolated Template flux')
    template_hdu.header['CONTENT']='Interpolated Template flux'
    
    hdul=pf.HDUList([empty_primary,flux_hdu,lam_hdu,template_hdu])
        
    hdul.writeto(output_file_name.split('.fits')[0]+'_norm.fits',overwrite=True)
    
    os.chdir(ccfspath)
    
    output_file_name=header['file']
    empty_primary=pf.PrimaryHDU(header=header)
    
    ccfs_hdu                  =pf.ImageHDU(matrices[3],name='CCFS' )
    ccfs_hdu.header['CONTENT']='CCFs'
    
    hdu2=pf.HDUList([empty_primary,ccfs_hdu])
    
    hdu2.writeto(output_file_name.split('.fits')[0]+'_ccfs.fits',overwrite=True) 
