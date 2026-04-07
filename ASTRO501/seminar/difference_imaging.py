from make_csv_table_ztf import *

df = pd.read_csv("ASTRO501/seminar/observations_df.csv")


imgs, psfs = [], []
for i, filefracday in enumerate(df["filefracday"].values):
    img, psf = get_image(filefracday)
    if i == 0:
        reference_img, reference_psf = img, psf
    else:
        imgs.append(img)
        psfs.append(psf)

# Img_diff = Img_science - (I_reference * K)
# Where K is the convolution kernel matching PSFs
def compute_kernel(psf_source, psf_target, eps=1e-6):
    """Function to get the convolution kernel"""
    F_source = np.fft.fft2(psf_source)
    F_target = np.fft.fft2(psf_target)

    kernel_fft = F_target / (F_source + eps)
    kernel = np.fft.ifft2(kernel_fft).real

    return np.fft.fftshift(kernel)

def difference_image(img, psf, ref_img, ref_psf):
    """Function to do difference imaging by convolving both images to the worse PSF before
    difference imaging."""
    # Normalize the PSFs
    psf = psf/np.sum(psf)
    ref_psf = ref_psf/np.sum(ref_psf)



science = fits.getdata("science_image.fits")
reference = fits.getdata("reference_image.fits")

psf_sci = fits.getdata("psf_science.fits")
psf_ref = fits.getdata("psf_reference.fits")