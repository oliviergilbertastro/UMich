import pandas as pd
from scipy.signal import fftconvolve
from reproject import reproject_interp
import numpy as np
import matplotlib.pyplot as plt
from make_csv_table_ztf import get_image
from astropy.visualization import ImageNormalize, LinearStretch
from astropy.stats import sigma_clipped_stats
from scipy.signal import correlate2d
from scipy.ndimage import shift

def get_vmin_vmax(data):
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    return median - 2*std, median + 10*std
# Img_diff = Img_science - (I_reference * K)
# Where K is the convolution kernel matching PSFs

def pad_to_same(psf1, psf2):
    ny = max(psf1.shape[0], psf2.shape[0])
    nx = max(psf1.shape[1], psf2.shape[1])
    def pad(psf):
        y, x = psf.shape
        pad_y = ny - y
        pad_x = nx - x
        return np.pad(
            psf,
            ((pad_y//2, pad_y - pad_y//2),
             (pad_x//2, pad_x - pad_x//2)),
            mode='constant'
        )
    return pad(psf1), pad(psf2)

def estimate_bkg(image, median=False):
    if median: return np.nanmedian(image)
    vals, bins = np.histogram(np.nan_to_num(image,nan=0), bins=np.linspace(0,np.quantile(np.nan_to_num(image,nan=0),0.8),100))
    max_index = np.argmax(vals)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    bkg_level = bin_centers[max_index]
    return bkg_level

def compute_kernel(psf_source, psf_target, eps=1e-6, smooth=False):
    """Function to get the convolution kernel"""
    F_source = np.fft.fft2(psf_source)
    F_target = np.fft.fft2(psf_target)
    
    # Suppress high-frequency noise
    kx = np.fft.fftfreq(psf_source.shape[1])
    ky = np.fft.fftfreq(psf_source.shape[0])
    kx, ky = np.meshgrid(kx, ky)
    k2 = kx**2 + ky**2

    taper = np.exp(-k2 / 1)  # tune this
    

    kernel_fft = F_target / (F_source + eps) * (taper if smooth else 1)
    kernel = np.fft.ifft2(kernel_fft).real

    return np.fft.fftshift(kernel)

def difference_image(img, psf, ref_img, ref_psf, mask=None, if_plot=False):
    """Function to do difference imaging by convolving both images to the worse PSF before
    difference imaging."""

    img_hdr, ref_hdr = img[0].header, ref_img[0].header
    img_data, ref_data = img[0].data, ref_img[0].data
    print(img_data)
    # Normalize the PSFs
    psf = psf[0].data
    ref_psf = ref_psf[0].data
    psf = psf/np.sum(psf)
    ref_psf = ref_psf/np.sum(ref_psf)

    ref_reproj, footprint = reproject_interp(
        (ref_data, ref_hdr),
        img_hdr
    )
    assert ref_reproj.shape == img_data.shape

    if False:
        # Try to correct for shift:
        corr = correlate2d(img_data, ref_reproj, mode='same')
        y0, x0 = np.unravel_index(np.argmax(corr), corr.shape)
        shift_y = y0 - img_data.shape[0]//2
        shift_x = x0 - img_data.shape[1]//2
        ref_reproj = shift(ref_reproj, shift=(shift_y, shift_x))

    ref_psf_padded, psf_padded = pad_to_same(ref_psf, psf)
    kernel = compute_kernel(ref_psf_padded, psf_padded)
    img_data -= estimate_bkg(img_data, median=True)
    ref_reproj -= estimate_bkg(ref_reproj, median=True)
    ref_reproj_clean = np.nan_to_num(ref_reproj, nan=0.0)
    ref_matched = fftconvolve(ref_reproj_clean, kernel, mode="same")
    ref_matched = fftconvolve(ref_reproj_clean, kernel, mode="same")
    diff = img_data - ref_matched

    if mask is not None:
        ref_reproj[mask] = np.nan
        img_data[mask] = np.nan
        diff[mask] = np.nan
    
    if if_plot:
        plt.figure(figsize=(15,5))
        vmin, vmax = get_vmin_vmax(ref_reproj)
        ax1 = plt.subplot(1,3,1)
        ax1.imshow(ref_reproj, origin="lower", cmap="Greys", vmin=vmin, vmax=vmax)
        ax1.set_title("Reference")

        vmin, vmax = get_vmin_vmax(img_data)
        ax2 = plt.subplot(1,3,2,sharex=ax1,sharey=ax1)
        ax2.imshow(img_data, origin="lower", cmap="Greys", vmin=vmin, vmax=vmax)
        ax2.set_title("Image")

        vmin, vmax = get_vmin_vmax(diff)
        ax3 = plt.subplot(1,3,3,sharex=ax1,sharey=ax1)
        ax3.imshow(diff, origin="lower", cmap="Greys", vmin=vmin, vmax=vmax)#, vmin=-v, vmax=v)
        ax3.set_title("Difference")

        plt.show()
    return diff


if __name__ == "__main__":
    df = pd.read_csv("ASTRO501/seminar/observations_df.csv")

    imgs, psfs = [], []
    for i, filefracday in enumerate(df["filefracday"].values):
        img, psf = get_image(filefracday, field=True)
        if i == 1:
            reference_img, reference_psf = img, psf
        else:
            imgs.append(img)
            psfs.append(psf)

    #mask = np.zeros_like(reference_img, dtype=bool)
    #mask[:,225:275] = np.ones_like(mask[0:460,225:275], dtype=bool)
    #reference_img[mask==1] = np.nan
    mask = np.zeros_like(reference_img[0].data, dtype=bool)
    mask[:, 225:275] = True # First star
    mask[80:260, 75:130] = True # Second star


    diffs = []
    for i in range(len(imgs)):
        #diffs.append(difference_image(imgs[i], psfs[i], reference_img, reference_psf, if_plot=True))
        diffs.append(difference_image(reference_img, reference_psf, imgs[i], psfs[i], if_plot=True, mask=mask))
    plt.figure(figsize=(17,7))
    for i in range(len(diffs)):
        if i==0:
            ax1 = plt.subplot(2,3,i+1)
            ax = ax1
        else:
            ax = plt.subplot(2,3,i+1,sharex=ax1,sharey=ax1)
        vmin, vmax = get_vmin_vmax(diffs[i])
        print(vmin, vmax)
        ax.imshow(diffs[i], origin="lower", cmap="Greys", vmin=vmin, vmax=vmax)
        ax_last = ax
    plt.tight_layout()
    plt.show()
    