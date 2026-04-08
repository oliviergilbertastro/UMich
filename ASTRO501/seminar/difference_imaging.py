import pandas as pd
from scipy.signal import fftconvolve
from reproject import reproject_interp
import numpy as np
import matplotlib.pyplot as plt
from make_csv_table_ztf import get_image

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

def estimate_bkg(image):
    vals, bins = np.histogram(np.nan_to_num(image,nan=0), bins=np.linspace(0,np.quantile(np.nan_to_num(image,nan=0),0.8),100))
    max_index = np.argmax(vals)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    bkg_level = bin_centers[max_index]
    return bkg_level

def compute_kernel(psf_source, psf_target, eps=1e-6):
    """Function to get the convolution kernel"""
    F_source = np.fft.fft2(psf_source)
    F_target = np.fft.fft2(psf_target)
    

    kernel_fft = F_target / (F_source + eps)
    kernel = np.fft.ifft2(kernel_fft).real

    return np.fft.fftshift(kernel)

def difference_image(img, psf, ref_img, ref_psf, if_plot=False):
    """Function to do difference imaging by convolving both images to the worse PSF before
    difference imaging."""

    img_hdr, ref_hdr = img[0].header, ref_img[0].header
    img_data, ref_data = img[0].data, ref_img[0].data

    # Normalize the PSFs
    psf = psf[0].data/np.sum(psf[0].data)
    ref_psf = ref_psf[0].data/np.sum(ref_psf[0].data)

    ref_reproj, footprint = reproject_interp(
        (ref_data, ref_hdr),
        img_hdr
    )
    assert ref_reproj.shape == img_data.shape

    ref_psf_padded, psf_padded = pad_to_same(ref_psf, psf)
    kernel = compute_kernel(ref_psf_padded, psf_padded)
    img_data -= estimate_bkg(img_data)
    ref_reproj -= estimate_bkg(ref_reproj)
    ref_reproj_clean = np.nan_to_num(ref_reproj, nan=0.0)
    ref_matched = fftconvolve(ref_reproj_clean, kernel, mode="same")
    ref_matched = fftconvolve(ref_reproj_clean, kernel, mode="same")
    diff = img_data - ref_matched

    v = 3 * np.std(diff)
    if if_plot:
        plt.figure(figsize=(15,5))

        plt.subplot(1,3,1)
        plt.imshow(ref_data, origin="lower", cmap="gray")
        plt.title("Reference")

        plt.subplot(1,3,2)
        plt.imshow(img_data, origin="lower", cmap="gray")
        plt.title("Image")

        plt.subplot(1,3,3)
        plt.imshow(diff, origin="lower", cmap="gray")#, vmin=-v, vmax=v)
        plt.title("Difference")

        plt.show()
    return diff


if __name__ == "__main__":
    df = pd.read_csv("ASTRO501/seminar/observations_df.csv")

    imgs, psfs = [], []
    for i, filefracday in enumerate(df["filefracday"].values):
        img, psf = get_image(filefracday)
        if i == 1:
            reference_img, reference_psf = img, psf
        else:
            imgs.append(img)
            psfs.append(psf)

    diffs = []
    for i in range(len(imgs)):
        diffs.append(difference_image(imgs[i], psfs[i], reference_img, reference_psf, if_plot=True))
    plt.figure(figsize=(25,6))
    for i in range(len(diffs)):
        ax = plt.subplot(1,len(diffs),i+1)
        ax.imshow(diffs[i], origin="lower", cmap="gray")
    plt.show()
    