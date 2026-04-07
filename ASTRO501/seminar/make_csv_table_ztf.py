import pandas as pd


df = pd.read_csv("ASTRO501/seminar/table_Science-Exposure.csv")
obslist = [2458963.7175579001,2458936.7646181001,2458900.8169328999,2459001.6712731002,2459037.7144559999,2458838.0099769002]
df_selected_images = df[df["obsjd"].isin(obslist)]
df_selected_images["mjd"] = df_selected_images["obsjd"]- 2400000.5
print(df_selected_images)

from astropy.io import fits
def get_image(filefracday):
    filefracday = str(filefracday)
    path = f"ASTRO501/seminar/sci_mrk421/{filefracday[:4]}/{filefracday[4:8]}/{filefracday[8:]}/ztf_{filefracday}_000713_zg_c03_o_q1_sciimg_ra166.1138_dec38.2088_asec29.fits"
    img = fits.open(path)
    psf_path = f"ASTRO501/seminar/sci_mrk421/{filefracday[:4]}/{filefracday[4:8]}/{filefracday[8:]}/ztf_{filefracday}_000713_zg_c03_o_q1_sciimgdaopsfcent.fits"
    psf = fits.open(psf_path)
    return img, psf

import matplotlib.pyplot as plt
import numpy as np



def net_flux(image, target_pos, noise_level, radius=40):
    '''
    Calculates the approximate flux of a source by substracting the noise of its nearby background

    Parameters:
    image: image data with the source (2D array)
    target_pos: list of position of the source (e.g. [100, 200])
    radius: approximate radius of the source object
    
    Returns the estimated flux of the source
    '''
    x0, y0 = target_pos[0], target_pos[1]
    pixels_inside = 0
    within_flux = 0
    for y, line in enumerate(image):
        for x, val in enumerate(line):
            if np.sqrt((x-x0)**2+(y-y0)**2)<=radius:
                within_flux+=val
                pixels_inside += 1

    flux_within = within_flux
    flux_source = flux_within - (noise_level)*pixels_inside
    if True:
        fig, ax1 = plt.subplots(1, 1)
        #fig.suptitle(f'Net flux of source: {flux_source}')

        circle3 = plt.Circle(tuple(target_pos), radius, color='red', fill=False)
        #ax1.set_title(label=f'Raw flux of source: {flux_within}\nNumber of pixels: {pixels_inside}')
        ax1.imshow(image, origin='lower', cmap="Greys")
        ax1.add_patch(circle3)
        plt.show()
    return flux_source

counts = []
for filefracday in df_selected_images["filefracday"].values:
    img, psf = get_image(filefracday)
    print(img.info())
    #plt.imshow(img[0].data, origin="lower", cmap="Greys")
    plt.figure()
    vals, bins = np.histogram(np.nan_to_num(img[0].data,nan=0), bins=np.linspace(0,np.quantile(img[0].data,0.8),100))
    max_index = np.argmax(vals)
    # Compute bin centers
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    # Get x-value corresponding to max bin
    noise_val = bin_centers[max_index]
    plt.stairs(vals, bins, fill=True, color="black")
    plt.axvline(noise_val, ls="--", color="red")
    print(noise_val)
    counts.append(net_flux(img[0].data, target_pos=(15,15), radius=12, noise_level=noise_val))
    plt.show()

plt.plot(df_selected_images["mjd"], counts, marker="o", ls="None")
plt.show()
plt.plot(df_selected_images["mjd"], df_selected_images["seeing"], marker="o", ls="None")
plt.show()