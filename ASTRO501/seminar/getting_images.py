from pulling_images import *

get_ztf_cutouts_from_zubercal("ASTRO501/seminar/153851661137291254.csv",
                              output_dir="ASTRO501/seminar/ztf_cutouts_153851661137291254/",
                              ps1id=153851661137291254,
                              radius_arcsec=2,
                              mjd_tolerance=0.01,
                              max_images_per_point=1)