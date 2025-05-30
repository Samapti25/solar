import os
from astropy.io import fits
import sunpy.map
from sunpy.coordinates import NorthOffsetFrame
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
import matplotlib.pyplot as plt
import matplotlib.colors
import glob
from scipy import ndimage
import numpy as np
import sunpy.timeseries
import scipy.ndimage as ndimage
from scipy.ndimage import label

folder_path = "/home/samapti-lakshan/carrington"
save_dir= "/home/samapti-lakshan/carrington/pixelimage"
files = os.listdir(folder_path)
fits_files = glob.glob(os.path.join(folder_path, '*.fits')) 
for file in fits_files:
	with fits.open(file)as f:
		header = f[0].header
		print (header)
		data=f[0].data
		fig = plt.figure()
		plt.imshow(data, origin="lower", cmap='gray')
		plt.colorbar()
		
		base_name = os.path.basename(file).replace('.fits', '')
		plt.title(f"{base_name}")
		image_name= f" {base_name}_pixel.jpeg"
		mag_image= os.path.join(save_dir, image_name)
		plt.savefig(mag_image, dpi=300)
		plt.show()
