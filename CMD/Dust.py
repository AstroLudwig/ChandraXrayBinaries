# -*- Copyright (c) 2020, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
	Dust Extinction Correction
PURPOSE:
	
Notes: 
	
"""
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import timeit
from astropy.io import fits
from astropy.wcs import WCS
#from astropy.wcs.utils import skycoord_to_pixel

class DustCorrection:
    def __init__(
        self,
        catalog_file,
        smc_hotav_file="Data/smc_hotav.fits",
        lmc_hotav_file="Data/lmc_hotav.fits",
        smc_alam_av_file="Data/SMC_Alambda_Av_Table.csv",
        lmc_alam_av_file="Data/LMC_Alambda_Av_Table.csv",
    ):

        print("Ignoring Milky Way values from IRSA table for now, check that this is ok.")
        
        print("Using Av/Alambda table generated in Foreground/DustRemoval.ipynb.")
        
        # Open Catalog  

        self.cat = pd.read_csv(catalog_file)
        
        # Make headers lower case
        
        self.cat.columns = map(str.lower,self.cat.columns)

        # Get Coordinates

        self.ra = self.cat.ra

        self.dec = self.cat.dec


        print("Assuming Coordinates are in degrees. No frame assumed.")

        self.coordinates = SkyCoord(self.ra,self.dec,unit=u.deg)    

        self.galaxy = self.which_galaxy(self.coordinates[0])

        print(f"Assuming all sources in catalog are within the {self.galaxy} galaxy.")   
    
        if self.galaxy == 'lmc':
            
            self.hotav_file = lmc_hotav_file
            
            self.alam_av_file = lmc_alam_av_file
            
        else:
            
            self.hotav_file = smc_hotav_file
            
            self.alam_av_file = smc_alam_av_file
            
            
        print("Assuming WCS orgin = 0.")
        
        self.hotav_x, self.hotav_y = self.coordinates.to_pixel(WCS(self.hotav_file),origin=0)
        
        # Make pixel coordinates indices.
        
        self.hotav_x = np.round(self.hotav_x).astype(int)
        
        self.hotav_y = np.round(self.hotav_y).astype(int)
        
        # Get Av values
        
        self.hotav = fits.open(self.hotav_file)[0].data[self.hotav_y,self.hotav_x]
        
        # Get DeRed columns
        
        self.alam_av = pd.read_csv(self.alam_av_file).set_index('Filter')
        
        self.alam_av.columns = ['Wav','Alam_Av']
        
        self.cat["dered_uvw2"] = self.cat.uvw2_mag - self.hotav * self.alam_av.loc['UVW2']['Alam_Av']
        self.cat["dered_uvm2"] = self.cat.uvm2_mag - self.hotav * self.alam_av.loc['UVM2']['Alam_Av']
        self.cat["dered_uvw1"] = self.cat.uvw1_mag - self.hotav * self.alam_av.loc['UVW1']['Alam_Av']
        self.cat["dered_u"] = self.cat.umag - self.hotav * self.alam_av.loc['U']['Alam_Av']
        self.cat["dered_b"] = self.cat.bmag - self.hotav * self.alam_av.loc['B']['Alam_Av']
        self.cat["dered_v"] = self.cat.vmag - self.hotav * self.alam_av.loc['V']['Alam_Av']
        self.cat["dered_i"] = self.cat.imag - self.hotav * self.alam_av.loc['I']['Alam_Av']
        
        print("Did not apply dereddening to errors. Check this. ")
        
        # Save New Catalog 
        
        save_name = catalog_file.split(".")[0]
        
        self.cat.to_csv(f"{save_name}_DeRed.csv",index=False)
    
    def which_galaxy(self,co):
        
        # Nearby galaxy coordinates, from simbad.
        
        smc = SkyCoord(13.15833333, -72.80027778,unit=u.deg).separation(co)
        
        lmc = SkyCoord(80.89416667, -69.75611111,unit=u.deg).separation(co)
        
        if lmc < smc :
            
            return 'lmc'
        
        return 'smc'
        
        
#############
# Graveyard #
#############        

#         self.MilkyWayTable = IrsaDust.get_extinction_table(self.Coordinate)
#         self.MilkyWayExtinction = pd.DataFrame(
#             {
#                 "Filter": self.MilkyWayTable["Filter_name"],
#                 "E_BV": (
#                     self.MilkyWayTable["A_over_E_B_V_SandF"]
#                     / self.MilkyWayTable["A_SandF"]
#                 )
#                 ** (-1),
#                 "A_lam": self.MilkyWayTable["A_SandF"],
#             },
#             columns=["Filter", "E_BV", "A_lam"],
#         )

