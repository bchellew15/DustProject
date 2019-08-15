#read the BOSS wavelength arrays and see what they look like
#and do some processing

#########
# WARNING: THIS CODE IS NOW OUTDATED
# THE DATA FILES IT USED WERE INCORRECT
# USE BOSS_TEST_UPDATE.PY
########


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys #for numpy threshold

hdulist = fits.open('../BrandtFiles/skyfibers_nativelam.fits')
flams = hdulist[0].data
ivars = hdulist[1].data
wavelengths = hdulist[2].data #remember the offset

#print("min: ", min(hdulist[2].data))
#np.set_printoptions(threshold=sys.maxsize)
#hdulist.info()

#take flam and pad with zeros so the columns align
#then save as a npy file
#use memory mapping to avoid filling up RAM
def create_flambda_array(flams):

    np.set_printoptions(threshold=sys.maxsize)
    
    #correct wavelength array
    lambda_0 = np.ones(flams.shape[0])
    lambda_0[2:] = hdulist[2].data
    lambda_0[:2] = lambda_0[2]
    
    #create empty array for padded lambdas
    num_intervals = int(round((max(hdulist[2].data) - min(hdulist[2].data)) / 0.0001))
    num_flams = num_intervals + flams[0].shape[0]
    padded_flams_mmap = np.memmap('padded_flams_memmap.dat', dtype=np.float32, mode='w+', shape=(flams.shape[0], num_flams))

    #before:
    #print(padded_flams_mmap[10])
    
    #fill up the array:
    print("filling array")
    for i in range(flams.shape[0]):
        #calculate indent
        num_spaces = int(round((lambda_0[i] - min(lambda_0)) / 0.0001))
        padded_flams_mmap[i, num_spaces:num_spaces+flams[0].shape[0]] = flams[i]

        #progress:
        if i % 1000 == 0:
            print("progress:", i)

    #after:
    #print(padded_flams_mmap[10])

    print("array is full")
        
    #convert to normal array:
    padded_flams = np.zeros(padded_flams_mmap.shape)
    padded_flams[:] = padded_flams_mmap[:]

    #save:
    np.save("/Volumes/TOSHIBA EXT/Dust_Overflow/padded_flams_boss.npy", padded_flams)

#create_flambda_array(flams)
#exit(0)

#create padded ivar array
def create_ivar_array(ivars):

    np.set_printoptions(threshold=sys.maxsize)
    
    #correct wavelength array
    lambda_0 = np.ones(ivars.shape[0])
    lambda_0[2:] = hdulist[2].data
    lambda_0[:2] = lambda_0[2]

    #create empty array for padded ivars
    num_intervals = int(round((max(hdulist[2].data) - min(hdulist[2].data)) / 0.0001))
    num_ivars = num_intervals + ivars[0].shape[0]
    padded_ivars_mmap = np.memmap('padded_ivars_memmap.dat', dtype=np.float32, mode='w+', shape=(ivars.shape[0], num_ivars))
    #fill with ones first:
    padded_ivars_mmap[:] = 1
    
    #fill up the array:
    print("filling array")
    for i in range(ivars.shape[0]):
        #calculate indent
        num_spaces = int(round((lambda_0[i] - min(lambda_0)) / 0.0001))
        padded_ivars_mmap[i, num_spaces:num_spaces+ivars[0].shape[0]] = ivars[i]

        #progress:
        if i % 1000 == 0:
            print("progress:", i)

    print("array is full")
    
    #convert to normal array:
    padded_ivars = np.zeros(padded_ivars_mmap.shape)
    padded_ivars[:] = padded_ivars_mmap[:]

    #save:
    np.save("/Volumes/TOSHIBA EXT/Dust_Overflow/padded_ivars_boss.npy", padded_ivars)
    
#create_ivar_array(ivars)
#exit(0)

def load_padded_arrs(flams, ivars, wavelengths):
    np.set_printoptions(threshold=sys.maxsize)

    print("loading arrays")
    
    ivar = np.load("/Volumes/TOSHIBA EXT/Dust_Overflow/padded_ivars_boss.npy", mmap_mode='r')

    print("ivar row 20 (padded)")
    print(ivar[20])
    print("ivar row 20 (original)")
    print(ivars[20])
    print("number of spaces")
    print(wavelengths[18]) #because of missing values
    print("")
    
    print("ivar row 100 (padded)")
    print(ivar[100])
    print("ivar row 100 (original)")
    print(ivars[100])
    print("number of spaces")
    print(wavelengths[98]) #because of missing values
    print("")

    flam = np.load("/Volumes/TOSHIBA EXT/Dust_Overflow/padded_flams_boss.npy", mmap_mode='r')

    print("flam row 20 (padded)")
    print(flam[20])
    print("flam row 20 (original)")
    print(flams[20])
    print("number of spaces")
    print(wavelengths[18]) #because of missing values
    print("")

    print("flam row 100 (padded)")
    print(flam[100])
    print("flam row 20 (original)")
    print(flams[100])
    print("number of spaces")
    print(wavelengths[100]) #because of missing values
    print("")

#load_padded_arrs(flams, ivars, wavelengths)
#exit(0)

#check whether fibers on same plate start with same lambda:

#lambda_0 = np.ones(flams.shape[0])
#lambda_0[2:] = hdulist[2].data
#lambda_0[:2] = lambda_0[2]

lambda_0 = hdulist[2].data
print(lambda_0.shape)

fiberinfo = np.genfromtxt('../BrandtFiles/sky_radec.dat')
plate = fiberinfo[:,0][:-2]
mjd = fiberinfo[:,1][:-2]
plate = 10000*plate + mjd%10000 #plate is now a unique identifier
print(plate.shape)
print(mjd.shape)

#for each unique plates, get all the lambdas, see if any don't match
unique_plates = np.unique(plate, return_index=True)[1] #2512 unique plates
for i in unique_plates:
    p = fiberinfo[:,0][i]
    d = mjd[i]
    l = lambda_0[(fiberinfo[:,0][:-2] == p) * (mjd == d)] #all the starting lambdas for the plate
    if len(np.unique(l)) != 1:
        print("fail")
        print(l)

exit(0)



#histograms of SFD i100 values for each set of locations:
'''
fiberinfo = np.loadtxt('/Users/blakechellew/Documents/DustProject/BrandtFiles/fiberinfo_halpha.dat')
i100_old = np.array(fiberinfo[:, 4]).astype('float32')  # 100 micron intensity (MJy/sr)
plt.hist(i100_old, bins=50, range=(0,10))
plt.xlabel("i100 value (SDSS)")
plt.ylabel("count")
plt.show()

i100_old_boss = np.loadtxt("/Users/blakechellew/Documents/DustProject/SFD_Maps/CodeC/SFD_i100_at_BOSS_locations.txt")[:,2].astype('float32')
plt.hist(i100_old_boss, bins=50, range=(0,10))
plt.xlabel("i100 value (BOSS)")
plt.ylabel("count")
plt.show()
'''

#histograms of ivar:
'''
ivar = np.load("/Volumes/TOSHIBA EXT/Dust_Overflow/padded_ivars_boss.npy", mmap_mode='r') #type: float64
print("check 1")
plt.hist(ivar[:,1037:1132], bins=50, range=(0, 30))
plt.show()

hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
ivar_boss = hdulist[3].data[:]
plt.hist(ivar[:,765:864], bins=50, range=(0, 30))
plt.show()
'''

#histograms of flambda:
flam = np.load("/Volumes/TOSHIBA EXT/Dust_Overflow/padded_flams_boss.npy", mmap_mode='r') #type: float64
plt.hist(flam[:,1037:1132], bins=50, range=(-3, 4))
plt.show()

hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
flam_boss = hdulist[2].data
plt.hist(flam_boss[:,765:864], bins=50, range=(-3, 4))
plt.show()



exit(0)


#figure out why the wavelength array has two less rows than flambda
#update: no need, brandt figured it out
'''
def wavelength_missing_elements():
    lams = hdulist[2].data
    print(lams.shape)
    fiberinfo = np.genfromtxt('./BrandtFiles/sky_radec.dat')
    plates = fiberinfo[:,0]
    print(plates.shape)

wavelength_missing_elements()
'''

    
# Flux densities Flam:
#flams = hdulist[0].data
# Inverse variances
#ivars = hdulist[1].data
lams = hdulist[2].data


flam = hdulist[0].data[0]
ivar = hdulist[1].data[0]
lam = (hdulist[2].data[0] + np.arange(flam.shape[0])*1e-4) #10**this to get actual lambdas)

print("length of arrays:", flam.shape[0])

print("flam shape:", flam.shape)
print("ivar shape:", ivar.shape)
#print(lam.shape)

print(lams)

#find out how many unique starting lambdas and how many of each
unique_lambdas, unique_counts = np.unique(lams, return_counts=True)
print("lambdas:")
print(unique_lambdas)
print("counts:")
print(unique_counts)

#for l in unique_lambdas:
#    plt.axvline(x=l)
#plt.show()

#look at the first lambda array, before converting away from log
print(lam)

#check skyfibers.dat to see if 1st column is plate number
plate_nums = np.genfromtxt('./BrandtFiles/skyfibers.dat', usecols=(2))
unique_plates = np.unique(plate_nums)
print("number of plates:")
print(len(unique_plates))
    

