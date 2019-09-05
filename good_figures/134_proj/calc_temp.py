#use images to calculate temperature
#given png screenshots taken from Brandt's paper

#also: code to help generage list of alphas given png of unbinned spectrum

import matplotlib.pyplot as plt
import numpy as np


from astropy.io import fits

alpha_img = plt.imread("unbinned_spectrum_screenshot.png")
plt.imshow(alpha_img)
plt.show()


#0 corresponds to alpha = 0.7
#949 corresponds to alpha = 0.05
pixel_values = np.array([745, 755, 784, 808, 878, 787, 757, 746, 774, 769, 722, 727, 727, 674, \
                         636, 666, 728, 774, 819, 823, 761, 768, 611, 236, 125, 300, \
                         558, 767, 744, 779, 776, 786, 739, 737, 765, 754, 670, 531, \
                         375, 438, 657, 721, 765, 777, 772, 778, 764, 770, 769, 727, \
                         764, 771, 762, 723, 739, 768, 739, 732, 750, 743, 741, 759, \
                         760, 764, 775, 770, 783, 764, 759, 734, 779, 751, 768, 772, \
                         778, 756, 750, 760, 755, 762, 765, 739, 781, 804, 781, 780, \
                         764, 750, 754, 766, 751, 749, 787, 762, 743, 751, 826, 771, \
                         760, 751, 766, 749, 760, 763, 762, 750, 754, 734, 725, 763, \
                         790, 770, 760, 766, 758, 776, 756, 727, 733, 778, 776, 752, \
                         762, 737, 761, 747, 741, 632, 472, 432, 549, 703, 767, 744, \
                         739, 755, 726, 633, 526, 520, 632, 732, 760, 749, 735, 770, \
                         752, 783, 761, 726, 717, 749, 759, 776, 779, 752, 744, 770, \
                         763, 769, 789, 755, 735, 751, 742])
alpha_values = pixel_values * (-.65/949) + 0.7
print(alpha_values)

hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
wavelength = np.array(hdulist[1].data)
plt.plot(wavelength[:len(alpha_values)], alpha_values, drawstyle='steps')
plt.show()


exit(0)



n2 = plt.imread("N2_img.png")
s2 = plt.imread("S2_img.png")
#x: 6500 to 8500 (2000 total, center 7500)
x_center_0 = 7500
#y: .52 to .86 (.34 total, center .69)
y_center_0 = .69

print(n2.shape)
print(s2.shape)
#482 by 890

#x: center 445, width 890
x_center_f = 445
#y: center 241, width 482
y_center_f = 241
#x difference = 7055
#y difference = 234.1


#original N
low_horiz = 0.60 
high_horiz = 0.74
left_vert = 6980
right_vert = 8250
#convert
low_horiz = 482 - (y_center_f + (low_horiz - y_center_0)*(482/.34))
high_horiz = 482 - (y_center_f + (high_horiz - y_center_0)*(482/.34)) #subtract because y-axis flipped in imshow
left_vert = x_center_f + (left_vert - x_center_0)*(890/2000)
right_vert = x_center_f + (right_vert - x_center_0)*(890/2000)
#show
plt.imshow(n2)
plt.axhline(y=low_horiz)
plt.axhline(y=high_horiz)
plt.axvline(x=left_vert)
plt.axvline(x=right_vert)
plt.show()

#tao N
low_horiz = 0.63 
high_horiz = 0.76
left_vert = 7080
right_vert = 8320
#convert
low_horiz = 482 - (y_center_f + (low_horiz - y_center_0)*(482/.34))
high_horiz = 482 - (y_center_f + (high_horiz - y_center_0)*(482/.34)) #subtract because y-axis flipped in imshow
left_vert = x_center_f + (left_vert - x_center_0)*(890/2000)
right_vert = x_center_f + (right_vert - x_center_0)*(890/2000)
#show
plt.imshow(n2)
plt.axhline(y=low_horiz)
plt.axhline(y=high_horiz)
plt.axvline(x=left_vert)
plt.axvline(x=right_vert)
plt.show()

#iris N
low_horiz = 0.62
high_horiz = 0.76
left_vert = 7040
right_vert = 8320
#convert
low_horiz = 482 - (y_center_f + (low_horiz - y_center_0)*(482/.34))
high_horiz = 482 - (y_center_f + (high_horiz - y_center_0)*(482/.34)) #subtract because y-axis flipped in imshow
left_vert = x_center_f + (left_vert - x_center_0)*(890/2000)
right_vert = x_center_f + (right_vert - x_center_0)*(890/2000)
#show
plt.imshow(n2)
plt.axhline(y=low_horiz)
plt.axhline(y=high_horiz)
plt.axvline(x=left_vert)
plt.axvline(x=right_vert)
plt.show()

#both N
low_horiz = 0.66 
high_horiz = 0.79
left_vert = 7170
right_vert = 8420
#convert
low_horiz = 482 - (y_center_f + (low_horiz - y_center_0)*(482/.34))
high_horiz = 482 - (y_center_f + (high_horiz - y_center_0)*(482/.34)) #subtract because y-axis flipped in imshow
left_vert = x_center_f + (left_vert - x_center_0)*(890/2000)
right_vert = x_center_f + (right_vert - x_center_0)*(890/2000)
#show
plt.imshow(n2)
plt.axhline(y=low_horiz)
plt.axhline(y=high_horiz)
plt.axvline(x=left_vert)
plt.axvline(x=right_vert)
plt.show()


#original S
low_horiz = 0.68 
high_horiz = 0.82
left_vert = 7280
right_vert = 8540
#convert
low_horiz = 482 - (y_center_f + (low_horiz - y_center_0)*(482/.34))
high_horiz = 482 - (y_center_f + (high_horiz - y_center_0)*(482/.34)) #subtract because y-axis flipped in imshow
left_vert = x_center_f + (left_vert - x_center_0)*(890/2000)
right_vert = x_center_f + (right_vert - x_center_0)*(890/2000)
#show
plt.imshow(s2)
plt.axhline(y=low_horiz)
plt.axhline(y=high_horiz)
plt.axvline(x=left_vert)
plt.axvline(x=right_vert)
plt.show()

#tao S
low_horiz = 0.68 
high_horiz = 0.82
left_vert = 7280
right_vert = 8540
#convert
low_horiz = 482 - (y_center_f + (low_horiz - y_center_0)*(482/.34))
high_horiz = 482 - (y_center_f + (high_horiz - y_center_0)*(482/.34)) #subtract because y-axis flipped in imshow
left_vert = x_center_f + (left_vert - x_center_0)*(890/2000)
right_vert = x_center_f + (right_vert - x_center_0)*(890/2000)
#show
plt.imshow(s2)
plt.axhline(y=low_horiz)
plt.axhline(y=high_horiz)
plt.axvline(x=left_vert)
plt.axvline(x=right_vert)
plt.show()

#iris S
low_horiz = 0.69
high_horiz = 0.83
left_vert = 7300
right_vert = 8560
#convert
low_horiz = 482 - (y_center_f + (low_horiz - y_center_0)*(482/.34))
high_horiz = 482 - (y_center_f + (high_horiz - y_center_0)*(482/.34)) #subtract because y-axis flipped in imshow
left_vert = x_center_f + (left_vert - x_center_0)*(890/2000)
right_vert = x_center_f + (right_vert - x_center_0)*(890/2000)
#show
plt.imshow(s2)
plt.axhline(y=low_horiz)
plt.axhline(y=high_horiz)
plt.axvline(x=left_vert)
plt.axvline(x=right_vert)
plt.show()

#both S
low_horiz = 0.69 
high_horiz = 0.83
left_vert = 7300
right_vert = 8560
#convert
low_horiz = 482 - (y_center_f + (low_horiz - y_center_0)*(482/.34))
high_horiz = 482 - (y_center_f + (high_horiz - y_center_0)*(482/.34)) #subtract because y-axis flipped in imshow
left_vert = x_center_f + (left_vert - x_center_0)*(890/2000)
right_vert = x_center_f + (right_vert - x_center_0)*(890/2000)
#show
plt.imshow(s2)
plt.axhline(y=low_horiz)
plt.axhline(y=high_horiz)
plt.axvline(x=left_vert)
plt.axvline(x=right_vert)
plt.show()


#plt.imshow("S2_img")
