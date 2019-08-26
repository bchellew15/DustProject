#use images to calculate temperature

import matplotlib.pyplot as plt

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
