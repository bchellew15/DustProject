from astropy.coordinates import SkyCoord

#convert to galactic coordinates
c = SkyCoord(244.8164, 55.33868, frame='fk5', unit="deg")

print(c.fk5.ra)
print(c.fk5.dec)
print(c.fk4.ra)
print(c.fk4.dec)
print(c.galactic.l)
print(c.galactic.b)


 
