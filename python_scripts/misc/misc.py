#this was in generate_plots.py
#histograms to verify that bootstrap dist is normal
'''

#sdss should be faster
data_hist = bootstrap_alphas_boss[0][1280]     #[int(bootstrap_alphas_sdss[0].shape[1]/2)]
data_hist = data_hist[(data_hist>-1)*(data_hist<2)] #cut off outliers

print("min and max")
print(min(data_hist))
print(max(data_hist))

range_min = 0
range_max = 1 #max(data_hist)

#avg_hist = np.mean(data_hist)
#var_hist = np.var(data_hist)
#pdf_x = np.linspace(range_min, range_max, 100)
#pdf_y = np.exp(-0.5*(pdf_x-avg_hist)**2/var_hist)
#removed prefactor 1.0/np.sqrt(2*np.pi*var_hist) from df_y

num_bins = 50
h, bin_edges, _ = plt.hist(data_hist, bins=num_bins, range=(range_min, range_max), density=True) #,normed=True)
print("hist:")
print(h)
bin_centers = (bin_edges[:-1] + bin_edges[1:])/2

#taken from stackoverflow
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

p0 = [1., 0., 1.]
coeff, var_matrix = curve_fit(gauss, bin_centers, h, p0=p0)
hist_fit = gauss(bin_centers, *coeff)

plt.plot(bin_centers, hist_fit)

print('Fitted mean = ', coeff[1])
print('Fitted standard deviation = ', coeff[2])

#bin_width = (range_max-range_min) / num_bins
#print("num samples:", data_hist.shape[0])
#a = data_hist.shape[0]*bin_width/(np.sqrt(2*np.pi)*np.sqrt(var_hist))
#print("amplitude:", a)

#plt.plot(pdf_x, a*pdf_y, 'k--')
#plt.show()

plt.show()
'''


'''
#vertical plot: 3 panels
(from generate_plots.py)
    
plt.figure(figsize=(6, 14))
    
plt.subplot(3, 1, 1)    
#compare original to tao
plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[1], c='r', drawstyle='steps', label=r'With $\tau$ Correction')
plt.plot(binned_lambdas, binned_stds[1], c='r', drawstyle='steps')
if bootstrap:
    plt.fill_between(binned_lambdas, bootstrap_binned_lower[0], bootstrap_binned_upper[0], linewidth=0.0, color='k', alpha=0.5, step='pre')
    plt.plot(binned_lambdas, bootstrap_binned_stds[0], c='m', drawstyle='steps')
    plt.fill_between(binned_lambdas, bootstrap_binned_lower[1], bootstrap_binned_upper[1], linewidth=0.0, color='r', alpha=0.5, step='pre')
    plt.plot(binned_lambdas, bootstrap_binned_stds[1], c='m', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(3, 1, 2)
#compare original to iris
plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[3], c='r', drawstyle='steps', label='With IRIS Data')
plt.plot(binned_lambdas, binned_stds[3], c='r', drawstyle='steps')
if bootstrap:
    plt.fill_between(binned_lambdas, bootstrap_binned_lower[0], bootstrap_binned_upper[0], linewidth=0.0, color='k', alpha=0.5, step='pre')
    plt.plot(binned_lambdas, bootstrap_binned_stds[0], c='m', drawstyle='steps')
    plt.fill_between(binned_lambdas, bootstrap_binned_lower[3], bootstrap_binned_upper[3], linewidth=0.0, color='r', alpha=0.5, step='pre')
    plt.plot(binned_lambdas, bootstrap_binned_stds[3], c='m', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(3, 1, 3)
#compare original to tao AND iris
#(or for SDSS, compare to BOSS)
plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[2], c='r', drawstyle='steps', label=r'With IRIS and $\tau$')
plt.plot(binned_lambdas, binned_stds[2], c='r', drawstyle='steps')
if bootstrap:
    plt.fill_between(binned_lambdas, bootstrap_binned_lower[0], bootstrap_binned_upper[0], linewidth=0.0, color='k', alpha=0.5, step='pre')
    plt.plot(binned_lambdas, bootstrap_binned_stds[0], c='m', drawstyle='steps')
    plt.fill_between(binned_lambdas, bootstrap_binned_lower[2], bootstrap_binned_upper[2], linewidth=0.0, color='r', alpha=0.5, step='pre')
    plt.plot(binned_lambdas, bootstrap_binned_stds[2], c='m', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.tight_layout()
if save:
    plt.savefig("../bootstrap_sdss_binned_102919.png")
plt.show()
'''

'''
#plot regions of sky (unbinned) same as original paper
#(if not bootstrap and not boss)

elif location == 3:
                ivar[b>35] = 0
                ivar[b<-35] = 0
            elif location == 4:
                ivar[(b>-35)*(b<35)] = 0
                ivar[b>50] = 0
                ivar[b<-50] = 0
            elif location == 5:
                ivar[(b>-50)*(b<50)] = 0
            elif location == 6:
                ivar[(l>60)*(l<300)] = 0
            elif location == 7:
                ivar[l<60] = 0
                ivar[(l>120)*(l<240)] = 0
                ivar[l>300] = 0
            elif location == 8:
                ivar[l<120] = 0
                ivar[l>240] = 0

plot_emissions([0, 4, 5, 6], ["Full Sky", "0<|b|<35", "35<|b|<50", "50<|b|<90"], ['k', 'b', 'r', 'g'])
plt.show()
plot_emissions([0, 7, 8, 9], ["Full Sky", "0<|l|<60", "60<|l|<120", "120<|l|<180"], ['k', 'b', 'r', 'g'])
plt.show()

if boss and not bootstrap:
plot_binned([2, 6, 7, 8], ['k', 'b', 'r', 'g'], ['Full Sky', '0<|b|<35', '35<|b|<50', '50<|b|<90'])
plt.show()
plot_binned([2, 9, 10, 11], ['k', 'b', 'r', 'g'], ['Full Sky', '0<|l|<60', '60<|l|<120', '120<|l|<180'])
plt.show()
'''
