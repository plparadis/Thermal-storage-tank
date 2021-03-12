# -*- coding: utf-8 -*-

from matplotlib import pyplot as plt


fig1, (ax1, ax11) = plt.subplots(2, 1, num='Montly Loads')
ax1.plot(paramRes2.y_num,np.flipud(Res2results.T[:,nb_t-1]), label='$\itt$={} h'.format(temps[nb_t-1]/3600))
ax1.plot(paramRes2.y_num,np.flipud(Res2results.T[:,round(nb_t/2)]), label='$\itt$={} h'.format(temps[round(nb_t/2)]/3600))  
ax1.plot(paramRes2.y_num,np.flipud(Res2results.T[:,round(nb_t/3.25)]), label='$\itt$={} h'.format(temps[round(nb_t/3.25)]/3600))
ax1.plot(paramRes2.y_num,np.flipud(Res2results.T[:,round(nb_t/5)]), label='$\itt$={} h'.format(temps[round(nb_t/5)]/3600))
ax1.plot(paramRes2.y_num,np.flipud(Res2results.T[:,round(nb_t/10)]), label='$\itt$={} h'.format(temps[round(nb_t/10)]/3600))
ax1.plot(paramRes2.y_num,np.flipud(Res2results.T[:,0]), label='$\itt$={} h'.format(temps[0]/3600)) 

ax1.set_xlabel('$\ity$ [m]', color='k', fontname='Lucida Bright',fontsize=9)
ax1.set_ylabel('$\itT_{SHW}$ [°C]', color='k', fontname='Lucida Bright',fontsize=9)

ax11.plot(temps/3600, Res2results.T[0,:], label='$\ity$={0:.2f} m'.format(paramRes2.y_num[paramRes2.nb_y-1]))        # Haut du réservoir
ax11.plot(temps/3600, Res2results.T[round(paramRes2.nb_y/2),:], label='$\ity$={0:.2f} m'.format(paramRes2.y_num[round(paramRes2.nb_y/2)])) # Mi hauteur
ax11.plot(temps/3600, Res2results.T[paramRes2.nb_y-1,:], label='$\ity$={0:.2f} m'.format(paramRes2.y_num[0]))    #Bas du réservoir
ax11.plot(temps/3600, HXefdresults.Tout[0,:], label='$\itDCW_{out}$')         # sortie DHW HX

ax11.set_xlabel('Time [h]', color='k', fontname='Lucida Bright',fontsize=9)
ax11.set_ylabel('$\itT_{SHW}$ or $\itT_{DCW_{out}}$ [°C]', color='k', fontname='Lucida Bright',fontsize=9)

ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=3, prop = {'family':'Lucida Bright', 'size':9})
ax11.legend(loc='upper center', bbox_to_anchor=(0.5, -0.20), ncol=2, prop = {'family':'Lucida Bright', 'size':9})

fig1.set_figheight(fig1.get_figheight()*1.5)
fig1.set_figwidth(fig1.get_figwidth())

fig1.tight_layout()  # otherwise the right y-label is slightly clipped


fig1.savefig('Tank Simulation results.jpg', format='JPG', dpi=300)

