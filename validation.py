import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

nT      = 0.5
N_Max   = 2000
delt_tt = 0.02
alpha   = 0.9
beta    = 0.9 

xticks = np.linspace( alpha*N_Max, N_Max , 5)
xrange = (xticks[0], xticks[-1])

df = pd.read_csv('data.txt',sep='\s+') 

index = df.index 

##create a figure of 2*2 subplot for each probe point
fig,ax=plt.subplots(4,1,constrained_layout=True)

ax[0].set_title('Newamrk Validation,nT={0}'.format(nT))

ax[0].plot( index,df['U1'].values       , label='newmark'  )
ax[0].plot( index,df['U1_exact'].values , label='analytical' )
ax[0].set_xticks([],[])
ax[0].set_xlim(xrange)
ax[0].set_ylabel('U1[m]')
ax[0].legend(loc=2)

##-------------------------------------------------------------
ax[1].plot(index, df['U2'] ,label='newmark'          )
ax[1].plot(index, df['U2_exact'] ,label='analytical' )
ax[1].set_xticks([],[])
ax[1].set_xlim(xrange)
ax[1].set_ylabel('U2[m]')
 

##-------------------------------------------------------------
ax[2].plot(index, df['U3']  ,label='newmark'         )
ax[2].plot(index, df['U3_exact'] ,label='analytical'  )
ax[2].set_xticks(xticks)
ax[2].set_xlim(xrange)
ax[2].set_ylabel('U3 [m]')
 

##-------------------------------------------------------------
xticks = np.arange( beta*N_Max, N_Max , 100)
xrange = (xticks[0], xticks[-1])

ax[3].plot(index, df['U3']  ,label='newmark'         )
ax[3].plot(index, df['U3_exact'] ,label='analytical' )
ax[3].set_xticks(xticks)
ax[3].set_xlim(xrange)

ax[3].set_xlabel('Iteration')
ax[3].set_ylabel('U3[m]')
 
##-------------------------------------------------------------
# save original data figs
fig_name='newmark validation nT={0}'.format(nT)+'.png'
fig.savefig( fig_name )

plt.close()
















