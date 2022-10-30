import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from cond_entropy_sym_k import *
from transfer_entropy_all_Hcond_sym_k import *

from pylab import plot, show, savefig, xlim, figure, \
                  ylim, legend, boxplot, setp, axes

m_embed = 1 # 1: data is already symbolized; 2,3,4: embedding dimension
dt = 0.2 # From paper 
k=1 # no. of time histories of target
l=1 # no. of time histories of source
flag_calc_rot = 0 # 1: calc rot ourselves from position data, 0: use encoded data of Valentini
flag_time_edit = 0 # 1: add the time steps missing as nan, 0: use ts as is ignoring the missing steps

minGen = 2; maxGen = 5
minRelease = 1; maxRelease = 12  # for each pair

# Modify string names for writing
fdir= './files/'; fstr = '_SS1_k'+str(k)+'_l'+str(l)+'.csv'

# Sampling 
isamp = 0
samptarr = np.array([0.2]) # make sure sampt/dt = integer 
sampN = int(samptarr[isamp]/dt) # skip these many,i.e. downsampling fac

str1 = './data/experimentalSelectedReleases'
str2 = './data/flightsExperimentalGen2to5'
str3 = './data/encoded/flightsSS1'

DF1 = pd.read_csv(str1 + '.csv', index_col=0)
within250 = DF1["Within250"]
print('No. of pairs with avg. linear dis > 250 m. : ',np.sum(within250)) 
if (flag_calc_rot==1):
    DF2 = pd.read_csv(str2 + '.csv')
else:
    if (dt>0.2): 
        print('Change SS number in str3 !!!!! ')
        time.sleep(20)
    DF3 = pd.read_csv(str3 + '.csv')
    

HcondE = np.nan*np.ones((maxGen-minGen+1)*(10*12)); HcondN = np.nan*np.ones((maxGen-minGen+1)*(10*12))
SDNcondE = np.nan*np.ones((maxGen-minGen+1)*(10*12)); SDEcondN = np.nan*np.ones((maxGen-minGen+1)*(10*12))
MI_EN = np.nan*np.ones((maxGen-minGen+1)*(10*12)); MI_NE = np.nan*np.ones((maxGen-minGen+1)*(10*12)) 
TE_EN = np.nan*np.ones((maxGen-minGen+1)*(10*12)); TE_NE = np.nan*np.ones((maxGen-minGen+1)*(10*12)) 
TE_star_EN = np.nan*np.ones((maxGen-minGen+1)*(10*12)); TE_star_NE = np.nan*np.ones((maxGen-minGen+1)*(10*12)) 
TE_star2_EN = np.nan*np.ones((maxGen-minGen+1)*(10*12)); TE_star2_NE = np.nan*np.ones((maxGen-minGen+1)*(10*12)) 
HcondNE = np.nan*np.ones((maxGen-minGen+1)*(10*12)); HcondEN = np.nan*np.ones((maxGen-minGen+1)*(10*12)); 
HcondNE2 = np.nan*np.ones((maxGen-minGen+1)*(10*12)); HcondEN2 = np.nan*np.ones((maxGen-minGen+1)*(10*12)); 

Gen = np.nan*np.ones((maxGen-minGen+1)*(10*12)); 
Pair = np.nan*np.ones((maxGen-minGen+1)*(10*12)); 
Release = np.nan*np.ones((maxGen-minGen+1)*(10*12))

icount = 0 # counter over all the pairs and all the releases considered in the analysis

for iGen in range(minGen,maxGen+1):
    df3i = DF3[DF3['Gen']==iGen]
    minPair = min(df3i['Pair']) # find range of pair nos. in this gen 
    maxPair = max(df3i['Pair'])

    for iPair in range(minPair,maxPair+1):
        print('Pair ', iPair)
        for iRelease in range(minRelease,maxRelease+1):
            print('Release ', iRelease)
            
            df1 = DF1[(DF1['Gen']==iGen) & (DF1['Pair']==iPair) & (DF1['Release']==iRelease)]
            df3 = DF3[(DF3['Gen']==iGen) & (DF3['Pair']==iPair) & (DF3['Release']==iRelease)]

            if ( (df1['Within250'].size>0) & (np.any(np.array(df1['Within250']))>0) ):
                if (np.array(df3['Er']).size==0):
                    print('Error: Within250>0 but Er,Nr not saved')
                    raise TypeError
                
                t = np.array(df3['Time']);  # Time
                Er0 = np.array(df3['Er']); Nr0 = np.array(df3['Nr']); # Rotation direction of E & N birds
                Nt = t.shape[0] # No. of time steps saved
            else:
                #print('Avg. diss between birds > 250 m. So, we ignore ... ')
                continue  #### skips the remaining lines when using in a for loop

            # Sampling 
            Er = Er0[0:Nt:sampN] 
            Nr = Nr0[0:Nt:sampN]

            # Remove nan's before calc entropy 
            Er = Er[~np.isnan(Er)]
            Nr = Nr[~np.isnan(Nr)]
            # print('Nt = ',Nt); time.sleep(0.5)

            # Calc. ER_E = HcondE = H(E^n+1|E^n), ER_N = HcondN = H(N^n+1|N^n)
            HcondE[icount] = cond_entropy_sym_k(Er,m_embed,k)
            HcondN[icount] = cond_entropy_sym_k(Nr,m_embed,k)

            # Calc. MI, TE, TE_C1, TE_C2
            TE_EN_all = transfer_entropy_all_Hcond_sym_k(Nr,Er,m_embed,k,l)
            TE_NE_all = transfer_entropy_all_Hcond_sym_k(Er,Nr,m_embed,k,l)

            MI_EN[icount] = TE_EN_all[0]; MI_NE[icount] = TE_NE_all[0]
            TE_EN[icount] = TE_EN_all[1]; TE_NE[icount] = TE_NE_all[1]
            TE_star_EN[icount] = TE_EN_all[2]; TE_star_NE[icount] = TE_NE_all[2]
            TE_star2_EN[icount] = TE_EN_all[3]; TE_star2_NE[icount] = TE_NE_all[3]
            HcondNE[icount] = TE_EN_all[4]; HcondEN[icount] = TE_NE_all[4]
            HcondNE2[icount] = TE_EN_all[5]; HcondEN2[icount] = TE_NE_all[5]

            Gen[icount] = iGen
            Pair[icount] = iPair
            Release[icount] = iRelease

            icount = icount + 1

    print('******* Gen ',iGen)

Gen = Gen[~np.isnan(Gen)]; Pair = Pair[~np.isnan(Pair)]; Release = Release[~np.isnan(Release)]
HcondE = HcondE[~np.isnan(HcondE)]; HcondN = HcondN[~np.isnan(HcondN)] 
MI_EN = MI_EN[~np.isnan(MI_EN)]; MI_NE = MI_NE[~np.isnan(MI_NE)] 
TE_EN = TE_EN[~np.isnan(TE_EN)]; TE_NE = TE_NE[~np.isnan(TE_NE)]
TE_star_EN = TE_star_EN[~np.isnan(TE_star_EN)]; TE_star_NE = TE_star_NE[~np.isnan(TE_star_NE)]
TE_star2_EN = TE_star2_EN[~np.isnan(TE_star2_EN)]; TE_star2_NE = TE_star2_NE[~np.isnan(TE_star2_NE)]
HcondEN = HcondEN[~np.isnan(HcondEN)]; HcondNE = HcondNE[~np.isnan(HcondNE)] 
HcondEN2 = HcondEN2[~np.isnan(HcondEN2)]; HcondNE2 = HcondNE2[~np.isnan(HcondNE2)] 

print('Total number of experiments = ', Gen.shape[0])

# Make Dataframe & Write in a csv file
df = pd.DataFrame()
df['Gen']=Gen
df['Pair']=Pair; df['Release']=Release
df['HcondE']=HcondE; df['HcondN']=HcondN
df['MI_EN']=MI_EN; df['MI_NE']=MI_NE
df['TE_EN']=TE_EN; df['TE_NE']=TE_NE
df['TE_star_EN']=TE_star_EN; df['TE_star_NE']=TE_star_NE
df['TE_star2_EN']=TE_star2_EN; df['TE_star2_NE']=TE_star2_NE
df['HcondEN']=HcondEN; df['HcondNE']=HcondNE
df['HcondEN2']=HcondEN2; df['HcondNE2']=HcondNE2

df.to_csv(fdir + 'TE_all_Hcond' +fstr)
