import numpy as np
import pandas as pd

def load_rake_data(datafile):

    PI = 1/180*np.pi

    data = pd.read_csv(datafile)

    u_max = np.max(data['ASl'])

    if u_max > 90:

        rk_max = np.max(data[data['ASl']<18.0]['rake'])
        rk_min = np.min(data[data['ASl']<18.0]['rake'])
        u_avg = np.mean(data[data['ASl']<18.0]['ASl'])
        rk_avg = np.mean(data[data['ASl']<18.0]['rake'])
        rk_50th = np.percentila(data[data['ASl']<18.0]['rake'],50)

    else:

        u_avg = np.mean(data[data['ASl']>1.0]['ASl'])
        
        rk_avg = np.mean(data[data['ASl']>1.0]['rake'])
        rk_max = np.max(data[data['ASl']>1.0]['rake'])
        rk_min = np.min(data[data['ASl']>1.0]['rake'])

        rk_50th = np.percentile(data[data['ASl']>1.0]['rake'],50)

    # rake = asin(-(Sld/sqrt(Sld^2+Sls^2)))/3.14*180+180

    return rk_max,rk_min, rk_avg, rk_50th, u_avg

def load_slip_data(datafile):

    PI = 1/180*np.pi

    data = pd.read_csv(datafile)
    u_max = np.max(data['ASl'])

    if u_max > 90:

        u_max = np.max(data[data['ASl']<18.0]['ASl'])
        u_avg = np.mean(data[data['ASl']<18.0]['ASl'])
        vr_avg = np.mean(data[data['ASl']<18.0]['Vr'])

    else:
        u_avg = np.mean(data[data['ASl']>1.0]['ASl'])
        vr_avg = np.mean(data[data['ASl']>1.0]['Vr'])

    # rake = np.cos(data['Sld']/data['Sls']*PI)

    return u_max, u_avg, vr_avg

def load_ss_data(ssfile):

    data = pd.read_csv(ssfile)

    # select Aal>0 will change cell selected in dataset, T=0 will be remove
    df = data[data["ASl"]>0.50]

    u_max = np.max(df['ASl'])


    if u_max > 90:

        u_max = np.max(df[df['ASl']<18.0]['ASl'])
        u_avg = np.mean(df[df['ASl']<18.0]['ASl'])

    else:
        u_avg = np.mean(df[df['ASl']>0.50]['ASl'])
    
    
    df1 = data[data["Time"]==0]
    df2 = data[data["TimeStep"].isin([80,154,160,100,82,68,63])]

    vr_mean = np.mean(df['Vr'])
    
    drop = df1['dstress'].to_numpy() - df2['dstress'].to_numpy()
    slp = df2['ASl'].to_numpy()

    # select where stress drop is positive

    index = np.where( (slp < 18.50) & (drop>10000.0) )

    drop_mean = np.sum(drop[index]*slp[index])/np.sum(slp[index])
    # drop_mean = np.mean(drop[index])

    drop_max = np.max(drop[index])
    drop_50th = np.percentile(drop[index], 50)
    drop_90th = np.percentile(drop[index], 80)
    drop_10th = np.percentile(drop[index], 20)
   

    return drop_mean, vr_mean, drop_max, drop_50th, drop_90th, drop_10th, u_max,u_avg

def load_ss_slp_data(ssfile):

    data = pd.read_csv(ssfile)

    # select Aal>0 will change cell selected in dataset, T=0 will be remove
    df = data[data["ASl"]>0.50]
    u_max = np.max(df['ASl'])

    if u_max > 90:

        u_max = np.max(df[df['ASl']<18.0]['ASl'])
        u_avg = np.mean(df[df['ASl']<18.0]['ASl'])

    else:
        u_avg = np.mean(df[df['ASl']>0.50]['ASl'])
    
    
    df1 = data[data["Time"]==0]
    df2 = data[data["TimeStep"].isin([80,154,160,100,82,68,63,76])]
    
    drop = df1['dstress'].to_numpy() - df2['dstress'].to_numpy()

    index = np.where( drop > 10000.0 )

    # setup critical slip distance 
    dc_mod = 0.5

    # area_all = np.sum(df2['Area'].to_numpy()[index])

    drop_mean = np.mean(drop[index])
    drop_max = np.max(drop[index])
    drop_50th = np.percentile(drop[index], 50)
    drop_90th = np.percentile(drop[index], 80)
    drop_10th = np.percentile(drop[index], 20)

    ss_ini = np.mean(df1['dstress'].to_numpy())
    ss_fin = np.mean(df2['dstress'].to_numpy())

    inf_eng = np.mean( 0.5* drop[index]  * df2['ASl'].to_numpy()[index] )
    rad_eng_mean = np.mean( 0.5* drop[index]  * df2['ASl'].to_numpy()[index]    - (0.5* drop[index] * dc_mod ))

    return drop_mean, drop_max, drop_50th, drop_90th, drop_10th, u_max,u_avg,ss_ini,ss_fin,inf_eng, rad_eng_mean

def load_ss_slp_area_radiate(ssfile):
    '''Load fault output variables including Time, TimeStep, ASl, Area, dstress (shear stress)
    calculate average slip, avearge stress drop, inferrable available energy, radiated energy, radiated energy per m^2, '''

    data = pd.read_csv(ssfile)

    df = data[data["ASl"]>0.50]
    u_max = np.max(df['ASl'])

    if u_max > 90:

        u_max = np.max(df[df['ASl']<18.0]['ASl'])
        u_avg = np.mean(df[df['ASl']<18.0]['ASl'])

    else:
        u_avg = np.mean(df[df['ASl']>0.50]['ASl'])
    
    
    df1 = data[data["Time"]==0]
    df2 = data[data["TimeStep"].isin([80,154,160,100,82,68,63,76])]
    
    drop = df1['dstress'].to_numpy() - df2['dstress'].to_numpy()

    # select where stress drop is positive
    # calcalute fracture energy and seismic efficiency


    index = np.where( drop > 1000.0 )

    # setup critical slip distance 
    dc_mod = 0.5

    
    # radiated energy Lambert et al. 2022, DeltaW- E_frac
    rad_eng = np.sum( 0.5* drop[index]  * df2['ASl'].to_numpy()[index] *df2['Area'].to_numpy()[index] )  -   np.sum(0.5* drop[index] * dc_mod * df2['Area'].to_numpy()[index]  )
    
    # rad_eng per Area
    rad_eng_mean = np.mean( 0.5* drop[index]  * df2['ASl'].to_numpy()[index]    - (0.5* drop[index] * dc_mod ))

    # area_all = np.sum(df2['Area'].to_numpy()[index])

    drop_mean = np.mean(drop[index])
    drop_max = np.max(drop[index])
    drop_50th = np.percentile(drop[index], 50)
    drop_90th = np.percentile(drop[index], 80)
    drop_10th = np.percentile(drop[index], 20)

    ss_ini = np.mean(df1['dstress'].to_numpy())
    ss_fin = np.mean(df2['dstress'].to_numpy())
   
    # inferable available energy DeltaW, Lambert et al. 2022, 1/2*average_stress_drop * average_slip
    
    inf_eng = np.mean( 0.5* drop[index]  * df2['ASl'].to_numpy()[index] )

    # inf_eng  = drop_mean * u_avg

    return drop_mean, drop_max, drop_50th, drop_90th, drop_10th, u_max,u_avg,ss_ini,ss_fin,inf_eng,rad_eng,rad_eng_mean



def load_ss_slp_area_data(ssfile):

    data = pd.read_csv(ssfile)

    # select Aal>0 will change cell selected in dataset, T=0 will be remove
    df = data[data["ASl"]>0.50]
    u_max = np.max(df['ASl'])

    if u_max > 90:

        u_max = np.max(df[df['ASl']<18.0]['ASl'])
        u_avg = np.mean(df[df['ASl']<18.0]['ASl'])

    else:
        u_avg = np.mean(df[df['ASl']>0.50]['ASl'])
    
    
    df1 = data[data["Time"]==0]
    df2 = data[data["TimeStep"].isin([80,154,160,100,82,68,63,76])]
    
    drop = df1['dstress'].to_numpy() - df2['dstress'].to_numpy()

    # select where stress drop is positive
    # calcalute fracture energy and seismic efficiency


    index = np.where( drop > 1000.0 )

    # setup critical slip distance 
    dc_mod = 0.5

    
    # radiated energy
    rad_eng = np.sum( 0.5* drop[index]  * df2['ASl'].to_numpy()[index] *df2['Area'].to_numpy()[index] )  -   np.sum(0.5* drop[index] * dc_mod * df2['Area'].to_numpy()[index]  )
    rad_eng_mean = np.mean( 0.5* drop[index]  * df2['ASl'].to_numpy()[index]    - (0.5* drop[index] * dc_mod ))

    # area_all = np.sum(df2['Area'].to_numpy()[index])

    drop_mean = np.mean(drop[index])
    drop_max = np.max(drop[index])
    drop_50th = np.percentile(drop[index], 50)
    drop_90th = np.percentile(drop[index], 80)
    drop_10th = np.percentile(drop[index], 20)

    ss_ini = np.mean(df1['dstress'].to_numpy())
    ss_fin = np.mean(df2['dstress'].to_numpy())
   
    inf_eng = np.mean( 0.5* drop[index]  * df2['ASl'].to_numpy()[index] )

    # inf_eng  = drop_mean * u_avg


    return drop_mean, drop_max, drop_50th, drop_90th, drop_10th, u_max,u_avg,ss_ini,ss_fin,inf_eng,rad_eng,rad_eng_mean