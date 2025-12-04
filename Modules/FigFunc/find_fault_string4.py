import geopandas as gpd
import pandas as pd

def find_fault_string(faultShp,strings):
    '''Find fault segments match 'string' and output sub table containing Fault_ID, Name, Dip_pref, SR_pref and Geometry'''
    gdf_xyz = gpd.read_file(faultShp,include_fields=["Fault_ID","Name", "Dip_pref","SR_pref",'Dip_dir',"Geometry"])
    gdf_out = gpd.GeoDataFrame(columns=gdf_xyz.keys())
    # print(gdf_out)

    for string in strings:
        # print(string)
        # print(gdf_xyz[gdf_xyz['Name'].str.contains(string)])
        gdf_out = pd.concat([gdf_out,gdf_xyz[gdf_xyz['Name'].str.contains(string)]])

    # gdf_out = pd.concat([gdf_out,gdf_xyz[gdf_xyz['Fault_ID']==570]])
    
    return gdf_out