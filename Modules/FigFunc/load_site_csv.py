import pandas as pd

def load_site_csv(siteFile):
    
    '''load station/site list, crop data, select attributes required and convert to XYZ coords for SeisSol input
    '''
    table1 = pd.read_csv(siteFile)
    print(table1.keys())

    site_lon = table1['lon']
    site_lat = table1['lat']
    site_elev = table1['elev']
    return site_lon, site_lat, site_elev
