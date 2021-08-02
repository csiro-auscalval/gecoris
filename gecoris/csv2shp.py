###############################################################################
#  csv2shp.py
#
#  Project:     G-MaPIT
#  Purpose:     Routine to convert g-mapit or sarproz generated InSAR PS CSV 
#               standard into OGR SHP file ready for import into QGIS
#  Author:      Richard Czikhardt, czikhardt.richard@gmail.com
#  Last edit:   20.02.2020
#  Python dependencies: pandas, geopandas, fiona
#
###############################################################################
def csv2shp(inCSV,outSHP):

    import pandas as pd
    import geopandas
    from fiona.crs import from_epsg
    # import multiprocessing as mp
    
    #inCSV = "/home/rc/sw/g-mapit/TS.csv"
    #outSHP = "/home/rc/CopernicusRUS/gmapit/Cunovo_DSC124.shp"
    
    # determine if gmapit/SRPRZ csv type:
    with open(inCSV) as f:
        fline = f.readline()
        if fline.startswith('ID'):
            # load CSV as pandas dataframe:
            data = pd.read_csv(inCSV)
        else:
            data= pd.read_csv(inCSV,skiprows=1)
            # remove duplicit dates:
            dupl = [s for s in data.columns if ".1" in s]
            data = data.drop(columns=dupl)
         
    # convert to geopandas dataframe:
    gdf = geopandas.GeoDataFrame(
        data, geometry=geopandas.points_from_xy(data.LON, data.LAT))    
    # fix CRS:
    gdf.crs = from_epsg(4326)    
    # save SHP:
    gdf.to_file(outSHP)
    
if __name__ == "__main__":
    import sys
    csv2shp(sys.argv[1],sys.argv[2])