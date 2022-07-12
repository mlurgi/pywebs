

import pyproj
import csv

def read_ico_data():
    records = csv.DictReader(open('../original.csv', 'rU'), delimiter=';')
    coords = []
    
    for row in records:
        current_record = dict()
        current_record['year'] = row['Any']
        current_record['month'] = row['Mes']
        current_record['day'] = row['dia']
        current_record['species'] = row['Nom_cient']
        current_record['lat'] = row['latitud'].split(' ')
        current_record['lon'] = row['longitud'].split(' ')
        coords.append(current_record)
        
    return coords

def write_ico_processed_data(dict_data):
    header_names = ['year','month','day','species','lat', 'lon', 'lat_utm','lon_utm']
    out = csv.DictWriter(open('../ico_preprocessed.csv', 'w'), header_names, delimiter=',')
    
    #we write the headers
    headers_dict = dict()
    for n in header_names:
        headers_dict[n] = n
    
    out.writerow(headers_dict)
    
    out_row = dict()
    for register in dict_data:
        out.writerow(register)
        
def transform_coords_to_utm(dict_data):
    for row in dict_data:
        longitude = row['lon']
        degs = None
        mins = None
        secs = None
        for j in longitude:
            if j != '':
                if degs == None:
                    degs = int(j)
                elif mins == None:
                    mins = int(j)
                elif secs == None:
                    secs = float(j)
        
        lon = float(degs) + (((float(mins)*60)+secs)/3600)
    
        latitude = row['lat']
        degs = None
        mins = None
        secs = None
        for j in latitude:
            if j != '':
                if degs == None:
                    degs = int(j)
                elif mins == None:
                    mins = int(j)
                elif secs == None:
                    secs = float(j)
        
        lat = float(degs) + (((float(mins)*60)+secs)/3600)
    
        utm = pyproj.Proj(init='epsg:23031')
        wgs84 = pyproj.Proj(init='epsg:4326')
        x,y = pyproj.transform(wgs84,utm,lon,lat)
        
        row['lon_utm'] = x
        row['lat_utm'] = y
    

if __name__ == '__main__':
    a = read_ico_data()
    transform_coords_to_utm(a)
    write_ico_processed_data(a)
    
    
    
    