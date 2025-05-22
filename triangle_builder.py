from itertools import combinations
import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord


catalogue = pd.read_csv('data/catalogue.csv')
catalogue = catalogue.set_index(['con', 'hip'])[['ra', 'dec']]
catalogue['ra'] *= 15


def to_cartesian(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ra_rad = np.radians(df['ra'].values)
    dec_rad = np.radians(df['dec'].values)
    cos_dec = np.cos(dec_rad)
    x = cos_dec * np.cos(ra_rad)
    y = cos_dec * np.sin(ra_rad)
    z = np.sin(dec_rad)
    return x, y, z


catalogue['x'], catalogue['y'], catalogue['z'] = to_cartesian(catalogue)


def get_angle_distance(first: pd.DataFrame, second: pd.DataFrame) -> float:
    first_coord = SkyCoord(ra=first['ra'].values[0] * u.deg,
                          dec=first['dec'].values[0] * u.deg)
    second_coord = SkyCoord(ra=second['ra'].values[0] * u.deg,
                           dec=second['dec'].values[0] * u.deg)
    separation = first_coord.separation(second_coord)
    return separation.deg


def build_triangle(first: pd.DataFrame, second: pd.DataFrame, third: pd.DataFrame, con: str) -> pd.DataFrame:
    hip1 = first.index.get_level_values('hip')[0]
    hip2 = second.index.get_level_values('hip')[0]
    hip3 = third.index.get_level_values('hip')[0]
    
    v1 = np.array([first['x'].values[0], first['y'].values[0], first['z'].values[0]])
    v2 = np.array([second['x'].values[0], second['y'].values[0], second['z'].values[0]])
    v3 = np.array([third['x'].values[0], third['y'].values[0], third['z'].values[0]])
    
    a = np.arccos(np.clip(np.dot(v2, v3), -1.0, 1.0))
    b = np.arccos(np.clip(np.dot(v1, v3), -1.0, 1.0))
    c = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
    
    cos_A = (np.cos(a) - np.cos(b) * np.cos(c)) / (np.sin(b) * np.sin(c))
    cos_B = (np.cos(b) - np.cos(a) * np.cos(c)) / (np.sin(a) * np.sin(c))
    cos_C = (np.cos(c) - np.cos(a) * np.cos(b)) / (np.sin(a) * np.sin(b))
    
    angles = np.degrees([
        np.arccos(np.clip(cos_A, -1.0, 1.0)),
        np.arccos(np.clip(cos_B, -1.0, 1.0)),
        np.arccos(np.clip(cos_C, -1.0, 1.0))
    ])

    angle_hip_pairs = [
        (angles[0], hip1),
        (angles[1], hip2),
        (angles[2], hip3)
    ]
    
    sorted_pairs = sorted(angle_hip_pairs, key=lambda x: x[0])
    sorted_angles = [pair[0] for pair in sorted_pairs]
    sorted_hips = [pair[1] for pair in sorted_pairs]

    return pd.DataFrame({
        'con': con,
        'hip1': sorted_hips[0],
        'hip2': sorted_hips[1],
        'hip3': sorted_hips[2],
        'angle1': sorted_angles[0],
        'angle2': sorted_angles[1],
        'angle3': sorted_angles[2]
    }, index=[0])
    

constellations = catalogue.index.get_level_values('con').unique()
rows = []
MAX_ANGLE = 20
iterator = 0

for con in constellations:
    single_con = catalogue.loc[con]
    for (i, j, k) in combinations(range(len(single_con)), 3):
        first = single_con.iloc[[i]]
        second = single_con.iloc[[j]]
        third = single_con.iloc[[k]]
        if (get_angle_distance(first, second) > MAX_ANGLE or
            get_angle_distance(first, third) > MAX_ANGLE or
            get_angle_distance(second, third) > MAX_ANGLE):
            continue
        row = build_triangle(first, second, third, con)
        rows.append(row)
        iterator += 1
        if iterator % 50 == 0:
            print(f"total triangles: {iterator}")
                
                
methrics = pd.concat(rows, ignore_index=True).dropna()
methrics.to_csv('data/triangles.csv', index=False)