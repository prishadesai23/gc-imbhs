#!/usr/bin/env python
# coding: utf-8

# In[4]:


from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os
import csv

directories = ['/work/10579/ag85537/red/spx15_reductions/reduce_ts23', 
               '/work/10579/ag85537/red/spx15_reductions/reduce_20250531', '/work/10579/ag85537/red/spx15_reductions/reduce_20250601', '/work/10579/ag85537/red/spx15_reductions/reduce_20250602', '/work/10579/ag85537/red/spx15_reductions/reduce_20250603', 
               '/work/10579/ag85537/red/spx15_reductions/reduce_20250912', '/work/10579/ag85537/red/spx15_reductions/reduce_20250913', '/work/10579/ag85537/red/spx15_reductions/reduce_20250915', 
               '/work/10579/ag85537/red/spx15_reductions/reduce_20260327', '/work/10579/ag85537/red/spx15_reductions/reduce_20260328', '/work/10579/ag85537/red/spx15_reductions/reduce_20260329', '/work/10579/ag85537/red/spx15_reductions/reduce_20260330']

# directories = ['old_bad_reductions/reduce_ts23', 
#                'old_bad_reductions/reduce_20250531', 'old_bad_reductions/reduce_20250601', 'old_bad_reductions/reduce_20250602', 'old_bad_reductions/reduce_20250603', 
#                'old_bad_reductions/reduce_20250912', 'old_bad_reductions/reduce_20250913', 'old_bad_reductions/reduce_20250915', 
#                'old_bad_reductions/reduce_20260327', 'old_bad_reductions/reduce_20260328', 'old_bad_reductions/reduce_20260329', 'old_bad_reductions/reduce_20260330']

name_map = {
    "m3": "M3",
    "m3star": "M3STAR" ,
    "m53 star": "M53STAR",
    "m53star": "M53STAR",
    "m53": "M53",
    "n2419": "NGC2419",
    "N2419": "NGC2419",
    "n2419_star": "NGC2419STAR",
    "N2419STAR": "NGC2419STAR",
    "n2419_star2": "NGC2419STAR2",
}

def normalize_name(name):
    return name_map.get(name.strip(), name.strip())  

with open("object_names.csv", mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Directory', 'Filename', 'Object'])
    for directory in directories:
        for filename in os.listdir(directory):
                    if filename.endswith('.fits'):
                        filepath = os.path.join(directory, filename)
                        try:
                            with fits.open(filepath) as hdul:
                                object_name = normalize_name(hdul[0].header.get('OBJECT', 'UNKNOWN'))
                                writer.writerow([directory, filename, object_name])
                        except Exception as e:
                            print(f'Error processing {filepath}: {e}')

object_stats = {}
for directory in directories:
    for filename in os.listdir(directory):
        if filename.endswith('.fits'):
            filepath = os.path.join(directory, filename)
            try:
                with fits.open(filepath) as hdul:
                    header = hdul[0].header
                    object_name = normalize_name(header.get('OBJECT', 'UNKNOWN'))
                    exptime = header.get('EXPTIME', 0)  # exposure time in seconds
                    hours = exptime / 3600.0
                    if object_name not in object_stats:
                        object_stats[object_name] = {'total_hours': 0.0, 'num_exposures': 0}
                    object_stats[object_name]['total_hours'] += hours
                    object_stats[object_name]['num_exposures'] += 1
            except Exception as e:
                print(f'Error processing {filepath}: {e}')

with open("object_stats.csv", mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Object', 'Total Hours', 'Total Exposures', 'Average Hours per Exposure'])
    for obj, stats in object_stats.items():
            avg_hours = stats['total_hours'] / stats['num_exposures'] if stats['num_exposures'] > 0 else 0
            writer.writerow([obj, stats['total_hours'], stats['num_exposures'], avg_hours])

