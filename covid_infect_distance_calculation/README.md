## Detection of SARS-CoV-2 infected area and calculation of distances


QuPath (v0.4.3) was used for detecting areas with SARS-CoV-2 infection. The image channel corresponding to SARS-CoV-2 nucleocapsid phosphoprotein was imported into QuPath, and a pixel classifier was trained to detect the CoV-2 infected areas. The classifier was applied on all images, and the resulting annotations were exported as geojson files.

These are in `qupath_outlines` folder

`001_calculate_distance_covid_infected_region.ipynb`: The annotations were imported as shapely(v2.1) polygon geometry in Python (3.10) and subsequently converted into MERSCOPE micron co-ordinate space. The distances between cells and CoV-2 infected area were computed using geopanda (v0.14.4).
