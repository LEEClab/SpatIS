# Metadata

## Files

Below we describe the data files and the columns they contain.

- `spatis_paper_bat_positions.csv`: movement data from _Sturnira lilium_ bats from 2009 and 2010.
  * Tag_ID: The ID of the radio tag attached to animals;
  * Animal_ID: An ID given to each individual bat, its identity;
  * x: x position (longitude);
  * y: y position (latitude);
  * date: date of the recorded position;
  * time: time of day of the recorded position;
  * position: whether the positions corresponded to the roost, to a horizontal position
  (flight event), or to a vertical position (animal standing). 
  
  The data was collected using UTM reference system, zone 23S, with datum and ellipsoid WGS84.  
  The proj4 code for that is: `"+proj=utm +datum=WGS84 +zone=23 +south +ellps=WGS84 +towgs84=0,0,0"`

- `spatis_paper_individual_bat_info.csv`: individual data from each bat.
  * Animal_ID: An ID given to each individual bat, its identity;
  * initial_timestamp: first timestamp (in format YYYY-MM-DD HH:MM:SS) recorded for each individual;
  * end_timestamp: last timestamp (in format YYYY-MM-DD HH:MM:SS) recorded for each individual;
  * sex: sex of each individual (M, F);
  * body_mass_g: body mass of each individual, in grams.
  
- `proportion_use_landcover_polygons.csv`: data on the number of position from each individual in each land cover polygon.
  * Animal_ID: An ID given to each individual bat, its identity;
  * polygon: ID of each land cover polygon visited by individuals;
  * success: number of positions of each animal to each land cover polygon;
  * sample.size: total number of positions by each individual;
  * failure: number of positions of each animal that were in other land cover polygons (not the focal one);
  * sex: sex of each individual (M, F);
  * body_mass_g: body mass of each individual, in grams.

- `spatis_paper_bat_positions.csv`: data on the number of position from each individual in each foraging area.
  * Animal_ID: An ID given to each individual bat, its identity;
  * area: ID of each foraging area visited by individuals;
  * success: number of positions of each animal to each foraging area;
  * sample.size: total number of positions by each individual;
  * failure: number of positions of each animal that were in other foraging areas (not the focal one);
  * sex: sex of each individual (M, F);
  * body_mass_g: body mass of each individual, in grams.

- `use_foraging_areas_density_fruits.csv`:
  * Foraging_area: ID of each foraging area visited by individuals;
  * Solanum: density of _Solanum_ fuits in each foraging area (number of individuals per 250 m<sup>2</sup>).
  * Piper: density of _Piper_ fuits in each foraging area (number of individuals per 250 m<sup>2</sup>).
  * Cecropia: density of _Cecropia_ fuits in each foraging area (number of individuals per 250 m<sup>2</sup>).

- `map_area` folder: folder with vector files with land cover maps and foraging areas.

- `files_movebank` folder: files uploaded to movebank (files with animal positions and information, 
excluding positions with no timestamp associated). 
