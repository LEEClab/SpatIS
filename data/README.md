# Metadata

## Files

Below we describe the data files and the columns they contain.

- spatis_paper_bat_positions.csv: movement data from _Sturnira lilium_ bats from 2009 and 2010.
  * Tag_ID: The ID of the radio tag attached to animals
  * Animal_ID: An ID given to each individual bat, its identity
  * x: x position (longitude)
  * y: y position (latitude)
  * date: date of the recorded position 
  * time: time of day of the recorded position 
  * position: whether the positions corresponded to the roost, to a horizontal position 
  (flight event), or to a vertical position (animal standing). 

- spatis_paper_individual_bat_info.csv: individual data from each bat.
  * Animal_ID: An ID given to each individual bat, its identity
  * initial_timestamp: first timestamp (in format YYYY-MM-DD HH:MM:SS) recorded for each individual
  * end_timestamp: last timestamp (in format YYYY-MM-DD HH:MM:SS) recorded for each individual
  * sex: sex of each individual
  * body_mass_g: body mass of each individual, in grams
