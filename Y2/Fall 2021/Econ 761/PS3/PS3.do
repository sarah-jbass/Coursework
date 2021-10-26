/*-------------------------------------------------------------------------------------------------------------------------------------------
Class: Econ 761
Task: IO Problem Set 3
Date created: 10/20/21
Programmer: Sarah Bass
Comments:
do "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Fall 2021\Econ 761\PS3\PS3.do"
----------------------------------------------------------------------------------------------------------------------------------------*/
clear all 

cd "C:\Users\19195\OneDrive\Documents\University of Wisconsin\Coursework\Y2\Fall 2021\Econ 761\PS3"

import excel "cereal_ps3.xls", firstrow case(lower)
tempfile cereal 
save `cereal'
clear

import excel "demog_ps3.xls", firstrow case(lower)
merge 1:m city year quarter using `cereal'

