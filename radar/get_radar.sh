#!/bin/bash

# Data location
data_loc="/g/data/rq0/rainfields3/"

# Radar site (see radar_site_list.csv)
radar_site="2"

# Date we want
year="2025"
month="10"
day="25"

filename="${data_loc}${radar_site}/${year}/prcp-crate/${radar_site}_${year}${month}${day}.prcp-crate.zip"

echo $filename
unzip $filename 



