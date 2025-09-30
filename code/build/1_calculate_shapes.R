library(here)
library(sf)
# This stores your repository path as a function "here()"
here::i_am('code/build/1_calculate_shapes.R')
# Depending on your computer, you may need to change the path to your Dropbox folder
# For example, on Prakash's computer:
if(grep('mishrap', here()) == 1) {
  dropbox_dir = '/Dropbox (Personal)/Protected Area Fragmentation/'
} else if(grep('ahaan', here()) == 1){
    # fill in the blank
    dropbox_dir = '...'
}

# Read in the shapefile of protected areas
pa <- st_read(paste0(dropbox_dir, 'Data/protected sites/protectedsites.shp'))
# ...