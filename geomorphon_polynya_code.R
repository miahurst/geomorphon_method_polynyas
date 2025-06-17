##### Code description ######

##this code will take tif files with longitude, laititude, time and sea ice concentration
##data and run them through the geomorphon pattern recognition algorithm
##the output will be tif files which have longitude, latitude, time, sea ice concentration,
##and also geomorphon landform type associated with each long, lat point

## Stepinski, T. F. and Jasiewicz, J.: Geomorphons – a new approach to classification of landforms,
## in: Proc. Geomorphometry 2011, 109–112, 2011.

## Jasiewicz, J. and Stepinski, T. F.: Geomorphons — a pattern recognition approach
## to classification and mapping of landforms, Geomorphology, 182, 147–156, 
## http://dx.doi.org/10.1016/j.geomorph.2012.11.005, 2013.

#### Clear memory and start afresh ####

rm(list = ls()) #clear R memory

##set the pathway to where you would like your tif geomorphon tif fies to save
getwd()
setwd("/Path/to where/you want/your output/tifs to/save")

#### Load in required packages ####

##Load in libraries

require(sp)

##couldn't find
require(rgdal)
#not available so need to load from website
install.packages("rgdal", repos="http://R-Forge.R-project.org", type="source")

##couldn't find
require(rgeos)
#not available so need to load from website
install.packages("rgeos", repos="http://R-Forge.R-project.org", type="source")

require(raster)

require(maptools)
install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")

require(rgrass)
library(rgrass)
#couldn't find this one - no package
#require(rgrass7) #to call GRASS GIS from R for geomorphons function

#### Setting up the paths for R to access ####
## and making sure CRS is set ##

##the main path should be set to where your data files are
##i.e., your sea ice concentration tifs
main_path<-"/Path/to where/your/tif/files/are"
##the grass path should be changed to where GRASS is installed on your machine
##and the version number should be changed accordingly
grass_path <- "/Applications/GRASS-8.3.app/Contents/Resources" 

#Define projection as UTM (Universal Transverse Mercator)
#a common map projection for small to medium-sized regions

utmCRS = CRS("+proj=utm +zone=23 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 
#specifies the UTM
##this should be changed to wherever the region is you are looking at
##in the case of our paper it defined:
#UTM zone 23 (corresponds to the zone of the Amundsen Sea)
#specifies the ellipsoid model used for the earth's shape - WGS84
#specifies the geodetic datum which aligns the ellipsoid to a specific location
#inidcates that coordinates are in meters (common for UTM)
#last part ensures no default parameters are applied

proj_LL  <- "+proj=longlat +ellps=WGS84 +datum=WGS84" 
#first part specifies a geographic cordinate system with lat and lon cordinates
#next part uses the WGS84 ellipsoid model for earth
#final bit uses the WGS84 geodetic datum

proj_UTM_m <- "+proj=utm +ellps=WGS84 +datum=WGS84 +zone=23 +south +units=m"
#this directly defines the same UTM projection as in utmCRS but without wrapping
#it in the CRS() function
#parameters are the same, but south indicates southern hemisphere
#coordinates measured in meters

#### Loading in sea ice concentration rasters and reprojecting them ####
## changing each tif file from epsg4326 to epsg3031 ##

##*if files are already in a coordinate reference system measuring units in meters
##*then this step can be missed out, but you will still need to set the target resolution

##because the geomorphon algorithm measures in meters the cordinate system needs
##to be changed from one which uses longitude and latitiude to one which uses
##meters
##in the case of this paper the tif files were in epsg4326 and needed to be changed
##to epsg3031

##set a target resolution for the tif files
##this needs to be set to whatever your source data resolution is
##the resolution is measured in meters
##the data used for this paper had a resolution of 0.083 degrees
##which was translated to 9.25km or 9250m
target_resolution <- 9250

##load the input raster (in this case they are in EPSG:4326)
input_raster <- raster("/path/to/raster/tif/file/")
##load each one in individually in their time stamps
##e.g., sea_ice_concentration_jan2020, sea_ice_concentration_feb2020, etc.
##and then run the below code for each tif file to create the geomorphon tif files

##define the target CRS (EPSG:3031) 
target_crs <- CRS("EPSG:3031")

##reproject the raster
output_raster <- projectRaster(input_raster, crs = target_crs, res=target_resolution)
##reprojects the raster using the target resolution specified and the target coordinate
##reference system

##save the reprojected raster
writeRaster(output_raster, "siconc_copernicus_epsg3031_2017oct9250m.tif",
            format = "GTiff", overwrite = TRUE)
##run through each time step file as above and save them as new reprojected files
##to be used in the next portion of code

###### running tif files through geomorphon algorithm ######

##load in the tif file you want to deal with
##one at a time as before, for each individual time step
siconc_oct2017 = terra::rast("Path/to/the/tif/file/siconc_copernicus_epsg3031_2017oct9250m.tif")

##this code sets up the environment for using GRASS GIS (a powerful open-source GIS software)

GRASS_INSTALLATION <- Sys.getenv("GRASS_INSTALLATION")
##this fetches the path to the GRASS GIS installation directory from the environment variables
##assumes you already have GRASS_INSTALLATION in your system's environment

file.info(GRASS_INSTALLATION)$isdir
##this checks if the fetched path is a directory (returns TRUE if it is)

##path to the GRASS GIS binaries
grass_bin <- file.path(grass_path, "bin")
##gass_path is assumed to be the base directory where GRASS GIS is installed
##file.path creates the full path to the bin directory inside the GRASS installation
##which contains the executable binaries for GRASS GIS

##set GRASS GIS environment variables
Sys.setenv(PATH = paste(Sys.getenv("PATH"), grass_bin, oct = ":"))
##this retrieves the current system PATH environment variables, which tells the
##operating system where to look executable programs

Sys.setenv(GRASS_INSTALLATION = grass_path)
##sets the GRASS_INSTALLATION environment variable to grass_path which is used
##internally by R and GRASS-related functions to locate the base GRASS GIS installation directory

Sys.setenv(GRASS_PROJSHARE = file.path(grass_path, "share", "proj"))
##GRASS GIS relies on projection definitions (e.g., EPSG codes) stored in a proj directory
##this sets the GRASS_PROJSHARE variable to point to the projection data files located 
##at grass_path/share/proj

Sys.setenv(GRASS_SH = "sh")
##specifies the shell that GRASS GIS commands should use
##in this case, it is set to sh which is the standard POSIX shell

##this code ^^ ensures that R is properly configured to locate and run GRASS GIS 
##commands, alongside the required dependencies, such as projection data

##initialize GRASS GIS
loc <- initGRASS(gisBase = Sys.getenv("GRASS_INSTALLATION"), home = tempdir(),
                 mapset = "PERMANENT", override = TRUE)
##initGRASS() is a function from the rgrass7 package that initialises a GRASS GIS
##session in R
##this code sets the base directory for the GRASS GIS session to a temporary directory
##created for this R session
##it also specifies the GRASS GIS mapset to be used - the permanenet mapset is a special
##default mapset where base data is usually stored
##the override = TRUE overwrites any existing GRASS session to start a new one

loc <-  initGRASS(grass_path, home=tempdir(), mapset = "PERMANENT", override = TRUE)
##this is a variant of the previous command where grass_path is directly provided
##instead of being retrieved from "GRASS INSTALLATION" as above
##both codes should work the same if correct

##use projection of your tif file (in this case was epsg3031)
execGRASS("g.proj", flags = c("c"), epsg = 3031)
##this code executes a GRASS GIS command from within R - this sets or displays projection
##information for the current mapset
##this sets the project's projection - which in this case is EPSG:3031 Antarctic Polar
##Steroegraphic projection

##writing the tif file as a raster to GRASS
write_RAST(x = siconc_oct2017,
           vname="siconc_oct2017",
           flags = c("overwrite")
)
##this function writes a raster layer to the GRASS GIS environment from R
##it allows you to save an R object (such as a raster)as a GRASS GIS raster

##setting the GRASS region to match the raster
execGRASS("g.region",
          parameters = list(raster = "siconc_oct2017",
                            res = as.character(res(siconc_oct2017)[1])
          )
)
##this command ensures that the region in GRASS is set according to the properties
##(extent and resolution) of the siconc_oct2017 raster

##displaying the current GRASS region
execGRASS("g.region", flags = c("p"))
##this command prints the current region settings in GRASS GIS, showing details like:
##the extent of the region (xmin, xmax, ymin, ymax)
##the resolution (cell size) of the region
##the coordinate system and projection used by the region

#### Calculate Geomorphons #####

##this is the section where we will have to change the search radius, 
##flat direction and distance threshold

##this code runs the GRASS GIS tool r.geomorphon
##which is used to classify terrain into geomorphons based on the elevation data
##provided - which in our case is sea ice conc data

##geomorphons are morphological features of the landscape that describe the 
##geometry of the terrain

##to find out more about the geomorphon pattern recognition algorithm
##please go to this link:
# https://grass.osgeo.org/grass-stable/manuals/r.geomorphon.html
##it provides in detail information about the paramters used in the algorithm
##and some of the customisable features which can be used and adapted

##run the geomorphon algorithm on the tif file loaded into the project
execGRASS("r.geomorphon",
          flags = c("overwrite", "verbose"),
          parameters = list(elevation = "siconc_oct2017",
                            forms = "geomorphons",
                            search = 20, 
                            flat = 0.0025,
                            dist = 1
          )
)

##search, flat, and dist are parameters:
##they are currently set at the limits which were found to be optimal

##search
##Search radius (in number of map cells)
##Determines length on the geodesic distances in all eight directions where line-of-sight is calculated
##To speed up calculation it determines only those cells whose centers fall within
##the distance... so one cell =9.25km in this map... search radius of 10 = 92.5km
#area of polynya (~27,000km2), meaning search radius is ~92km

##flat
##The difference (in degrees) between zenith and nadir line-of-sight which
##indicate flat direction 
##If higher threshold will produce more flat (homogeneous) maps
##If resolution of the map is low (more than 1 km per cell) threshold should be very
##small (much smaller than 1 degree) because on such distance 1 degree of difference
##means several meters of high difference
##e.g., a flatness threshold of 0.25 degrees means that if the difference in elevation
##is less than 0.25 degrees, the area will be considered flat

##dist
##flat distance - this is additional parameter defining the distance above which
##the threshold starts to octrease to avoid problems with pseudo-flat line-of-sights
##if real elevation difference appears on the distance where its value is higher

##you can read more about the parameters here:
# https://grass.osgeo.org/grass-stable/manuals/r.geomorphon.html

#### Plotting geomorphon landforms ####
##to visualise the geomorphon landforms we have jsut created you can plot them

##reading in the raster data and setting up the plot layout

##reading raster data
geom_siconc <- raster::raster(read_RAST(vname = "geomorphons"))
##this reads the raster data for geomorphons and then converts it to a raster object

##sets up the plot layout
par(mfrow=c(1,1)) 

##plotting the geomorphons raster
plot(geom_siconc,
     col = c(
       "grey",
       "black",
       "red",
       "orange",
       "yellow3",
       "yellow",
       "chartreuse",
       "cyan",
       "blue",
       "blue4"),
     legend = FALSE, main="")

par(xpd = TRUE)

##adding a legend
legend(
  "bottomleft",
  legend = c("Flat", "Summit", "Ridge", "Shoulder", "Spur","Slope",
             "Hollow","Footslope","Valley","Depression"),
  fill = c(
    "grey",
    "black",
    "red",
    "orange",
    "yellow3",
    "yellow",
    "chartreuse",
    "cyan",
    "blue",
    "blue4"),
  horiz = FALSE,
  #inset = -0.175
)

##adding a title with the parameters used
title(main="siconc_2017oct ...")

###### saving the geomorphon tif files to be used elsewhere #####
##then to save the geomorphon tif files you have created you use the following code
##these tif files can then be used to analyse the geomorphon output to explore
##the extent and distribution of polynyas identified by the algorithm

##if you want to save the produced geomorphon raster as a tif for QGIS
##save as GeoTIFF
writeRaster(geom_siconc, "label_for_your_geomorphon_file.tif", format = "GTiff", overwrite = TRUE)



