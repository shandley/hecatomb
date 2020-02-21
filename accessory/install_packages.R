# Install any required packages
# This code came directly from https://stackoverflow.com/questions/38928326/is-there-something-like-requirements-txt-for-r
# You should refer to that site. 

# You should probably also update the URL!

pkgLoad <- function( packages = "favourites" ) {

	if( length( packages ) == 1L && packages == "favourites" ) {
		packages <- c(
			      "tidyverse"
			      )
	}

	packagecheck <- match( packages, utils::installed.packages()[,1] )

	packagestoinstall <- packages[ is.na( packagecheck ) ]

	closestURL <- "https://cloud.r-project.org/"

	if( length( packagestoinstall ) > 0L ) {
		utils::install.packages( packagestoinstall,
					repos = closestURL
					)
	} else {
		print( "All requested packages already installed" )
	}

	for( package in packages ) {
		suppressPackageStartupMessages(
					       library( package, character.only = TRUE, quietly = TRUE )
					       )
	}

}

pkgLoad()
