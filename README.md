# GEYSERS
Getting Exciting Young Supernova Experiment Recent Supernovae (GEYSERS) - version 0.1

Python script for grabbing recent transients from YSE-PZ with properties you're interested in (i.e. of a specific colour, rise time, host environment).

This script requires python 3 and astropy, alongside other classical libraries like numpy, scipy and pandas.

To run the script, you must have an output .csv from the YSE-PZ SQL Explorer: 
https://ziggy.ucolick.org/yse/explorer/ 
There is a pre-written query, 'DARK_YSE_Recent_Transients' which grabs data from all transients within the last 20 days (we can change this if requested). Run this query and download the output .csv file from. Save it into the GEYERS-master folder.  
The main code is the 'GEYSERS.py' script and can be easily launched from the shell as such:

python GEYSERS.py DARK_YSE_Recent_Transients.csv

There are several options to select candidates. You can see these with the following command:

python GEYSERS.py DARK_YSE_Recent_Transients.csv --help

The selection options are as follows. All colour and rise properties are limits (i.e. the minimum and maximum colour/rise time you wish to search for):

  --grmin       The desired g-r colour minimum limit. Default set to infinite.
  --grmax       The desired g-r colour maximum limit. Default set to infinite.
  --gimin       The desired g-i colour minimum limit. Default set to infinite.
  --gimax       The desired g-i colour maximum limit. Default set to infinite.
  --gzmin       The desired g-z colour minimum limit. Default set to infinite.
  --gzmax       The desired g-z colour maximum limit. Default set to infinite.
  --rimin       The desired r-i colour minimum limit. Default set to infinite.
  --rimax       The desired r-i colour maximum limit. Default set to infinite.
  --rzmin       The desired r-z colour minimum limit. Default set to infinite.
  --rzmax       The desired r-z colour maximum limit. Default set to infinite.
  --izmin       The desired i-z colour minimum limit. Default set to infinite.
  --izmax       The desired i-z colour maximum limit. Default set to infinite.
  --risemin     The rise rate (minimum). Use negative values for rising, positive for declining. Must specify band.
  --risemax     The rise rate (maximum). Use negative values for rising, positive for declining. Must specify band.
  --riseband    The rise band (g,r,i,z). This is required if you provide any rise limit.
  --host        The apparent host environment. Options are "host" / "nuclear" / "orphan"
  --maglim      The limit for faintest objects of interest (all bands). Default is set to 20.0 unless specified.

Add these options to the end of the command line. E.g., if you were only interested in slowly evolving red transients in host galaxies, you might run:

python GEYSERS.py DARK_YSE_Recent_Transients.csv --grmin 1.0 --gimin 1.2 --gzmin 1.3 --host host --risemin 0.00001 --riseband r

The script will print out a list of object IDs to your terminal which match the requirements. It will also generate a txt file of YSE-PZ URLs for candidates of interest ("GEYERS_output.txt").

This github contains a test SQL Query output (in the event we haven't observed for a few days....).

Please, please, play around with this script and send me your feedback! Let me know what you think is missing/what could be done to improve it!

Also...if it breaks, let me know that too!

I acknowledge support from Luca Izzo (for his amazing PS1_query script which helps select the host environments) and Christa Gall (for motivating me to do this).
