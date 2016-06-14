To install the stereogene to a local galaxy, https://wiki.galaxyproject.org/Admin/Tools/AddToolTutorial is very useful. 
Add Stereogene to path.
Create correlation/ folders in tools/ folder of galaxy. 
Put stereogene.xml to the folder.
Take the tool_conf.xml.sample as sample, copy it to tool_conf.xml and add the content of our tool_conf_add.xml (do not break the xml structure!). 
Edit to refer to tools/correlation folder you prepred a line ago.

You can upload the chromosome length file from our local folder or get it from GetData -> UCSC Main -> Human, hg38, AllTables: chromInfo table.
Modify the creature name according to what you need.

Data files - I did not succeed in downloading full tracks from UCSC using 'Get data' interface. 
Galaxy can download files by ftp (Fetch Data button in Upload File tool) if you have a bed/wigg URL.

UCSC provides BigWig for ftp, you are to convert it to Wig or Bed and then upload to Galaxy.

You can download data files from say Epigenomic roadmap directly to Galaxy using the data URL's. Examples:

ftp://ftp.genboree.org/EpigenomeAtlas/Current-Release/sample-experiment/CD8_Primary_Cells/Histone_H3K27ac/UW.CD8_Primary_Cells.H3K27ac.RO_01679.Histone.DS21779.wig.gz
ftp://ftp.genboree.org/EpigenomeAtlas/Current-Release/sample-experiment/CD8_Primary_Cells/Histone_H3K4me3/UW.CD8_Primary_Cells.H3K4me3.RO_01679.Histone.DS21778.wig.gz
ftp://ftp.genboree.org/EpigenomeAtlas/Current-Release/sample-experiment/CD8_Primary_Cells/Chromatin_Accessibility/UW.CD8_Primary_Cells.ChromatinAccessibility.RO_01701.DS17885.wig.gz

Eventually, we use sh script as wrapper, and the thing started to work.

Critical: galaxy is to know the file format, do not use 'unknown'.

Galaxy can give full name of a data set and its type (extension): ${query.name} ${query.ext}
