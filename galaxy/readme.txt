To install the stereogene to a local galaxy, https://wiki.galaxyproject.org/Admin/Tools/AddToolTutorial is very useful. 
Add Stereogene to path.
Create correlation/ folders in tools/ folder of galaxy. 
Put stereogene.xml to the folder.
Take the tool_conf.xml.sample as sample, copy it to tool_conf.xml and add the content of our tool_conf_add.xml (do not break the xml structure!). 
Edit to refer to tools/correlation folder you prepred a line ago.

You can upload the chromosome length file from our local folder or get it from GetData -> UCSC Main -> Human, hg38, AllTables: chromInfo table.
Modify the creature name according to what you need.

Data files - I did not succeed in downloading full tracks from UCSC using 'Get data' interface. 
Happily, Galaxy can download files by ftp (Fetch Data button in Upload File tool).
UCSC provides BigWig for ftp, you are to convert it to Wig or Bed and then upload to Galaxy.
But, you con download data files from say Epigenomic roadmap directly to Galaxy using the data URL's.


