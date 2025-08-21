For processing steps, see the paper (Methods).

Data: 

large files (.tif rasters) come from "ESA WorldCover" and "Copernicus CLC+".  
R scripts assume the default download filenames (do not rename the tiles).
ESA WorldCover:
- Data access: [https://esa-worldcover.org/en/data-access](https://esa-worldcover.org/en/data-access)
CLC+ Backbone:
- Data access: [https://land.copernicus.eu/en/products/clc-backbone/clc-backbone-2021](https://land.copernicus.eu/en/products/clc-backbone/clc-backbone-2018)

Data tables：

"NENP_data_extraction_spreadsheet.xlsx" : the complete data extraction spreadsheet，including metadata.
"data7.csv" : analysis-ready extracted data used by the R scripts.

Scripts overview：


"5.18.regular expression.R" – Data extraction

"5.23.ShinyDigitise.R", "6.10.data_1028.R", "6.13.data_1153.R", "6.13.data_1815.R", "6.3.data_2613.R", "6.3.data_516.R", "6.9.data_1448.R" – Effect size data extraction

"6.26.urban_snh_ESA.R", "6.28.rural_snh_ESA.R", "7.19.rural_snh_CLC.R", "7.19.urban_snh_CLC.R" – Moderator variable extraction

"7.6.effect_size.R" – Effect size calculation

"7.6.model_total.R", "7.13.model.bumb.R" – Meta-analysis models
