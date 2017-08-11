# Mouse-Diet-Sequencing-Data

File Guide and Explanation:

"MouseDietDADA2.R" is the original file that readies the data for creation of a phyloseq object for analysis. Generation of heatmap compares all diets and time points. Random forest modeling compares all diets and time points. ***This file must be run before any of the following files.***

"Mouse_Diet_noAnti.R" is the same as "MouseDietDADA2.R" but with all antibiotic time points removed.

"Diet_group_analysis.R" compares only the three diet groups (dN, dPD, dZD). Generates a heatmap and random forest model comparing only the diets. Antibiotic time points removed.

"Time_group_analysis.R" compares only the time points (d0, d5, d8, d12, d15). Generates a random forest model comparing only the time points. Antibiotic time points removed.

"Protein_Def_Analysis.R" compares only the time points of the protein deficient (dPD) samples. Generates a heatmap comparing only the protein deficient samples. Antibiotic time points removed.

"Zinc_Def_Analysis.R" compares only the time points of te zinc deficient (dZD) samples. Generates a heatmap comparing only the zinc deficient samples.

"Heatmaps_Each_Time_Point.R" compares the diets over each separate time point. Generates 4 separate heatmaps that compares the diets at each separate time point (d5, d8, d12, d15).

The "Heatmaps" folder contains heatmaps generated from the files above. "All_wAnti" was generated from "MouseDietDADA2.R". "All_noAnti" was generated from "Mouse_Diet_noAnti.R". The "Diet_Day5_Heatmap" through "Diet_Day15_Heatmap" were generated from "Heatmaps_Each_Time_Point". "Diet_Heatmap" was generated from "Diet_group_analysis.R". "dPD_heatmap" and dZD_heatmap" were generated from "Protein_Def_Analysis.R" and "Zinc_Def_Analysis.R" respectively. 

The "RF_analysis" folder contains the random forest modeling results from the above files and the ranked most important OTU's for classification. It also contains results from an alpha diversity analysis. 
