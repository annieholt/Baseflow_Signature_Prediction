# Scripts to calculate hydrologic signatures, running TOSSH Toolbox functions, and to predict hydrologic signatures using random forest modeling.

* main.py: calls signature calculation functions, and some PCA analysis/plotting that was done as part of coursework
  
* signatures/calculate_sigs_camels: functions to run TOSSH toolbox and retrieve/organize outputs
* r/random_forest_R: new/primary random forest models, predicting hydrologic signatures based on catchment attributes (currently not organized into functions)
* r/rf_results_plotting: script to read in and generate plots based on random forest outputs, mostly related to variable importance and model performance
* r/assemble_catchment_attributes: assembling attributes datasets, from CARAVAN, CAMELS, and new attributes on wetland and geologic age. rely on downloaded data stored locally. includes some correlation analysis

* camels_gw_sigs_pi_v2.ojb: pickle object storing the signature data for CAMELS catchments, now including BFI_90

# Old scripts or analysis for coursework
* signatures/random_forest: Python workflow to generate random forest models and write outputs, predicting hydrologic signatures based on catchment attributes. Not convinced the variable importance is correct.
* signatures/pca and signatures/viz were scripts related to PCA clustering analysis for coursework, but not used in thesis work
