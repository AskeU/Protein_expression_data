# Protein_expression_data repository description
## Script description 
This is the initial code for the article: "Application of long-read sequencing to establish a cost-efficient high-throughput automated pipeline for generation of organism specific protein-expression data" 
The ```Lead_script.py``` runs the ```Merger_blast.py```, ```Double_barcode_Demultiplexer_vectorized.py``` and ```Demultiplexed_To_Stuctured_With_Lum_data.py``` in order to transform Nanopore FastQ files and Excel-sheets from out Hamilton robot into structured data and plots. 

Fw/rev primers are the primers containing barcodes and the independent scripts ```Hamilton_Data_Extractor(Individual).py``` and ```Clustered_luminescence.py``` visualizes  the data from the Hamilton alone and from the output of the ```Lead_script.py```. 

The scripts can still be improved but they work consistently with 16 mtplates at a time. They will be improved later.
