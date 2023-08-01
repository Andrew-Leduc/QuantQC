## data path
data_path <- '/Users/andrewleduc/Library/CloudStorage/GoogleDrive-leduc.an@husky.neu.edu/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/cellenONE/Trachea_3_6_pSCoPE/Trachea_3_6_pSCoPE_part1/evidence.txt'

## Links MS runs to plate inject wells (list runs in order they ran on mass spec)
linker <- read.csv('/Users/andrewleduc/Library/CloudStorage/GoogleDrive-leduc.an@husky.neu.edu/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/cellenONE/Trachea_3_6_pSCoPE/linker_age.csv')

## Read in cell isolation files from CellenONE and assign cell type
one <-read.delim("/Users/andrewleduc/Library/CloudStorage/GoogleDrive-leduc.an@husky.neu.edu/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/cellenONE/Trachea_3_6_pSCoPE/Young_isolated.xls")
two <- read.delim("/Users/andrewleduc/Library/CloudStorage/GoogleDrive-leduc.an@husky.neu.edu/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/cellenONE/Trachea_3_6_pSCoPE/old_isolated.xls")
all_cells <- rbind(one,two)




#Read and format raw data, link MS runs to injection wells
Trachea_3_7_prot <- MQ_to_nPOP(data_path,linker,PIF = .75,PEP = .04)


# Normalize single cell runs to reference channel, 
# filter out data points over twice reference
# Generate turn long table format to Peptide X Cell matrix
Trachea_3_7_prot <- TMT_Reference_channel_norm(Trachea_3_7_prot)


## mapping cellenONE meta data to raw data
Trachea_3_7_prot <- link_cellenONE_Raw(Trachea_3_7_prot,all_cells)



# Calculate statistics to look for decrease in LC/MS performance
# over time
Trachea_3_7_prot <- Calculate_run_order_statistics(Trachea_3_7_prot)


# Plot MS Intensity Statistics
PlotIntensityDrift(Trachea_3_7_prot)


#Plot Retention time Statistics
PlotRTDrift(Trachea_3_7_prot)




# Test negative controls, i.e. samples with no cell
Trachea_3_7_prot <- EvaluateNegativeControls(Trachea_3_7_prot,.43)


PlotNegCtrl(Trachea_3_7_prot,.43)


Trachea_3_7_prot <- FilterBadCells(Trachea_3_7_prot,.43)


# Compute the size relative to the carrier of each cell

PlotSCtoCarrierRatio(Trachea_3_7_prot)


# Collapse and give statistics

Trachea_3_7_prot <- CollapseToProtein(Trachea_3_7_prot,1)


PlotProtAndPep(Trachea_3_7_prot)


#

Trachea_3_7_prot <- SharedPeptideCor(Trachea_3_7_prot)

PlotPepCor(Trachea_3_7_prot)


# Impute 

Trachea_3_7_prot <- KNN_impute(Trachea_3_7_prot)





