# SparCC | Pvals and Corelation to GEPHI Edge File (Edit Required After File is Made)

karan_read = function(file_location){
  
  # read everything
  raw =  read.csv(file_location, stringsAsFactors = F, sep = '\t', header = F)
  
  # row & col names, appended with "1. <x>" to ensure uniqueness
  species = raw[1,]
  species = species[-1]
  species = apply(rbind((1:length(species)), species), 2, function(x) paste(x, sep = '', collapse = '. '))
  
  raw = matrix(as.numeric(unlist(raw[-1,-1])), ncol = length(species)) #raw is now only data
  rownames(raw) = species
  colnames(raw) = species
  
  return(raw)
  
}

karan_filter = function(corrMat, pMat, cor = 0.3, p = 0.05){
  
  n = nrow(corrMat) #NOTE: they're all square, so i'm using nrow as n
  species = rownames(corrMat)
  species = sapply(species, function(x) paste(unlist(strsplit(x, ' '))[-1], sep = '', collapse = ''))
  
  
  filtered = corrMat * (abs(corrMat)>cor) * (pMat < p)
  indices = which(filtered != 0)
  
  i_row = ((indices-1) %% n) + 1
  i_col = ceiling(indices/n)
  
  OTU.1 = species[i_row]
  OTU.2 = species[i_col]
  corr = corrMat[indices]
  pval = pMat[indices]
    
  return(data.frame(OTU.1, OTU.2, corr, pval))
  
}

## Using it on my personal four file (Please change the name here)

genus_20 = karan_read('Matrixes/cov_mat_SparCC_Genus_Top_20.out')
genus_30 = karan_read('Matrixes/cov_mat_SparCC_Genus_Top_30.out')
species_30 = karan_read('Matrixes/cov_mat_SparCC_Species_Top_30.out')
species_50 = karan_read('Matrixes/cov_mat_SparCC_Species_Top_50.out')

p_genus_20 = karan_read('P_Vals/pvals_two_sided_Genus_Top_20.txt')
p_genus_30 = karan_read('P_Vals/pvals_two_sided_Genus_Top_30.txt')
p_species_30 = karan_read('P_Vals/pvals_two_sided_Species_Top_30.txt')
p_species_50 = karan_read('P_Vals/pvals_two_sided_Species_Top_50.txt')

filter_genus_20 = karan_filter(genus_20, p_genus_20)
filter_genus_30 = karan_filter(genus_30, p_genus_30)
filter_species_30 = karan_filter(species_30, p_species_30)
filter_species_50 = karan_filter(species_50, p_species_50)

write.table(filter_genus_20, 'filter_genus_20')
write.table(filter_genus_30, 'filter_genus_30')
write.table(filter_species_30, 'filter_species_30')
write.table(filter_species_50, 'filter_species_50')
