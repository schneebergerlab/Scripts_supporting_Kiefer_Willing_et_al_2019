library(ape)
library(caper)

args = commandArgs(trailingOnly=TRUE)

file_tree = args[1] 
file_geno = args[2] 
file_pheno = args[3]
species_num = as.numeric(args[4])
cnv_flag = as.numeric(args[5])

#######################

# load species tree
pls.tree <- read.tree(file_tree)

########################

# select phenotype: comment out other phenotype
# load phenotype: Flowering behaviour
m <- read.table(file_pheno)
pheno <- m$V2
names(pheno) <- m$V1

#######################

# load genotypes
d <- read.table(file_geno, header=T)
geno <- d[,c(4:(species_num+3))]
rownames(geno) <- d[,1]

# subset genotypes for species with phenotype
asso_matrix <- geno[,names(pheno)]
asso_matrix[asso_matrix > 1] <- 1 
asso_matrix[asso_matrix == -1] <- NA                   
t_asso_matrix <- t(asso_matrix)

if (!is.factor(pheno)) {
	t_asso_matrix <- t_asso_matrix[pheno >= 0,]
}

t_asso_matrix_uniq =  t_asso_matrix

# iterate individually
total_groups <- 0
count_groups <- 0
results <- list()

for(og_id in colnames(t_asso_matrix_uniq)){
	total_groups <- total_groups + 1
	#print(paste(total_groups, og_id, "(", length(colnames(t_asso_matrix_uniq)), ")", sep=" "))

	# build up OG specific data 
	c_pheno = c()
	c_og = c()
	if (!is.factor(pheno)) {
		c_pheno <- pheno[pheno >= 0 & is.na(t_asso_matrix_uniq[,og_id])==FALSE]
		c_og <- t_asso_matrix_uniq[pheno >= 0 & is.na(t_asso_matrix_uniq[,og_id])==FALSE, og_id]
	} else {
		c_pheno <- pheno[is.na(t_asso_matrix_uniq[,og_id])==FALSE]
		c_og <- t_asso_matrix_uniq[is.na(t_asso_matrix_uniq[,og_id])==FALSE, og_id]
	}
	n_og <- c_og
	n_og[n_og == 0] <- "No"
	n_og[n_og >= 1 & n_og != "No"] <- "Yes"
	n_og <- as.factor(n_og)
	c_species <- names(c_pheno)

	if(length(n_og[n_og == "Yes"]) > 3 && length(n_og[n_og == "No"]) > 3 && length(n_og) > length(t_asso_matrix_uniq[, 1])*0.80){
		count_groups <- count_groups + 1
		#print(count_groups)

		 # Prepare data for caper
                pp.data <- as.data.frame(c_species)
                pp.data <- cbind(pp.data, c_pheno)
                pp.data <- cbind(pp.data, n_og)
                pp.data <- cbind(pp.data, c_og)
                pp <- comparative.data(pls.tree, pp.data, c_species, na.omit=FALSE)

		if (!is.factor(pheno)) {
			try({
				if (cnv_flag == 1) {
				    	mod1 <- anova(crunch(c_pheno ~ c_og, pp, equal.branch.length=TRUE)) #NEEDED FOR ASTRAL TREE
    				}
				else {			
	    				mod1 <- anova(brunch(c_pheno ~ n_og, pp, equal.branch.length=TRUE))
				}

			    	results[[og_id]] <- mod1
			})
		}
		else {
			#try({
			#    mod1 <- anova(brunch(c_og ~ c_pheno, pp))
			#    results[[og_id]] = mod1
			#})
		}
	}
}


# Extract pvalues from anova outputs
pval = list()
for (i in 1:length(results)) {
	for (name in names(results[i])) { 
		pval[name] = results[i][[name]]$"Pr(>F)"[1]
	}
}

pval_sig = pval[pval<0.05]
pval_hsig = pval[pval<0.01]

# short report:
#print(paste("\n", "num pval:", length(pval), "num pval sig:", length(pval_sig), "num pval hsig:", length(pval_hsig), "\n", sep= " "))

# write results
write.table(as.matrix(pval_sig), file = paste(file_pheno, ".sig.out", sep=""), row.names=TRUE, col.names=FALSE, append=FALSE)
write.table(as.matrix(pval_hsig), file = paste(file_pheno, ".hsig.out", sep=""), row.names=TRUE, col.names=FALSE, append=FALSE)
