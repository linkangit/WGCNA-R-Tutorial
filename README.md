# WGCNA Tutorial in R: From Raw Data to Gene Networks

*A step-by-step guide to building weighted gene co-expression networks for any RNA-seq dataset*

## What is WGCNA and Why Should You Care?

Weighted Gene Co-expression Network Analysis (WGCNA) is like finding friend groups in your high school - except instead of teenagers, we're looking at genes that "hang out together" (get expressed together) under different conditions. 

Think of it this way: if you notice that every time Student A is at a party, Students B, C, and D are also there, you might conclude they're part of the same social group. WGCNA does the same thing with genes - it finds groups of genes that consistently increase or decrease together across your experimental conditions.

**Why is this useful?**
- **Disease research**: Find groups of genes that work together in cancer, diabetes, or neurological disorders
- **Drug discovery**: Identify gene networks that drugs might target
- **Plant/crop research**: Understand which genes coordinate responses to drought, disease, or growth
- **Development studies**: See how gene networks change as organisms develop

## What You'll Need Before Starting

**Your Data Requirements:**
- RNA-seq count data (raw counts, not normalized)
- At least 15-20 samples (more is better for network stability)
- Multiple experimental conditions or time points
- Sample metadata (which samples belong to which groups)

**Software Requirements:**
- R (version 4.0 or higher)
- RStudio (recommended for easier analysis)

**Time Investment:**
- Initial setup: 30 minutes
- Analysis: 2-4 hours depending on dataset size
- Interpretation: Variable based on your research goals

## Step 1: Setting Up Your Environment

First, let's install all the packages you'll need. Run these commands in R:

```r
# Install packages if you don't have them
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "WGCNA"))
install.packages(c("tidyverse", "magrittr"))

# Load libraries
library(tidyverse)     # For data manipulation and plotting
library(magrittr)      # For the %>% pipe operator
library(WGCNA)         # The star of our show
library(DESeq2)        # For normalizing count data
```

## Step 2: Understanding Your Data Structure

Before diving in, you need to understand how your data should be organized. WGCNA expects your data in a specific format:

**What your raw data probably looks like:**
```
Gene_ID    Sample1  Sample2  Sample3  Sample4
Gene_A     150      200      175      180
Gene_B     50       75       60       65
Gene_C     1000     1200     1100     1150
```

**What WGCNA wants:**
- Rows = samples/treatments
- Columns = genes
- No missing values
- Normalized expression data

Let's load and examine your data:

```r
# Load your count data
# Replace "your_data.txt" with your actual file name
data <- readr::read_delim("your_data.txt", delim = "\t")

# Look at the first few rows and columns
data[1:5, 1:10]

# Check the structure
str(data)
dim(data)  # Shows [number of genes, number of samples + 1]
```

**Common data issues and fixes:**
- **Problem**: Gene names in first column aren't named properly
  - **Fix**: `names(data)[1] = "Gene_ID"`
- **Problem**: Sample names have spaces or special characters
  - **Fix**: `names(data) <- gsub("[^A-Za-z0-9_]", "_", names(data))`

## Step 3: Exploratory Data Analysis - Know Your Data

Before building networks, you need to understand your data's quality and characteristics. Think of this as getting to know your data personally.

```r
# Create a metadata dataframe for your samples
# Modify this based on your experimental design
sample_names <- names(data)[-1]  # All columns except gene names

# Example: if your samples are named like "Treatment1_Rep1", "Treatment1_Rep2", etc.
metadata <- data.frame(
  Sample = sample_names,
  Treatment = gsub("_.*", "", sample_names),  # Extract treatment name
  Replicate = gsub(".*_", "", sample_names)   # Extract replicate number
)

# View your metadata
print(metadata)
```

Now let's visualize the data to spot any obvious problems:

```r
# Convert to long format for plotting
data_long <- data %>%
  pivot_longer(cols = -1, names_to = "Sample", values_to = "Count") %>%
  left_join(metadata, by = "Sample")

# Plot expression distributions by sample
data_long %>%
  ggplot(aes(x = Sample, y = Count)) +
  geom_violin(aes(fill = Treatment)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_y_log10() +  # Log scale because RNA-seq data is highly skewed
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Distribution of Gene Expression by Sample",
       y = "Gene Expression (log10 scale)")
```

**What to look for:**
- **Outlier samples**: Samples with very different distributions might be technical failures
- **Batch effects**: Samples processed on different days might cluster together
- **Low count samples**: Samples with overall low counts might need to be removed

## Step 4: Data Normalization - Making Genes Comparable

Raw RNA-seq counts aren't directly comparable between samples due to differences in library size and composition. We'll use DESeq2 to normalize the data.

```r
# Prepare data for DESeq2
count_matrix <- as.matrix(data[,-1])  # Remove gene names column
rownames(count_matrix) <- data[[1]]   # Set gene names as row names

# Make sure counts are integers (DESeq2 requirement)
count_matrix <- round(count_matrix)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ Treatment  # Adjust based on your experimental design
)

# Run DESeq2 normalization
dds <- DESeq(dds)

# Get variance-stabilized data (better for correlation analysis than raw normalized counts)
vsd <- varianceStabilizingTransformation(dds)
normalized_data <- assay(vsd)

print(paste("Data dimensions:", nrow(normalized_data), "genes,", ncol(normalized_data), "samples"))
```

## Step 5: Gene Filtering - Focus on the Signal

Not all genes are informative for network analysis. We'll keep only genes that show meaningful variation across samples.

```r
# Calculate variance for each gene
gene_variances <- apply(normalized_data, 1, var)
summary(gene_variances)

# Keep top 25% most variable genes (you can adjust this)
variance_threshold <- quantile(gene_variances, 0.75)
filtered_data <- normalized_data[gene_variances > variance_threshold, ]

print(paste("After filtering:", nrow(filtered_data), "genes remain"))

# Visualize the effect of filtering
data.frame(
  Variance = gene_variances,
  Kept = gene_variances > variance_threshold
) %>%
  ggplot(aes(x = Variance, fill = Kept)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  geom_vline(xintercept = variance_threshold, color = "red", linetype = "dashed") +
  scale_x_log10() +
  theme_bw() +
  labs(title = "Gene Filtering Based on Variance",
       x = "Gene Variance (log10 scale)")
```

**Filtering considerations:**
- **Too strict**: You might lose biologically relevant genes
- **Too loose**: Analysis will be slower and might include noise
- **Sweet spot**: Usually keeping 20-50% of genes works well

## Step 6: Prepare Data for WGCNA

WGCNA expects samples as rows and genes as columns (opposite of typical RNA-seq format).

```r
# Transpose the data
wgcna_input <- t(filtered_data)

# Check the new dimensions
print(paste("WGCNA input:", nrow(wgcna_input), "samples,", ncol(wgcna_input), "genes"))

# Check for missing values (WGCNA can't handle them)
any(is.na(wgcna_input))

# If you have missing values, you can:
# wgcna_input[is.na(wgcna_input)] <- 0  # Replace with zeros (not ideal)
# Or remove genes/samples with missing values
```

## Step 7: Choose Soft Thresholding Power - Finding the Goldilocks Zone

This is one of the most critical steps. We need to find a "soft threshold" that creates a network that's not too dense (everything connected to everything) or too sparse (nothing connected to anything).

```r
# Enable multi-threading for faster computation
allowWGCNAThreads()

# Test different soft threshold powers
powers <- c(1:10, seq(12, 20, by = 2))

# This step takes a few minutes - perfect time for coffee!
sft <- pickSoftThreshold(
  wgcna_input,
  powerVector = powers,
  verbose = 5
)

# Plot the results
par(mfrow = c(1, 2))

# Plot 1: Scale-free topology fit
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit (R²)",
     type = "n", 
     main = "Scale Independence")

text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, 
     cex = 0.9, 
     col = "red")

# Add reference line at R² = 0.80
abline(h = 0.80, col = "red")

# Plot 2: Mean connectivity
plot(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     type = "n",
     main = "Mean Connectivity")

text(sft$fitIndices[,1], 
     sft$fitIndices[,5], 
     labels = powers, 
     cex = 0.9, 
     col = "red")
```

**How to choose the right power:**
1. Look for the lowest power where the R² curve flattens out above 0.80
2. Don't go too high (mean connectivity shouldn't be below 100)
3. Common values are between 6-20 for RNA-seq data

```r
# Choose your power based on the plots
chosen_power <- 12  # Adjust based on your plots
print(paste("Selected soft threshold power:", chosen_power))
```

## Step 8: Build the Network - The Main Event

Now comes the exciting part - actually building the co-expression network and finding gene modules (clusters).

```r
# This is the computationally intensive step
# For large datasets (>5000 genes), this might take 30-60 minutes
network <- blockwiseModules(
  wgcna_input,
  
  # Power parameter
  power = chosen_power,
  networkType = "signed",  # Distinguishes positive and negative correlations
  
  # Module detection parameters
  deepSplit = 2,           # How finely to split modules (0-4, higher = more modules)
  minModuleSize = 30,      # Minimum genes per module
  
  # Module merging
  mergeCutHeight = 0.25,   # Merge similar modules (lower = more merging)
  
  # Performance options
  maxBlockSize = 5000,     # Max genes to process at once
  
  # Save intermediate results
  saveTOMs = TRUE,
  saveTOMFileBase = "NetworkTOM",
  
  # Output options
  numericLabels = TRUE,
  verbose = 3
)

# Check how many modules were found
module_colors <- labels2colors(network$colors)
table(module_colors)
```

**Understanding the parameters:**
- **networkType = "signed"**: Genes that correlate negatively are treated differently from those that correlate positively
- **deepSplit**: Higher values create more, smaller modules
- **minModuleSize**: Smaller modules might be noise; larger modules are more robust
- **mergeCutHeight**: Lower values merge more similar modules together

## Step 9: Visualize Your Network Structure

Let's see what modules we discovered:

```r
# Plot the gene dendrogram with module colors
plotDendroAndColors(
  network$dendrograms[[1]], 
  module_colors[network$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, 
  hang = 0.03,
  addGuide = TRUE, 
  guideHang = 0.05
)

# Create a summary table
module_summary <- data.frame(
  Gene_ID = names(network$colors),
  Module = module_colors
)

# Count genes per module
module_counts <- table(module_summary$Module)
print("Genes per module:")
print(sort(module_counts, decreasing = TRUE))
```

## Step 10: Relate Modules to Your Experimental Conditions

This is where biology meets statistics. We want to know which gene modules are associated with your experimental treatments.

```r
# Calculate module eigengenes (the "average" gene in each module)
module_eigengenes <- moduleEigengenes(wgcna_input, module_colors)$eigengenes

# Reorder modules by similarity
module_eigengenes <- orderMEs(module_eigengenes)

# Add sample information
module_eigengenes$Sample <- rownames(module_eigengenes)
module_eigengenes <- module_eigengenes %>%
  left_join(metadata, by = "Sample")

# Create correlation heatmap
me_long <- module_eigengenes %>%
  select(-Sample, -Replicate) %>%  # Adjust based on your metadata columns
  pivot_longer(cols = starts_with("ME"), 
               names_to = "Module", 
               values_to = "Eigengene") %>%
  mutate(Module = gsub("ME", "", Module))

# Plot module-trait relationships
me_long %>%
  ggplot(aes(x = Treatment, y = Module, fill = Eigengene)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Module-Treatment Relationships",
       fill = "Module\nEigengene")
```

**Interpreting the heatmap:**
- **Red**: Module genes are highly expressed in that treatment
- **Blue**: Module genes are lowly expressed in that treatment  
- **White**: No strong association

## Step 11: Export Results for Further Analysis

Let's save our results in formats useful for downstream analysis:

```r
# Save module assignments
write_tsv(module_summary, "gene_modules.txt")

# Save module eigengenes
write_tsv(module_eigengenes, "module_eigengenes.txt")

# For network visualization in Cytoscape, create an edge list
# (This focuses on modules of interest to keep file size manageable)
interesting_modules <- c("blue", "turquoise", "brown")  # Adjust based on your results
interesting_genes <- module_summary %>%
  filter(Module %in% interesting_modules) %>%
  pull(Gene_ID)

# Calculate network connections for interesting genes only
interesting_data <- filtered_data[interesting_genes, ]
TOM <- TOMsimilarityFromExpr(t(interesting_data), power = chosen_power)

# Create edge list
rownames(TOM) <- colnames(TOM) <- interesting_genes
edge_list <- TOM %>%
  as.data.frame() %>%
  rownames_to_column("Gene1") %>%
  pivot_longer(-Gene1, names_to = "Gene2", values_to = "Weight") %>%
  filter(Gene1 < Gene2, Weight > 0.1) %>%  # Remove weak connections and duplicates
  arrange(desc(Weight))

write_tsv(edge_list, "network_edges.txt")
```

## Step 12: Biological Interpretation - Making Sense of Your Results

Now comes the most important part - what do your results mean biologically?

```r
# Identify hub genes (highly connected genes in each module)
hub_genes <- chooseTopHubInEachModule(
  wgcna_input, 
  module_colors,
  omitColors = "grey"  # Grey module contains unassigned genes
)

print("Hub genes for each module:")
print(hub_genes)

# Get expression profiles for modules of interest
plot_modules <- c("blue", "turquoise")  # Adjust based on your results

expression_profiles <- filtered_data[module_summary$Gene_ID[module_summary$Module %in% plot_modules], ] %>%
  as.data.frame() %>%
  rownames_to_column("Gene_ID") %>%
  left_join(module_summary, by = "Gene_ID") %>%
  pivot_longer(cols = -c(Gene_ID, Module), names_to = "Sample", values_to = "Expression") %>%
  left_join(metadata, by = "Sample")

# Plot expression profiles
expression_profiles %>%
  ggplot(aes(x = Treatment, y = Expression, group = Gene_ID)) +
  geom_line(alpha = 0.3, aes(color = Module)) +
  stat_summary(aes(group = Module), fun = mean, geom = "line", size = 2, color = "black") +
  facet_wrap(~Module) +
  theme_bw() +
  labs(title = "Gene Expression Profiles by Module",
       y = "Normalized Expression")
```

## Troubleshooting Common Issues

**Problem: "Not enough samples" error**
- **Solution**: WGCNA needs at least 15-20 samples. Consider combining conditions or using a different analysis method.

**Problem: No modules found or everything in one module**
- **Solution**: Try different `deepSplit` values (0-4) or adjust `minModuleSize`.

**Problem: Too many small modules**
- **Solution**: Increase `minModuleSize` or decrease `mergeCutHeight`.

**Problem: Analysis is very slow**
- **Solution**: Filter more stringently (keep fewer genes) or use `maxBlockSize` parameter.

**Problem: Soft threshold plots don't show clear cutoff**
- **Solution**: Your data might not follow scale-free topology. This can still work, just pick a reasonable power (6-15).

## Next Steps: Where to Go From Here

1. **Gene Ontology Analysis**: Use tools like DAVID or clusterProfiler to find what biological processes your modules represent.

2. **Network Visualization**: Import your edge list into Cytoscape for beautiful network visualizations.

3. **Experimental Validation**: Design experiments to test if your hub genes really control the biological processes you're studying.

4. **Integration with Other Data**: Combine with protein-protein interaction networks, ChIP-seq data, or clinical outcomes.

5. **Time Course Analysis**: If you have time series data, look at how modules change over time.

## Key Takeaways

- WGCNA finds groups of genes that behave similarly across your conditions
- The "soft threshold" is critical - spend time getting this right
- Modules should make biological sense - if they don't, adjust your parameters
- Focus on modules that correlate with your experimental conditions
- Hub genes within interesting modules are great candidates for follow-up experiments

Remember: WGCNA is an exploratory tool. It generates hypotheses about gene relationships that need experimental validation. The real biological insights come from combining these computational predictions with your domain expertise and further experiments.

## Final Checklist

Before considering your WGCNA analysis complete:

- [ ] Soft threshold achieves reasonable scale-free topology (R² > 0.8)
- [ ] Module sizes are reasonable (30-1000 genes typically)
- [ ] Modules correlate with experimental conditions in biologically sensible ways
- [ ] Hub genes are known players in your biological system (if known)
- [ ] Results are saved in formats suitable for downstream analysis
- [ ] You have a plan for experimental validation of key findings

Good luck with your gene network analysis!
