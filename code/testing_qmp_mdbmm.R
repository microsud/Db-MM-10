#### Get seq with copies
```{r}
# read in the function Doris Vandeputte obtained from github
# CopyNumber
copy_num
# otu table
counts_tables <- meta(ps)
counts_tables$SampleID <- rownames(counts_tables)
counts_tables <- counts_tables[,c("SampleID", "CopyNumber")]
counts_tables$total_copies <- as.numeric(counts_tables$CopyNumber)
#rownames(counts_tables) <- counts_tables[,1]


otu_qmp <- as.data.frame(t(abundances(ps))) # extract abundance table form copy number corrected pseq object

otu_qmp_2 <- otu_qmp / rowSums(otu_qmp)
otu_qmp_2$SampleID <- rownames(otu_qmp_2)
counts_tables$SampleID <- counts_tables$SampleID

otu_qmp_ldf <- reshape2::melt(otu_qmp_2)
colnames(otu_qmp_ldf) <- c("SampleID", "taxa", "reabcounts")

#The relative abundances obtained based on the sequencing reads were translated into total abundances by multiplying the relative abundance of each taxon by the total bacterial abundance (qPCR) in the sample. 
otu_qmp_pcr <- otu_qmp_ldf %>% 
  left_join(counts_tables, by="SampleID") %>% 
  mutate(counts = ceiling(reabcounts*total_copies))
otu_qmp_mat <- reshape2::dcast(otu_qmp_pcr, SampleID ~ taxa, value.var="counts")
rownames(otu_qmp_mat) <- otu_qmp_mat$SampleID
otu_qmp_mat <- otu_qmp_mat[,-1]

ps1.qmp <- ps
otu_qmp_A <- otu_table(t(otu_qmp_mat), taxa_are_rows = T)
otu_table(ps1.qmp) <- otu_qmp_A
ps1.qmp

```

**Copy number correction**  
  
  The copy numbers were obtained from barnap genome sequences.  

```{r eval=FALSE}
copy_num <- read.csv("data/rrna_MMDb/dbmm_copynumbers.csv", 
                     header = T, 
                     row.names = 1, 
                     stringsAsFactors = F)


DT::datatable(copy_num)
```

The copy numbers for individial strains was estimated using barrnap.  

```{r eval=FALSE}
# Extract strain counts from phyloseq object
#counts <- as.data.frame(abundances(ps))
#copy_num_file <- copy_num
#column_with_ids = "Genus_species"

ps1 <- copy_number_correction(ps= ps1.qmp, column_with_ids = "Genus_species", copy_num_tab=copy_num)
ps1
taxa_names(ps1) <- gsub("_", " ", taxa_names(ps1))
taxa_names(ps1)
qmp_df <- psmelt(ps1)


# remove the innoculum
ps.df.qmp <- subset(qmp_df, Timepoint != "X")
#ps.df.qmp$Species <- gsub("_", " ", ps.df.qmp$Species)


ps.df.qmp$TimepointNum <- ps.df.qmp$Timepoint
ps.df.qmp$TimepointNum <- as.numeric(gsub("T", "",ps.df.qmp$TimepointNum))
comp_plot_qmp <- ggline(ps.df.qmp, 
                        x="TimepointNum", 
                        y="Abundance", 
                        add= "mean_se",
                        color = "OTU", 
                        palette = strain.colors,
                        #linetype = 2,
                        #shape= 1,
                        #point.size=2,
                        size= 1,
                        show.line.label = F,
                        #scales = "free",
                        add.params = list(group="Species"),
                        ylab = "Abundance (log10)",
                        xlab = "Time (hours)") + theme_bw(base_size = 12) 
comp_plot_qmp <-  comp_plot_qmp +  geom_vline(aes(xintercept=1.5), linetype="dashed", alpha = 0.5) + geom_vline(aes(xintercept=2.5), linetype="dashed", alpha = 0.7) 

comp_plot_qmp <- ggpar(comp_plot_qmp, legend = "right", font.legend = c(10, "italic")
                       #, yscale = "log10"
)
comp_plot_qmp <- comp_plot_qmp + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                               labels = trans_format("log10", math_format(10^.x)))
comp_plot_qmp
