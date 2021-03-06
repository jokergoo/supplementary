
Supplementary S1. Making Enhanced OncoPrint
==============================================

**Author**: Zuguang Gu ( z.gu@dkfz.de )

**Date**: 2016-05-13

----------------------------------------



To successfully run this example, version of **ComplexHeatmap** must be >= 1.10.1.
The newest version can be obtained by:


```r
library(devtools)
install_github("jokergoo/Complexheatmap")
```

Load the package.


```r
library(ComplexHeatmap)
```


An [OncoPrint](http://www.cbioportal.org/faq.jsp#what-are-oncoprints) is a way to visualize various genomic alterations
in a set of genes over multiple patients. The `oncoPrint()` function implemented in the **ComplexHeatmap** package provides great 
enhancement for the original [cBioPortal](http://www.cbioportal.org/index.do) style OncoPrint.
With the functionality in the **ComplexHeatmap** package, the OncoPrint can be highly flexibly customized, e.g. by
splitting genes into several groups by categorical variables,
appending more annotations or other heatmaps to associate other layers of information.

In this supplementary, we demonstrate the power of **ComplexHeatmap** to generate enhanced OncoPrints by visualizing 
mutations in the Lung Adenocarcinoma cohort from TCGA.
Mutation data is obtained from [cBioPortal](http://www.cbioportal.org/index.do) through the
following steps:

1. Select "**Lung Adenocarcinoma (TCGA, Provisoinal)**";
2. Click "**select From Recurrently Mutated Genes**" and select all the genes.
3. Since cBioPortal only allows visualizing 100 genes maximum, you need to split
   the gene list into two parts and merge the results afterwards.
4. Submit the gene list, and once the OncoPrint is loaded, click "**Download**" tab
   and copy the text below "**Type of Genetic alterations across all cases: (Alterations are summarized as MUT, Gain, HetLoss, etc.)**"

The processed mutation data is already stored in `lung_adenocarcinoma_TCGA_provisional_MutSig.rds`. Please note that rows should correspond to genes and columns to patients. 


```r
mat = readRDS("lung_adenocarcinoma_TCGA_provisional_MutSig.rds")
mat[1:4, 1:4]
```

```
##       TCGA-05-4384-01 TCGA-05-4390-01 TCGA-05-4425-01 TCGA-38-4631-01
## STK11 "  "            "  "            "  "            "  "           
## TP53  "MUT: Y205C"    "  "            "  "            "  "           
## KRAS  "  "            "MUT: G12V"     "  "            "  "           
## KEAP1 "MUT: G524C"    "  "            "  "            "MUT: F139L"
```

```r
dim(mat)
```

```
## [1] 173 230
```

In order to clearly visualize the data, genes for which less than 10% patients have mutations and
patients which have mutations in less than 5% of the genes are removed.


```r
l1 = apply(mat, 1, function(x) sum(!grepl("^\\s*$", x))/length(x) > 0.1)
l2 = apply(mat, 2, function(x) sum(!grepl("^\\s*$", x))/length(x) > 0.05)
mat = mat[l1, l2]
```

The biological functions of the recurrently mutated genes can give insights in the molecular basis of the cancer.
This information can be attached to the OncoPrint as an additional heatmap.

The biological functions annotated by Gene Ontology are obtained from [MSigDB](http://software.broadinstitute.org/gsea/msigdb).
The file `c5.bp.v5.0.symbols.gmt` is used.


```r
gene_set = strsplit(readLines("c5.bp.v5.0.symbols.gmt"), "\t")
names(gene_set) = sapply(gene_set, "[", 1)
gene_set = lapply(gene_set, "[", -(1:2))
```

`mat_gs` is a binary matrix which represents whether the gene has the corresponding biological function.
In order to reduce the amount of GO terms, only these terms having more than 3 genes annotated are kept.


```r
mat_gs = matrix(nrow = nrow(mat), ncol = length(gene_set))
colnames(mat_gs) = names(gene_set)
rownames(mat_gs) = rownames(mat)
for(i in seq_along(gene_set)) {
	mat_gs[, i] = rownames(mat) %in% gene_set[[i]] + 0
}
mat_gs = mat_gs[, colSums(mat_gs) > 3, drop = FALSE]
```

Clinical information can be added to the OncoPrint as column annotations. The [**TCGA2STAT** package](https://cran.r-project.org/web/packages/TCGA2STAT/index.html) is used
to retrieve clinical information directly from TCGA.


```r
# Following code didn't run when building the document.
library(TCGA2STAT)
Sys.setenv("TAR" = Sys.which("tar")) # just in case TAR environment variable is not se
mut = getTCGA(disease = "LUAD", data.type = "Mutation", type = "all", clinical=TRUE)
anno_df = mut$clinical[gsub("-01$", "", colnames(mat)), ]
```

Alternatively, the clinical data is already provided in `clinical_data.rds`. One thing that should be noted
is that the table returned by **TCGA2STAT** is a character matrix, so you may need to convert it to a data frame 
and convert numeric columns to real numbers.


```r
anno_df = readRDS("clinical_data.rds")
```

The following code defines how to visualize different alterations. The style is the same as 
[cBioPortal](http://www.cbioportal.org/index.do). Here `HOMDEL`, `AMP` and `MUT` are alteration
types encoded in `mat`.


```r
col = c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")
alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    HOMDEL = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["HOMDEL"], col = NA))
    },
    AMP = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["AMP"], col = NA))
    },
    MUT = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["MUT"], col = NA))
    }
)
```

If `alter_fun` is specifyed as a list of functions, graphics are added in a layer-by-layer mode. Graphics can
also be added in a grid-by-grid mode by specifying `alter_fun` as a single function. Now `alter_fun` accepts 
a fifth argument `v` which is a logical vector that shows whether the corresponding alteration exists in the current grid.
E.g:


```
## HOMDEL    AMP    MUT 
##   TRUE  FALSE  FALSE
```

Then graphics can be defined according to the existance of corresponding alterations:


```r
# NOTE: following code is not used when generating the heatmaps, it is just for demonstration
alter_fun = function(x, y, w, h, v) {
    # background
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    # graphics for three alterations
    if(v["HOMDEL"]) grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["HOMDEL"], col = NA))
    if(v["AMP"])    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["AMP"], col = NA))
    if(v["MUT"])    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["MUT"], col = NA))
}
```

Genes are split into two groups based on the amplification rate across patients. The variable
`amp` contains labels for genes indicating which group they are in. `amp` is converted to a 
factor to control the order of two row-slices on the plot.


```r
amp = ifelse(apply(mat, 1, function(x) sum(grepl("AMP", x))/length(x) > 0.1), "high AMP events", "low AMP events")
amp = factor(amp, levels = c("low AMP events", "high AMP events"))
```

Column annotations which contain clinical data are defined by `ha`. There are two simple annotations which are gender and stage
and one complex age annotation which is represented as points.


```r
gender = anno_df[, "gender"]
yearstobirth = as.numeric(anno_df[, "yearstobirth"])
pathologicstage = anno_df[, "pathologicstage"]
ha = HeatmapAnnotation(gender = gender, stage = pathologicstage,
	age = anno_points(yearstobirth, ylim = c(0, max(yearstobirth, na.rm = TRUE)), axis = TRUE),
	col = list(gender = c("male" = "red", "female" = "blue"),
		       stage = c("stage i" = "#FF0000", "stage ia" = "#FF6060", "stage ib" = "#FFB0B0", 
		       	         "stage iia" = "#60FF60", "stage iib" = "#B0FFB0",
		       	         "stage iiia" = "#6060FF", "stage iiib" = "#B0B0FF",
		       	         "stage iv" = "#FFFF00")),
	annotation_height = unit(c(5, 5, 15), "mm"),
    annotation_legend_param = list(gender = list(title = "Gender"),
                                   stage = list(title = "Stage"))
)
```

When everything for the OncoPrint is ready, the whole complex OncoPrint can be made.

Let's look back what is stored in `mat`.


```r
mat[1:4, 1:4]
```

```
##       TCGA-05-4390-01 TCGA-38-4631-01 TCGA-38-4632-01 TCGA-44-6145-01
## STK11 "  "            "  "            "  "            "  "           
## TP53  "  "            "  "            "MUT: A159P"    "  "           
## KRAS  "MUT: G12V"     "  "            "  "            "MUT: G12V"    
## KEAP1 "  "            "MUT: F139L"    "  "            "  "
```

The `get_type` argument in `oncoPrint()` expects a self-defined function which will be applied 
to extract all alteration types encoded in `mat`.
So for `get_type` defined in the following code, `"MUT: A159P; AMP"` will be converted to `c("MUT", "AMP")` internally.


```r
ht = oncoPrint(mat, get_type = function(x) gsub(":.*$", "", strsplit(x, ";")[[1]]),
    alter_fun = alter_fun, col = col, 
    column_title = "OncoPrint for recurrently mutated genes in Lung Adenocarcinoma",
    heatmap_legend_param = list(title = "Alterations", at = c("AMP", "HOMDEL", "MUT"), 
        labels = c("Amplification", "Deep deletion", "Mutation")), split = amp,
    bottom_annotation = ha)
```

For the matrix of biological functions, column names are added as a column annotation
using text rotated by 45 degrees.


```r
ha_cn = HeatmapAnnotation(cn = anno_text(colnames(mat_gs), rot = -45, just = "left", 
	offset = unit(1, "npc") - unit(1, "mm"), gp = gpar(fontsize = 8)), annotation_height = unit(6, "cm"))
```

Now the heatmap for the binary matrix can be added to the OncoPrint.


```r
ht_list = ht + Heatmap(mat_gs, col = c("0" = "white", "1" = "purple"), 
	rect_gp = gpar(col = "grey"), show_row_names = FALSE, cluster_columns = TRUE, 
    show_column_dend = FALSE, bottom_annotation = ha_cn, show_column_names = FALSE, 
    show_heatmap_legend = FALSE, width = unit(6, "cm"), column_title = "Map to Gene Ontology (BP)")
```

Finally the whole plot is drawn and customizations are applied to add the labels for the annotations afterwards.


```r
draw(ht_list, row_sub_title_side = "left")
decorate_annotation("gender", {
	grid.text("Gender", x = unit(-2, "mm"), just = "right")
})
decorate_annotation("stage", {
	grid.text("Stage", x = unit(-2, "mm"), just = "right")
})
decorate_annotation("age", {
	grid.text("Age", x = unit(-10, "mm"), just = "right")
})
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)

By default, rows are ordered by the mutation rate among patients and columns are sorted to show mutual exclusivity across patients.

## Session info


```r
sessionInfo()
```

```
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.11.4 (El Capitan)
## 
## locale:
## [1] C/en_US.UTF-8/C/C/C/C
## 
## attached base packages:
## [1] methods   grid      stats     graphics  grDevices utils     datasets 
## [8] base     
## 
## other attached packages:
## [1] ComplexHeatmap_1.10.1
## 
## loaded via a namespace (and not attached):
##  [1] circlize_0.3.7       dendextend_1.1.8     formatR_1.4         
##  [4] magrittr_1.5         evaluate_0.9         stringi_1.0-1       
##  [7] GlobalOptions_0.0.10 whisker_0.3-2        GetoptLong_0.1.3    
## [10] RColorBrewer_1.1-2   rjson_0.2.15         tools_3.2.3         
## [13] stringr_1.0.0        colorspace_1.2-6     shape_1.4.2         
## [16] knitr_1.13
```

