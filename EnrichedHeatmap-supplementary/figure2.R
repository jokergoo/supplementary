
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(GetoptLong)
library(circlize)
library(RColorBrewer)

load("roadmap_normalized_matrices.RData")

expr_mean = rowMeans(expr[, SAMPLE$subgroup == "subgroup1"]) - rowMeans(expr[, SAMPLE$subgroup == "subgroup2"])
expr_split = ifelse(expr_mean > 0, "high", "low")
expr_split = factor(expr_split, levels = c("high", "low"))

set.seed(123)
upstream_index = length(attr(meth_mat_mean, "upstream_index"))
meth_split = kmeans(meth_mat_mean[, seq(round(upstream_index*0.8), round(upstream_index*1.4))], centers = 2)$cluster
x = tapply(rowMeans(meth_mat_mean[, seq(round(upstream_index*0.8), round(upstream_index*1.4))]), meth_split, mean)
od = structure(order(x), names = names(x))
meth_split = paste0("cluster", od[as.character(meth_split)])

combined_split = paste(meth_split, expr_split, sep = "|")

l = combined_split != "cluster2|high"
tss = tss[l]
expr = expr[l, ]
hist_mat_corr_list = lapply(hist_mat_corr_list, function(x) x[l, ])
hist_mat_mean_list = lapply(hist_mat_mean_list, function(x) x[l, ])
hist_mat_diff_list = lapply(hist_mat_diff_list, function(x) x[l, ])
mat_neg_cr = mat_neg_cr[l, ]
mat_cgi = mat_cgi[l, ]
meth_mat_corr = meth_mat_corr[l, ]
meth_mat_mean = meth_mat_mean[l, ]
meth_mat_diff = meth_mat_diff[l, ]
expr_split = expr_split[l]
meth_split = meth_split[l]
combined_split = combined_split[l]
n_row_cluster = length(unique(combined_split))

merge_row_order = function(l_list) {
	do.call("c", lapply(l_list, function(l) {
		if(sum(l) == 0) return(integer(0))
		if(sum(l) == 1) return(which(l))
		dend1 = as.dendrogram(hclust(dist_by_closeness(mat_neg_cr[l, ])))
		dend1 = reorder(dend1, wts = -enriched_score(mat_neg_cr[l, ]))
		od = order.dendrogram(dend1)
		which(l)[od]
	}))
}

row_order = merge_row_order(list(
	combined_split == "cluster1|high",
	combined_split == "cluster1|low",
	combined_split == "cluster2|low"
))

dend1 = as.dendrogram(hclust(dist(t(expr[, SAMPLE$subgroup == "subgroup1"]))))
hc1 = as.hclust(reorder(dend1, colMeans(expr[, SAMPLE$subgroup == "subgroup1"])))
expr_col_od1 = hc1$order
dend2 = as.dendrogram(hclust(dist(t(expr[, SAMPLE$subgroup == "subgroup2"]))))
hc2 = as.hclust(reorder(dend2, colMeans(expr[, SAMPLE$subgroup == "subgroup2"])))
expr_col_od2 = hc2$order
expr_col_od = c(which(SAMPLE$subgroup == "subgroup1")[expr_col_od1], 
	            which(SAMPLE$subgroup == "subgroup2")[expr_col_od2])

ht_list = Heatmap(expr, name = "expr", show_row_names = FALSE,
	show_column_names = FALSE, width = unit(4, "cm"), show_column_dend = FALSE, 
	cluster_columns = FALSE, column_order = expr_col_od,
	top_annotation = HeatmapAnnotation(df = SAMPLE[, -1], col = COLOR, 
		show_annotation_name = TRUE, annotation_name_side = "left"),
	column_title = "Expression", column_title_gp = gpar(fontsize = 12),
	show_row_dend = FALSE, use_raster = TRUE, raster_quality = 2)

library(genefilter)
df = rowttests(expr, factor(SAMPLE$subgroup))
top_genes = rownames(df[order(df$p.value)[1:20], ])

index =  which(rownames(expr) %in% top_genes)
labels = gene_symbol[rownames(expr)[index]]
ht_list = rowAnnotation(sig_gene = row_anno_link(at = index, labels = labels,
		side = "left", labels_gp = gpar(fontsize = 8), link_width = unit(5, "mm"), 
		padding = 0.8, extend = unit(c(2, 0), "cm")), 
	width = max_text_width(labels, gp = gpar(fontsize = 8)) + unit(5, "mm")) + ht_list

gl = width(gene[names(tss)])
gl[gl > quantile(gl, 0.95)] = quantile(gl, 0.95)
ht_list = ht_list + rowAnnotation(gene_len = row_anno_points(gl, size = unit(1, "mm"), gp = gpar(col = "#00000040")), 
	width = unit(1.5, "cm"))

axis_name = c("-5kb", "TSS", "10kb")
ht_list = ht_list + EnrichedHeatmap(mat_cgi, col = c("white", "darkorange"), name = "CGI",
	column_title_gp = gpar(fontsize = 12),
	top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = "darkorange", 
		lty = 1:n_row_cluster), yaxis_facing = "left")), 
	top_annotation_height = unit(2, "cm"), column_title = "CGI", axis_name = axis_name,
	axis_name_gp = gpar(fontsize = 8), use_raster = TRUE, raster_quality = 2) 

cor_col_fun = colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "red"))
ht_list = ht_list + EnrichedHeatmap(meth_mat_corr, col = cor_col_fun, name = "meth_corr", 
	top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(pos_col = "red", 
		neg_col = "darkgreen", lty = 1:n_row_cluster), yaxis_facing = "left")), 
	top_annotation_height = unit(2, "cm"), column_title = "meth_corr", 
	axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 12),
	use_raster = TRUE, raster_quality = 2)

meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
ht_list = ht_list + EnrichedHeatmap(meth_mat_mean, col = meth_col_fun, name = "meth_mean", 
	column_title = "meth_mean", column_title_gp = gpar(fontsize = 12), axis_name = axis_name,
	top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = "red", 
		lty = 1:n_row_cluster), yaxis_facing = "left")),
	axis_name_gp = gpar(fontsize = 8), use_raster = TRUE, raster_quality = 2)

generate_diff_color_fun = function(x) {
	q = quantile(x, c(0.05, 0.95))
	max_q = max(abs(q))
	colorRamp2(c(-max_q, 0, max_q), c("#3794bf", "#FFFFFF", "#df8640"))
}

meth_diff_col_fun = generate_diff_color_fun(meth_mat_diff)
ht_list = ht_list + EnrichedHeatmap(meth_mat_diff, col = meth_diff_col_fun,
	name = "meth_diff", column_title = "meth_diff", column_title_gp = gpar(fontsize = 12), axis_name = axis_name,
	top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(pos_col = "#df8640", 
		neg_col = "#3794bf", lty = 1:n_row_cluster), yaxis_facing = "left")),
	axis_name_gp = gpar(fontsize = 8), use_raster = TRUE, raster_quality = 2)

ht_list2 = NULL
ht_list1 = NULL
mark_name = names(hist_mat_corr_list)
for(i in seq_along(hist_mat_corr_list)) {
	if(i == 2) {
		ht_list1 = ht_list
		ht_list = NULL
	}

	ht_list = ht_list + EnrichedHeatmap(hist_mat_corr_list[[i]], col = cor_col_fun, 
		name = qq("@{mark_name[i]}_corr"),
		top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(pos_col = "red",
			neg_col = "darkgreen", lty = 1:n_row_cluster), yaxis_facing = "left")), 
		top_annotation_height = unit(2, "cm"), column_title = qq("@{mark_name[i]}_corr"), 
		column_title_gp = gpar(fontsize = 12), axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), 
		use_raster = TRUE, raster_quality = 2)

    ht_list = ht_list + EnrichedHeatmap(hist_mat_mean_list[[i]], 
    	col = colorRamp2(c(0, quantile(hist_mat_mean_list[[i]], 0.95)), c("white", "purple")), 
    	name = qq("@{mark_name[i]}_mean"), column_title = qq("@{mark_name[i]}_mean"), 
    	column_title_gp = gpar(fontsize = 12), axis_name = axis_name,
		top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = "purple", 
			lty = 1:n_row_cluster), yaxis_facing = "left")),
		axis_name_gp = gpar(fontsize = 8), use_raster = TRUE, raster_quality = 2)

	ht_list = ht_list + EnrichedHeatmap(hist_mat_diff_list[[i]], name = qq("@{mark_name[i]}_diff"), 
		col = generate_diff_color_fun(hist_mat_diff_list[[i]]), column_title = qq("@{mark_name[i]}_diff"), 
		column_title_gp = gpar(fontsize = 12), axis_name = axis_name,
		top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(pos_col = "#df8640", 
			neg_col = "#3794bf", lty = 1:n_row_cluster), yaxis_facing = "left")),
		axis_name_gp = gpar(fontsize = 8), use_raster = TRUE, raster_quality = 2)

}
ht_list2 = ht_list

split = as.vector(combined_split)
split[combined_split == "cluster1|high"] = "cluster1"
split[combined_split == "cluster1|low"] = "cluster2"
split[combined_split == "cluster2|low"] = "cluster3"

ht_list = Heatmap(expr_split, show_row_names = FALSE, name = "expr_diff", 
	col = c("high" = "red", "low" = "darkgreen"), 
	show_column_names = FALSE, width = unit(2, "mm")) + ht_list1

gb1 = grid.grabExpr({
ht_list = draw(ht_list, 
	cluster_rows = FALSE, row_order = row_order, show_row_dend = FALSE,
	split = split, heatmap_legend_side = "bottom", gap = unit(2, "mm"),
	show_heatmap_legend = FALSE, show_annotation_legend = FALSE)

add_boxplot_of_gene_length = function(ht_list) {
	
	row_order_list = row_order(ht_list)
	lt = lapply(row_order_list, function(ind) gl[ind])
	bx = boxplot(lt, plot = FALSE)$stats
	n = length(row_order_list)
	x_ind = (seq_len(n) - 0.5)/n
	w = 1/n*0.5
	decorate_annotation("gene_len", slice = 1, {
		rg = range(bx)
		rg[1] = rg[1] - (rg[2] - rg[1])*0.1
		rg[2] = rg[2] + (rg[2] - rg[1])*0.1
		pushViewport(viewport(y = unit(1, "npc") + unit(1, "mm"), just = "bottom", height = unit(2, "cm"), yscale = rg))
		grid.rect(gp = gpar(col = "black"))
		grid.segments(x_ind - w/2, bx[5, ], x_ind + w/2, bx[5, ], default.units = "native", gp = gpar(lty = 1:n))
		grid.segments(x_ind - w/2, bx[1, ], x_ind + w/2, bx[1, ], default.units = "native", gp = gpar(lty = 1:n))
		grid.segments(x_ind, bx[1, ], x_ind, bx[5, ], default.units = "native", gp = gpar(lty = 1:n))
		grid.rect(x_ind, colMeans(bx[c(4, 2), ]), width = w, height = bx[4, ] - bx[2, ], default.units = "native", 
			gp = gpar(fill = c("red", "darkgreen", "darkgreen"), lty = 1:n))
		grid.segments(x_ind - w/2, bx[3, ], x_ind + w/2, bx[3, ], default.units = "native", gp = gpar(lty = 1:n))
		grid.text("Gene length", y = unit(1, "npc") + unit(2.5, "mm"), gp = gpar(fontsize = 12), just = "bottom")
		grid.segments(unit(1, "npc") - unit(1, "mm"), c(0, 100000, 200000), unit(1, "npc"), c(0, 100000, 200000), default.units = "native")
		grid.text("200kb", unit(1, "npc") - unit(2, "mm"), unit(200000, "native") + unit(2, "mm"), default.units = "native", rot = 90, just = c("right", "bottom"), 
			gp = gpar(fontsize = 8))
		upViewport()
	})
}

add_boxplot_of_gene_length(ht_list)
i = 0
for(f in names(ht_list1@ht_list)) {
	if(grepl("meth|H3K4me1|H3K4me3|H3K27ac|H3K27me3", f)) {
		decorate_column_title(f, {
			grid.rect(height = unit(0.8, "npc"), gp = gpar(fill = brewer.pal(8, "Set2")[as.integer(i/3)+1], col = NA))
			grid.text(ht_list1@ht_list[[f]]@column_title, gp = gpar(fontsize = 12))
		})
		i = i + 1
	}
}
decorate_annotation("gene_len", slice = n_row_cluster, {
	grid.segments(c(0, 200000), unit(0, "npc"), c(0, 200000), unit(-1, "mm"), default.units = "native")
	grid.text("0kb", unit(0, "native") - unit(2, "mm"), unit(-2, "mm"), gp = gpar(fontsize = 8), just = c("left", "top"))
	grid.text("200kb", unit(200000, "native") + unit(1, "mm"), unit(-2, "mm"), gp = gpar(fontsize = 8), 
		just = c("right", "top"))
}) 
})

ht_list = Heatmap(expr_split, show_row_names = FALSE, name = "expr_diff", 
	col = c("high" = "red", "low" = "darkgreen"), 
	show_column_names = FALSE, width = unit(2, "mm")) + ht_list2

gb2 = grid.grabExpr({
draw(ht_list,
	cluster_rows = FALSE, row_order = row_order, show_row_dend = FALSE,
	split = split, heatmap_legend_side = "bottom", gap = unit(2, "mm"),
	show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
for(f in names(ht_list2@ht_list)) {
	decorate_column_title(f, {
		grid.rect(height = unit(0.8, "npc"), gp = gpar(fill = brewer.pal(8, "Set2")[as.integer(i/3)+1], col = NA))
		grid.text(ht_list2@ht_list[[f]]@column_title, gp = gpar(fontsize = 12))
	})
	i = i + 1
}
})


pdf("figure2.pdf", width = 16, height = 10)

grid.newpage()

lgd_list = packLegend(
	Legend(at = names(COLOR$group), title = "Group", legend_gp = gpar(fill = COLOR$group)),
	Legend(at = names(COLOR$sample_type), title = "Sample type", legend_gp = gpar(fill = COLOR$sample_type)),
	Legend(at = names(COLOR$subgroup), labels = c("Embryonic cell", "Mature cell"), title = "Subgroup", legend_gp = gpar(fill = COLOR$subgroup)),
	color_mapping_legend(ht_list1@ht_list[["expr"]]@matrix_color_mapping, plot = FALSE, title = "Expression", at = c(0, 4, 8)),
	color_mapping_legend(ht_list1@ht_list[["CGI"]]@matrix_color_mapping, plot = FALSE, title = "Overlap to CGI"),
	Legend(at = c(-1, 0, 1), title = "Correlation", col_fun = cor_col_fun),
	Legend(at = c(0, 0.5, 1), title = "Methylation", col_fun = meth_col_fun),
	color_mapping_legend(ht_list1@ht_list[["meth_diff"]]@matrix_color_mapping, plot = FALSE, title = "Methylation difference", at = c(-0.4, 0, 0.4)),
	Legend(at = c(0, 1), labels = c("no intensity", "high intensity"), title = "Histome modification", col_fun = colorRamp2(c(0, 1), c("white", "purple"))),
	Legend(at = c(0, 0.5, 1), labels = c("high in mature cells", "no difference", "high in embryonic cells"), title = "Difference", 
		col_fun = colorRamp2(c(0, 0.5, 1), c("#3794bf", "#FFFFFF", "#df8640"))),
	Legend(at = c("cluster1", "cluster2", "cluster3", "positive correlation", "negative correlation", "positive difference", "negative difference"), 
		title = "Enrichment lines", legend_gp = gpar(lty = c(1:n_row_cluster, rep(1, 4)), col = c(rep("black", n_row_cluster), "red", "darkgreen", "#df8640", "#3794bf")), type = "lines")
)

lgd_width = grobWidth(lgd_list) + unit(4, "mm")

pushViewport(viewport(x = 0, y = 0.75, width = unit(1, "npc") - lgd_width, height = 0.5, just = "left"))
grid.draw(gb1)
upViewport()

pushViewport(viewport(x = 0, y = 0.25, width = unit(1, "npc") - lgd_width, height = 0.5, just = "left"))
grid.draw(gb2)
upViewport()

pushViewport(viewport(x = unit(1, "npc") - unit(2, "mm"), width = lgd_width - unit(2, "mm"), just = "right"))
grid.draw(lgd_list)
upViewport()
dev.off()
