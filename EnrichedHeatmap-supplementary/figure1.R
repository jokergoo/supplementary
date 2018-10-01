
library(EnrichedHeatmap)
library(circlize)

# figure A
gd_figure_A = grid.grabExpr({
	grid.lines(c(0.1, 0.9), c(0.3, 0.3), gp = gpar(col = "red", lwd = 4))
	grid.text("window", 0.5, unit(0.3, "npc") - unit(4, "mm"), just = "top", gp = gpar(col = "red"))
	grid.lines(c(0, 0.2), c(0.6, 0.6))
	grid.lines(c(0.3, 0.5), c(0.6, 0.6))
	grid.lines(c(0.7, 1), c(0.6, 0.6))
	grid.lines(c(0.1, 0.2), c(0.6, 0.6), gp = gpar(lwd = 4))
	grid.text("x1", 0.15, unit(0.6, "npc") - unit(4, "mm"))
	grid.lines(c(0.3, 0.5), c(0.6, 0.6), gp = gpar(lwd = 4))
	grid.text("x3", 0.4, unit(0.6, "npc") - unit(4, "mm"))
	grid.lines(c(0.7, 0.9), c(0.6, 0.6), gp = gpar(lwd = 4))
	grid.text("x5", 0.8, unit(0.6, "npc") - unit(4, "mm"))
	grid.lines(c(0.15, 0.35), c(0.45, 0.45), gp = gpar(lwd = 4))
	grid.text("x2", 0.25, unit(0.45, "npc") - unit(4, "mm"))
	grid.lines(c(0.6, 0.75), c(0.45, 0.45), gp = gpar(lwd = 4))
	grid.text("x4", 0.675, unit(0.45, "npc") - unit(4, "mm"))
	grid.lines(c(0.1, 0.1), c(0, 1), gp = gpar(lty = 2, col = "grey"))
	grid.lines(c(0.9, 0.9), c(0, 1), gp = gpar(lty = 2, col = "grey"))
	grid.text("genomic signal regions", 0.5, unit(0.6, "npc") + unit(4, "mm"), just = "bottom")
})


axis_name = c("-5kb", "TSS", "5kb")
# figure B
load(system.file("extdata", "chr21_test_data.RData", package = "EnrichedHeatmap"))
tss = promoters(genes, upstream = 0, downstream = 1)
meth = readRDS("roadmap_lung_wgbs_chr21_meth.rds")
mat2 = normalizeToMatrix(meth, tss, value_column = "meth", mean_mode = "absolute",
    extend = 5000, w = 50, background = NA)
meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
gd_figure_B1 = grid.grabExpr(draw(EnrichedHeatmap(mat2, col = meth_col_fun, name = "Methylation", 
	top_annotation = HeatmapAnnotation(enriched = anno_enriched()),
	row_title = "original methylation", use_raster = TRUE, heatmap_legend_param = list(legend_direction = "horizontal"),
		axis_name = axis_name), 
		heatmap_legend_side = "bottom"))

mat2 = normalizeToMatrix(meth, tss, value_column = "meth", mean_mode = "absolute",
    extend = 5000, w = 50, background = NA, smooth = TRUE)
gd_figure_B2 = grid.grabExpr(draw(EnrichedHeatmap(mat2, col = meth_col_fun, name = "Methylation", 
		top_annotation = HeatmapAnnotation(enriched = anno_enriched()),
		row_title = "smoothed methylation", use_raster = TRUE, heatmap_legend_param = list(legend_direction = "horizontal"),
		axis_name = axis_name), 
		heatmap_legend_side = "bottom"))

# figure C

load("roadmap_normalized_matrices.RData")

gb_figure_C1 = grid.grabExpr(draw(EnrichedHeatmap(mat_neg_cr, name = "negCR", col = c("white", "darkgreen"),
	top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = "darkgreen"))),
	row_title = "by default enriched scores", use_raster = TRUE, heatmap_legend_param = list(legend_direction = "horizontal", nrow = 1, border = "black"),
		axis_name = axis_name),
	heatmap_legend_side = "bottom"))
gb_figure_C2 = grid.grabExpr(draw(EnrichedHeatmap(mat_neg_cr, name = "negCR", col = c("white", "darkgreen"),
	top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = "darkgreen"),)),
	cluster_rows = TRUE, row_title = "by hierarchcal clustering + Euclidean distance", use_raster = TRUE, 
	heatmap_legend_param = list(legend_direction = "horizontal", nrow = 1, border = "black"),
		axis_name = axis_name),
	heatmap_legend_side = "bottom"))
gb_figure_C3 = grid.grabExpr(draw(EnrichedHeatmap(mat_neg_cr, name = "negCR", col = c("white", "darkgreen"),
	top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = "darkgreen"))),
	cluster_rows = TRUE, clustering_distance_rows = dist_by_closeness, 
	row_title = "by hierarchcal clustering + closeness distance", use_raster = TRUE, 
	heatmap_legend_param = list(legend_direction = "horizontal", nrow = 1, border = "black"),
		axis_name = axis_name),
	heatmap_legend_side = "bottom"))



pdf("figure1.pdf", width = 14, height = 8)
grid.newpage()
pushViewport(viewport(y = unit(0, "npc"), height = unit(1, "npc") - unit(2, "cm"), just = "bottom",
	layout = grid.layout(nr = 1, nc = 9,
		width = unit.c(unit(1, "null"), unit(1, "cm"), unit(c(1, 1), "null"), 
			unit(1, "cm"), unit(c(1.1, 1.1, 1.1), "null"), unit(1, "cm")))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
pushViewport(viewport(height = 0.3, y = 1, just = "top"))
grid.draw(gd_figure_A)
popViewport()
grid.text("A", y = unit(1, "npc") + unit(1, "cm"), gp = gpar(fontsize = 20))
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
grid.draw(gd_figure_B1)
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
grid.draw(gd_figure_B2)
grid.text("B", x = 0, y = unit(1, "npc") + unit(1, "cm"), gp = gpar(fontsize = 20))
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 6))
grid.draw(gb_figure_C1)
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 7))
grid.draw(gb_figure_C2)
grid.text("C", y = unit(1, "npc") + unit(1, "cm"), gp = gpar(fontsize = 20))
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 8))
grid.draw(gb_figure_C3)
popViewport()

dev.off()

