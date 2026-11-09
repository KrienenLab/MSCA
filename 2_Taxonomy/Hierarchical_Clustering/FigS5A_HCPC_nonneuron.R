## Subcortex group-level taxonomy tree of the non-neuronal groups (incl. Pineal gland) with cell / cluster
## count bars and a marker gene dot plot aligned to the tree.
## Run from the Hierarchical_Clustering folder.
##
## Inputs (all in this folder):
##   Group_metacell_log1p_hvg.csv    group-level metacells (log1p, HVG); also used by Fig1B_HCPC_neuron.R
##   nonneuron_tree_input.csv        one row per non-neuronal group: annotation + color
##                                   (derived from Supplementary Table S5)
##   nonneuron_marker_dotplot.tsv    mean expression and % expressing cells of marker genes per group
##                                   (columns '<group>:mean' and '<group>:percent_expressed')

library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
library(FactoMineR)
library(tidytree)
library(patchwork)

metacell_path <- "Group_metacell_log1p_hvg.csv"
group_info_path <- "nonneuron_tree_input.csv"
marker_path <- "nonneuron_marker_dotplot.tsv"
out_dir <- "."

## marker genes left out of the dot plot
drop_genes <- c("NFIA", "POU3F2", "NR2F1")

### ============================================================================ ###
###                                Functions                                     ###
### ============================================================================ ###
## PCA -> hierarchical clustering on PCs (HCPC); returns the tree as a tidytree tibble
do_hcpc <- function(metacell) {
  res.pca  <- PCA(metacell, ncp = 100, graph = FALSE)
  res.hcpc <- HCPC(res.pca, graph = FALSE, nb.clust = -1)
  as_tibble(ape::as.phylo(res.hcpc$call$t$tree))
}

## pad labels with "·" so that all labels reach the same width (aligned dotted leaders)
label_pad_align <- function(label, fontsize_pt = 8, extra_mm = 0.6, overshoot_dots = 1) {
  w <- function(x) grid::convertWidth(
    grid::grobWidth(grid::textGrob(x, gp = grid::gpar(fontfamily = "Arial", fontsize = fontsize_pt))),
    "mm", valueOnly = TRUE)
  w_lab  <- vapply(label, w, numeric(1))
  n_dots <- pmax(0, ceiling((max(w_lab) + extra_mm - w_lab) / w("·")) + overshoot_dots)
  paste0(label, vapply(n_dots, strrep, "", x = "·"))
}

## horizontal bars drawn in the same y-space as the tree tips
hbar_panel <- function(df, x_col, x_lab, ylim, log10_axis = TRUE, bar_h = 0.6) {
  x <- df[[x_col]]
  if (log10_axis) x <- pmax(x, 1)  # zeros shown as 1 on the log scale
  df_rect <- mutate(df, xmin = if (log10_axis) 1 else 0, xmax = x,
                    ymin = y - bar_h / 2, ymax = y + bar_h / 2)

  p <- ggplot(df_rect) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = label)) +
    scale_fill_manual(values = group_color_pal, na.value = "grey80", guide = "none") +
    coord_cartesian(ylim = ylim, expand = FALSE, clip = "off") +
    labs(x = x_lab, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid   = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
          plot.margin  = margin(2, 0, 2, 0),
          axis.ticks.x = element_line(linewidth = 0.25),
          axis.ticks.length = unit(1.2, "mm"),
          axis.line.x  = element_line(linewidth = 0.25, color = "black"))

  if (log10_axis) {
    p + scale_x_continuous(trans = "log10", breaks = scales::log_breaks(n = 4),
                           labels = scales::trans_format("log10", scales::math_format(10^.x)),
                           expand = expansion(mult = c(0, 0.03)))
  } else {
    p + scale_x_continuous(expand = expansion(mult = c(0, 0.03)))
  }
}

### ============================================================================ ###
###                    Load group table and metacells                            ###
### ============================================================================ ###
taxonomy_group_class <- read.csv(group_info_path)

## group metacells (rows = groups, columns = genes); keep expressed genes, then the non-neuronal groups
Group_hvg_metacell <- read.csv(metacell_path, row.names = 1)
Group_hvg_metacell <- Group_hvg_metacell[, colSums(Group_hvg_metacell) > 10]
Group_hvg_metacell <- Group_hvg_metacell[rownames(Group_hvg_metacell) %in% taxonomy_group_class$Subcortex_Group_v4, ]

if (!setequal(rownames(Group_hvg_metacell), taxonomy_group_class$Subcortex_Group_v4)) {
  stop("Group metacells and group table do not match: ",
       paste(c(setdiff(rownames(Group_hvg_metacell), taxonomy_group_class$Subcortex_Group_v4),
               setdiff(taxonomy_group_class$Subcortex_Group_v4, rownames(Group_hvg_metacell))),
             collapse = ", "))
}

### ============================================================================ ###
###                     Hierarchical clustering -> tree                          ###
### ============================================================================ ###
tree <- do_hcpc(Group_hvg_metacell)

p_base_tree <- ggtree(as.treedata(tree), branch.length = "none", size = 0.3, layout = "rectangular") +
  ggtree::theme_tree()

### ============================================================================ ###
###                       Tip-level annotation table                             ###
### ============================================================================ ###
group.df <- p_base_tree$data %>%
  filter(isTip) %>%
  select(label, x, y) %>%
  mutate(Subcortex_Group_v4 = gsub("_", " ", label),
         new_label = label_pad_align(Subcortex_Group_v4)) %>%
  left_join(taxonomy_group_class, by = "Subcortex_Group_v4")

## attach padded labels to the tree
p_base_tree$data <- p_base_tree$data %>%
  mutate(label = gsub("_", " ", label)) %>%
  left_join(group.df[, c("Subcortex_Group_v4", "new_label")], by = c("label" = "Subcortex_Group_v4"))

p_tree <- p_base_tree +
  geom_tiplab(aes(label = new_label), size = 3, offset = 0.02, color = "black", family = "arial", hjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(2, 75, 2, 2))  # room on the right for the tip labels

### ============================================================================ ###
###                                Bar plots                                     ###
### ============================================================================ ###
ylim <- range(group.df$y) + c(-0.5, 0.5)

group_color_pal <- setNames(group.df$Group_color_hex, group.df$Subcortex_Group_v4)
plot_df <- transmute(group.df, label = Subcortex_Group_v4, y, n_cells, n_clusters)

g_cells   <- hbar_panel(plot_df, "n_cells", "Number of\nCells", ylim, log10_axis = TRUE)
g_cluster <- hbar_panel(plot_df, "n_clusters", "Number of\nClusters", ylim, log10_axis = FALSE)

### ============================================================================ ###
###                     Marker gene dot plot (aligned to tree)                   ###
### ============================================================================ ###
marker_wide <- read.delim(marker_path, check.names = FALSE, stringsAsFactors = FALSE)
colnames(marker_wide)[1] <- "gene"
marker_wide <- filter(marker_wide, !gene %in% drop_genes)

## columns "<group>:mean" and "<group>:percent_expressed" -> one row per gene and group
marker_long <- marker_wide %>%
  pivot_longer(cols = -gene, names_to = c("label", "metric"), names_sep = ":", values_to = "value") %>%
  pivot_wider(names_from = metric, values_from = value)

dot_df <- group.df %>%
  transmute(label = Subcortex_Group_v4, y) %>%
  left_join(marker_long, by = "label") %>%
  mutate(gene = factor(gene, levels = marker_wide$gene),
         percent_expressed = pmax(percent_expressed, 0)) %>%
  filter(!is.na(mean), mean > 0)

## order the genes by the expression-weighted position along the tree, giving a diagonal pattern
gene_order <- dot_df %>%
  group_by(gene) %>%
  summarise(y_cm = weighted.mean(y, w = (mean + 1e-6) * (percent_expressed + 1e-6)), .groups = "drop") %>%
  arrange(y_cm) %>%
  pull(gene) %>%
  as.character()
dot_df <- mutate(dot_df, gene = factor(gene, levels = gene_order))

x_pad <- 0.40  # the frame extends this far beyond the first / last gene

g_markers <- ggplot(dot_df, aes(x = gene, y = y)) +
  geom_point(aes(size = percent_expressed, color = mean), stroke = 0.25) +
  scale_size_continuous(range = c(0.2, 3.2), name = "% expressed") +
  scale_color_distiller(palette = "Reds", direction = 1, name = "Mean expr") +
  scale_x_discrete(expand = expansion(add = 0.25)) +
  scale_y_continuous(limits = ylim, expand = c(0, 0)) +
  coord_cartesian(ylim = ylim, expand = FALSE, clip = "off") +
  annotate("rect", xmin = 0.5 - x_pad, xmax = nlevels(dot_df$gene) + 0.5 + x_pad,
           ymin = ylim[1], ymax = ylim[2], fill = NA, color = "black", linewidth = 0.4) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks.x = element_line(linewidth = 0.25),
        axis.ticks.length = unit(1.2, "mm"),
        axis.line.x  = element_line(linewidth = 0.25, color = "black"),
        panel.grid   = element_blank(),
        panel.border = element_blank(),
        plot.margin  = margin(2, 0, 2, 0),
        legend.position = "right",
        legend.key.width  = unit(3, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.spacing.y  = unit(1, "mm"))

### ============================================================================ ###
###                              Combine and save                                ###
### ============================================================================ ###
p_all <- p_tree + g_cells + g_cluster + g_markers +
  plot_layout(widths = c(0.15, 0.05, 0.05, 0.16), guides = "collect")

ggsave(file.path(out_dir, "nonneuron_taxonomy_tree.svg"), p_all, width = 5.75, height = 3, bg = "transparent")
