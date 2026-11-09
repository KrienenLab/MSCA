## Subcortex group-level taxonomy tree (neurons) with Class / Region / NT tiles and
## cell / cluster count bars.
## Run from the Hierarchical_Clustering folder.
##
## Inputs (both in this folder):
##   Group_metacell_log1p_hvg.csv   group-level metacells (log1p, HVG); builds the tree
##   taxonomy_tree_input.csv        one row per neuronal group: annotation + colors
##                                                   (derived from Supplementary Table S5)

library(ggtree)
library(ggplot2)
library(dplyr)
library(ggnewscale)
library(ggtreeExtra)
library(FactoMineR)
library(tidytree)
library(patchwork)

metacell_path <- "Group_metacell_log1p_hvg.csv"
group_info_path <- "taxonomy_tree_input.csv"
out_dir <- "."

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
label_pad_align <- function(label, fontsize_pt = 8, extra_mm = 0.6, overshoot_dots = 3) {
  w <- function(x) grid::convertWidth(
    grid::grobWidth(grid::textGrob(x, gp = grid::gpar(fontfamily = "Arial", fontsize = fontsize_pt))),
    "mm", valueOnly = TRUE)
  w_lab  <- vapply(label, w, numeric(1))
  n_dots <- pmax(0, ceiling((max(w_lab) + extra_mm - w_lab) / w("·")) + overshoot_dots)
  paste0(label, vapply(n_dots, strrep, "", x = "·"))
}

## horizontal bars drawn in the same y-space as the tree tips
hbar_panel <- function(df, x_col, x_lab, ylim, log10_axis = TRUE, bar_h = 0.6, breaks = NULL) {
  x <- df[[x_col]]
  if (log10_axis) x <- pmax(x, 1)  # zeros shown as 1 on the log scale
  df_rect <- mutate(df, xmin = if (log10_axis) 1 else 0, xmax = x,
                    ymin = y - bar_h / 2, ymax = y + bar_h / 2)

  p <- ggplot(df_rect) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = label)) +
    coord_cartesian(ylim = ylim, expand = FALSE, clip = "off") +
    scale_fill_manual(values = group_color_pal, guide = "none") +
    labs(x = x_lab, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          panel.grid   = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
          plot.margin  = margin(t = 0, r = 4, b = 0, l = 4),
          axis.ticks.x = element_line(linewidth = 0.25, color = "black"),
          axis.ticks.length = grid::unit(2, "pt"),
          axis.line.x  = element_line(linewidth = 0.25, color = "black"))

  if (log10_axis) {
    p + scale_x_continuous(trans = "log10", breaks = breaks,
                           labels = function(x) parse(text = paste0("10^", round(log10(x)))))
  } else {
    p + scale_x_continuous(breaks = breaks)
  }
}

### ============================================================================ ###
###                    Load group table and metacells                            ###
### ============================================================================ ###
taxonomy_group_class <- read.csv(group_info_path, na.strings = "")

## group metacells (rows = groups, columns = genes); keep expressed genes, then neurons only
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
###                              Color palettes                                  ###
### ============================================================================ ###
## named palettes from the "<level>" / "<level>_color_hex" columns; levels without a color are left out (drawn as NA grey)
get_pal <- function(col, color_col) {
  d <- distinct(taxonomy_group_class, .data[[col]], .data[[color_col]]) %>% filter(!is.na(.data[[color_col]]))
  setNames(d[[color_col]], d[[col]])
}
class_color_pal  <- get_pal("Subcortex_Class_v4", "Class_color_hex")
region_color_pal <- get_pal("Group_Dominant_Region", "Region_color_hex")
NT_color_pal     <- get_pal("Group_NT", "NT_color_hex")

## legend order of the neurotransmitter tile
NT_levels <- c("Glut", "GABA", "Glut-GABA", "Glut-Chol", "Chol", "Dopa", "Hist", "Sero")

### ============================================================================ ###
###                     Hierarchical clustering -> tree                          ###
### ============================================================================ ###
tree <- do_hcpc(Group_hvg_metacell)

p_base_tree <- ggtree(as.treedata(tree), branch.length = "none", size = 0.2, layout = "rectangular")

## flip y so the display order is reversed (top-to-bottom)
p_base_tree$data$y <- max(p_base_tree$data$y, na.rm = TRUE) + 1 - p_base_tree$data$y

### ============================================================================ ###
###                       Tip-level annotation table                             ###
### ============================================================================ ###
label_every <- 3  # print a label on every 3rd tip only

group.df <- p_base_tree$data %>%
  filter(isTip) %>%
  select(x, y, label) %>%
  mutate(Subcortex_Group_v4 = gsub("_", " ", label)) %>%
  left_join(taxonomy_group_class, by = "Subcortex_Group_v4") %>%
  arrange(y) %>%
  mutate(
    ## "CS20250428_GROUP_0057: name" -> "57: name"
    uid_int = as.integer(sub(".*_", "", Group_UID)),
    new_label = label_pad_align(Subcortex_Group_v4),
    new_label_uid = ifelse(is.na(uid_int), new_label, paste0(uid_int, ": ", new_label)),
    new_label_sparse = ifelse(row_number() %% label_every == 1, new_label_uid, ""),
    Group_NT = factor(Group_NT, levels = NT_levels)
  )

## attach annotations to the tree
p_base_tree$data <- p_base_tree$data %>%
  mutate(label = gsub("_", " ", label)) %>%
  left_join(group.df[, c("Subcortex_Group_v4", "new_label_sparse", "Subcortex_Class_v4",
                         "Group_Dominant_Region", "Group_NT")],
            by = c("label" = "Subcortex_Group_v4"))

### ============================================================================ ###
###                  Tree with labels, BG marker and annotation tiles            ###
### ============================================================================ ###
p2 <- p_base_tree +
  geom_treescale(linesize = 0.0005, offset = 10, width = 70) +
  geom_tiplab(aes(label = new_label_sparse), size = 2.12, offset = 2, color = "black",
              family = "arial", hjust = 0) +
  geom_point(data = group.df, aes(color = Subcortex_Class_v4_BG),
             position = position_nudge(x = 1), size = 0.8) +
  scale_color_manual(values = c("BG" = "grey", "Non-BG Subcortex" = "#4DBBD5FF")) +
  guides(color = guide_legend(title = "", override.aes = list(size = 6))) +
  ## Neurotransmitter tile
  new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping = aes(fill = Group_NT), width = 2.5, offset = 3.5) +
  scale_fill_manual(name = "Neurotransmitter", values = NT_color_pal) +
  ## Dominant region tile
  new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping = aes(fill = Group_Dominant_Region), width = 2.5, offset = 0.2) +
  scale_fill_manual(name = "Dominant Region", values = region_color_pal) +
  ## Class tile
  new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping = aes(fill = Subcortex_Class_v4), width = 2.5, offset = 0.2) +
  scale_fill_manual(name = "Class", values = class_color_pal)

### ============================================================================ ###
###                                Bar plots                                     ###
### ============================================================================ ###
ylim <- range(group.df$y) + c(-0.5, 0.5)

## bars colored by class (grey if class has no color)
group_color_pal <- class_color_pal[group.df$Subcortex_Class_v4]
group_color_pal[is.na(group_color_pal)] <- "grey"
names(group_color_pal) <- group.df$Subcortex_Group_v4

plot_df <- group.df %>%
  transmute(label = Subcortex_Group_v4, y, n_cells, n_clusters)

g_cells   <- hbar_panel(plot_df, "n_cells", "Number of\nCells", ylim,
                        log10_axis = TRUE, breaks = c(10, 100, 1000, 10000, 100000))
g_cluster <- hbar_panel(plot_df, "n_clusters", "Number of\nClusters", ylim,
                        log10_axis = FALSE, breaks = scales::pretty_breaks(n = 4)(plot_df$n_clusters))

p_tree <- p2 + coord_cartesian(ylim = ylim, expand = FALSE, clip = "off")

p_all <- (p_tree + g_cells + g_cluster +
            plot_layout(widths = c(1, 0.1, 0.1), guides = "collect")) &
  theme(plot.margin = margin(0, 0, 0, 0), legend.position = "right", text = element_text(size = 7))

ggsave(file.path(out_dir, "neuron_taxonomy_tree.svg"),
       p_all, width = 5, height = 8, device = svglite::svglite)
