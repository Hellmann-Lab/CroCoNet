setwd("data-raw/intermediate_files/")

# download JASPAR2024 vertebrate core
download.file("https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt", "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt")

# download JASPAR2024 vertebrate unvalidated
download.file("https://jaspar.elixir.no/download/data/2024/collections/JASPAR2024_UNVALIDATED_non-redundant_pfms_jaspar.txt", "JASPAR2024_UNVALIDATED_non-redundant_pfms_jaspar.txt")

# download IMAGE motifs
download.file("http://bioinformatik.sdu.dk/solexa/webshare/IMAGE/IMAGE_v1.1.tar.gz", "IMAGE_v1.1.tar.gz")
system("tar -xf IMAGE_v1.1.tar.gz")
system("mv IMAGE/utils/Genename_Motif.txt IMAGE_motifs.txt")
system("rm -r IMAGE IMAGE_v1.1.tar.gz")

# list of TRs
jaspar_core_TRs <- data.table::fread("grep '>' JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
                                     col.names = c("motif_id","SYMBOL")) %>%
  separate_rows(SYMBOL, sep="::") %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort() %>%
  toupper()

jaspar_unvalidated_TRs <- data.table::fread("grep '>' JASPAR2024_UNVALIDATED_non-redundant_pfms_jaspar.txt",
                                            col.names = c("motif_id","SYMBOL")) %>%
  separate_rows(SYMBOL, sep="::") %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort() %>%
  toupper()

image_TRs <- data.table::fread("IMAGE_motifs.txt",
                               col.names = c("SYMBOL", "motif_id", "Evidence")) %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort()

# color scales
eigengene_colors <- rev(c("#650015", "#84001e", "#940022","#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#DAF9D9", "#AFEEEE", "#84C7C8", "#5AA1A3", "#307B7D", "#065558", "#044B4E", "#024244", "#013436"))
pseudotime_colors <- c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF")
species_color_ramp <- scales::colour_ramp(c("#07beb8", "steelblue3", "#2E4172", "#256F5C", "forestgreen", "#83AA3E", "#FFB600",  "#E28100"))
cell_type_color_ramp <- scales::colour_ramp(c("indianred2", "#B24141", "maroon", "palevioletred", "grey65", "grey40", "peachpuff4", "tan"))
dist_colors <- c("#006837", "#1A9850", "#66BD63", "#A6D96A", "#D9EF8B", "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D73027", "#A50026")
module_size_colors <- c("#F1BB7B", "#FD6467", "#5B1A18")
edge_div_colors <- rev(c("#73001a", "#84001e", "#940022","#A50026", "#D73027", "#F46D43", "#E1A97C", "#E0C789", "#E2E2AB", "#C4D9DE", "#9EC0CC", "#7FA2BA", "#4575B4" ,"#313695", "#2c3086", "#272b77", "#222568"))
residual_colors <- c("darkgreen","#2B823A","olivedrab3", "gold", "salmon1", "#AA4139", "red4")

# valid graph layouts
graph_layouts <- c("tree", "sugiyama", "star", "circle", "nicely", "dh", "gem", "graphopt", "grid", "mds", "sphere", "fr", "kk", "drl", "lgl", "auto",  "eigen", "fabric", "linear", "matrix", "stress", "unrooted")

# save as internal data
usethis::use_data(jaspar_core_TRs, jaspar_unvalidated_TRs, image_TRs, eigengene_colors, pseudotime_colors, species_color_ramp, cell_type_color_ramp, dist_colors, module_size_colors, edge_div_colors, residual_colors, graph_layouts, internal = TRUE, overwrite = TRUE)
