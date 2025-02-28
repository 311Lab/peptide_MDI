options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))  


required_packages <- c("yaml", "ggplot2","ggsci","dplyr","fs","pheatmap","curl","httr","plotly")

# Iterate through each package, checking if it is installed, and installing it if not
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

invisible(lapply(required_packages, library, character.only = TRUE))

# Read configuration file
config <- yaml::read_yaml("config.yaml")

# Constructing the input file path
input_file <- file.path(config$step9_output_folder, "prodigy_peptide_summary.csv")

output_folder <- file.path(config$output_folder, config$plot_output_folder)
if (!dir_exists(output_folder)) {
  dir_create(output_folder)
}

cat("Generated input file path:", input_file, "\n")

# Check if the file exists
if (!file.exists(input_file)) {
  stop("Error: File not found at path:", input_file)
}

# Reading CSV files
data <- read.csv(input_file, stringsAsFactors = FALSE)

# Heat map data processing
data1 <- data[, c(1,5:10)]
row.names(data1) <- data1$sequence
data1 <- data1[, -1]

# Draw heat maps and save as PDF
heatmap_pdf <- file.path(output_folder, "heatmap.pdf")
pdf(heatmap_pdf, width = 10, height = 8)
pheatmap(data1,
         cellwidth = 30,
         show_rownames = FALSE)
dev.off()  # Turning off the graphics device

# Scatterplot data processing
data2 <- data[, c(1, 13)]
data2 <- data2 %>% 
  mutate(label = ifelse(Binding_Affinity < -9, "Yes", "No"))

# Generate static ggplot scatterplot
p <- ggplot(data2, aes(x = Filename, y = Binding_Affinity)) +
  geom_point(aes(fill = label), size = 2, shape = 21) +
  scale_fill_aaas() +
  labs(x = "") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.line.x = element_blank(),
        panel.border = element_blank()) +
  scale_y_continuous(expand = c(0.05, 0)) +
  scale_x_discrete(expand = c(0.05, 1))

# Save static scatterplot as PDF
scatter_pdf <- file.path(output_folder, "scatter_plot.pdf")
pdf(scatter_pdf, width = 10, height = 8)
print(p)  # Printing ggplot objects in PDF
dev.off()

# Tips for saving interactive scatterplots as HTML (interactive charts cannot be exported to PDF)
interactive_plot <- ggplotly(p, tooltip = c("x", "y"))
interactive_html <- file.path(output_folder, "interactive_plot.html")
htmlwidgets::saveWidget(interactive_plot, interactive_html)

# Output Save Path Hint
cat("All results have been saved to the folder:", output_folder, "\n",
    "heat map: heatmap.pdf\n",
    "static scatterplot: scatter_plot.pdf\n",
    "Interactive Scatterplot: interactive_plot.html\n")


