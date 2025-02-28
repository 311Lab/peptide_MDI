options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))  


required_packages <- c("yaml", "ggplot2","ggsci","cowplot","fs")

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
input_file <- file.path(config$step10_output_folder, "output.csv")

output_folder <- file.path(config$output_folder, config$plot_output_folder)

print(paste("Output folder:", output_folder))


if (is.null(output_folder) || output_folder == "") {
    stop("Error: output_folder is not defined in the configuration file.")
}

if (!dir_exists(output_folder)) {
     dir_create(output_folder)
    print(paste("Directory created:", output_folder))
}



# Print the path to the generated file
cat("Generated input file path:", input_file, "\n")

# Check if the file exists
if (!file.exists(input_file)) {
  stop("Error: File not found at path:", input_file)
}

# Reading CSV files
data <- read.csv(input_file, stringsAsFactors = FALSE)

# Calculate the number of amino acids per peptide
data$AAs <- nchar(data$sequence)
data$AAs = as.factor(data$AAs)

# View Results
print(data)


p1 <- ggplot(data, aes(x=AAs, y=binding_energy)) +
            geom_jitter(width = 0.32,height = 0.01, aes(fill = AAs),alpha=0.9, size=3, shape=21, color="black") +
            scale_fill_d3() +
            labs(y="Docking energy (kcal/mol)", x="") +
            theme_bw()+
            theme(legend.position = "none")

p2 <- ggplot(data, aes(x=binding_energy)) +
            geom_density(alpha=0.5,fill = "pink") +
            labs(x="", y="") +
            theme_bw() +
            theme(legend.position = "none",
                  axis.title.y = element_blank(),
                  axis.text.y = element_blank())+
            coord_flip()

p <- plot_grid(p1, p2, ncol=2, rel_widths=c(5, 1),axis = "tblr",align = "hv")

output_file <- file.path(output_folder, "Docking_plot.pdf")
ggsave(output_file, plot = p, width = 8, height = 6)
cat("The image has been successfully saved to the：", output_file, "\n")

