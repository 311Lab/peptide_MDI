options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))  # Using Tsinghua Mirror


required_packages <- c("yaml", "ggplot2","ggsci","dplyr","fs","ggrepel")

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
input_file <- file.path(config$step8_output_folder, "interface_analysis_results.csv")

output_folder <- file.path(config$output_folder, config$plot_output_folder)
if (!dir_exists(output_folder)) {
  dir_create(output_folder)
}

output_folder <- file.path(config$output_folder, config$plot_output_folder)
if (!dir_exists(output_folder)) {
  dir_create(output_folder)
}

cat("Generated input file path:", input_file, "\n")

# Check if the file exists
if (!file.exists(input_file)) {
  stop("Error: File not found at path:", input_file)
}

# 检查文件是否存在
data <- read.csv(input_file, stringsAsFactors = FALSE)

data = data %>%
        mutate(group = case_when(
                                 interface_dG <= -5 & dSASA_int >= 1000 ~ "A1",
                                 interface_dG <= -5 & dSASA_int < 1000  ~ "A2",
                                 interface_dG > -5 & dSASA_int > 1000   ~ "A3",
                                 TRUE                                   ~ "A4"
                                  ))

p = ggplot(data,aes(x=interface_dG,y=dSASA_int))+
              geom_point(aes(fill = group),shape = 21, color = "black",size=3,alpha=0.8)+
              geom_hline(yintercept = 1000,color="red",linetype=2,linewidth=0.5)+
              geom_vline(xintercept = -5,color="red",linetype=2,linewidth=0.5)+
              geom_text_repel(data = subset(data,group == "A1"),
                             aes(label = description),
                             size= 3,
                  )+
              scale_fill_d3()+
              theme_bw()+
              theme(legend.position = "none",
                    axis.text.x = element_text(size=10,color="black"),
                    axis.title.x = element_text(size=15,color="black"),
                    axis.text.y = element_text(size=10,color="black"),
                    axis.title.y = element_text(size=15,color="black"))


output_file <- file.path(output_folder, "Interface_plot.pdf")
ggsave(output_file, plot = p, width = 10, height = 10)
cat("The image has been successfully saved to the：", output_file, "\n")

