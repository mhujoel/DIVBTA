# Load required libraries
library(haven)        # For reading SAS files
library(fitdistrplus) # For fitting gamma distribution
library(sn)           # Loaded as per provided code, though not used here

# Get the current Windows username
username <- Sys.getenv("USERNAME")

# Construct the dynamic file path for input
path <- file.path("C:", "Users", username, "Dropbox", "Carlisle", "R-diabetes", "SAS datafiles")
data_file <- file.path(path, "DPPT.sas7bdat")

# Read the dataset
A1c_data <- read_sas(data_file)

# Fit gamma distribution
set.seed(123)
fit_gamma <- fitdist(A1c_data$hba1c_6, "gamma", method = "mle")

# Extract shape and scale
summary(fit_gamma)
shape_est <- fit_gamma$estimate["shape"]
scale_est <- 1 / fit_gamma$estimate["rate"]
rate_est <- 1 / scale_est
cat("Estimated shape:", shape_est, "\n")
cat("Estimated scale:", scale_est, "\n")
cat("Estimated rate:", rate_est, "\n")
# Define output directory (as specified)
output_dir <- "C:/Users/hujoe/Dropbox/Carlisle/R-diabetes/Submission_manuscripts/LnCVR-SAS/LNCVR-publication_SAS/Revision/export Excel Files for Github"

# Ensure directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define JPEG file path
jpeg_file <- file.path(output_dir, "HbA1c_Gamma_Fit_Histogram.jpg")

# Open JPEG device
jpeg(filename = jpeg_file, width = 800, height = 600, units = "px", res = 120)

# Create histogram with NO main title
hist(A1c_data$hba1c_6, breaks = 30, prob = TRUE, 
     xlab = "Six-month HbA1c measure (%)", 
     main = "",  # Explicitly suppress default title
     col = "lightblue", border = "black",
     cex.lab = 1.2, cex.axis = 1.1)

# Add fitted gamma curve
curve(dgamma(x, shape = shape_est, scale = scale_est), 
      from = 0, to = max(A1c_data$hba1c_6, na.rm = TRUE) * 1.2, 
      add = TRUE, col = "red", lwd = 2.5)

# Add two custom titles using mtext
mtext("Diabetes and Periodontal Therapy Trial (JAMA, 2013)", 
      side = 3, line = 3, cex = 1.2, font = 1)
mtext("Fit of the Gamma Distribution to 6-month HbA1c measures", 
      side = 3, line = 1.6, cex = 1.2, font = 1)


# Close device
dev.off()

