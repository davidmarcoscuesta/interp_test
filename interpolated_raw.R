library(iterators)
library(foreach)
library(Rmpi)
library(doMPI)
library(data.table)
library(pracma)

# Initialize the MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Load the dataset (adjust the file path as needed)
merged_df <- fread("data/vrg_raw_cleaned.csv")  

# Check for required columns
if (!all(c("participant_id", "session", "x", "y", "dt") %in% colnames(merged_df))) {
  stop("Required columns are missing from the dataset")
}

# Filter out rows with missing participant_id or session
merged_df <- merged_df %>% filter(!is.na(participant_id), !is.na(session))

# Define target delta time for 250 Hz
target_dt <- 1 / 250

# Split the data by participant/session
data_chunks <- split(merged_df, list(merged_df$participant_id, merged_df$session))

# Parallel interpolation using MPI
results <- foreach(chunk = data_chunks, .combine = rbind, .packages = c("data.table", "pracma")) %dopar% {
  tryCatch({
    cur <- as.data.table(chunk)
    if (nrow(cur) < 2) return(NULL)  # Skip small chunks
    
    cur[, t := cumsum(dt)]
    if (any(is.na(t))) return(NULL)  # Skip if cumulative time has issues
    
    interp_points <- seq(from = min(t, na.rm = TRUE), to = max(t, na.rm = TRUE), by = target_dt)
    if (length(interp_points) < 2) return(NULL)  # Skip invalid interpolation
    
    interp_x <- pracma::interp1(x = t, y = x, xi = interp_points, method = "linear")
    interp_y <- pracma::interp1(x = t, y = y, xi = interp_points, method = "linear")
    
    data.table(
      participant_id = unique(cur$participant_id),
      session = unique(cur$session),
      t = interp_points,
      x = interp_x,
      y = interp_y
    )
  }, error = function(e) {
    cat("Error in chunk: ", conditionMessage(e), "\n")
    NULL
  })
}

# Combine all results into one data.table
interpolated_df <- rbindlist(results)

# Save the interpolated data
fwrite(interpolated_df, "interpolated_merged_df.csv")

# Shut down the MPI cluster
closeCluster(cl)
mpi.quit()

cat("Interpolation complete! Results saved to 'interpolated_merged_df.csv'\n")
