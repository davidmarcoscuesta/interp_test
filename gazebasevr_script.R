# Load necessary libraries
library(data.table)
library(dplyr)
library(DFA)
library(ggplot2)
library(foreach)
library(doMPI)

# --- Step 1: Set up MPI Cluster ---
cl <- startMPIcluster()
registerDoMPI(cl)

# --- Step 2: Filter and Merge Only VRG Files ---
data_path <- "data"  # Directory containing your data files
output_path <- "."   # Save outputs in the same directory as the script

# List all files and filter only VRG tasks
file_list <- list.files(data_path, pattern = "VRG\\.csv$", full.names = TRUE)

# Function to calculate angular and Euclidean velocities
calculate_velocities <- function(data) {
  data <- data %>%
    mutate(
      dt = c(NA, diff(n)),          # Time difference (in ms)
      x_diff = c(NA, diff(x)),      # Change in horizontal gaze direction
      y_diff = c(NA, diff(y)),      # Change in vertical gaze direction
      
      # Angular velocity (degrees per time step)
      ang_vel = atan2(y_diff, x_diff),
      
      # Euclidean velocity (distance per time step)
      euc_vel = sqrt(x_diff^2 + y_diff^2) / (dt / 1000)  # Convert ms to seconds
    )
  
  return(data)
}

# Function to process individual VRG files
process_file <- function(file) {
  # Read the file efficiently
  data <- fread(file)
  
  # Extract metadata from the filename
  file_info <- str_match(basename(file), "S_(\\d+)_S(\\d+)_(\\d+)_(\\w+)\\.csv")
  participant_id <- file_info[2]
  session <- file_info[3]
  task <- file_info[5]
  
  # Add metadata columns to the data
  data[, `:=`(participant_id = participant_id, session = session, task = task)]
  
  # Calculate velocities
  data <- calculate_velocities(data)
  
  return(data)
}

# Process all VRG files and merge into a single dataset
merged_data <- rbindlist(lapply(file_list, process_file))
saveRDS(merged_data, file.path(output_path, "merged_vrg_dataset_with_velocities.rds"))

# --- Step 3: Downsampling and Interpolation ---
# Downsample to half the original frequency (e.g., from 120 Hz to 60 Hz)
downsampled_data <- merged_data %>%
  filter(row_number() %% 2 == 0)

saveRDS(downsampled_data, file.path(output_path, "vrg_downsampled_dataset.rds"))

# Interpolate to restore to the original frequency (from 60 Hz back to 120 Hz)
interpolated_data <- downsampled_data %>%
  mutate(
    x = approx(n, x, n)$y,
    y = approx(n, y, n)$y,
    ang_vel = approx(n, ang_vel, n)$y,
    euc_vel = approx(n, euc_vel, n)$y
  )

saveRDS(interpolated_data, file.path(output_path, "vrg_interpolated_dataset.rds"))

# --- Step 4: RQA and DFA Analyses ---
# Function for RQA analysis
perform_rqa <- function(data, prefix) {
  winsize <- 360
  winstep <- winsize * 0.25
  
  rqa_output <- foreach(trial = unique(data$uniqueID), .combine = 'rbind', .packages = c("dplyr", "wincrqa")) %dopar% {
    tmp <- data %>% filter(uniqueID == trial)
    
    res_win <- wincrqa(
      ts1 = tmp$ang_vel[2:nrow(tmp)], 
      ts2 = tmp$ang_vel[2:nrow(tmp)], 
      windowstep = winstep,
      windowsize = winsize,
      delay = 11, 
      embed = 4, 
      radius = 0.005, 
      rescale = 2, 
      normalize = 0, 
      mindiagline = 2, 
      minvertline = 2, 
      tw = 0, 
      whiteline = FALSE, 
      recpt = FALSE, 
      side = "both", 
      method = "rqa", 
      metric = "euclidean", 
      datatype = "continuous", 
      trend = FALSE
    )
    
    res <- data.frame(
      entr = res_win$ENTR, win_num = res_win$win,
      uniqueID = trial
    )
    return(res)
  }
  
  saveRDS(rqa_output, file.path(output_path, paste0(prefix, "_rqa_output.rds")))
}

# Function for DFA analysis
perform_dfa <- function(data, prefix) {
  dfa_results <- foreach(subject_id = unique(data$participant_id), .combine = 'c', .packages = c("DFA", "dplyr")) %:%
    foreach(subject_trial = unique(data$trial), .combine = 'rbind') %dopar% {
      subject_data <- data %>% 
        filter(participant_id == subject_id & trial == subject_trial)
      
      if (nrow(subject_data) == 0) return(NULL)
      
      mean_ang_vel <- mean(subject_data$ang_vel, na.rm = TRUE)
      subject_data <- subject_data[-1, ] %>%
        mutate(mean_centered_ang_vel = ang_vel - mean_ang_vel,
               cumulative_sum = cumsum(mean_centered_ang_vel))
      
      y <- subject_data$cumulative_sum
      n <- length(y)
      window_size <- 100
      step_size <- 10
      
      subject_trial_results <- list()
      
      for (start in seq(1, n - window_size, by = step_size)) {
        end <- start + window_size - 1
        if (end > n) break
        window_data <- subject_data[start:end, ]
        y_window <- window_data$cumulative_sum
        dfa_result <- DFA(y_window, scale = 2^(1/8), box_size = 4, m = 1)
        box_sizes <- dfa_result[, "box"]
        dfa_values <- dfa_result[, "DFA"]
        
        log_box_sizes <- log(box_sizes)
        log_dfa_values <- log(dfa_values)
        fit <- lm(log_dfa_values ~ log_box_sizes)
        dfa_alpha <- coef(fit)[2]
        
        if (!is.na(dfa_alpha)) {
          subject_trial_results[[length(subject_trial_results) + 1]] <- data.frame(
            subject_id = subject_id,
            trial = subject_trial,
            start_time = start,
            dfa_alpha = dfa_alpha
          )
        }
      }
      
      if (length(subject_trial_results) > 0) {
        do.call(rbind, subject_trial_results)
      } else {
        NULL
      }
    }
  
  saveRDS(dfa_results, file.path(output_path, paste0(prefix, "_dfa_results.rds")))
}

# Perform RQA and DFA on original data
perform_rqa(merged_data, "vrg_original")
perform_dfa(merged_data, "vrg_original")

# Perform RQA and DFA on interpolated data
perform_rqa(interpolated_data, "vrg_interpolated")
perform_dfa(interpolated_data, "vrg_interpolated")

# --- Step 5: Stop MPI Cluster ---
closeCluster(cl)
mpi.quit()
