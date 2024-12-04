increment_and_delete <- function(file_path) {
  # Check if the file exists
  if (!file.exists(file_path)) {
    message("File does not exist: ", file_path)
    return(FALSE)
  }
  
  # Extract the directory, base name, and extension
  file_dir <- dirname(file_path)
  file_base <- tools::file_path_sans_ext(basename(file_path))
  file_ext <- tools::file_ext(file_path)
  
  # Match the number at the end of the file base name
  match <- regexpr("\\d+$", file_base)
  
  if (match > 0) {
    # If a number exists, increment it
    number <- as.numeric(substr(file_base, match, nchar(file_base))) + 1
    new_file_base <- paste0(substr(file_base, 1, match - 1), number)
  } else {
    # If no number exists, append 1
    new_file_base <- paste0(file_base, "1")
  }
  
  # Construct the new file name
  new_file_path <- file.path(file_dir, paste0(new_file_base, ifelse(file_ext == "", "", paste0(".", file_ext))))
  
  # Delete the original file
  file.remove(file_path)
  message("Deleted file: ", file_path)
  
  # Return the new file name
  return(new_file_path)
}

delete_if_exists <- function(file_path) {
  # Check if the file exists
  if (file.exists(file_path)) {
    # Delete the file
    if (file.remove(file_path)) {
      message("Deleted file: ", file_path)
      return(TRUE)
    } else {
      warning("Failed to delete file: ", file_path)
      return(FALSE)
    }
  } else {
    message("File does not exist: ", file_path)
    return(FALSE)
  }
}
