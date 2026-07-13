# ============================================================
# Tn and T combinatorial subtraction from the precursor ions (MGF files)
#
# Diagnostic-dependent decision tree:
#
# 204 present, 366 absent:
#   Apply only Tn subtraction series
#
# 204 present, 366 present:
#   Apply T and mixed Tn/T subtraction series
#
# 204 absent:
#   Do not process and exclude the corresponding MSMS spectrum
# ============================================================


# ============================================================
# 1. USER SETTINGS
# ============================================================
# Select input MGF interactively
input_file <- file.choose()

# Diagnostic sugar oxonium fragment-ion m/z values
diagnostic_204_mz <- 204.0867
diagnostic_366_mz <- 366.1395


# Diagnostic-ion tolerances in Da (+/- 0.02 Da)

diagnostic_204_tolerance <- 0.02
diagnostic_366_tolerance <- 0.02

# Minimum remaining precursor mass
# Precursor_mz * charge - sugar_mass >= 300

minimum_remaining_mass <- 300


# Decimal places used for modified PEPMASS
pepmass_digits <- 8

# Add MGF files header
add_mass_header <- TRUE


# ============================================================
# 2. SUGAR DELTA MASS TABLE
# ============================================================

# Each row produces one output MGF file.
#
# diagnostic_group determines file group and sugar increment


sugar_deltas <- data.frame(
  
  label = c(
    "(1G)_1xTn",
    "(2G)_2xTn",
    "(3G)_3xTn",
    "(1G)_1xT",
    "(2G)_1xTn_1xT",
    "(3G)_2xTn_1xT",
    "(2G)_2xT",
    "(3G)_1xTn_2xT",
    "(3G)_3xT",
    "(4G)_4xTn",
    "(4G)_4xT",
    "(5G)_5xTn",
    "(5G)_5xT"
  ),
  
  delta_mass = c(
    203.07937,    # 1 x Tn
    406.15874,    # 2 x Tn
    609.23811,    # 3 x Tn
    365.13219,    # 1 x T
    568.21157,    # 1 x Tn + 1 x T
    771.29094,    # 2 x Tn + 1 x T
    730.26440,    # 2 x T
    933.34377,    # 1 x Tn + 2 x T
    1095.39660,   # 3 x T
    812.31748,    # 4 x Tn
    1460.52880,   # 4 x T
    1015.39685,   # 5 x Tn
    1825.66095    # 5 x T
  ),
  
  diagnostic_group = c(
    "204_only",       # 1 x Tn
    "204_only",       # 2 x Tn
    "204_only",       # 3 x Tn
    "204_and_366",    # 1 x T
    "204_and_366",    # 1 x Tn + 1 x T
    "204_and_366",    # 2 x Tn + 1 x T
    "204_and_366",    # 2 x T
    "204_and_366",    # 1 x Tn + 2 x T
    "204_and_366",    # 3 x T
    "204_only",       # 4 x Tn
    "204_and_366",    # 4 x T
    "204_only",       # 5 x Tn
    "204_and_366"     # 5 x T
  ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# 3. FUNCTIONS
# ============================================================


# Split the MGF file into separate spectra.
#
# Each spectrum starts with BEGIN IONS and ends with END IONS.

split_mgf_spectra <- function(lines) {
  
  begin_positions <- grep(
    "^\\s*BEGIN IONS\\s*$",
    lines,
    ignore.case = TRUE
  )
  
  end_positions <- grep(
    "^\\s*END IONS\\s*$",
    lines,
    ignore.case = TRUE
  )
  
  if (length(begin_positions) == 0) {
    stop("No BEGIN IONS entries were found.")
  }
  
  if (length(begin_positions) != length(end_positions)) {
    stop(
      paste0(
        "MGF structure error: found ",
        length(begin_positions),
        " BEGIN IONS lines and ",
        length(end_positions),
        " END IONS lines."
      )
    )
  }
  
  if (any(end_positions < begin_positions)) {
    stop("MGF structure error: END IONS occurs before BEGIN IONS.")
  }
  
  Map(
    function(start_position, end_position) {
      lines[start_position:end_position]
    },
    begin_positions,
    end_positions
  )
}


# Extract precursor m/z from the PEPMASS line.
#
# Handles:
#
# PEPMASS=m/z Intensity

get_pepmass <- function(spectrum) {
  
  position <- grep(
    "^\\s*PEPMASS\\s*=",
    spectrum,
    ignore.case = TRUE
  )
  
  if (length(position) == 0) {
    return(NA_real_)
  }
  
  line <- spectrum[position[1]]
  
  value_part <- sub(
    "^\\s*PEPMASS\\s*=\\s*",
    "",
    line,
    ignore.case = TRUE
  )
  
  fields <- strsplit(
    trimws(value_part),
    "\\s+"
  )[[1]]
  
  suppressWarnings(
    as.numeric(fields[1])
  )
}


# Extract charge from the CHARGE line.
#
# Handles examples such as:
#
# CHARGE=3+
# CHARGE=3
# CHARGE=+3

get_charge <- function(spectrum) {
  
  position <- grep(
    "^\\s*CHARGE\\s*=",
    spectrum,
    ignore.case = TRUE
  )
  
  if (length(position) == 0) {
    return(NA_integer_)
  }
  
  line <- spectrum[position[1]]
  
  value_part <- sub(
    "^\\s*CHARGE\\s*=\\s*",
    "",
    line,
    ignore.case = TRUE
  )
  
  charge_match <- regexpr(
    "[0-9]+",
    value_part
  )
  
  if (charge_match[1] == -1) {
    return(NA_integer_)
  }
  
  charge_text <- regmatches(
    value_part,
    charge_match
  )
  
  suppressWarnings(
    as.integer(charge_text)
  )
}


# Extract all fragment-ion m/z values from a spectrum.
#
# Only lines containing two numerical columns are treated as
# fragment peaks.

get_fragment_mz <- function(spectrum) {
  
  peak_pattern <- paste0(
    "^\\s*",
    "[0-9]+(?:\\.[0-9]+)?",
    "(?:[eE][+-]?[0-9]+)?",
    "\\s+",
    "[0-9]+(?:\\.[0-9]+)?",
    "(?:[eE][+-]?[0-9]+)?",
    "\\s*$"
  )
  
  peak_lines <- grep(
    peak_pattern,
    spectrum,
    value = TRUE,
    perl = TRUE
  )
  
  if (length(peak_lines) == 0) {
    return(numeric(0))
  }
  
  peak_fields <- strsplit(
    trimws(peak_lines),
    "\\s+"
  )
  
  mz_values <- vapply(
    peak_fields,
    function(fields) {
      suppressWarnings(
        as.numeric(fields[1])
      )
    },
    numeric(1)
  )
  
  mz_values[is.finite(mz_values)]
}


# Test whether a diagnostic ion is present within the defined
# absolute tolerance.

contains_diagnostic_ion <- function(
    fragment_mz,
    target_mz,
    tolerance
) {
  
  if (length(fragment_mz) == 0) {
    return(FALSE)
  }
  
  any(
    abs(fragment_mz - target_mz) <= tolerance,
    na.rm = TRUE
  )
}


# Determine whether 204 and 366 are present in one spectrum.

get_diagnostic_status <- function(
    spectrum,
    mz_204,
    tolerance_204,
    mz_366,
    tolerance_366
) {
  
  fragment_mz <- get_fragment_mz(spectrum)
  
  has_204 <- contains_diagnostic_ion(
    fragment_mz = fragment_mz,
    target_mz = mz_204,
    tolerance = tolerance_204
  )
  
  has_366 <- contains_diagnostic_ion(
    fragment_mz = fragment_mz,
    target_mz = mz_366,
    tolerance = tolerance_366
  )
  
  c(
    has_204 = has_204,
    has_366 = has_366
  )
}


# Replace only the precursor m/z in the PEPMASS line.
#
# Any precursor intensity after the m/z value is retained.

replace_pepmass <- function(
    spectrum,
    new_pepmass,
    digits = 8
) {
  
  position <- grep(
    "^\\s*PEPMASS\\s*=",
    spectrum,
    ignore.case = TRUE
  )
  
  if (length(position) == 0) {
    return(spectrum)
  }
  
  position <- position[1]
  original_line <- spectrum[position]
  
  value_part <- sub(
    "^\\s*PEPMASS\\s*=\\s*",
    "",
    original_line,
    ignore.case = TRUE
  )
  
  fields <- strsplit(
    trimws(value_part),
    "\\s+"
  )[[1]]
  
  new_pepmass_text <- formatC(
    new_pepmass,
    format = "f",
    digits = digits
  )
  
  # Remove unnecessary terminal zeros
  new_pepmass_text <- sub(
    "0+$",
    "",
    new_pepmass_text
  )
  
  new_pepmass_text <- sub(
    "\\.$",
    "",
    new_pepmass_text
  )
  
  if (length(fields) >= 2) {
    
    remaining_fields <- paste(
      fields[-1],
      collapse = " "
    )
    
    new_line <- paste0(
      "PEPMASS=",
      new_pepmass_text,
      " ",
      remaining_fields
    )
    
  } else {
    
    new_line <- paste0(
      "PEPMASS=",
      new_pepmass_text
    )
  }
  
  spectrum[position] <- new_line
  
  spectrum
}


# Write spectra to one MGF output file.

write_mgf_file <- function(
    spectra,
    output_file,
    add_mass_header = TRUE
) {
  
  connection <- file(
    output_file,
    open = "wt"
  )
  
  on.exit(
    close(connection),
    add = TRUE
  )
  
  if (add_mass_header) {
    writeLines(
      "MASS=Monoisotopic",
      connection
    )
  }
  
  for (spectrum in spectra) {
    
    writeLines(
      spectrum,
      connection
    )
    
    writeLines(
      "",
      connection
    )
  }
}


# ============================================================
# 4. MAIN PROCESSING FUNCTION
# ============================================================

process_mgf_sugar_subtractions <- function(
    input_file,
    sugar_deltas,
    diagnostic_204_mz = 204.0867,
    diagnostic_204_tolerance = 0.02,
    diagnostic_366_mz = 366.1395,
    diagnostic_366_tolerance = 0.02,
    minimum_remaining_mass = 300,
    pepmass_digits = 8,
    add_mass_header = TRUE
) {
  
  if (!file.exists(input_file)) {
    stop(
      paste0(
        "Input file does not exist: ",
        input_file
      )
    )
  }
  
  required_columns <- c(
    "label",
    "delta_mass",
    "diagnostic_group"
  )
  
  if (!all(required_columns %in% names(sugar_deltas))) {
    stop(
      paste0(
        "sugar_deltas must contain the columns: ",
        "label, delta_mass and diagnostic_group."
      )
    )
  }
  
  allowed_groups <- c(
    "204_only",
    "204_and_366"
  )
  
  if (
    any(
      !sugar_deltas$diagnostic_group %in% allowed_groups
    )
  ) {
    stop(
      paste0(
        "diagnostic_group must be either ",
        "'204_only' or '204_and_366'."
      )
    )
  }
  
  if (anyDuplicated(sugar_deltas$label)) {
    stop("The sugar-delta labels must be unique.")
  }
  
  if (any(!is.finite(sugar_deltas$delta_mass))) {
    stop("All delta masses must be finite numeric values.")
  }
  
  
  message("Reading input MGF: ", input_file)
  
  all_lines <- readLines(
    input_file,
    warn = FALSE
  )
  
  spectra <- split_mgf_spectra(
    all_lines
  )
  
  number_of_spectra <- length(spectra)
  
  message(
    "Number of spectra found: ",
    number_of_spectra
  )
  
  
  # ----------------------------------------------------------
  # Detect 204 and 366 separately in every spectrum
  # ----------------------------------------------------------
  
  diagnostic_matrix <- t(
    vapply(
      spectra,
      get_diagnostic_status,
      logical(2),
      mz_204 = diagnostic_204_mz,
      tolerance_204 = diagnostic_204_tolerance,
      mz_366 = diagnostic_366_mz,
      tolerance_366 = diagnostic_366_tolerance
    )
  )
  
  has_204 <- diagnostic_matrix[, "has_204"]
  has_366 <- diagnostic_matrix[, "has_366"]
  
  
  # Spectra used for Tn-only subtraction
  spectra_204_only <- has_204 & !has_366
  
  # Spectra used for T and mixed Tn/T subtraction
  spectra_204_and_366 <- has_204 & has_366
  
  
  message(
    "Spectra with 204 but without 366: ",
    sum(spectra_204_only)
  )
  
  message(
    "Spectra with both 204 and 366: ",
    sum(spectra_204_and_366)
  )
  
  message(
    "Spectra without 204 and therefore ignored: ",
    sum(!has_204)
  )
  
  
  # Read precursor m/z and charge once for every spectrum.
  
  precursor_mz <- vapply(
    spectra,
    get_pepmass,
    numeric(1)
  )
  
  precursor_charge <- vapply(
    spectra,
    get_charge,
    integer(1)
  )
  
  
  valid_metadata <- (
    is.finite(precursor_mz) &
      !is.na(precursor_charge) &
      precursor_charge > 0
  )
  
  
  relevant_invalid_metadata <- (
    has_204 & !valid_metadata
  )
  
  if (any(relevant_invalid_metadata)) {
    warning(
      sum(relevant_invalid_metadata),
      paste0(
        " spectra containing 204 have missing or invalid ",
        "PEPMASS/CHARGE and will be skipped."
      )
    )
  }
  
  
  # Determine output folder and base filename.
  
  input_directory <- dirname(
    normalizePath(
      input_file,
      winslash = "/",
      mustWork = TRUE
    )
  )
  
  input_filename <- basename(
    input_file
  )
  
  input_base_name <- tools::file_path_sans_ext(
    input_filename
  )
  
  
  processing_summary <- data.frame(
    label = sugar_deltas$label,
    delta_mass = sugar_deltas$delta_mass,
    diagnostic_group = sugar_deltas$diagnostic_group,
    spectra_written = integer(nrow(sugar_deltas)),
    output_file = character(nrow(sugar_deltas)),
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Create one output file for every sugar delta
  # ----------------------------------------------------------
  
  for (delta_index in seq_len(nrow(sugar_deltas))) {
    
    label <- sugar_deltas$label[delta_index]
    delta_mass <- sugar_deltas$delta_mass[delta_index]
    diagnostic_group <-
      sugar_deltas$diagnostic_group[delta_index]
    
    
    output_file <- file.path(
      input_directory,
      paste0(
        input_base_name,
        "_",
        label,
        ".mgf"
      )
    )
    
    
    # Select the appropriate diagnostic class.
    
    if (diagnostic_group == "204_only") {
      
      diagnostic_filter <- spectra_204_only
      
    } else if (diagnostic_group == "204_and_366") {
      
      diagnostic_filter <- spectra_204_and_366
      
    } else {
      
      stop(
        paste0(
          "Unknown diagnostic group: ",
          diagnostic_group
        )
      )
    }
    
    
    # Precursor * charge - delta_mass >= 300
    
    remaining_precursor_mass <- (
      precursor_mz * precursor_charge
    ) - delta_mass
    
    
    spectrum_is_eligible <- (
      diagnostic_filter &
        valid_metadata &
        is.finite(remaining_precursor_mass) &
        remaining_precursor_mass >= minimum_remaining_mass
    )
    
    
    eligible_indices <- which(
      spectrum_is_eligible
    )
    
    
    modified_spectra <- lapply(
      eligible_indices,
      function(spectrum_index) {
        
        spectrum <- spectra[[spectrum_index]]
        
        new_precursor_mz <- (
          precursor_mz[spectrum_index] -
            delta_mass / precursor_charge[spectrum_index]
        )
        
        replace_pepmass(
          spectrum = spectrum,
          new_pepmass = new_precursor_mz,
          digits = pepmass_digits
        )
      }
    )
    
    
    write_mgf_file(
      spectra = modified_spectra,
      output_file = output_file,
      add_mass_header = add_mass_header
    )
    
    
    processing_summary$spectra_written[delta_index] <-
      length(modified_spectra)
    
    processing_summary$output_file[delta_index] <-
      output_file
    
    
    message(
      "Created: ",
      basename(output_file),
      " | diagnostic group: ",
      diagnostic_group,
      " | spectra written: ",
      length(modified_spectra)
    )
  }
  
  
  message("Processing completed.")
  
  processing_summary
}


# ============================================================
# 5. RUN THE SCRIPT
# ============================================================

result_summary <- process_mgf_sugar_subtractions(
  input_file = input_file,
  sugar_deltas = sugar_deltas,
  diagnostic_204_mz = diagnostic_204_mz,
  diagnostic_204_tolerance = diagnostic_204_tolerance,
  diagnostic_366_mz = diagnostic_366_mz,
  diagnostic_366_tolerance = diagnostic_366_tolerance,
  minimum_remaining_mass = minimum_remaining_mass,
  pepmass_digits = pepmass_digits,
  add_mass_header = add_mass_header
)


# Display output summary
print(result_summary)