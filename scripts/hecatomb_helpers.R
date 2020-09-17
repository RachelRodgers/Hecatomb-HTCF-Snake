
#----- Function to load or install require packages -----#

TryInstall <- function(reqPackage) {
  # Attempts to install required package to a path in .libPaths.
  # Will continue trying paths until one works or all fail.
  #
  # Args:
  #   reqPackage: A required package
  #
  # Returns:
  #   Loads package to continue or stops with exception.
  
  systemLibPaths <- .libPaths()
  
  if (!require(reqPackage, character.only = TRUE, quietly = TRUE)) { # if not installed
    for (i in 1:length(systemLibPaths)) { # try installing in each path till it works
      currPath <- systemLibPaths[i]
      
      tryCatch({ # install + load
        cat("\nAttempting to install", reqPackage, "to", currPath, "\n")
        install.packages(reqPackage, lib = currPath, verbose = FALSE, quiet = TRUE)
        require(reqPackage, character.only = TRUE, quietly = TRUE)
        
      }, warning = function(w) {
        if (i < length(systemLibPaths)) {
          message("\nProblem with package installation. Attempting installation again...\n")
          
        } else {
          stop(cat("\nFATAL: Package", reqPackage, "could not be automatically installed into R. Please install manually.\n\n"),
               call. = FALSE)
        }
      }, error = function(e) "error")
    }
  } else { # loads package if available
    require(reqPackage, character.only = TRUE)}
}
