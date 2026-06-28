# Run this script in RStudio or R console to initialize dependency tracking
install.packages("renv")
renv::init(bare = TRUE)
renv::snapshot()
