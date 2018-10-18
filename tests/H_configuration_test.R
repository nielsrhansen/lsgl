library(lsgl)
# warnings = errors
options(warn=2)

c_config <- lsgl.c.config()

if(c_config$timing) stop()

if(c_config$debugging) stop()
