##Install R Language
```shell
sudo apt-get install r-base r-base-core
```

##Installing packages

###Option 1: package install in R (for fast and parallel versions)

Default package manager:
```r
install.packages(c("snow", "data.table", "readr"))
```

Bioconductor package manager (may be neccessary for `readr`)
```r
source("https://bioconductor.org/biocLite.R"
biocLite(c("snow", "data.table", "readr"))
```

###Option 2: pacakge install via terminal
```shell
sudo apt-get install r-cran-snow r-cran-data.table r-cran-readr
```
