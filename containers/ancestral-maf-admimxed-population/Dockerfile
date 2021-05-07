FROM uwgac/topmed-master:2.8.1
USER root
RUN Rscript -e 'install.packages("gridExtra", repos = "https://cran.rstudio.com/")'
RUN Rscript -e 'install.packages("Epi", repos = "https://cran.rstudio.com/")'
RUN Rscript -e 'install.packages("splitstackshape", repos = "https://cran.rstudio.com/")'
RUN Rscript -e 'install.packages("rtracklayer", repos = "https://cran.rstudio.com/")'
RUN Rscript -e 'BiocManager::install("rtracklayer")'
RUN Rscript -e 'install.packages("tidyverse", repos = "https://cran.rstudio.com/")'
