## build command: docker build -f gg-tidyr.Dockerfile -t biolockjdevteam/gg-tidyr:v2 .

FROM rocker/tidyverse:3.6.3

#1.) set shell to bash
SHELL ["/bin/bash", "-c"]
ARG DEBIAN_FRONTEND=noninteractive


#5.) Install more R Packages
RUN Rscript -e "install.packages('tidyr', dependencies=c('Depends', 'Imports') )"
RUN Rscript -e "install.packages('ggpubr', dependencies=c('Depends', 'Imports') )"
RUN Rscript -e "install.packages('vegan', dependencies=c('Depends', 'Imports') )"
#RUN Rscript -e "install.packages('reshape2', dependencies=c('Depends', 'Imports') )"



#6.) check that packages installed
RUN Rscript -e "library('ggpubr'); library('tidyr');"
RUN Rscript -e "library('vegan'); "
RUN Rscript -e "library('reshape2'); "
RUN Rscript -e "library('ggplot2'); "
RUN Rscript -e "library('ggrepel'); "


#7.) Cleanup
RUN	apt-get clean && \
	find / -name *python* | xargs rm -rf && \
	rm -rf /tmp/* && \
	rm -rf /usr/share/* && \
	rm -rf /var/cache/* && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/log/*