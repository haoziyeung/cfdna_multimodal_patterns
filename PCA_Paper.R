library(ggfortify)
library(ggplot2)

my_data <- read.csv("precar_pca_baseline.txt",sep="\t",header=T)

pca_res <- prcomp(my_data[3:10], scale. = T)

autoplot(pca_res, data = my_data, colour = 'class',frame.type='norm',frame=T,size=3) +
  theme_bw()+
  theme_classic()+
  theme(legend.position="right",
        axis.text.x=element_text(colour="black",family="Times",size=21),
        axis.text.y=element_text(colour="black",family="Times",size=21,face="plain"),
        axis.title.y=element_text(family="Times",size =21,face="plain"),
        axis.title.x=element_text(family="Times",size =21,face="plain"),
        plot.title = element_text(family="Times",size=24,face="bold",hjust = 0.5),
        axis.line.x=element_line(linetype=1,color="black",size=1.5),
        axis.line.y=element_line(linetype=1,color="black",size=1.5)
  )+
  theme(plot.title = element_text(face = "bold"))
