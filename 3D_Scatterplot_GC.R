### Creates 3D scatterplot of phage GC content vs host GC content vs phage length

library(scatterplot3d)
# create column indicating point color
a$pcolor[a$Phylum=='Proteobacteria'] <- "red"
a$pcolor[a$Phylum=='Firmicutes'] <- "blue"
a$pcolor[a$Phylum=='Actinobacteria'] <- "darkgreen"
a$pcolor[a$Phylum=='Spirochaetes'] <- "cyan"
a$pcolor[a$Phylum=='Bacteroidetes'] <- "gold"
a$pcolor[a$Phylum=='Other'] <- "black"
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

with(a, {
  s3d <- scatterplot3d(Host_GC, PP_GC, PP_Length,        # x y and z axis
                       color=pcolor,pch='.',cex.symbols=4,        # circle color indicates no. of cylinders
                       main="GC Content",
                       xlab="Host GC%",
                       ylab="Prophage GC%",
                       zlab="Prophage Length",
                       col.grid = 'blue',
                       col.axis = 'blue'
              
                       )
  #Regression line
  my.lm <- lm(PP_Length ~ Host_GC + PP_GC)
  s3d$plane3d(my.lm)
  
  # add the legend
  legend("topleft", inset=.05,      # location and inset
         bty="n", cex=.5,              # suppress legend box, shrink text 50%
         title="Phylum",
         c("Proteobacteria", "Firmicutes", "Actinobacteria", "Spirochaetes", "Bacteriodetes", "Other"), fill=c("red", "blue", "darkgreen", "cyan", "gold", "black"))
})
# 3. Add grids
addgrids3d(a[,1:3],grid = c("xy", "xz", "yz"))

