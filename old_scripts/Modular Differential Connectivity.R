# Modular Differential Connectivity (MDC): 
# Compute ratio of the mean intra-module connectivity in one network to that of the same set of genes in another network
# This code uses functions from: https://github.com/songw01/melanoma_network_codes
# Author(s): Bin Zhang Cell 153(3):707-720. PMID: 23622250 
#
#
#    Input:
#         1) inputfnameB       -- mRNA data (tab-delimited text file) for one disease state B
#         2) inputfnameBModule -- module assignment file as output from WINA or WGCNA
#         3) inputfname        -- mRNA data (tab-delimited text file) for another disease state A
#         4) shortnames        -- shortnames for B and A, respectively. Order is important
#         5) headCol           -- # of gene information columns in inputfname
#         6) headCol2          -- # of gene information columns in inputfnameB
#         7) corrlpower        -- the exponent of the power function, determined by WINA or WGCNA
#
#   Output:
#         1) "DiffConnectvty\*_wFDR.xls"             -- MDC with false discovery rates
#         2) "DiffConnectvty\*.png"                  -- MDC plot
#         3) "DiffConnectvty\Heatmaps\*_iHMP0_*.png" -- MDC heatmaps of mulitple modules 
#         3) "DiffConnectvty\Heatmaps\*_iHMP-*.png" -- MDC heatmaps of individual modules
#
#
################################################################################################
#memory.limit(size = 35000)
rm(list = ls())
library(lattice) # require is design for use inside functions 
library(plotrix)
#library(sma) # this is needed for plot.mat below
source("R-functions_MDC.R")
source("colPalette_200.txt")


##### 
inputfnameB       = "adj_mat_young.Rds"; headCol2 =2;
inputfnameBModule = "Clusters_table_young.txt"
inputfname        = "adj_mat_old.Rds"; headCol =2;
shortnames        = c("young", "old")
corrlpower        = 6

#
# -----------------------------End of Parameters to be changed --------------------------------------

no.perms = 50
imgwid=600; imghei=400

#heatmapColorRG = rev( rgcolors.func(50) )
heatmapColor = heat.colors(50)

# specify the directory for holding analysis results, you need make such a sub-directory 
#    under your working directory if it doesn't exist
outputDir1  ="DiffConnectvty_young_to_old_new/"
dir.create(outputDir1)

outputDir  ="DiffConnectvty_young_to_old_new/Random/"
dir.create(outputDir)

outputDirTom = "DiffConnectvty_young_to_old_new/Heatmaps/"
dir.create(outputDirTom)

# image type
imagetype="png" #"pdf"


fname      =getFileNameNopath(inputfname)
fnameB     =getFileNameNopath(inputfnameBModule)
extname    =getFileExtension(inputfname)

fname = paste(fnameB, "_VS_", shortnames[2], "_MDC", sep="")
flog  = paste(outputDir1, fname, "_Randoms", ".xls",       sep='')
flog2 = paste(outputDir1, fname, "_wFDR", ".xls",       sep='')
fimgRt= paste(outputDir1, fname, ".png",sep='')  #png only

#------------------------- second file --------------------------------------

datExpr2 = readRDS(inputfnameB)
dim(datExpr2)

no.genes <- dim(datExpr2)[2]
genesInfor2 <- colnames(datExpr2)

# ------ modules from 2nd file --------------------------------
allMatrix2modinfo <- read.delim(inputfnameBModule,sep="\t", header=T)
allMatrix2modinfo$Gene.Symbol <- paste(allMatrix2modinfo$Tissue,"_",allMatrix2modinfo$Gene.Symbol,sep="")
dim(allMatrix2modinfo)
ncols2m = dim(allMatrix2modinfo)[2]

allMatrix2modinfo <- as.matrix(allMatrix2modinfo)[,c(1,ncols2m)]
genes2Idx = cbind(as.matrix(genesInfor2)[,1], c(1:no.genes) )
merged = merge(genes2Idx, allMatrix2modinfo, by.x=1, by.y=2, all.x=T)
merged = as.matrix(merged)
dim(merged)
morder = order(as.integer(merged[,2]))
allMatrix2modinfoXX = merged[morder,] # aligned to the other gene infos
modulescolorB = as.character(merged[morder,3])
modulescolorB = ifelse(is.na(modulescolorB), "grey", modulescolorB)

modtb = table(modulescolorB)
modulenames = names(modtb)
modtb = modtb[modulenames!="grey"]

if (length(modulenames) > length(col.names))
{
  col.names <- colors();
  col.names <- col.names[-grep("grey|gray|white",col.names)]
  if (length(modulenames) > length(col.names))
  {
    col.names <- rep(col.names,ceiling(length(modulenames)/length(col.names)))
  }
}


modulenames = names(modtb)

umodulesizes = union(sort(as.integer(modtb)), NULL)

no.modules  = length(modulenames)
no.nets     = 2

#------------------------- first file --------------------------------------

datExpr <- readRDS(inputfname)
dim(datExpr)

#corhelp<- abs(corhelp)

#*-------------------------------------------------------------------------------------
#* initilization

uno = length(umodulesizes)

meanPermodule   = matrix(0, no.modules, no.nets)
sdPermodule     = matrix(0, no.modules, no.nets)
ks.pvalues      = rep(1, no.modules)
linkCor         = matrix(-1, no.modules, no.nets)

meanPermoduleTop   = matrix(0, no.modules, no.nets)
sdPermoduleTop     = matrix(0, no.modules, no.nets)
ks.pvaluesTop      = rep(1, no.modules)

heatmaps.list      = as.list(rep(NA, no.modules))
modsizes           = rep(NA, no.modules)



#############################################################################
#-------------------- permute genes ---------------------------------------
#
no.perms = 50
GlobalyRandomMDC = NULL
for( y in c(1:no.perms) ) {
  if(y%%5==0) {
    print(paste("********* random genes", y, " ........"))
  }
  
  GRmeanPermodule   = matrix(0, no.modules, no.nets)
  for ( x in c(1:no.modules) ) {
    
    xsel = sample(c(1:no.genes), as.integer(modtb[x]), replace=F)
    corhelp <- datExpr[xsel, xsel]
    diag(corhelp) <- 0
    links   <- apply(abs(corhelp), 1, sum, na.rm=T)
    
    lpanel  <- lower.tri(corhelp)
    corhelp <- abs(corhelp[lpanel])
    isel    <- !is.na(corhelp)
    no.corrls = length(corhelp)
    
    corhelp <- corhelp[isel]
    
    # ---------------------------------------------------------------
    corhelp2 <-datExpr2[xsel, xsel]
    diag(corhelp2) <- 0
    links2   <- apply(abs(corhelp2), 1, sum, na.rm=T)
    
    corhelp2 <- abs(corhelp2[lpanel])
    isel     <- !is.na(corhelp2)
    corhelp2 <- corhelp2[isel]
    
    GRmeanPermodule[x, 1] = mean(corhelp2, na.rm=T)
    GRmeanPermodule[x, 2] = mean(corhelp, na.rm=T)
    
    rm(corhelp2)
    rm(corhelp)
    collect_garbage()
  }
  GlobalyRandomMDC = cbind(GlobalyRandomMDC, GRmeanPermodule[,1]/GRmeanPermodule[,2])
}

GRMDC_mean = apply(GlobalyRandomMDC, 1, mean)
GRMDC_sd = apply(GlobalyRandomMDC, 1, sd)
GR_Dat = cbind(yfold, GRMDC_mean, GRMDC_sd)
GR_FDR_N = apply(GR_Dat, 1, MDC_FDR_by_normal_distr)

GR_FDR = apply(cbind(yfold, GlobalyRandomMDC), 1, MDC_FDR)

xfinal = cbind(modulenames, meanFoldChange, GlobalyRandomMDC)

colnames(xfinal) <- c("module", paste("MDC_random_genes_", c(1:no.perms), sep="") )

write.table(xfinal, flog, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

#--------------- compute FDR --------------------
#
rmean = apply(meanFoldChange, 1, mean)
rsd   = apply(meanFoldChange, 1, sd)
mdcDat = cbind(yfold, rmean, rsd)
mdcFDR_N = apply(mdcDat, 1, MDC_FDR_by_normal_distr) # FDR by distribution

mdcFDR = apply(cbind(yfold, meanFoldChange), 1, MDC_FDR)

comFDR = apply(cbind(GR_FDR, mdcFDR), 1, max)

final = cbind(modulenames, round(yfold,2),  signif(comFDR,4), signif(mdcFDR,4), signif(GR_FDR,4))
colnames(final) = c("module", "MDC",  "FDR", "FDR_random_samples", "FDR_random_genes")
final

write.table(final, flog2, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)


#################################################################################
#---------------------- plot ratio of means -----------------------------------
#
rylab = paste("MDC: k_", shortnames[1]," / k_", shortnames[2], sep="")

imgwid=1400; imghei=600;
ignore_small_modules=F; module_sizecut=100

hlines = c(2,1.5,1)
if (max(yfold)>20){ hlines = c(20, 10, 5) }

rylab = paste("MDC(", shortnames[1],", ", shortnames[2], ")", sep="")
od = order(-yfold)
dispDat = cbind(yfold[od]); dispNames=modulenames[od];
rownames(dispDat) <- dispNames; dispModuleSize = modtb[od]
dispSign = ifelse(comFDR[od]<= 0.10, "*", "")

yymax = closest_integer_bound(max(yfold))
barplot_by_oneGroup(datMatrix=dispDat, maskMatrix=cbind(dispSign), 
                    rowcol=col.names[1:length(dispNames)], vlines=hlines,
                    title="", #"Modular Differential Connectivity", 
                    ylab=rylab, 
                    imgfile=fimgRt, imgwid=2000, imghei=800, 
                    mcex=0.5, xlab_cex=0.65, xaxisLabel_rotate=45, 
                    legendcols=1, legendxDelta=0, legendCex=1, 
                    ipointsize=12, iunits="px", ires=300, icompression="lzw",
                    show_value=FALSE, showLegend=FALSE, show_xaxisLabel=TRUE, show_mask=TRUE,
                    xmargin=3, ymargin=5, zmargin=0.5, xlabel_orientation=0, ValueXdelta=0, 
                    ValueYdelta=2*yymax/100, valuecex=1, yMAX=yymax, yMIN=0)

#************************************************************************************
#
#  Plot Heatmaps
#
#************************************************************************************

orderedColor = modulescolorB

moduleHeatmaps = as.list(rep(NA, no.modules))

yhei=1000; ydelta = 0.2
if(no.modules<25) {
  rowfigs = as.integer(no.modules^0.5); subfigs = as.integer(no.modules/rowfigs+0.5)*rowfigs;
  if(subfigs < no.modules) {subfigs =subfigs + rowfigs}
  #subfigs = 12; rowfigs = 4
  print(paste(subfigs, rowfigs)) 
} else {
  subfigs = 24; rowfigs = 4
}

ywid = as.integer(yhei*(rowfigs^2)/subfigs);
tcexbase = 2.2 #1.5-unix #3.5, windows

#zrange=c(0, 1-max(max(dist1), max(dist1)) )
zrange=c(0,1)

module_order_by_size = FALSE

# ---------------- combined plot heatmap ----------------------
#
for (z in c(1:no.modules) ) {
  
  #image(heatmaps.list[[z]], xlim=c(0,1), ylim=c(0,1+ydelta), zlim=zrange, axes=F, col=heatmapColor)
  
  if(z%%subfigs==1) {
    end = subfigs+z-1; if(end>no.modules){end=no.modules}
    if(module_order_by_size) {
      imgHeatMap  =paste(outputDirTom, fname, "_iHMP0_", z, "-", end, ".png",sep='')  #png only
    } else{
      imgHeatMap  =paste(outputDirTom, fname, "_iHMP0_", modulenames[z], "-", modulenames[end], ".png",sep='')  #png only
    }
    
    openImgDev(imgHeatMap, iwidth =ywid, iheight = yhei)
    par(mfrow=c(subfigs/rowfigs,rowfigs), mar=c(0, 0, 0, 0),cex=1)
    #par(mfrow=c(2,2), mar=c(0, 0, 0, 0), cex=1)
    print(imgHeatMap)
  }
  
  par(mar=c(0.2, 0, 0.2, 0)+0.4)
  image((heatmaps.list[[z]])^2, xlim=c(0,1), ylim=c(0,1+ydelta), zlim=zrange, axes=F, col=heatmapColor)
  gradient.rect(0,1.05,1,1+ydelta,col=col.names[modulenames[z]],border=F)
  
  # font size
  colorchars = unlist(strsplit(modulenames[z], ""))
  if (length(colorchars)>12) {
    tcex = tcexbase#*(10/length(colorchars))^0.4
  } else{tcex=tcexbase}
  
  tcolor = SetTextContrastColor(col.names[z])
  text(0.5,1+ydelta*0.7,labels=modulenames[z],cex=tcex, col=tcolor)
  
  if(z%%subfigs==0) {
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    dev.off();
  }       
}

if(z%%subfigs>0) {
  
  scales = rev(seq(0,1, 0.01)); nscales = length(scales)
  scalesM= repmat(scales, 1, 100)
  textsel= (scales*100)%%25 ==0;
  scalesLabel = scales[textsel];
  #scalesLabel = as.integer(100*(zrange[1]+(zrange[2]-zrange[1])*scales[textsel])/100;
  ydel = 1.0
  
  par(mar=c(0, 0, 0, 0),cex=1)
  image(scalesM, xlim=c(-0.05,1.05), ylim=c(-ydel,1+1), zlim=zrange, axes=F, col=heatmapColor)
  # draw axis and ticks
  linecol = "brown"; ypos = -ydel/8; tickhei=ydel/4
  lines(x=c(0,1),y=c(ypos, ypos), col=linecol, lwd=1)
  for(iv in scalesLabel){
    lines(x=c(iv, iv), y=c(ypos,ypos-tickhei), col=linecol, lwd=1)
  }
  text(x=scalesLabel,y=rep(-ydel*0.7), labels=as.character(scalesLabel),cex=tcexbase*0.7, col="brown")
  
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
  dev.off();
  
}

#------------ plot individual heatmap ---------------------------------
#
yywid = 200
for (z in c(1:no.modules) ) {
  
  print(paste("Plot iHMP: ", modulenames[z]))
  
  imgHeatMap  =paste(outputDirTom, "/", fname, "_iHMP-", modulenames[z],".png",sep='')  #png only
  openImgDev(imgHeatMap, iwidth =yywid, iheight =yywid)
  par(mfrow=c(1,1), mar=c(0, 0, 0, 0),cex=1)
  
  image((heatmaps.list[[z]])^2, xlim=c(0,1), ylim=c(0,1), zlim=zrange, axes=F, col=heatmapColor)
  
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
  dev.off();
}

#quit(save = "no",status = 0)
rm(heatmaps.list)

