library(RColorBrewer)
library(reshape2)
library(gplots)

rm(list=ls())


# Load outout
load('determinePower/output/con_eval_v2auto.RData')


# R^2 powers plot for each tissue

tissues = c("Adipose", "Muscle", "Brain")

tissue_col = brewer.pal(9, "Set1")[c(1:5, 7:9)]

# Separate tissue codes
tissue_pairs = lapply(
  names(con_eval_pairs),
  function(s) strsplit(s, "-")[[1]]
)


pdf("determinePower/plots/beta_series.pdf", height=9, width=5)
par(mfrow=c(4, 2))
for (i in 1:length(con_eval)) {
  tissue = names(con_eval)[i]
  
  # Plot within-tissue R^2 curve
  plot(opts$powers,
       # con_eval[[i]][1,],
       -sign(con_eval[[i]]$fitIndices[,3])*con_eval[[i]]$fitIndices[,2],
       type="l",
       main=tissue,
       ylim=c(0, 1),
       lwd=2,
       ylab=expression(R^2),
       xlab=expression(beta)
  )
  
  # find pairwise ids involving ith tissue
  idx = which(sapply(tissue_pairs, function(x) tissue %in% x))
  
  for (j in idx) {
    to_tis = tissue_pairs[[j]][tissue_pairs[[j]] != tissue]
    
    lines(opts$powers,
          # con_eval_pairs[[j]][1,],
          -sign(con_eval_pairs[[j]]$fitIndices[,3])*con_eval_pairs[[j]]$fitIndices[,2],
          col=tissue_col[match(to_tis, tissues)],
          lwd=1.5
    )
    tissue_pairs
  }
  
  abline(v=2.7, col="grey", lty=2)
  abline(v=5.2, col="grey", lty=1)
}

plot(0, 0, type="n")
legend("bottomleft", legend=tissues, col=tissue_col, pch=15, cex=1.0)
dev.off()





# Slope and mean connectivity comparisons
# -----------------------------------------

for (i in 1:length(con_eval)) {
  con_eval[[i]]$fitIndices$signed_R2 = -sign(con_eval[[i]]$fitIndices[,3])*con_eval[[i]]$fitIndices[,2]
}


pdf("determinePower/plots/beta_series_tissue_slope.pdf", height=16, width=5)

par(mfrow=c(4, 1))
yrange = range(sapply(con_eval, function(x) x$fitIndices$signed_R2))
plot(con_eval[[1]]$fitIndices$Power, con_eval[[1]]$fitIndices$signed_R2, type="b", col=tissue_col[i], ylim=yrange, pch=16, xlab=expression(beta), ylab=expression(R^2))
for (i in 2:length(con_eval)) {
  points(con_eval[[i]]$fitIndices$Power, con_eval[[i]]$fitIndices$signed_R2, col=tissue_col[i], type="b", pch=16)
}
abline(v=5.2, col="grey", lty=1)



yrange = range(sapply(con_eval, function(x) x$fitIndices$slope))
yrange = range(c(yrange, -1))
plot(con_eval[[1]]$fitIndices$Power, con_eval[[1]]$fitIndices$slope, type="b", col=tissue_col[i], ylim=yrange, pch=16, xlab=expression(beta), ylab="Slope")
for (i in 2:length(con_eval)) {
  points(con_eval[[i]]$fitIndices$Power, con_eval[[i]]$fitIndices$slope, col=tissue_col[i], type="b", pch=16)
}
abline(v=5.2, col="grey", lty=1)



yrange = range(sapply(con_eval, function(x) x$fitIndices$mean.k.))
plot(con_eval[[1]]$fitIndices$Power, con_eval[[1]]$fitIndices$mean.k., type="b", col=tissue_col[i], ylim=yrange, pch=16, xlab=expression(beta), ylab="Mean connectivity")
for (i in 2:length(con_eval)) {
  points(con_eval[[i]]$fitIndices$Power, con_eval[[i]]$fitIndices$mean.k., col=tissue_col[i], type="b", pch=16)
}
abline(v=5.2, col="grey", lty=1)


yrange = range(sapply(con_eval, function(x) log10(x$fitIndices$mean.k.)))
plot(con_eval[[1]]$fitIndices$Power, log10(con_eval[[1]]$fitIndices$mean.k.), type="b", col=tissue_col[i], ylim=yrange, pch=16, xlab=expression(beta), ylab="Mean connectivity (log10)")
for (i in 2:length(con_eval)) {
  points(con_eval[[i]]$fitIndices$Power, log10(con_eval[[i]]$fitIndices$mean.k.), col=tissue_col[i], type="b", pch=16)
}
abline(v=5.2, col="grey", lty=1)


legend("topright", legend=tissues, col=tissue_col, pch=16, cex=1.0)
dev.off()


