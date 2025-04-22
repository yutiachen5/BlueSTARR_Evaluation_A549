library(dplyr)
# on log2 scale
read_df = function(path,x){
  df = read.table(path, header = FALSE, sep = "\t")
  print(nrow(df))
  colnames(df) = c('location','DMSO.refScore','Dex.refScore','sequence')
  df$fc = df$Dex.refScore - df$DMSO.refScore
  # df$distance = rep(seq(23, 243, by = 20), 50)
  df$distance = rep(x, nrow(df)/length(x))
  df_fc <- df %>%
    group_by(distance) %>%
    summarise(
      mean_fc = mean(fc),
      sd_fc = sd(fc)
    )
  
  return (df_fc)
}

############################# full model #######################################
# model trained by all chromosomes
x = seq(23, 563, by = 20)
# Dex-DMSO
# df_fc1 = read_df('/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/600-bases/chr1-b-ap1-b-gr-fc.txt',x)
# df_fc2 = read_df('/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/600-bases/chr1-b-ap1-f-gr-fc.txt',x)
# df_fc3 = read_df('/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/600-bases/chr1-f-ap1-b-gr-fc.txt',x)
# df_fc4 = read_df('/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/600-bases/chr1-f-ap1-f-gr-fc.txt',x)

df_fc1 = read_df('/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/600-bases/tf+pos/chr1-b-ap1-b-gr-fc.txt',x)
df_fc2 = read_df('/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/600-bases/tf+pos/chr1-b-ap1-f-gr-fc.txt',x)
df_fc3 = read_df('/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/600-bases/tf+pos/chr1-f-ap1-b-gr-fc.txt',x)
df_fc4 = read_df('/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/600-bases/tf+pos/chr1-f-ap1-f-gr-fc.txt',x)

pdf(file = '/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/600-bases/tf+pos/tf+pos.pdf', width = 12, height = 8)

par(mfrow=c(2,2), mar=c(3, 3, 1, 1), mgp=c(2, 1, 0)) 

plot(x, df_fc1$mean_fc,type = "b", pch = 19,
     main = "chr1 ap1-/gr-", xaxt = "n", ylim = c(0,2),
     xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
arrows(x, df_fc1$mean_fc - df_fc1$sd_fc, x, df_fc1$mean_fc + df_fc1$sd_fc,
       angle = 90, code = 3, length = 0.1, lwd = 1.5)
axis(side = 1, at = x, labels = x, cex.axis = 0.6)


plot(x, df_fc2$mean_fc,type = "b", pch = 19,
     main = "chr1 ap1-/gr+", xaxt = "n", ylim = c(0,2),
     xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
arrows(x, df_fc2$mean_fc - df_fc2$sd_fc, x, df_fc2$mean_fc + df_fc2$sd_fc,
       angle = 90, code = 3, length = 0.1, lwd = 1.5)
axis(side = 1, at = x, labels = x, cex.axis = 0.6)

plot(x, df_fc3$mean_fc,type = "b", pch = 19,
     main = "chr1 ap1+/gr-", xaxt = "n", ylim = c(0,2),
     xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
arrows(x, df_fc3$mean_fc - df_fc3$sd_fc, x, df_fc3$mean_fc + df_fc3$sd_fc,
       angle = 90, code = 3, length = 0.1, lwd = 1.5)
axis(side = 1, at = x, labels = x, cex.axis = 0.6)

plot(x, df_fc4$mean_fc,type = "b", pch = 19,
     main = "chr1 ap1+/gr+", xaxt = "n", ylim = c(0,2),
     xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
arrows(x, df_fc4$mean_fc - df_fc4$sd_fc, x, df_fc4$mean_fc + df_fc4$sd_fc,
       angle = 90, code = 3, length = 0.1, lwd = 1.5)
axis(side = 1, at = x, labels = x, cex.axis = 0.6)

dev.off()

# negative control K562-DMSO
# df_nc1 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/K562-DMSO/chr1-b-ap1-b-gr-fc.txt')
# df_nc2 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/K562-DMSO/chr1-b-ap1-f-gr-fc.txt')
# df_nc3 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/K562-DMSO/chr1-f-ap1-b-gr-fc.txt')
# df_nc4 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/K562-DMSO/chr1-f-ap1-f-gr-fc.txt')

# par(mfrow=c(2,2), mar=c(3, 3, 1, 1), mgp=c(2, 1, 0)) 
# x = seq(23, 243, by = 20)
# plot(x, df_nc1$mean_fc,type = "b", pch = 19,
#      main = "chr1 ap1-/gr-", xaxt = "n", ylim = c(-0.8,-0.2),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (K562/DMSO)")
# arrows(x, df_nc1$mean_fc - df_nc1$sd_fc, x, df_nc1$mean_fc + df_nc1$sd_fc,
#        angle = 90, code = 3, length = 0.1, lwd = 1.5)
# axis(side = 1, at = x, labels = x, cex.axis = 0.6)

# plot(x, df_nc2$mean_fc,type = "b", pch = 19,
#      main = "chr1 ap1-/gr+", xaxt = "n", ylim = c(-0.8,-0.2),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (K562/DMSO)")
# arrows(x, df_nc2$mean_fc - df_nc2$sd_fc, x, df_nc2$mean_fc + df_nc2$sd_fc,
#        angle = 90, code = 3, length = 0.1, lwd = 1.5)
# axis(side = 1, at = x, labels = x, cex.axis = 0.6)

# plot(x, df_nc3$mean_fc,type = "b", pch = 19,
#      main = "chr1 ap1+/gr-", xaxt = "n", ylim = c(-0.8,-0.2),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (K562/DMSO)")
# arrows(x, df_nc3$mean_fc - df_nc3$sd_fc, x, df_nc3$mean_fc + df_nc3$sd_fc,
#        angle = 90, code = 3, length = 0.1, lwd = 1.5)
# axis(side = 1, at = x, labels = x, cex.axis = 0.6)

# plot(x, df_nc4$mean_fc,type = "b", pch = 19,
#      main = "chr1 ap1+/gr+", xaxt = "n", ylim = c(-0.8,-0.2),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (K562/DMSO)")
# arrows(x, df_nc4$mean_fc - df_nc4$sd_fc, x, df_nc4$mean_fc + df_nc4$sd_fc,
#        angle = 90, code = 3, length = 0.1, lwd = 1.5)
# axis(side = 1, at = x, labels = x, cex.axis = 0.6)




