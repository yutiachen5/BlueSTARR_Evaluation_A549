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
x = seq(23, 243, by = 20)
# Dex-DMSO, without max pooling layer
# df_fc1 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/chr1-b-ap1-b-gr-fc.txt',x)
# df_fc2 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/chr1-b-ap1-f-gr-fc.txt',x)
# df_fc3 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/chr1-f-ap1-b-gr-fc.txt',x)
# df_fc4 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/chr1-f-ap1-f-gr-fc.txt',x)

# # Dex-DMSO, max-pooling layer added before the last convolution layer
# x = seq(23, 243, by = 20)
# df_fc1 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/max-pooling/chr1-b-ap1-b-gr-fc.txt',x)
# df_fc2 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/max-pooling/chr1-b-ap1-f-gr-fc.txt',x)
# df_fc3 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/max-pooling/chr1-f-ap1-b-gr-fc.txt',x)
# df_fc4 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/max-pooling/chr1-f-ap1-f-gr-fc.txt',x)

# # Dex-DMSO, last convolution layer was replaced by 3 attention layers
# x = seq(23, 243, by = 20)
# df_fc1 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/attention/chr1-b-ap1-b-gr-fc.txt',x)
# df_fc2 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/attention/chr1-b-ap1-f-gr-fc.txt',x)
# df_fc3 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/attention/chr1-f-ap1-b-gr-fc.txt',x)
# df_fc4 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/attention/chr1-f-ap1-f-gr-fc.txt',x)

# Dex-DMSO, transformer encoder layer and positional encoding layer
x = seq(23, 243, by = 20)
df_fc1 = read_df('/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/300-bases/attention-pos/chr1-b-ap1-b-gr-fc.txt',x)
df_fc2 = read_df('/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/300-bases/attention-pos/chr1-b-ap1-f-gr-fc.txt',x)
df_fc3 = read_df('/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/300-bases/attention-pos/chr1-f-ap1-b-gr-fc.txt',x)
df_fc4 = read_df('/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/300-bases/attention-pos/chr1-f-ap1-f-gr-fc.txt',x)


pdf("/datacommons/igvf-pm/A549/GR-AP1/enhancer-seq/fc-full/Dex-DMSO/300-bases/attention-pos/GR-AP1.pdf")

par(mfrow=c(2,2), mar=c(3, 3, 1, 1), mgp=c(2, 1, 0)) 

plot(x, df_fc1$mean_fc,type = "b", pch = 19,
     main = "chr1 ap1-/gr-", xaxt = "n", ylim = c(0,1.7),
     xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
arrows(x, df_fc1$mean_fc - df_fc1$sd_fc, x, df_fc1$mean_fc + df_fc1$sd_fc,
       angle = 90, code = 3, length = 0.1, lwd = 1.5)
axis(side = 1, at = x, labels = x, cex.axis = 0.6)

plot(x, df_fc2$mean_fc,type = "b", pch = 19,
     main = "chr1 ap1-/gr+", xaxt = "n", ylim = c(0,1.7),
     xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
arrows(x, df_fc2$mean_fc - df_fc2$sd_fc, x, df_fc2$mean_fc + df_fc2$sd_fc,
       angle = 90, code = 3, length = 0.1, lwd = 1.5)
axis(side = 1, at = x, labels = x, cex.axis = 0.6)

plot(x, df_fc3$mean_fc,type = "b", pch = 19,
     main = "chr1 ap1+/gr-", xaxt = "n", ylim = c(0,1.7),
     xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
arrows(x, df_fc3$mean_fc - df_fc3$sd_fc, x, df_fc3$mean_fc + df_fc3$sd_fc,
       angle = 90, code = 3, length = 0.1, lwd = 1.5)
axis(side = 1, at = x, labels = x, cex.axis = 0.6)

plot(x, df_fc4$mean_fc,type = "b", pch = 19,
     main = "chr1 ap1+/gr+", xaxt = "n", ylim = c(0,1.7),
     xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
arrows(x, df_fc4$mean_fc - df_fc4$sd_fc, x, df_fc4$mean_fc + df_fc4$sd_fc,
       angle = 90, code = 3, length = 0.1, lwd = 1.5)
axis(side = 1, at = x, labels = x, cex.axis = 0.6)

# Dex-DMSO, single-based
# x = seq(8, 277, by = 1)
# df_fc1 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/single-based/chr1-b-ap1-b-gr-fc.txt',x)
# df_fc2 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/single-based/chr1-b-ap1-f-gr-fc.txt',x)
# df_fc3 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/single-based/chr1-f-ap1-b-gr-fc.txt',x)
# df_fc4 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/single-based/chr1-f-ap1-f-gr-fc.txt',x)

# # Dex-DMSO, single-based with max-pooling layer
# x = seq(8, 277, by = 1)
# df_fc1 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/single-based/max-pooling/chr1-b-ap1-b-gr-fc.txt',x)
# df_fc2 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/single-based/max-pooling/chr1-b-ap1-f-gr-fc.txt',x)
# df_fc3 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/single-based/max-pooling/chr1-f-ap1-b-gr-fc.txt',x)
# df_fc4 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-full/Dex-DMSO/single-based/max-pooling/chr1-f-ap1-f-gr-fc.txt',x)

# par(mfrow=c(2,2), mar=c(3, 3, 1, 1), mgp=c(2, 1, 0)) 

# plot(x, df_fc1$mean_fc, cex=0.2, type = "l", pch = 19,
#      main = "chr1 ap1-/gr-",  xaxt = "n", ylim = c(0,1.7),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
# axis(side = 1, at = c(8,seq(23, 243, by = 20),277), 
#      labels = c(8,seq(23, 243, by = 20),277), cex.axis = 0.5)

# plot(x, df_fc2$mean_fc,type = "l", pch = 19,cex=0.2,
#      main = "chr1 ap1-/gr+", xaxt = "n", ylim = c(0,1.7),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
# axis(side = 1, at = c(8,seq(23, 243, by = 20),277), 
#      labels = c(8,seq(23, 243, by = 20),277), cex.axis = 0.5)


# plot(x, df_fc3$mean_fc,type = "l", pch = 19,cex=0.2,
#      main = "chr1 ap1+/gr-", xaxt = "n", ylim = c(0,1.7),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
# axis(side = 1, at = c(8,seq(23, 243, by = 20),277), 
#      labels = c(8,seq(23, 243, by = 20),277), cex.axis = 0.5)


# plot(x, df_fc4$mean_fc,type = "l", pch = 19,cex=0.2,
#      main = "chr1 ap1+/gr+", xaxt = "n", ylim = c(0,1.7),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
# axis(side = 1, at = c(8,seq(23, 243, by = 20),277), 
#      labels = c(8,seq(23, 243, by = 20),277), cex.axis = 0.5)






# # AP1 and GR were put at the longest possible distance
# read_df_long = function(path){
#   df = read.table(path, header = FALSE, sep = "\t")
#   colnames(df) = c('location','DMSO.refScore','Dex.refScore','sequence')
#   df$fc = df$Dex.refScore - df$DMSO.refScore
#   df$distance = rep(287, 50)
#   df_fc <- df %>%
#     group_by(distance) %>%
#     summarise(
#       mean_fc = mean(fc),
#       sd_fc = sd(fc)
#     )
  
#   return (df_fc)
# }

# df_long1 = read_df_long('/hpc/home/yc583/BlueSTARR/GR-AP1/Dex-DMSO/longest-distance/chr1-b-ap1-b-gr-fc-longest-dis.txt')
# df_long2 = read_df_long('/hpc/home/yc583/BlueSTARR/GR-AP1/Dex-DMSO/longest-distance/chr1-b-ap1-f-gr-fc-longest-dis.txt')
# df_long3 = read_df_long('/hpc/home/yc583/BlueSTARR/GR-AP1/Dex-DMSO/longest-distance/chr1-f-ap1-b-gr-fc-longest-dis.txt')
# df_long4 = read_df_long('/hpc/home/yc583/BlueSTARR/GR-AP1/Dex-DMSO/longest-distance/chr1-f-ap1-f-gr-fc-longest-dis.txt')

# # concat
# df1 = rbind(df_fc1, df_long1)
# df2 = rbind(df_fc2, df_long2)
# df3 = rbind(df_fc3, df_long3)
# df4 = rbind(df_fc4, df_long4)

# par(mfrow=c(2,2), mar=c(3, 3, 1, 1), mgp=c(2, 1, 0)) 
# x = c(seq(23, 243, by = 20), 287)
# plot(x, df1$mean_fc,type = "b", pch = 19,
#      main = "chr1 ap1-/gr-", xaxt = "n", ylim = c(0,1.7),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
# arrows(x, df1$mean_fc - df1$sd_fc, x, df1$mean_fc + df1$sd_fc,
#        angle = 90, code = 3, length = 0.1, lwd = 1.5)
# axis(side = 1, at = x, labels = x, cex.axis = 0.6)

# plot(x, df2$mean_fc,type = "b", pch = 19,
#      main = "chr1 ap1-/gr+", xaxt = "n", ylim = c(0,1.7),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
# arrows(x, df2$mean_fc - df2$sd_fc, x, df2$mean_fc + df2$sd_fc,
#        angle = 90, code = 3, length = 0.1, lwd = 1.5)
# axis(side = 1, at = x, labels = x, cex.axis = 0.6)

# plot(x, df3$mean_fc,type = "b", pch = 19,
#      main = "chr1 ap1+/gr-", xaxt = "n", ylim = c(0,1.7),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
# arrows(x, df3$mean_fc - df3$sd_fc, x, df3$mean_fc + df3$sd_fc,
#        angle = 90, code = 3, length = 0.1, lwd = 1.5)
# axis(side = 1, at = x, labels = x, cex.axis = 0.6)

# plot(x, df4$mean_fc,type = "b", pch = 19,
#      main = "chr1 ap1+/gr+", xaxt = "n", ylim = c(0,1.7),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
# arrows(x, df4$mean_fc - df4$sd_fc, x, df4$mean_fc + df4$sd_fc,
#        angle = 90, code = 3, length = 0.1, lwd = 1.5)
# axis(side = 1, at = x, labels = x, cex.axis = 0.6)




# # negative control K562-DMSO
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


# # sequence used in the original paper
# df_fc_original = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/Dex-DMSO/original-seq-fc.txt')
# x = seq(23, 243, by = 20)
# plot(x, df_fc_original$mean_fc,type = "b", pch = 19,
#     main = "sequence used in original paper", xaxt = "n",
#     xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
# arrows(x, df_fc_original$mean_fc - df_fc_original$sd_fc, x, df_fc_original$mean_fc + df_fc_original$sd_fc,
#       angle = 90, code = 3, length = 0.1, lwd = 1.5)
# axis(side = 1, at = x, labels = x, cex.axis = 0.6)

# ############################# chr1 model #######################################
# # model not trained on chr1
# df_fc1 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-chr1/Dex-DMSO/chr1-b-ap1-b-gr-fc.txt')
# df_fc2 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-chr1/Dex-DMSO/chr1-b-ap1-f-gr-fc.txt')
# df_fc3 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-chr1/Dex-DMSO/chr1-f-ap1-b-gr-fc.txt')
# df_fc4 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-chr1/Dex-DMSO/chr1-f-ap1-f-gr-fc.txt')


# par(mfrow=c(2,2), mar=c(3, 3, 1, 1), mgp=c(2, 1, 0)) 
# x = seq(23, 243, by = 20)
# plot(x, df_fc1$mean_fc,type = "b", pch = 19,
#      main = "chr1 ap1-/gr-", xaxt = "n", ylim = c(-0.5,2),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
# arrows(x, df_fc1$mean_fc - df_fc1$sd_fc, x, df_fc1$mean_fc + df_fc1$sd_fc,
#        angle = 90, code = 3, length = 0.1, lwd = 1.5)
# axis(side = 1, at = x, labels = x, cex.axis = 0.6)

# plot(x, df_fc2$mean_fc,type = "b", pch = 19,
#      main = "chr1 ap1-/gr+", xaxt = "n", ylim = c(-0.5,2),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
# arrows(x, df_fc2$mean_fc - df_fc2$sd_fc, x, df_fc2$mean_fc + df_fc2$sd_fc,
#        angle = 90, code = 3, length = 0.1, lwd = 1.5)
# axis(side = 1, at = x, labels = x, cex.axis = 0.6)

# plot(x, df_fc3$mean_fc,type = "b", pch = 19,
#      main = "chr1 ap1+/gr-", xaxt = "n", ylim = c(-0.5,2),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
# arrows(x, df_fc3$mean_fc - df_fc3$sd_fc, x, df_fc3$mean_fc + df_fc3$sd_fc,
#        angle = 90, code = 3, length = 0.1, lwd = 1.5)
# axis(side = 1, at = x, labels = x, cex.axis = 0.6)

# plot(x, df_fc4$mean_fc,type = "b", pch = 19,
#      main = "chr1 ap1+/gr+", xaxt = "n", ylim = c(-0.5,2),
#      xlab = "Base Pairs between AP1 and GR", ylab = "log2FC (Dex/DMSO)")
# arrows(x, df_fc4$mean_fc - df_fc4$sd_fc, x, df_fc4$mean_fc + df_fc4$sd_fc,
#        angle = 90, code = 3, length = 0.1, lwd = 1.5)
# axis(side = 1, at = x, labels = x, cex.axis = 0.6)



# # negative control
# df_nc1 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-chr1/K562-DMSO/chr1-b-ap1-b-gr-fc.txt')
# df_nc2 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-chr1/K562-DMSO/chr1-b-ap1-f-gr-fc.txt')
# df_nc3 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-chr1/K562-DMSO/chr1-f-ap1-b-gr-fc.txt')
# df_nc4 = read_df('/hpc/home/yc583/BlueSTARR/GR-AP1/fc-chr1/K562-DMSO/chr1-f-ap1-f-gr-fc.txt')

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
