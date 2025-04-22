scatter_plots = function(df, outfile){
  # pdf(outfile, width=7, height=5)  # ✅ Open PDF inside function
  
  dev.new(width=7, height=5)
  
  par(mfrow=c(1,3), mar=c(3, 3, 1, 1), mgp=c(2, 1, 0))  
  plot(df$rep1, df$rep2, main=paste("rho=", rho_cal(df$rep1, df$rep2)), xlab="Replicate1", ylab="Replicate2")
  abline(a=0, b=1, col="red")
  
  plot(df$rep1, df$rep3, main=paste("rho=", rho_cal(df$rep1, df$rep3)), xlab="Replicate1", ylab="Replicate3")
  abline(a=0, b=1, col="red")
  
  plot(df$rep2, df$rep3, main=paste("rho=", rho_cal(df$rep2, df$rep3)), xlab="Replicate2", ylab="Replicate3")
  abline(a=0, b=1, col="red")
}

density_plots = function(df) {
  dev.new(width=7, height=5)
  par(mfrow=c(1,1))  # Reset layout for a single plot
  dens = density(df$log2_theta)  # Compute density
  plot(dens, main='', xlab="log2 theta", ylab="Density", col="lightblue", lwd=2)
}

rho_cal = function(x,y){
  print(round(cor(x,y,method = "spearman"),2))
  return (round(cor(x,y,method = "spearman"),2))
}

read_data = function(dir){
  # header = c('DNA1','DNA2','DNA3','RNA1','RNA2','RNA3')
  df = read.table(gzfile(dir), sep="\t", header=TRUE)
  print(head(df))
  num_dna = 3
  num_rna = 3
  
  #DNA
  input_rep1 = 123384250
  input_rep2 = 130680200
  input_rep3 = 134539576

  #RNA
  output_rep1 = 116270841
  output_rep2	= 217503220
  output_rep3	= 219086330

  # normalize counts
  c = 10**6
  df$DNA1 = df$input_rep1/input_rep1
  df$DNA2 = df$input_rep2/input_rep2
  df$DNA3 = df$input_rep3/input_rep3
  df$RNA1 = df$output_rep1/output_rep1*c
  df$RNA2 = df$output_rep2/output_rep2*c
  df$RNA3 = df$output_rep3/output_rep3*c
  
  
  df$DNA_aver = (df$DNA1+df$DNA2+df$DNA3)/num_dna
  df$RNA_aver = (df$RNA1+df$RNA2+df$RNA3)/num_rna
  df$rep1 = df$RNA1/df$DNA_aver
  df$rep2 = df$RNA2/df$DNA_aver
  df$rep3 = df$RNA3/df$DNA_aver
  # df$log2_theta = log2(df$RNA_aver/df$DNA_aver)
  # subset = df[1:10000,]
  # return (subset)
  return (df)
}



dir = "/work/igvf-pm/K562/full-set/600bp/K562.combined.input_and_output.w600s50.gt_200.thres_dna_only.log2FC.sequenc.na_removed.txt.gz"
df = read_data(dir)
subset = df[1:100000,]
scatter_plots(subset, outfile)
density_plots(df)


