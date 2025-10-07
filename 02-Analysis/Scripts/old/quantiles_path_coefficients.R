a <- read.table("C:/Users/mbeneat/Desktop/path coefficients.csv", row.names = c(1), header = c(2), sep=",", dec=".")


quantile_jackall_notstd <- quantile(as.numeric(a[-c(1:2),1]))
quantile_jackRMR0_notstd <- quantile(as.numeric(a[-c(1:2),2]))
quantile_phylosem_notstd <- quantile(as.numeric(a[-c(1:2),3]))
quantile_jackall_std <- quantile(as.numeric(a[-c(1:2),4]))
quantile_jackRMR0_std <- quantile(as.numeric(a[-c(1:2),5]))
# quantile_phylosem_std <- quantile(as.numeric(a[-c(1:2),1]))


quantile_jackall_notstd <- quantile(as.numeric(a[grep("c_m",rownames( a)),1]))
quantile_jackRMR0_notstd <- quantile(as.numeric(a[grep("c_m",rownames( a)),2]))
quantile_phylosem_notstd <- quantile(as.numeric(a[grep("c_m",rownames( a)),3]))
quantile_jackall_std <- quantile(as.numeric(a[grep("c_m",rownames( a)),4]))
quantile_jackRMR0_std <- quantile(as.numeric(a[grep("c_m",rownames( a)),5]))
