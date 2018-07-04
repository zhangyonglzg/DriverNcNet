tumor_type <- "GBM";
dir_main <- "/pub5/xiaoyun/Jobs/J19/";
setwd(dir_main);

###1. 获取拷贝数变异谱
deal_GISTIC.file <- function(path_file, name_file){
	temp <- read.table(paste(path_file, name_file,sep="/"), header = T, sep = "\t", quote = "");
	pos1 <- grep("^ENSG00", temp$Gene.Symbol); 	#lncRNA
	pos2 <- grep("^hsa-", temp$Gene.Symbol);	#hsa-mir
	pos3 <- grep("^LOC", temp$Gene.Symbol);	#  去掉没有刻画的转录本
	temp <- temp[c(-pos1,-pos2,-pos3),];
	pos <- which(!duplicated(temp$Locus.ID));
	samples <- sub("TCGA.\\w+.(\\w+).\\d+", "\\1", colnames(temp))[-c(1:3)];
	result <- temp[pos, -c(1:3)];
	rownames(result) <- temp$Locus.ID[pos];
	colnames(result) <- samples;
	return(result);
}

path_file <- paste0(dir_main,tumor_type,"/result/GISTIC");
cnv_M1 <- deal_GISTIC.file(path_file, name_file="all_thresholded.by_genes.txt");
gene_seg <- deal_GISTIC.file(path_file, name_file="all_data_by_genes.txt");
save(cnv_M1,file=paste(dir_main, "/result2/cnv_M1_gene.RData", sep=tumor_type));
save(gene_seg,file=paste(dir_main, "/result2/gene_seg_gene.RData", sep=tumor_type));

###2.去噪
bin_denoising <- function(cnv_M, qvalue){
		#拷贝数高水平变化
		cnv_M[abs(cnv_M) == 1] = 0;
		A <- apply(cnv_M == 2, 1, sum);
		D <- apply(cnv_M == -2, 1, sum);
		index <- which(abs(A - D) >= qbinom(qvalue, A + D, prob=0.5, lower.tail = F));
		result <- cnv_M[index,];
		amp <- apply(result == 2, 1, sum);
		del <- apply(result == -2, 1, sum);
		result1 <- result[(amp - del) > 0,];
		result1[result[(amp - del) > 0,] == -2] <- 0;
		result2 <- result[(amp - del) < 0,];
		result2[result[(amp - del) < 0,] == 2] <- 0;
		result <- rbind(result1,result2);
		result <- sign(result);
		return(result);
}

###3. 离散拷贝数与表达相关
find_cor_wil <- function(gene_id, cnv_M, gene_M, cnv_amp, cnv_del){
	pvalue <- sapply(gene_id, function(x){
				##秩和检验
				group <- as.numeric(cnv_M[x,]);
					pos <- which(group != 0);
					exp_x <- as.numeric(gene_M[x,-pos]);
					exp_y <- as.numeric(gene_M[x,pos]);
					#单侧秩和检验，注意分开考虑
				if(x %in% cnv_amp){	#检验扩增
					w1 <- wilcox.test(exp_y, exp_x, "greater");	
				}
				if(x %in% cnv_del){	#检验缺失
					w1 <- wilcox.test(exp_y, exp_x, "less");		
				}
				p_w <- w1$p.value;
				
			});
	return(pvalue);
}

###2-3 作为一个重要筛选步骤
candidate_M<-function(gene_M, cnv_M, fdr, qvalue){
###通过改变频率和对表达的影响筛选候选基因
###gene_M全基因组表达谱
###cnv_M全基因组的遗传改变谱
###percent最小改变频率阈值
###fdr阈值
		 cnv_M1 <- bin_denoising(cnv_M, qvalue);   ##二项分布去噪
		 cnv_M_temp <- cnv_M1;
		 gene_id <- intersect(rownames(gene_M),rownames(cnv_M1));
		 common_sample<-intersect(colnames(gene_M),colnames(cnv_M1));
		 gene_M1 <- gene_M[gene_id,common_sample];
		 cnv_M2 <- cnv_M1[gene_id,common_sample];
		
		nums <- rowSums(cnv_M2);
		amp <- which(nums > 0);
		cnv_amp <- rownames(cnv_M2)[amp];
		del <- which(nums < 0);
		cnv_del <- rownames(cnv_M2)[del];
		
		gene_id <- c(cnv_amp,cnv_del);
		pvalue <- find_cor_wil(gene_id, cnv_M2, gene_M1, cnv_amp, cnv_del);
		
		pos <- which(pvalue <= fdr);
		gene_id <- gene_id[pos];
		return(cnv_M_temp[gene_id,]);
} 
###############################################################


load(paste(dir_main, "/result2/gene_case_normal_exp.RData", sep=tumor_type));
cnv_M_coding <- candidate_M(gene_M,cnv_M=cnv_M1, fdr=0.05,qvalue=0.05);
save(cnv_M_coding,file=paste(dir_main, "/result2/cnv_M_coding.RData", sep=tumor_type));

###4 拷贝数与突变合并
cnv_bind_mut<-function(cnv_M,mut_M){
	cnv_M <- abs(cnv_M);
   sample_name<-intersect(colnames(cnv_M),colnames(mut_M))
   cnv_gene_ID<-rownames(cnv_M)
   mut_gene_ID<-rownames(mut_M)
   
 ##提取相同样本和相同基因的标签   
   gene_ID<-intersect(cnv_gene_ID,mut_gene_ID)
   cnv_uniq<-setdiff(cnv_gene_ID,gene_ID)
   mut_uniq<-setdiff(mut_gene_ID,gene_ID)
   
   M1 <- cnv_M[gene_ID,sample_name] + mut_M[gene_ID,sample_name];
##将和大于1全部赋值为一  
   M1[M1==2]<-1;
   M2<-cnv_M[cnv_uniq,sample_name]
   M3<-mut_M[mut_uniq,sample_name]
#将三个矩阵并在一起   
   cnv_mut_M<-rbind(M1,rbind(M2,M3))
   pos<-which(rowSums(cnv_mut_M)==0)
   if(length(pos)){
       cnv_mut_M<-cnv_mut_M[-pos,]
   }
  return(cnv_mut_M);
}

load(paste(dir_main, "/result/mut_M.RData", sep=tumor_type)); 	#gene的突变谱 mut_M
colnames(mut_M)<-sapply(colnames(mut_M),function(x){substr(x,9,12)})

cnv_mut_M <- cnv_bind_mut(cnv_M_coding, mut_M);		
pos<-which(rowSums(cnv_mut_M)/ncol(cnv_mut_M)>=0.1)
cnv_mut_M<-cnv_mut_M[pos,]

##5 与其他组学数据交集
find_tf_M <- function(tf_gene, gene_M){
	m<-unique(as.character(tf_gene[,1]))
	n<-intersect(m,rownames(gene_M))
	tf_M<-gene_M[n,]
	return(tf_M);
}

sample_intersect <- function(gene_M, copy_M, meth_M, mir_M, cnv_mut_M){
	con_sample <- intersect(colnames(gene_M), colnames(copy_M));
	con_sample <- intersect(con_sample, colnames(meth_M));
	con_sample <- intersect(con_sample, colnames(mir_M));
	con_sample <- intersect(con_sample, colnames(cnv_mut_M));	#候选gene
	return(con_sample);
}##样本交集

load(paste(dir_main, "/result/M450_promoter_meth_M.RData", sep=tumor_type)); 	
meth_M <- promoter_meth_M;

load(paste(dir_main, "/result/mir_M.RData", sep=tumor_type));

load(paste(dir_main, "/result/tf_gene.RData", sep=tumor_type)); #tf_gene的调控关系
tf_M <- find_tf_M(tf_gene, gene_M);

save(gene_M, gene_seg, meth_M, tf_M, mir_M, mut_M,
	file=paste(dir_main, "/result2/all_data_factors_coding.RData", sep=tumor_type));
	##后面需要的各因子数据
	
con_sample <- sample_intersect(gene_M, gene_seg, meth_M, mir_M, cnv_mut_M);

gene_M <- gene_M[,con_sample];
copy_M <- gene_seg[,con_sample];
meth_M <- meth_M[,con_sample];
mir_M <- mir_M[,con_sample];
tf_M <- tf_M[,con_sample];
cnv_mut_M <- cnv_mut_M[,con_sample];	#候选gene


###6 要求最小变异频数
find_alter_gene<-function(mutation_cnv_M,percent){
 ##通过改变频率筛选基因
 ##mutation_cnv_M全基因组的遗传改变谱
 ##percent最小改变频率阈值
         sample_num<-dim(mutation_cnv_M)[2]
         num_thresh<-floor(sample_num*percent)##最小改变频数
         temp_sum<-abs(rowSums(mutation_cnv_M));
         pos<-which(temp_sum>=num_thresh)
         mutation_cnv_M<-mutation_cnv_M[pos,]
         return(mutation_cnv_M)       
}

num_sample <- 10;	
percent <- num_sample/length(con_sample);
cnv_mut_M <- find_alter_gene(cnv_mut_M,percent);		##样本交集之后，基于变异频率进一步筛选候选基因


####7 要求这些基因必须出现在CGC中
load(paste(dir_main, "/result/cancerGenes.RData", sep=tumor_type))
cancerGenes.CGC<-as.character(cancerGenes.CGC)
jiao_gene<-intersect(cancerGenes.CGC,rownames(cnv_mut_M))
index<-match(jiao_gene,rownames(cnv_mut_M))
cnv_mut_M<-cnv_mut_M[index,]
save(cnv_mut_M,file=paste(dir_main, "/result2/cnv_mut_M.RData", sep=tumor_type));


load(paste(dir_main, "/result/gene_tf_list.RData", sep=tumor_type)); #tf_gene的调控关系
load(paste(dir_main, "/result/gene_mir_list.RData", sep=tumor_type)); #mir_gene的调控关系
save(cnv_mut_M,gene_M,copy_M,meth_M,tf_M,mir_M,gene_tf_list,gene_mir_list, 
	file=paste(dir_main, "/result2/all_data_coding.RData", sep=tumor_type));	


