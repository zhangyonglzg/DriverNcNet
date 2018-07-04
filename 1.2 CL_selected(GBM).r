tumor_type <- "GBM";
dir_main <- "/pub5/xiaoyun/Jobs/J19/";
setwd(dir_main);

###1.读取copy number数据
deal_GISTIC.file <- function(path_file, name_file){
	temp <- read.table(paste(path_file, name_file,sep="/"), header = T, sep = "\t", quote = "");
	pos <- grep("^ENSG00", temp$Gene.Symbol);
	samples <- sub("TCGA.\\w+.(\\w+).\\d+", "\\1", colnames(temp))[-c(1:3)];
	result <- temp[pos, -c(1:3)];
	rownames(result) <- temp$Gene.Symbol[pos];
	colnames(result) <- samples;
	return(result);
}

path_file <- paste0(dir_main,tumor_type,"/result/GISTIC");
cnv_M1 <- deal_GISTIC.file(path_file, name_file="all_thresholded.by_genes.txt");
gene_seg1 <- deal_GISTIC.file(path_file, name_file="all_data_by_genes.txt");
save(cnv_M1,file=paste(dir_main, "/result2/cnv_M1_lnc.RData", sep=tumor_type));
save(gene_seg1,file=paste(dir_main, "/result2/gene_seg_lnc.RData", sep=tumor_type));

###2. 对离散拷贝数数据进行去噪
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

###3.去除低表达的lncRNA
lncRNA_expr_deal <- function(file_path, file_name, thr){
###在TANRIC数据库中下载的lncRNA表达谱,并将30%的样本低于0.3的lncRNA剔除
	setwd(file_path);
	expr <- read.table(file_name, header = T, sep = "\t", quote = "\"'");
	rownames(expr) <- sub("\\.\\d+", "", expr[,1], fixed=FALSE);	#去掉lncRNA的版本号
	colnames(expr) <- sub(".*TCGA.\\w+.(\\w+)", "\\1", colnames(expr));
	expr <- expr[,-1];
	expr_thr <- apply(expr >= 0.3, 1, sum);
	index <- which(expr_thr/dim(expr)[2] >= thr);

	expr <- expr[index,];
	expr <- as.matrix(expr);
	return(expr);
}

###4. 表达与离散copy number做相关
find_cor_wil <- function(gene_id, cnv_M, gene_M, cnv_amp, cnv_del){
	pvalue <- sapply(gene_id, function(x){
			group <- as.numeric(cnv_M[x,]);
			pos <- which(group != 0);
			exp_x <- as.numeric(gene_M[x,-pos]);
			exp_y <- as.numeric(gene_M[x,pos]);
		if(x %in% cnv_amp){	#检验扩增
			w1 <- wilcox.test(exp_y, exp_x, "greater");	 #秩和检验
		}
		if(x %in% cnv_del){	#检验缺失
			w1 <- wilcox.test(exp_y, exp_x, "less");		#秩和检验
		}
		p_w <- w1$p.value;
						
	});
	return(pvalue);
}
	
###2-4	三步一起作为前期筛选lnc的重要步骤
candidate_M <- function(file_path,file_name,cnv_M,thr,qvalue,fdr){
###gene_M全基因组表达谱
###cnv_M全基因组的遗传改变谱, mutation_cnv_M
###cnv_value为gene拷贝数连续值
###percent最小改变频率阈值
###fdr阈值
	cnv_M1 <- bin_denoising(cnv_M, qvalue);   ##二项分布去噪
	cnv_M_temp <- abs(cnv_M1);
	gene_M <- lncRNA_expr_deal(file_path, file_name, thr=0.3);  ##将表达均值低于0.3的lncRNA剔除
	gene_id <- intersect(rownames(gene_M),rownames(cnv_M1));
	common_sample <- intersect(colnames(gene_M),colnames(cnv_M1));
	cnv_M2 <- cnv_M1[gene_id,common_sample];
	gene_M1 <- gene_M[gene_id,common_sample];
	
	nums <- rowSums(cnv_M2);
	amp <- which(nums > 0);
	cnv_amp <- rownames(cnv_M2)[amp];
	del <- which(nums < 0);
	cnv_del <- rownames(cnv_M2)[del];
	
	gene_id <- c(cnv_amp,cnv_del);
	pvalue <- find_cor_wil(gene_id, cnv_M2, gene_M1, cnv_amp, cnv_del);

	pos <- which(pvalue <= fdr)
	gene_id <- gene_id[pos];
	return(cnv_M_temp[gene_id,]);
}  

file_path <- paste(dir_main, "/result/", sep=tumor_type)
file_name <- paste("TCGA-", "-rnaexpr.tsv", sep=tumor_type);
lncRNA_expr <- lncRNA_expr_deal(file_path, file_name, thr=0.3);		#lnc的表达
save(lncRNA_expr,file=paste(dir_main, "/result2/lncRNA_expr.RData", sep=tumor_type));
cnv_M_lnc <- candidate_M(file_path, file_name, cnv_M1, thr=0.3, qvalue=0.05, fdr=0.05);
###############################################################
pos<-which(rowSums(cnv_M_lnc)/ncol(cnv_M_lnc)>=0.1)
cnv_M_lnc<-cnv_M_lnc[pos,]

###5.与其他组学数据交集
find_tf_M <- function(tf_gene, gene_M){
	m<-unique(as.character(tf_gene[,1]))
	n<-intersect(m,rownames(gene_M))
	tf_M<-gene_M[n,]
	return(tf_M);
}

##样本交集
sample_intersect <- function(gene_M, copy_M, meth_M, mir_M, cnv_M_lnc){
	con_sample <- intersect(colnames(gene_M), colnames(copy_M));
	con_sample <- intersect(con_sample, colnames(meth_M));
	con_sample <- intersect(con_sample, colnames(mir_M));
	con_sample <- intersect(con_sample, colnames(cnv_M_lnc));	#候选gene
	return(con_sample);
}##样本交集

load(paste(dir_main, "/result/gene_case_normal_exp.RData", sep=tumor_type));
load(paste(dir_main, "/result/M450_promoter_meth_M.RData", sep=tumor_type)); 	
meth_M <- promoter_meth_M;

load(paste(dir_main, "/result/mir_M.RData", sep=tumor_type));	#需要miRNA的表达谱

load(paste(dir_main, "/result/tf_gene.RData", sep=tumor_type)); #tf_gene的调控关系
tf_M <- find_tf_M(tf_gene, gene_M);

####这里是基因的copy number，不是lnc本身的
load(paste(dir_main, "/result2/gene_seg_gene.RData", sep=tumor_type));

save(gene_M, gene_seg, meth_M, tf_M, mir_M, lncRNA_expr,
		file=paste(dir_main, "/result2/all_data_factors_lnc.RData", sep=tumor_type));
	##后面需要的各因子数据
	
con_sample <- sample_intersect(gene_M, gene_seg, mir_M, meth_M, cnv_M_lnc);

gene_M <- gene_M[,con_sample];
copy_M <- gene_seg[,con_sample];
meth_M <- meth_M[,con_sample];
mir_M <- mir_M[,con_sample];
tf_M <- tf_M[,con_sample];
cnv_M_lnc <- cnv_M_lnc[,con_sample];	#候选lnc


###6 要求最小变异频数
find_alter_gene<-function(mutation_cnv_M,percent){
 ##通过改变频率筛选基因(先去噪，再算根据最小改变频率)
 ##mutation_cnv_M全基因组的遗传改变谱
 ##percent最小改变频率阈值
         sample_num<-dim(mutation_cnv_M)[2]
         num_thresh<-floor(sample_num*percent)##最小改变频数
         temp_sum<-abs(rowSums(mutation_cnv_M));
         pos<-which(temp_sum>=num_thresh)
         mutation_cnv_M<-mutation_cnv_M[pos,]
         return(mutation_cnv_M)       
}

num_sample <- 10;	#变异样本数
percent <- num_sample/length(con_sample);
cnv_M_lnc <- find_alter_gene(cnv_M_lnc,percent);		##样本交集之后，基于变异频率进一步筛选候选基因	 
save(cnv_M_lnc, file=paste(dir_main, "/result2/cnv_M_lnc.RData", sep=tumor_type));

load(paste(dir_main, "/result/gene_tf_list.RData", sep=tumor_type)); #tf_gene的调控关系
load(paste(dir_main, "/result/gene_mir_list.RData", sep=tumor_type)); #mir_gene的调控关系
save(cnv_M_lnc,gene_M,copy_M,meth_M,tf_M,mir_M,gene_tf_list,gene_mir_list, 
		file=paste(dir_main, "/result2/all_data_lnc.RData", sep=tumor_type));		  