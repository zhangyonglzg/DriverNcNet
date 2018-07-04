
get_plrs_result<-function(x){
	for(i in x){
	#################################################################################
	###针对某一个gene来说，提取变异样本和野生型样本
	###@gene_factor_profile ：某个gene的因子谱，是一个矩阵，行是样本，列是因子
	###@case_list ：变异样本
	###@normal_list ：野生型样本
	###结果为长度为2的list，分别为变异样本和野生样本
	case_normal_profile<-function(gene_factor_profile,case_list,normal_list){

			   result<-list();
			   result[[1]]<-gene_factor_profile[case_list,];
			   result[[2]]<-gene_factor_profile[normal_list,];
			
			   names(result)<-c("case","normal");
			   return(result);
	}
	#################################################################################
	###对所有gene，提取变异样本和野生样本，根据候选gene拷贝数变异与否将样本分为两组，变异组与野生组
	###@gene_id ：候选geneID(可以是所有的)
	###@ cnv_mutation_profile：候选gene的拷贝数变异普
	###@ gene_factor_more_list ：这是一个以全基因组所有gene的ID命名的list，每个元素为一个矩阵，
	###行是样本，列代表该gene的因子
	###% 结果是一个以全基因组所有gene的ID命名的list，每个元素为一个长度为2的list，即该gene在两组样本的因子普
	###行是样本，列代表该gene的因子
	altered_gene_case_normal_profile<-function(gene_id,cnv_mutation_profile,gene_factor_more_list){
		
		  normal_pos<-which(cnv_mutation_profile[gene_id,]==0);
		  normal_case<-list();
		  normal_case[[1]]<-colnames(cnv_mutation_profile[,-normal_pos]);
		  normal_case[[2]]<-colnames(cnv_mutation_profile[,normal_pos]);
		  result<-lapply(gene_factor_more_list,case_normal_profile,case_list=normal_case[[1]],normal_list=normal_case[[2]]);
		  return(result);
	}
									
	###########################################################################
	#############################################################################

	##2.计算因子对基因表达的贡献度 
	################################################################################
	###!  建立最小二乘回归，提取出下游gene在回归中的回归系数及其显著性p值
	###@gene_factor_profile：某一个gene的因子普，行是样本，列是该基因的表达及因子
	###@pthr：p的阈值
	###% 返回结果是长度为2的list，分别为该gene的回归系数及p值
	library("pls");
	plrs_find_factor<-function(gene_factor_profile,pthr){
		 num_d<-dim(gene_factor_profile)[2]
		 a<-data.frame(y=gene_factor_profile[,1]);
		 a$x<-as.matrix(gene_factor_profile[,2:num_d]);
		 
		 mod<- plsr(y ~ x, data = a, validation = "CV", jackknife = TRUE, scale=T);
			##计算最小二乘回归
		 ncomp<-mod$ncomp 	##提取所有的主成分个数
		 adjcv<-(RMSEP(mod)$val)[2,1,2:(ncomp+1)]; ##通过计算每个主成分的adjcv值，确定最佳主成分个数，adjcv值越小越好
		 pos<-which(adjcv=="NaN")
		 if(length(pos)){
		 adjcv[pos]<-max(adjcv[-pos])
		 }
		 min_adjcv<-min(adjcv)
		 ncomp<-which(adjcv==min_adjcv)

		 jtest<-jack.test(mod,ncomp = ncomp);
			#使用jack.test计算系数以及p值
		 variable_coef<-(jtest$coefficients)[1:(num_d-1),1,1];
		 pvalue<-(jtest$pvalue)[1:(num_d-1),1,1];
		 pos<-which(pvalue<=pthr);

		 if(length(pos)){
			result<-list(variable_coef[pos],pvalue[pos]);
			names(result)<-c("coefficients","pvalue")
			 return(result)
		 }else{
			 result<-list(variable_coef,pvalue);
			 names(result)<-c("coefficients","pvalue")
			 return(result)
		 }
	}
	####################################################################

	#################################################################################
	###! 对某个gene，构建回归模型，提取出其分别在两组样本中的回归系数及p值
	###@case_normal_gene_factor_profile ：长度为2的list，代表某个gene在两组样本的因子普
	###@pthr ：阈值
	###%返回结果是一个长度为2的list，即某gene在两组样本中的回归系数及p值(2层list)
	gene_case_normal_find_factor<-function(case_normal_gene_factor_profile,pthr){

			result<-lapply(case_normal_gene_factor_profile,plrs_find_factor,pthr=pthr)
			names(result)<-c("case","normal")
			return(result)

	}

	#########################################################################
	###!all_gene_case_normal_find_factor ：对所有gene，构建回归模型，提取出其分别在两组样本中的回归系数及p值
	###@case_normal_allgene_factor_profile ：所有gene的因子普
	###%返回结果是一个长度为所有gene的list，即某gene在两组样本中的回归系数及p值(3层list)
	all_gene_case_normal_find_factor<-function(case_normal_allgene_factor_profile,pthr){

			result<-lapply(case_normal_allgene_factor_profile,gene_case_normal_find_factor,pthr=pthr)
			names(result)<-names(case_normal_allgene_factor_profile)
			return(result)
	}
											

	##########################################################################################
	### convert the list to matrix	
	### @a_list只有一层
	list_p_c<-function(a_list){ ##convert single component such as normal
		p_c<-data.frame(factors=names(a_list[[1]]),coefs=a_list$coefficients,pvalue=a_list$pvalue)
		return(p_c)
	}
	#list_p_c1 <- list_p_c(a_list=all_gene_case_normal_find_factor1[[1]][[1]]);
	####################################
	###convert a list containing case and normal
	###两层list，只有一个基因内层，第一层是一个基因的case与normal
	case_control_p_c<-function(a_list){

		 c_n_c_p_M<-c();
		 c_c_p_c<-lapply(a_list,list_p_c)
		 c_n_c_p_M<-cbind(c_c_p_c[[1]],c_c_p_c[[2]][,c(2,3)])
		 colnames(c_n_c_p_M)<-c("factors","c_c","c_p","n_c","n_p")
		 return(t(c_n_c_p_M))

	}
	#case_control_p_c1 <- case_control_p_c(a_list=all_gene_case_normal_find_factor1[[1]]); 

	#####################################
	coln_n<-function(M){## calculate the number of factors
		  return(dim(M)[2])
	}
	#####################################
	###convert the list containing multiple genes to matrix
	###三层list，所一个gene外层
	list2M<-function(P_C_result){
		F_P_C<-lapply(P_C_result,case_control_p_c)
		l<-unlist(lapply(F_P_C,coln_n))
		gene_name<-rep(names(F_P_C),l)
		result<-t(rbind(gene_name,matrix(unlist(F_P_C),nrow=5)))
		colnames(result)<-c("gene_id","factors","c_c","c_p","n_c","n_p")
		return(result)
	}
	############################################################################


	#############################################################################
	###一个候选gene的ID
	###lnc的表达谱
	##cnv_mutation_profile ：候选gene的拷贝数变异普
	##gene_factor_list ：所有gene的因子普
	## add_flag为TRUE时，表示加入lnc的表达值作为其中一个factor
	plsr_c_p<-function(gene_id,cnv_mutation_profile,gene_factor_list){
	  ##form the gene_factor_profile for all genes in the gene_id case samples and gene_id normal sample
		 case_control_profile<-altered_gene_case_normal_profile(gene_id,
												cnv_mutation_profile=cnv_mutation_profile,
												gene_factor_more_list=gene_factor_list);
 
		 P_C_result<-all_gene_case_normal_find_factor(case_control_profile,1);
	   ##convert the list result in to matrix result   
		 C_N_C_P<-list2M(P_C_result)
		 return(C_N_C_P)
	}

	lnc_id<-rownames(cnv_M_lnc)[i]
	C_N_C_P_temp<- plsr_c_p(gene_id=lnc_id,cnv_mutation_profile=cnv_M_lnc,
							gene_factor_list=gene_factor_list2)
	save_file<-paste0("/pub5/xiaoyun/Jobs/J19/GBM_revise/lnc_result/true_result/",i,".RData");
	save(C_N_C_P_temp, file=save_file);
	}
}

