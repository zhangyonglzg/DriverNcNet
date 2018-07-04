##4 �ܲ���
get_plrs_result<-function(x){#x��ʾҪ����Ļ������
	for(i in x){
##1.����������ĺ�ѡdriver gene��Ž������ֳɱ���������Ұ��������
	case_normal_profile<-function(gene_factor_profile,case_list,normal_list){
		result<-list();
		result[[1]]<-gene_factor_profile[case_list,];
		result[[2]]<-gene_factor_profile[normal_list,];
			
		names(result)<-c("case","normal");
		return(result);
	}
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


##2.�������ӶԻ�����Ĺ��׶� 
	################################################################################
	###!  ������С���˻ع飬��ȡ������gene�ڻع��еĻع�ϵ������������pֵ
	###@gene_factor_profile��ĳһ��gene�������գ��������������Ǹû���ı�Ｐ����
	###@pthr��p����ֵ
	###% ���ؽ���ǳ���Ϊ2��list���ֱ�Ϊ��gene�Ļع�ϵ����pֵ
	library("pls");
	plrs_find_factor<-function(gene_factor_profile,pthr){
		 num_d<-dim(gene_factor_profile)[2]
		 a<-data.frame(y=gene_factor_profile[,1]);
		 a$x<-as.matrix(gene_factor_profile[,2:num_d]);
		 
		 mod<- plsr(y ~ x, data = a, validation = "CV", jackknife = TRUE, scale=T);
			##������С���˻ع�
		 ncomp<-mod$ncomp 	##��ȡ���е����ɷָ���
		 adjcv<-(RMSEP(mod)$val)[2,1,2:(ncomp+1)]; ##ͨ������ÿ�����ɷֵ�adjcvֵ��ȷ��������ɷָ�����adjcvֵԽСԽ��
		 pos<-which(adjcv=="NaN")
		 if(length(pos)){
		 adjcv[pos]<-max(adjcv[-pos])
		 }
		 min_adjcv<-min(adjcv)
		 ncomp<-which(adjcv==min_adjcv)

		 jtest<-jack.test(mod,ncomp = ncomp);
			#ʹ��jack.test����ϵ���Լ�pֵ
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
	################################################################################

	#################################################################################
	###! ��ĳ��gene�������ع�ģ�ͣ���ȡ����ֱ������������еĻع�ϵ����pֵ
	###@case_normal_gene_factor_profile ������Ϊ2��list������ĳ��gene������������������
	###@pthr ����ֵ
	###%���ؽ����һ������Ϊ2��list����ĳgene�����������еĻع�ϵ����pֵ(2��list)
	gene_case_normal_find_factor<-function(case_normal_gene_factor_profile,pthr){

			result<-lapply(case_normal_gene_factor_profile,plrs_find_factor,pthr=pthr)
			names(result)<-c("case","normal")
			return(result)

	}
	
	#########################################################################
	###!all_gene_case_normal_find_factor ��������gene�������ع�ģ�ͣ���ȡ����ֱ������������еĻع�ϵ����pֵ
	###@case_normal_allgene_factor_profile ������gene��������
	###%���ؽ����һ������Ϊ����gene��list����ĳgene�����������еĻع�ϵ����pֵ(3��list)
	all_gene_case_normal_find_factor<-function(case_normal_allgene_factor_profile,pthr){

			result<-lapply(case_normal_allgene_factor_profile,gene_case_normal_find_factor,pthr=pthr)
			names(result)<-names(case_normal_allgene_factor_profile)
			return(result)
	}
	################################################################################
																									

	##########################################################################################
	### convert the list to matrix	
	### @a_listֻ��һ��
	list_p_c<-function(a_list){ ##convert single component such as normal
		p_c<-data.frame(factors=names(a_list[[1]]),coefs=a_list$coefficients,pvalue=a_list$pvalue)
		return(p_c)
	}
	#list_p_c1 <- list_p_c(a_list=all_gene_case_normal_find_factor1[[1]][[1]]);

	case_control_p_c<-function(a_list){

		 c_n_c_p_M<-c();
		 c_c_p_c<-lapply(a_list,list_p_c)
		 c_n_c_p_M<-cbind(c_c_p_c[[1]],c_c_p_c[[2]][,c(2,3)])
		 colnames(c_n_c_p_M)<-c("factors","c_c","c_p","n_c","n_p")
		 return(t(c_n_c_p_M))

	}

	#####################################
	coln_n<-function(M){## calculate the number of factors
		  return(dim(M)[2])
	}
	#####################################
	###convert the list containing multiple genes to matrix
	###����list����һ��gene���
	list2M<-function(P_C_result){
		F_P_C<-lapply(P_C_result,case_control_p_c)
		l<-unlist(lapply(F_P_C,coln_n))
		gene_name<-rep(names(F_P_C),l)
		result<-t(rbind(gene_name,matrix(unlist(F_P_C),nrow=5)))
		colnames(result)<-c("gene_id","factors","c_c","c_p","n_c","n_p")
		return(result)
	}

	###�����ع鲢�ҽ������table��ʽչʾ
	############################################################################
	plsr_c_p<-function(gene_id,gene_M,cnv_mutation_profile,gene_factor_list){
		 case_control_profile<-altered_gene_case_normal_profile(gene_id,
												cnv_mutation_profile=cnv_mutation_profile,
												gene_factor_more_list=gene_factor_list);
	 
		 P_C_result<-all_gene_case_normal_find_factor(case_control_profile,1); 
		 C_N_C_P<-list2M(P_C_result)
		 return(C_N_C_P)
	}
	#############################################################################
	gene_id<-rownames(cnv_mut_M)[i]
	C_N_C_P_temp<-plsr_c_p(gene_id=gene_id,cnv_mutation_profile=cnv_mut_M,gene_factor_list=gene_factor_list2);
	save_file<-paste0("/pub5/xiaoyun/Jobs/J19/GBM_revise/gene_result/true_result/",i,".RData");
    save(C_N_C_P_temp,file=save_file);
	}
}

