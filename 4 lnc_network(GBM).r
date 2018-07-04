main_path<-"/pub5/xiaoyun/Jobs/J19/"
cancer <- "GBM"; 
path1<-paste0(main_path,cancer,"/gene_result")
path2<-paste0(main_path,cancer,"/lnc_result")
path3<-paste0(main_path,cancer,"/analysis")
setwd(path2)
load("lnc_downstream_effectors_new.RData")

#######
load(paste0(main_path,cancer,"/result/mir_lnc_gene.RData")
	#ceRNA
load(paste0(main_path,cancer,"/result2/all_data_lnc.RData"));
load(paste0(main_path,cancer,"/result2/lncRNA_expr.RData"));
	#各因子谱的数据

lnc_gene_con <- intersect(colnames(lncRNA_expr), colnames(gene_M));
lnc_mir_con <- intersect(colnames(lncRNA_expr), colnames(mir_M));
lnc_meth_con <- intersect(colnames(lncRNA_expr), colnames(meth_M));
##

####第三步  识别候选driver lncRNA 的下游功能效应子

lnc_dfes <- lapply(1:length(lnc_downstream_effectors_new), function(x){
	lncRNA <- names(lnc_downstream_effectors_new)[x];
	des<-lnc_downstream_effectors_new[[x]]
	pos1<-grep("hsa",des)
	pos2<-grep("meth",des)
	pos3<-setdiff(1:length(des),union(pos1,pos2))
	###1.对miRNA进行挑选
	mir<-des[pos1]
	mir_new<-intersect(mir,names(mir_lnc_gene))
	ceRNA <- mir_lnc_gene[mir_new];
	flag<-sapply(ceRNA,function(y){
		(lncRNA%in%y[[1]])
	})
	if(any(flag)){
		index<-which(flag)
		mir_new_final<-abs(sapply(mir_new[index],function(y){
			cor(lncRNA_expr[lncRNA,lnc_mir_con], mir_M[y,lnc_mir_con]);
		}))
	}else{
	mir_new_final<-NULL
	}
	###2 对甲基化基因进行挑选
	meth_gene<-des[pos2]
	meth_gene_new<-sapply(meth_gene,function(y){unlist(strsplit(y,"meth"))})
	names(meth_gene_new)<-NULL
	lnc_meth <- t(sapply(meth_gene_new,function(y){
		a<-cor.test(lncRNA_expr[lncRNA,lnc_meth_con], meth_M[y,lnc_meth_con]);
		c(a$estimate,a$p.value)
	}))
	flag_meth<-apply(lnc_meth,1,function(y){
	     f1<-(y[2]<=0.05)
		 f2<-(abs(y[1])>=0)
		 all(c(f1,f2))
	})
	pos_meth<-which(flag_meth)
	if(length(pos_meth)>0){
		meth_gene_final<-abs(sapply(meth_gene_new[pos_meth],function(y){
				cor(lncRNA_expr[lncRNA,lnc_meth_con], meth_M[y,lnc_meth_con]);
		}))
	}else{
		meth_gene_final<-NULL
	}
	
	###3 对基因和tf进行挑选
	ord_gene<-des[pos3]
	lnc_gene <- t(sapply(ord_gene,function(y){
		a<-cor.test(lncRNA_expr[lncRNA,lnc_gene_con], gene_M[y,lnc_gene_con]);
		c(a$estimate,a$p.value)
	}))
	flag_gene<-apply(lnc_gene,1,function(y){
	     f1<-(y[2]<=0.05)
		 f2<-(abs(y[1])>=0)
		 all(c(f1,f2))
	})
	pos_gene<-which(flag_gene)
	if(length(pos_gene)>0){
		ord_gene_final<-abs(sapply(ord_gene[pos_gene],function(y){
				cor(lncRNA_expr[lncRNA,lnc_gene_con], gene_M[y,lnc_gene_con]);
		}))
		}else{
		ord_gene_final<-NULL
	}
	return(c(mir_new_final,meth_gene_final,ord_gene_final))
})
	
names(lnc_dfes)<-names(lnc_downstream_effectors_new)
setwd(path3)
save(lnc_dfes,file="lnc_dfes.RData");


####第四步  使用网络相似测度识别driver lncRNA 
load(paste0(path3,"/response_genes_diffK.RData"))
response_genes<-response_genes_diffK
index<-which(sapply(response_genes,length)>=3)
response.genes<-response_genes[index]
library(igraph)
load(paste0(main_path,cancer,"/PPI_string.Rdata"))
PPI<-as.data.frame(PPI,stringsAsFactors=F)
PPI[,3]<-as.numeric(PPI[,3])/1000
g<-graph.data.frame(PPI[,1:2],directed=FALSE)
g1<-set.graph.attribute(g, "weight", PPI[,3])
E(g1)$"weight"<-PPI[,3]
g2<-set.graph.attribute(g, "weight", 1-PPI[,3])

###获取邻接矩阵
M<-get.adjacency(g1, attr="weight")
#sp_matrix_all<-shortest.paths(g2,weights=1-PPI[,3])
###获取随机所需要的背景的节点
gene_ppi<-V(g)$name  ####可能有偏倚
load(paste0(path1,"/gene_factor_list2.RData"))
gene_all<-intersect(names(gene_factor_list2),gene_ppi)
###不进行模块选择，直接使用lnc效应子
lnc_genes<-sapply(lnc_dfes,function(x){
	intersect(names(x),gene_ppi)
	})
index2<-which(sapply(lnc_genes,length)>0)
lnc_genes<-	lnc_genes[index2]
##4.1 
###随机游走，使用lncRNA的下游基因作为种子，基因与lnc表达相关性作为种子基因权重，而driver基因的下游基因作为响应基因，随机游走之后，响应基因的概率和作为统计量与随机比较
get_random_walk_p<-function(k,random_time=1000){

	seed.genes<-lnc_dfes[[k]]
	index<-which(names(seed.genes)%in%lnc_genes[[k]])
	seed.genes<-seed.genes[index]

	goal<-0.000001
	r<-0.3
	RSA(seed.genes,gene_ppi,M,goal,r)
	true_result<-RSE(IP("RandomWalk",seed.genes,gene_ppi,M,goal,r))
	true_result_new<-rev(sort(unlist(true_result)))
	
	true_score<-lapply(response.genes,function(temp_gene){
		pos<-na.omit(match(temp_gene,names(true_result_new)))
		# index<-which(pos<=500)
		# name_id[names(true_result_new)[pos][index],1]
		rank_op_sum<-sum(1/pos)
	})
	
	####随机情况下响应基因的概率与真实情况下比较产生p值
	random_score_all<-sapply(1:random_time,function(x){
		names(seed.genes)<-sample(gene_all,length(seed.genes))
		RSA(seed.genes,gene_ppi,M,goal,r)
		random_result<-RSE(IP("RandomWalk",seed.genes,gene_ppi,M,goal,r))
		random_result_new<-rev(sort(unlist(random_result)))
		temp_score<-sapply(response.genes,function(temp_gene){
			pos<-na.omit(match(temp_gene,names(random_result_new)))
			rank_op_sum<-sum(1/pos)
		})
	})
	p_value<-sapply(1:length(true_score),function(x){
		length(which(random_score_all[x,]>=true_score[x]))/random_time
	})
	names(p_value)<-names(response.genes)
	p_value
}

random_walk_p_diffK<-t(sapply(1:length(lnc_genes),get_random_walk_p,random_time=100))
rownames(random_walk_p_diffK)<-names(lnc_genes)
setwd(path3)
save(random_walk_p_diffK,file="random_walk_p_diffK.RData")
###


####4.2
##将driver 基因的下游基因集合标记为gene1，driver lnc的下游基因集合标记为gene2
###在网络中计算任意一个lnc下游基因形成模块与基因的下游基因形成模块的疾病模块overlap分值
library(igraph)
load(paste0(main_path,cancer,"/sp_matrix_all.RData"))
m<-sp_matrix_all
get_overlap_score<-function(gene1,gene2){ 
	###将driver 基因的下游基因集合标记为gene1，driver lnc的下游基因集合标记为gene2
	pos1<-na.omit(match(gene1,rownames(m)))
	pos2<-na.omit(match(gene2,rownames(m)))

	Daa_ori<-sapply(pos1,function(x){
		temp<-m[x,setdiff(pos1,x)]
		return(min(temp))
	})
	Daa<-mean(Daa_ori[Daa_ori!="Inf"])
	
	Dbb_ori<-sapply(pos2,function(x){
		temp<-m[x,setdiff(pos2,x)]
		return(min(temp))
	})
	Dbb<-mean(Dbb_ori[Dbb_ori!="Inf"])
	
	
	Dab1_ori<-sapply(pos1,function(x){
		temp<-m[x,pos2]
		return(min(temp))
	})
	Dab1<-Dab1_ori[Dab1_ori!="Inf"]
	
	Dab2_ori<-sapply(pos2,function(x){
		temp<-m[x,pos1]
		return(min(temp))
	})
	Dab2<-Dab2_ori[Dab2_ori!="Inf"]
	
	Dab<-mean(c(Dab1,Dab2))
	Sab<-Dab-(Daa+Dbb)/2
	Sab
}

get_overlap_score_p<-function(gene1,gene2,random_time=1000){ 
	#gene1<-response.genes[[1]]
	#gene2<-lnc_genes[[2]]
	true_score<-get_overlap_score(gene1,gene2)
	random_score_all<-sapply(1:random_time,function(x){
		gene2_rand<-sample(gene_all,length(gene2))
		temp_score<-get_overlap_score(gene1,gene2_rand)	
	})
	p_value<-length(which(random_score_all<=true_score))/random_time###此处方向别弄反了
	return(p_value)
}

get_overlap_score(response.genes[[1]],lnc_genes[[2]])
get_overlap_score_p(response.genes[[1]],lnc_genes[[2]],random_time=10)
###如果得分为负或者说越小，认为两者功能联系比较紧密
###由于0作为阈值，太理想化，使用随机的方式来确定两者关系是否关联紧密

module_associated_p_diffK<-sapply(response.genes,function(x){
	sapply(lnc_genes,function(y){
			get_overlap_score_p(x,y,random_time=1000)
	})
})
setwd(path3)
save(module_associated_p_diffK,file="module_associated_p_diffK.RData")



####第五步：获取最终的driver lncRNA

setwd(path3)
load("lnc_dfes.RData")
load("random_walk_p_diffK.RData")
load("module_associated_p_diffK.RData")
dim(module_associated_p_diffK)
dim(random_walk_p_diffK)


####两种方法得到结果的并集
result1<-apply(random_walk_p_diffK,1,function(x){ names(which(x<=0.05))})
result2<-apply(module_associated_p_diffK,1,function(x){names(which(x<=0.05))})

final_result<-sapply(1:nrow(module_associated_p_diffK),function(i){union(result1[[i]],result2[[i]])})
names(final_result)<-rownames(module_associated_p_diffK)


driver_lnc<-names(final_result)[which(sapply(final_result,length)>0)]
driver_lnc_dfes<-sapply(lnc_dfes[names(driver_lnc)],names)

save(driver_lnc,file="driver_lnc.RData")
save(driver_lnc_dfes,file="driver_lnc_dfes.RData")