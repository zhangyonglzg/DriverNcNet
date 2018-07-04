main_path<-"/pub5/xiaoyun/Jobs/J19/"
cancer<-"GBM"  
load(paste0(main_path,cancer,"/result2/all_data_lnc.RData"));
path0<-paste0(main_path,cancer,"/")
path1<-paste0(main_path,cancer,"/gene_result")
path2<-paste0(main_path,cancer,"/lnc_result")


####1 构建并行随机所需要的数据池
setwd(path1)
load("gene_factor_list2.RData")
##如果所有基因的所有因子一起校正，随机需要用全部的基因
lnc_diff_gene<-lapply(rownames(cnv_M_lnc),function(x){names(gene_factor_list2)})
names(lnc_diff_gene)<-rownames(cnv_M_lnc)

# ##如果单基因内fdr校正，随机下游基因可以挑选
#setwd(path2)
# load("lnc_diff_relation_list_new1.RData")
# lnc_diff_gene<-sapply(lnc_diff_relation_list_new,function(x){as.character(unique(x[,1]))})

###构建随机的基因组变异谱（针对候选的driver基因）
library(BiRewire)
cnv_M_lnc_random<-lapply(1:1000,function(x){
	birewire.rewire.bipartite(cnv_M_lnc,verbose=FALSE)
})
setwd(path2)
save(cnv_M_lnc_random,file="cnv_M_lnc_random.RData")
###构建并行所需要的索引表
big_table<-cbind(cbind(rep(names(lnc_diff_gene),each=1000),rep(1:nrow(cnv_M_lnc),each=1000)),rep(1:1000,nrow(cnv_M_lnc)))
setwd(path2)
save(big_table,file="big_table.RData")

##一共11个数据，都不能少，注意加上 lnc_diff_gene,cnv_M_lnc_random
save(cnv_M_lnc,copy_M,gene_M,meth_M,gene_mir_list,mir_M,gene_tf_list,tf_M,gene_factor_list2,lnc_diff_gene,cnv_M_lnc_random,file="random_data_required1.RData")
###################################################################################


##2. 随机并行开始
##############################################################################
source(paste0(path0,"/lnc_parallel_random.r"))
load(paste0(path2,"/big_table.RData"))
RSA(big_table);
RSA(get_plrs_result);
interIndex <-lapply(1:(nrow(big_table)/10),function(x){
	start_pos<-(x-1)*10+1;
	end_pos<-x*10;
	c(start_pos:end_pos)
})
RSA(interIndex) ;
file.RData<-paste0(path2,"random_data_required1.RData")
RSA(file.RData) ;
RSE(IPP(interIndex,get_plrs_result,fire="weak", file.RData=file.RData, numThred.weak=64))
###############################################################
	

####3. 获取并行随机并进行识别显著差异调控边
##############################################################	
load(paste0(path2,"/lnc_diff_relation_list_new1.RData"))

setwd(paste0(path2,"/random_result"));
my_file<-dir()
my_file_new<-sapply(my_file,function(x){unlist(strsplit(x,"\\.RData"))})
my_file_new_order<-sapply(sort(as.numeric(my_file_new)),function(x){paste0(x,".RData")})


##### 使用点差异网络的方法
fdr<-function(p,fdr_flag){ ##calculate the fdr
   p<-p.adjust(p,method=fdr_flag)
   return(p)
}

get_diffK_score<-function(C_N_C_P,fdr_flag,thr){
	   ## C_N_C_P=C_N_C_P_lnc_true[[7]]
	   ###先对调控p值进行校正，不显著的调控边设为0
	   fac<-factor(C_N_C_P[,1],level=unique(C_N_C_P[,1]))
       p1<-unlist(tapply(C_N_C_P[,4],fac,function(x){
	   fdr(as.numeric(x),fdr_flag="fdr")
	   }))
	   #p1<-fdr(as.numeric(C_N_C_P[,4]),fdr_flag="fdr")
       pos1<-which(p1<=thr)
      if(length(pos1)){
          C_N_C_P[-pos1,3]<-0
       }else{
          C_N_C_P[,3]<-0
       }
	   

       p2<-unlist(tapply(C_N_C_P[,6],fac,function(x){
	   fdr(as.numeric(x),fdr_flag="fdr")
	   }))
	   #p2<-fdr(as.numeric(C_N_C_P[,6]),fdr_flag="fdr")
       pos2<-which(p2<=thr)
       if(length(pos2)){
          C_N_C_P[-pos2,5]<-0
       }else{
          C_N_C_P[,5]<-0
       }
	   ###将meth，cnv添上自己的名称
	   index1<-which("meth"==C_N_C_P[,2])
	   index2<-which("cnv"==C_N_C_P[,2])
	   C_N_C_P[index1,2]<-paste0(C_N_C_P[index1,1],"meth")
	   C_N_C_P[index2,2]<-paste0(C_N_C_P[index2,1],"cnv")
	   ###对两种情况分别构建网络
	   pos<-union(pos1,pos2)
       C_N_C_P1<-as.data.frame(C_N_C_P[pos,c(1,2,3)],stringsAsFactors=F)
	   C_N_C_P1[,3]<-abs(as.numeric(C_N_C_P1[,3]))
	   C_N_C_P2<-as.data.frame(C_N_C_P[pos,c(1,2,5)],stringsAsFactors=F)
	   C_N_C_P2[,3]<-abs(as.numeric(C_N_C_P2[,3]))
	   g1<-graph.data.frame(C_N_C_P1,directed=F)
	   g2<-graph.data.frame(C_N_C_P2,directed=F)
	   ###设置边的权重
	   g1<-set.edge.attribute(g1, name="weight", index=E(g1), value=C_N_C_P1[,3])
	   g2<-set.edge.attribute(g2, name="weight", index=E(g2), value=C_N_C_P2[,3])
	   degree1<-graph.strength(g1,vids = V(g1))
	   names(degree1)<-names(V(g1))
	   degree2<-graph.strength(g2,vids = V(g2))
	   names(degree2)<-names(V(g2))
	   ####  DiffK (加权的度)
	   DiffK_score<-abs(degree1/max(degree1)-degree2/max(degree2))
	   DiffK_score
	  
}

get_downstream_effectors<-function(k){##i表示真实gene下标
	#k=7
	###获取随机回归的结果
	file_pos<-((k-1)*100+1):(k*100);
	zy_file<-paste0(((file_pos-1)*10+1),".RData")
	diffK_score_random_list<-NULL
	for(my_file in zy_file){
		file_name_temp<-paste0(path2,"/random_result/",my_file)
		load(file_name_temp)
		temp_score<-lapply(1:10,function(k){
		 get_diffK_score(C_N_C_P=C_N_C_P_temp[[k]],fdr_flag="fdr",thr=0.05)
		})
		diffK_score_random_list<-c(diffK_score_random_list,temp_score)
	}
	

	###计算差异的显著性
	diffK_score_list_temp<-lnc_diffK_score_list[[k]]
	p_value_all<-sapply(1:length(diffK_score_list_temp),function(i){
		gene_name<-names(diffK_score_list_temp)[i]
		a<-sapply(1:1000,function(j){
			pos<-which(gene_name==names(diffK_score_random_list[[j]]))
			if(length(pos)==0)
			return(F)
			else{
			flag<-diffK_score_random_list[[j]][pos]>diffK_score_list_temp[i]
			}
		})
		p_value<-length(which(a==T))/1000
		p_value
	})
	names(p_value_all)<-names(diffK_score_list_temp)
	downstream_effectors<-p_value_all
	save.file<-paste0(path2,"/temp/",k,"aa.RData")
	save(downstream_effectors,file=save.file)
	downstream_effectors
}

lnc_downstream_effectors<-sapply(1:length(lnc_diffK_score_list),get_downstream_effectors)
names(lnc_downstream_effectors)<-names(lnc_diffK_score_list)
setwd(path2)
save(lnc_downstream_effectors,file="lnc_downstream_effectors.RData")


lnc_downstream_effectors_new<-sapply(lnc_downstream_effectors,function(x){
	# x=lnc_downstream_effectors[[1]]
	index1<-grep("cnv",names(x))
	x<-x[-index1]

	pos<-which(x<=0.05)
	names(pos)
})

save(lnc_downstream_effectors_new,file="lnc_downstream_effectors_new.RData")


