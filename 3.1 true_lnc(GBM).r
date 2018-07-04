main_path<-"/pub5/xiaoyun/Jobs/J19/"
cancer<-"GBM"
load(paste0(main_path,cancer,"/result2/gene_case_normal_exp.RData"));
gene_case<-gene_M
load(paste0(main_path,cancer,"/result2/all_data_lnc.RData"));
path0<-paste0(main_path,cancer,"/")
path1<-paste0(main_path,cancer,"/gene_result")
path2<-paste0(main_path,cancer,"/lnc_result")

load(paste0(path1,"/gene_factor_list2.RData"))
##############################################################


####4.识别真实情况下的候选调控关系
#########################################
source(paste0(path0,"/lnc_parallel_true.r"))
get_plrs_result(1:nrow(cnv_M_lnc))
#########################################

####5.获取所有关系对
############################################
setwd(paste0(path2,"/true_result"))
my_file<-dir()
my_file_new<-sapply(my_file,function(x){unlist(strsplit(x,"\\.RData"))})
my_file_new_order<-sapply(sort(as.numeric(my_file_new)),function(x){paste0(x,".RData")})

C_N_C_P_lnc_true<-NULL
for(my_file in my_file_new_order){
load(my_file)
C_N_C_P_lnc_true<-c(C_N_C_P_lnc_true,list(C_N_C_P_temp))
}
names(C_N_C_P_lnc_true)<-rownames(cnv_M_lnc)
setwd(path2)
save(C_N_C_P_lnc_true,file="C_N_C_P_lnc_true.RData")
###################################################################

####使用差异网络（点）的方法求差异关系
###6 删除两种状态下都不显著的关系
#########################################################
fdr<-function(p,fdr_flag){ ##calculate the fdr
   p<-p.adjust(p,method=fdr_flag)
   return(p)
}
library(igraph)
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
	   ####方法  DiffK 
	   DiffK_score<-abs(degree1/max(degree1)-degree2/max(degree2))
	   DiffK_score
	  
}

lnc_diffK_score_list<-lapply(1:length(C_N_C_P_lnc_true),function(x){
	temp<-C_N_C_P_lnc_true[[x]]
	temp_score<- get_diffK_score(C_N_C_P=temp,fdr_flag="fdr",thr=0.05)
})

names(lnc_diffK_score_list)<-rownames(cnv_M_lnc)
setwd(path2)
save(lnc_diffK_score_list,file="lnc_diffK_score_list.RData")
###################################################################
