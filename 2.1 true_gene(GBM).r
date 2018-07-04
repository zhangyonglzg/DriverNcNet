main_path<-"/pub5/xiaoyun/Jobs/J19/"
cancer<-"GBM"
load(paste0(main_path,cancer,"/result2/gene_case_normal_exp.RData"));
gene_case<-gene_M
load(paste0(main_path,cancer,"/result2/all_data_coding.RData"));
path0<-paste0(main_path,cancer,"/")
path1<-paste0(main_path,cancer,"/gene_result")

find_var_thresh<-function(gene_M,per_thresh){
   ###使用标准差和相对阈值
   gene_cv<-apply(gene_M,1,function(x){sd(2^x)})
   gene_cv_new<-rev(sort(gene_cv))[1:floor(length(gene_cv)*per_thresh)]
   return(names(gene_cv_new))
}
library(siggenes)
dif_sam2_fold_up<-function(expr1,expr2,fold_thr,thr){

        l1<-dim(expr1)[2];##the number of normal sample
        l2<-dim(expr2)[2];##the nunber of case sample
        coln<-c(rep(0,l1),rep(1,l2));##construct the class label for sample
        expr<-cbind(expr1,expr2);##combind the case and normal expression

        result<-list(); #construct the frame of result
        length(result)<-3;
        names(result)<-c("dif_genes","up_genes","down_genes");

        sam.out <- sam(expr, coln, method=d.stat,rand=123)##perform the SAM
        fc<-sam.out@fold; ##extract the fold changes
        delta <- findDelta(sam.out, fdr=thr);

        if(length(delta)){
              if(length(dim(delta))){
                  genes <- list.siggenes(sam.out, delta[1,1]);

                  fc<-fc[genes];
                  result[[1]]<-names(which(fc>fold_thr|fc<(1/fold_thr)));

                  gene_up<-names(which(fc>fold_thr));
                  result[[2]]<-gene_up;

                  gene_down<-names(which(fc<(1/fold_thr)));
                  result[[3]]<-gene_down;
                  return(result);
               }else{
                   genes <- list.siggenes(sam.out, delta[1]);
                   fc<-fc[genes];
                   result[[1]]<-names(which(fc>fold_thr|fc<(1/fold_thr)));
                   gene_up<-names(which(fc>fold_thr));
                   result[[2]]<-gene_up;
                   gene_down<-names(which(fc<(1/fold_thr)));
                   result[[3]]<-gene_down;
                   return(result);
               }

        }  else{
             return(NULL)
      }
}

find_var_gene <- find_var_thresh(gene_M=gene_case,per_thresh=0.8);
diff_gene <- dif_sam2_fold_up(expr1=gene_normal,expr2=gene_case,fold_thr=2,thr=0.05);


do_normalization<-function(gene_matrix){
	t(apply(gene_matrix,1,function(x){
		(x-mean(x))/sd(x)
	}))
}
gene_case<-do_normalization(gene_case)
gene_normal<-do_normalization(gene_normal)
gene_M<-do_normalization(gene_M)
mir_M<-do_normalization(mir_M)
tf_M<-do_normalization(tf_M)
copy_M<-do_normalization(copy_M)
meth_M<-do_normalization(meth_M)



##1.构建基因—因子谱 
###################################################################
gene_factor_M<-function(gene_name,gene_M,copy_M,meth_M,mir_M,tf_M,gene_mir_list,gene_tf_list){
##提取一个基因层面数据，组成gene-factor谱
#分析gene
	temp_M<-c();
	temp_name<-c();
	temp_M<-data.frame(gene_M[gene_name,]);##表达谱
	temp_name<-c(temp_name,"exp");

	if(length(which(rownames(copy_M)==gene_name))){##拷贝数谱
		temp_M<-cbind(temp_M,as.numeric(copy_M[gene_name,]));
		temp_name<-c(temp_name,"cnv");
	}

	if(length(which(rownames(meth_M)==gene_name))){##甲基化谱

		temp_M<-cbind(temp_M,as.numeric(meth_M[gene_name,]));
		temp_name<-c(temp_name,"meth");
	}

	if(length(which(names(gene_mir_list)==gene_name))){##调控的miRNA表达谱

		temp_mir<-gene_mir_list[[gene_name]];
		temp_mir<-intersect(rownames(mir_M),temp_mir);

		if(length(temp_mir)>1){
			temp_M<-cbind(temp_M,t(mir_M[temp_mir,]));
			temp_name<-c(temp_name,temp_mir);
		}else if(length(temp_mir)==1){
			temp_M<-cbind(temp_M,mir_M[temp_mir,]);
			temp_name<-c(temp_name,temp_mir);
		}
	}

	if(length(which(names(gene_tf_list)==gene_name))){##调控的转录因子谱

		temp_tf<-as.character(gene_tf_list[[gene_name]]);
		temp_tf<-intersect(rownames(tf_M),temp_tf);

		if(length(temp_tf)>1){
			temp_M<-cbind(temp_M,t(tf_M[temp_tf,]));
			temp_name<-c(temp_name,temp_tf);
		}else if(length(temp_tf)==1){

			temp_M<-cbind(temp_M,tf_M[temp_tf,]);
			temp_name<-c(temp_name,temp_tf);
		}

	}

	colnames(temp_M)<-temp_name;
	return(temp_M)
}

gene_name <- rownames(gene_M);
gene_factor_list  <- lapply(gene_name,gene_factor_M,gene_M,copy_M,meth_M,mir_M,tf_M,gene_mir_list,gene_tf_list);
names(gene_factor_list) <- gene_name;
######################################################################


##2.删除影响因子小于2个的基因
#####################################################################
factor_num<-function(M){
## 计算影响因子数目
      return((dim(M)[2]-1)) 
}
delete_lessF_genes<-function(gene_factor_list,num_thresh){
    factor_nums<-unlist(lapply(gene_factor_list,factor_num))
    pos<-which(factor_nums>=num_thresh)##
    gene_factor_list<-gene_factor_list[pos]
    return(gene_factor_list)
}
gene_factor_list1 <- delete_lessF_genes(gene_factor_list,num_thresh=2);
###################################################################


##3.使用case_control差异表达并在case样本中广泛变化的基因（提前）过滤
##################################################################

gene_factor_list_deal <- function(find_var_gene, diff_gene, gene_factor_list){
	pos1 <- match(find_var_gene, names(gene_factor_list), nomatch=0);
	result <- gene_factor_list[pos1];
	pos2 <- match(diff_gene$dif_genes, names(result), nomatch=0);
	result <- result[pos2];
	return(result);
}

gene_factor_list2 <- gene_factor_list_deal(find_var_gene=find_var_gene, diff_gene=diff_gene, gene_factor_list=gene_factor_list1);
setwd(path1)
save(gene_factor_list2,file="gene_factor_list2.RData")
##############################################################


##4.识别真实情况下的候选调控关系
###############################################
source(paste0(path0,"/gene_parallel_true.r"))
get_plrs_result(1:nrow(cnv_mut_M))
################################################


##5 获取所有关系对
####################################################
setwd(paste0(path1,"/true_result"))
my_file<-dir()
my_file_new<-sapply(my_file,function(x){unlist(strsplit(x,"\\.RData"))})
my_file_new_order<-sapply(sort(as.numeric(my_file_new)),function(x){paste0(x,".RData")})

C_N_C_P_gene_true<-NULL
for(my_file in my_file_new_order){
load(my_file)
C_N_C_P_gene_true<-c(C_N_C_P_gene_true,list(C_N_C_P_temp))
}
names(C_N_C_P_gene_true)<-rownames(cnv_mut_M)
setwd(path1)
save(C_N_C_P_gene_true,file="C_N_C_P_gene_true.RData")
######################################################


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

gene_diffK_score_list<-lapply(1:length(C_N_C_P_gene_true),function(x){
	temp<-C_N_C_P_gene_true[[x]]
	temp_score<- get_diffK_score(C_N_C_P=temp,fdr_flag="fdr",thr=0.05)
})

names(gene_diffK_score_list)<-rownames(cnv_mut_M)
setwd(path1)
save(gene_diffK_score_list,file="gene_diffK_score_list.RData")
###################################################################
