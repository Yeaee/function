#function
#######准备工作：
#######在R打开后新建文件R script后，要保存并关闭，然后再从文件夹内打开Rscript，才能更改目录。






###################
#function1：GSEnumber_to_visualized_raw_exprSet
#输入：GSE号，纯GSE+数字
#输出：return该GSE的原始表达集raw_exprSet。
#输出：本地save下txt_gz_file_eSet.Rdata和GSE_raw_exprSet文件。
GSEnumber_to_visualize_raw_exprSet <- function(GSE_number){
  library(GEOquery)
  GSE_txt_gz_flie <- getGEO(GSE_number,destdir = '.',
                            AnnotGPL = F,
                            getGPL = F)
  save(GSE_txt_gz_flie,file=paste0(GSE_number,"_eSet.Rdata"))
  group_list = c(rep('control',3),rep('case',3)) 
  
  the_first_object <- GSE_txt_gz_flie[[1]]
  new_name = array(paste0(GSE_number,'_raw_exprSet'))
  raw_exprSet <- exprs(the_first_object)
  save(raw_exprSet,group_list,
       file = paste0(GSE_number,'_raw_exprSet.Rdata'))
  return(raw_exprSet)}
#使用案例：
GSE194261_raw_exprSet<-GSEnumber_to_visualize_raw_exprSet('GSE194261')




###################
#function2:GSEnumber_to_visualized_excel
#输入：GSE号，纯GSE+数字
#输出：默认当前文件夹下该GSE的可视化excel表格
GSEnumber_to_visualized_excel <- function(studyID,destdir = '.'){
  
  library(GEOquery)
  eSet <- getGEO(studyID,destdir = destdir,getGPL = F)
  
  exprSet = exprs(eSet[[1]])
  pData = pData(eSet[[1]])
  
  write.csv(exprSet,paste0(studyID,'_exprSet.csv'))
  write.csv(pData,paste0(studyID,'_metadata.csv'))
}
#使用案例（在该R脚本的默认文件夹下可以看到csv文件）：
GSEnumber_to_visualized_excel('GSE194261')
####注释：元数据（Metadata）文件，是指关于该文件或对象的一切信息都是元数据。



###################
#function3:make_GSEs_ids
#人工：搜索GSE号在GEO首页找到platform一栏，找到其对应platform为GPL****或其他.
#人工: 在网站http://www.bio-info-trainee.com/1399.html寻找该GPL对应的bioc_package
#输入：该GSE_number芯片平台对应的bioc_package的.db与SYMBOL形式。
#输出：制作好的ids
make_GSEs_ids <- function(symbol){
  ####if (!require("BiocManager", quietly = TRUE))
  ####  install.packages("BiocManager")
  #library(BiocManager)
  #BiocManager::install(array(paste0(bioc_package,".db")))
  #这里报错,怎么搞都导不进去，算了，这步就手动挡导入库吧，放在最前面就行
  #library(db)
  ids = toTable(symbol)
  return(ids)
}
#使用案例：
##搜索GSE号在GEO首页找到platform一栏，找到其对应platform为GPL****或其他.
##在网站http://www.bio-info-trainee.com/1399.html
##寻找该GPL对应的bioc_package为hugene10sttranscriptcluster；
##手动输入.db格式和SYMBOL格式。
#BiocManager::install("hugene10sttranscriptcluster.db")
library("hugene10sttranscriptcluster.db")
GSE_42872_ids <- make_GSEs_ids(hugene10sttranscriptclusterSYMBOL)



################
#分支function3.1：注释平台文件内找不到对应的R包的处理（https://www.jianshu.com/p/8ffa5b602b9c）
################
##失败案例
#这里试一下这个网站提供的idmap2
#尝试一：
#BiocManager::install('idmap2') #这个装都装不了，cran也试过不行。
#尝试二：
#BiocManager::install('AnnoProbe') #这个能装，但不知道怎么用
#插件：
#先更新一下biocmanager吧。
#BiocManager::install(version = "devel")
#尝试三：
#BiocManager::install('BiocInstaller')
#BiocInstaller::biocLite('idmap2')
#install.packages('BiocInstaller')
#屌用没有。试一下方法二。
##################
#成功案例：
library(GEOquery)
gpl <- getGEO('GPL17077', destdir=".")
colnames(Table(gpl)) 
head(Table(gpl)[,c(1,6,7)])  ## you need to check this , which column do you need 
write.csv(Table(gpl)[,c(1,6,7)],"GPL17077.csv")
ids=read.csv('GPL17077.csv')
GSE_194261_ids = ids
#做完ids后查看一下，只要有这两列，一列是id，一列是gene_symbol就行。
#可以再去raw_exprset里面查看一下，随便找一个id复制，去ids里面ctrlF一下，看看有没有。
#成功，还是这个方法好用，下次都可以直接用这个方法了。

################
#分支function3.2：raw_exprset做id转换，拆分步骤。
rownames(GSE194261_raw_exprSet) <- ids$GENE_SYMBOL


###################
#function4：raw_exprSet_to_final_exprSet
#输入：function1得到的raw_exprSet以及function3得到的ids。
#输出：final_exprSet
#函数简介：这个函数滤了两次,根据symbol分了个类,最后把数字换成了基因名。
raw_exprSet_to_final_exprSet <- function(exprSet,ids){
  table(rownames(exprSet)%in%ids$probe_id)
  exprSet = exprSet[rownames(exprSet)%in%ids$probe_id,]
  ids = ids[match(rownames(exprSet),ids$probe_id),]
  tmp = by(exprSet,
           ids$symbol,
           function(x) rownames(x)[which.max(rownames(x))]) 
  probes = as.character(tmp)
  new_exprSet = exprSet[rownames(exprSet)%in% probes,]
  rownames(new_exprSet)=ids[match(rownames(new_exprSet),ids$probe_id),2]
  return(new_exprSet)
}
#使用范例：
GSE194261_final_exprSet <- raw_exprSet_to_final_exprSet(GSE194261_raw_exprSet,GSE_194261_ids)


###################
#function5:make_group_list_of_GSE_exprSet
#人工：去该GSE主页看下sample一栏中的样本情况。
#输入：控制组个数，对照组个数。
#输出：根据样本情况写出group_list。
make_group_list_of_GSE_exprSet <- function(number_of_control,number_of_case){
  group_list = c(rep('control',number_of_control),rep('case',number_of_case))
  return(group_list)
}
#使用案例：
GSE_42872_group_list <- make_group_list_of_GSE_exprSet(3,3)


###################
#function6:check_the_analysis_value_of_final_exprSets
#输入：GSE*****_final_exprSet,GSE_*****_group_list
#输出：箱线图
#解释：就是画箱线图（检测图里的线是不是差不多在一条线上，是则代表矩阵有分析价值）
check_the_analysis_value_of_final_exprSets <- function(new_exprSet,group_list){
  library(reshape2)
  library(ggplot2)
  exprSet_L=melt(new_exprSet)
  colnames(exprSet_L)=c('probe','sample','value')
  exprSet_L$group=rep(group_list,each=nrow(new_exprSet))
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank(),panel.background = element_rect(fill = 'lightgrey'))
  print(p)
}
#案例：
check_the_analysis_value_of_final_exprSets(GSE42872_final_exprSet,GSE_42872_group_list)



###################
#function7:check_if_difference_exsits_between_control_and_case
#输入：GSE*****_final_exprSet,GSE_*****_group_list,sample_size(为control加case数)
#输出：hclust_plot，一种树状图，显示数据的分层结构和聚类关系。
#解释：看图里的case和control分的开不开，分的开就是好，说明存在差异。
check_if_difference_exsits_between_control_and_case <- function(new_exprSet,group_list,sample_size){
  colnames(new_exprSet)=paste(group_list,1:sample_size,sep='')
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                  cex = 0.7, col = "blue")
  hc=hclust(dist(t(new_exprSet)))
  par(mar=c(5,5,5,10)) 
  plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)
}
#案例：
check_if_difference_exsits_between_control_and_case(GSE42872_final_exprSet,GSE_42872_group_list,6)



###################
#function8:make_sure_the_final_exprSet_can_be_used
#输入：GSE*****_final_exprSet,GSE_*****_group_list。
#输出：PCA_plot，用来绘制主成分分析（PCA）的结果，也就是一种散点图。
#解释：PCA图：也是看分地开不开，三个图已经够了，可以确定矩阵没什么问题。
make_sure_the_final_exprSet_can_be_used <- function(new_exprSet,group_list){
  library(ggfortify)
  df=as.data.frame(t(new_exprSet))
  df$group=group_list 
  autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')
}
#案例：
make_sure_the_final_exprSet_can_be_used(GSE42872_final_exprSet,GSE_42872_group_list)

###################
#function9:make_design_matrix
#输入：GSE*****_final_exprSet,GSE_*****_group_list。
#输出：design矩阵。
make_design_matrix <- function(new_exprSet,group_list){
  design <- model.matrix(~0+factor(group_list))
  colnames(design) = levels(factor(group_list))
  rownames(design) = colnames(new_exprSet)
  return(design)
}
#案例：
GSE_42872_design_matrix <- make_design_matrix(GSE42872_final_exprSet,GSE_42872_group_list)



###################
#function10:make_contrast_matrix
#输入：GSE_*****_group_list,GSE_*****_design_matrix
#输出：contrast矩阵。
#解释：把实验组作为1，把控制组作为-1.
make_contrast_matrix <- function(group_list,design){
  library(limma)
  contrast.matrix<-makeContrasts("case-control",
                                 levels = design)
  return(contrast.matrix)
}
#案例：
GSE_42872_contrast_matrix = make_contrast_matrix(GSE_42872_group_list,GSE_42872_design_matrix)



###################
#function11:make_nrDEGs_by_differential_analysis
#输入：GSE*****_final_exprSet,GSE_*****_design_matrix,GSE_*****_contract_matrix.
#输出：非冗余差异表达基因，就是有显著差异的基因片段。nrDEG(non-redundant differentially expressed genes)
#解释：把三个矩阵喂进去，得到差异分析结果,注意这个版本的nrDEG已经按logFC的绝对值排好序了，其中logFC有正有负。(DEG:Differential Expressed Genes)
make_nrDEGs_by_differential_analysis<- function(new_exprSet,design,contrast_matrix){
  library(limma)
  fit <- lmFit(new_exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast_matrix) ##这一步很重要，大家可以自行看看效果
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  nrDEG <- nrDEG[order(abs(nrDEG[,1]), decreasing=TRUE),]
  return(nrDEG)
}
GSE_42872_nrDEG = make_nrDEGs_by_differential_analysis(GSE42872_final_exprSet,GSE_42872_design_matrix,GSE_42872_contrast_matrix)



###################
#function12:visualize_nrDEG_by_volcano_plot
#输入：nrDEG
#输出：可视化的火山图，把logFC>2,p<0.05的显著差异基因用重点色标注。
visualize_nrDEG_by_volcano_plot<-function(DEG){
  logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)) )
  DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                                ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                      '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                      '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
  )
  this_tile
  head(DEG)
  g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)
  print(g)
}
visualize_nrDEG_by_volcano_plot(GSE_42872_nrDEG)


###################
#function12.1:visualize_nrDEG_by_volcano_plot_v2
#解释：可以在图里显示几个显著差异基因名称。
#解释：ifelse语句：#test：条件语句 #yes：条件为True时执行 #no：条件为False时执行 ifelse(test, yes, no)
visualize_nrDEG_by_volcano_plot_v2 <- function(deg){
  library(ggpubr)
  nrDEG=deg
  head(nrDEG)
  attach(nrDEG)
  plot(logFC,-log10(P.Value))
  library(ggpubr)
  df=nrDEG
  df$v= -log10(P.Value) #df新增加一列'v',值为-log10(P.Value)
  ggscatter(df, x = "logFC", y = "v",size=0.5)
  
  df$g=ifelse(df$P.Value>0.01,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
              ifelse( df$logFC >2,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
                      ifelse( df$logFC < -2,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
  )
  table(df$g)
  df$name=rownames(df)
  head(df)
  ggscatter(df, x = "logFC", y = "v",size=0.5,color = 'g')
  ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
            label = "name", repel = T,
            #label.select = rownames(df)[df$g != 'stable'] ,
            label.select = head(rownames(deg)), #挑选一些基因在图中显示出来
            palette = c("#00AFBB", "#E7B800", "#FC4E07") )
}
visualize_nrDEG_by_volcano_plot_v2(GSE_42872_nrDEG)




###################
#function13:visualize_the_top25_Differential_Expressed_Genes_by_pheatmap
#输入：GSE*****_final_exprSet,nrDEG
#输出：top25差异基因热图
visualize_top25_Differential_Expressed_Genes_by_pheatmap <- function(new_exprSet,nrDEG){
  library(pheatmap)
  # 选取头部的25个差异最显著的基因，做一个小表达矩阵
  choose_gene=head(rownames(nrDEG),25)
  choose_matrix=new_exprSet[choose_gene,]
  choose_matrix=t(scale(t(choose_matrix)))
  pheatmap(choose_matrix)
  #颜色深浅就是表达量的高低
}
visualize_top25_Differential_Expressed_Genes_by_pheatmap(GSE42872_final_exprSet,GSE_42872_nrDEG)



###################
#function14:visualize_all_gene_expression_level_by_pheatmap
#输入：GSE*****_final_exprSet,nrDEG
#输出：全部基因片段表达量的热图，颜色深浅即为表达量的高低。
#备注：这个没什么可读性，可以不用，把rowname去掉就可以用了(已去除）。
visualize_all_gene_expression_level_by_pheatmap <- function(new_exprSet,nrDEG){
  library(pheatmap)
  choose_gene=rownames(nrDEG)
  choose_matrix=new_exprSet[choose_gene,]
  choose_matrix=t(scale(t(choose_matrix)))
  pheatmap(choose_matrix,show_rownames = FALSE)
  #颜色深浅就是表达量的高低
}
visualize_all_gene_expression_level_by_pheatmap(GSE42872_final_exprSet,GSE_42872_nrDEG)



###################
#function16:pick_significant_DEGs_by_fixed_number
#输入：目标选择的基因数量，nrDEG矩阵。
#输出：前多少个显著基因的矩阵,格式为三列。
#解释：直接把原来的矩阵排序后，取前面固定个就行了。
pick_significant_DEGs_by_fixed_number<- function(fixed_number,nrDEG){
  a = head(rownames(nrDEG),fixed_number)
  a.df <- bitr(gene, fromType = "SYMBOL",
               toType = c("ENSEMBL", "ENTREZID"),
               OrgDb = org.Hs.eg.db)
  return(a.df)
}
GSE_42872_significant_DEGs.df<- pick_significant_DEGs_by_fixed_number(25,GSE_42872_nrDEG)



###################
#function17:pick_significant_DEGs_by_logFC
#输入：想要设定的阈值，即为logFC_t,以及原先的nrDEG矩阵。
#输出：阈值筛选后的significant_DEG矩阵，格式为三列。
pick_significant_DEGs_by_logFC <- function(logFC_t,nrDEG){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  gene = rownames(nrDEG)
  data(geneList, package="DOSE")
  
  gene <- names(geneList)[abs(geneList) > logFC_t]
  gene.df <- bitr(gene, fromType = "ENTREZID",
                  toType = c("ENSEMBL", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
  return(gene.df)
}
#范例：
GSE_42872_significant_DEGs.df <-pick_significant_DEGs_by_logFC(1.5,GSE_42872_nrDEG)



###################
#function18：make_KEGG_pathway_enrichment_analysis_for_picked_significant_DEGs
#直接使用现成的函数：
kk <- enrichKEGG( gene = GSE_42872_significant_DEGs.df$ENTREZID,
                  organism ='hsa',
                  pvalueCutoff = 0.05
)
head(kk)[,1:6]







###################
#目前还没用到
pick_significant_DEGs <- function(logFC_t,deg){
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  deg$g=ifelse(deg$P.Value>0.05,'stable',
               ifelse( deg$logFC > logFC_t,'UP',
                       ifelse( deg$logFC < -logFC_t,'DOWN','stable'))
  )
  DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
  deg$symbol=rownames(deg)
  df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  DEG=deg
  DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
  save(DEG,file = 'anno_DEG.Rdata')
  gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
  gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
  gene_diff=c(gene_up,gene_down)
  return(gene_diff)
}
a <- pick_significant_DEGs(1.5,GSE_42872_nrDEG)
