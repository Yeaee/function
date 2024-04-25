#第一步：用getGEO下载数据。
#说明：先安装BiocManager，再用BiocManager安装GEOquery。
library(GEOquery)
gset = getGEO('GSE194261',destdir = '.',
              getGPL = FALSE,
              AnnotGPL = FALSE
              )
#提示：用str(getGEO)来查看该函数内所需的对象。


#第二步：查看该获取数据的格式，养成习惯。
class(gset)
#说明：得出数据格式为list，list为对象的集合,是二维甚至多维的数据类型。
#拓展：向量是一维的只包含数组，而列表内可包括向量、矩阵、数据框、甚至是其他列表。
#特别提醒：向量中的元素可以通过索引进行访问，而列表中的元素需要使用双方括号[[]]进行访问。


#第三步：把获取的数据保存到本地,必须用引号设立文件名。
#注意：gset文件惯例命名为GSE号，文件后缀一定要写，后缀文件属性为.Rdata，为list文件。
#额外尝试：好像不加后缀也能被load函数提取出来，但是加了后缀更方便理解文件属性，
#不加后缀虽然也没什么不好，但加上Rdata对我来说更加地海阔天空嘛。
save(gset,file = 'GSE194261.Rdata')


#第四步：清空列表后尝试读取刚刚存下的文件，文件名目前基本都要加引号。
rm(list = ls())
#说明：ls() 函数返回当前环境中的所有对象的名称。list为要删除对象的名称。
#说明：rm(objectname)为删除特定变量名称的对象。
load('GSE194261.Rdata')


#第五步：提取该列表中的第一个元素。
#注意：在R语言中，使用单方括号[]访问多维数组时，返回的是一个子数组，而不是一个元素。
#注意：因此，为了访问列表中的元素，需要使用双方括号[[]]来指定要访问的元素的位置。
a=gset[[1]]


#第六步：观察a的属性
class(a)
#类别为"ExpressionSet"，在R语言中，ExpressionSet是一个用于高通量测序数据的类。


#第七步：使用exprs()函数获取表达矩阵。
#补充：exprs()是Bioconductor中Biobase包中的一个函数，该函数需要对属性为ExpressionSet的对象使用。
dat=exprs(a)


#第八步：使用pData函数来获取临床信息。
pd=pData(a)
#说明：后续需要根据临床信息表内的title来制定分组，该表结构也比较清晰易懂，可以查看。
#扩展：但其实在GEO主页看sunnmary更加易懂。


#第九步：检查一下获取的pd临床信息是否和GEO官方库里对这个基因集的介绍一样。
View(pd)
#分析：每个样本的基本信息应该是没有差别的，偷个懒就不去官网对照了。


#第十步：获取ids转换文件所需对应的包。
#说明：该步骤一般是找特定平台对应的官方包，如果没有官方包再用GPL。
#前人经验：getGPL获得平台的注释信息，但下载速度会慢很多，且注释文件格式大多不如bioconductor包好用。
#特殊标记：针对该数据库，这一步存在较大改进空间，因为找不到官方注释库，用的是泛id库。
gpl <- getGEO('GPL17077', destdir=".")
#细节：这里要开梯子，不然五十多个mb有点难下。



#第十一步：查找包内id转换需要的列。
colnames(Table(gpl)) 
#分岔路口：这步也可也直接点UI右侧的gpl，一层一层从dataTable，columns里查看。


#第十二步：挑选目的列，存贮为ids文件。
#需求：做完ids后查看一下，只要有这两列，一列是id，一列是gene_symbol就行。
#源代码内自带注释：you need to check this , which column do you need
head(Table(gpl)[,c(1,6,7,8)]) 
write.csv(Table(gpl)[,c(1,6,7,8)],"GPL17077.csv")
ids=read.csv('GPL17077.csv')
GSE_194261_ids = ids
#验证：可以再去raw_exprset里面查看一下，随便找一个id复制，去ids里面ctrlF一下，看看有没有。
#针对性结果：成功，还是这个方法好用，下次都可以直接用这个方法了。



#第十三步：清除ids中无法找到对应id的行。
dat <- dat[rownames(dat)%in% ids$ID,]
#bug：这一步有个bug，怎么所有的id都能找到啊，太不正常了。而且条数都是50739。
#古夫，观察并思考：观察发现数据库里的x.**是按规律排得，估计这里就是代指找不到的id了。
#抽样搜索一下，抽了一个能搜出来是基因。

#第十四步：替换ids。
#内核：这一步其实就是把表里的行名从平台id换成gene_symbol。
rownames(dat)=ids[match(rownames(dat),ids$ID),4]
save(dat,file = 'GSE194261_idtransed_exprSet')





#阶段总结：第一部分准备工作完成。

###########################################
##############第二部分#####################
rm(list = ls())


#第一步：读取记忆。
load('GSE194261_idtransed_exprSet')
#右边点开检查一下，是不是已经转换过ID了。


#第二步：确定分组，制作group_list。
#去GSE首页看下，应该怎么分组。就拿2-5分第成control组，6-9分为case组吧。
#要点：当前理解下，第二阶段的核心就是列名的替换，把GSM带数字的列名替换成control与case，和ID转换本质上是一样的。
group_list <- c(rep('control', 4), rep('case', 4),rep('other',12))


#第三步：检验一下前八个样本内的表达分布
boxplot(dat[,1:8])


#第四步：检验内参表达量，通过管家基因表达水平判断数据是否有效。
dat["GAPDH",]#检查内参表达量
dat["ACTB",]#检查内参表达量
#基准：一般gapdh在8-9是正常的，参考https://zhuanlan.zhihu.com/p/344426350。


#第五步：group_list制作完之后再SL大法一下，熟练掌握一下R语言里的断续处理。
save(dat,group_list,
     file = "./stage2_middle.Rdata")


#第六步：把列名进行替换为group_list，完成分组。
rm(list = ls())
load('./stage2_middle.Rdata')
colnames(dat) <- paste(group_list,1:20)


#第七步：查看替换完的聚类是否按照group_list分组。
#结果分析：大部分聚在一起，除了case5。
plot(hclust(dist(t(dat))))


#第七步的平行替换：图不一样，目的一样，看聚类。
# BiocManager::install('ggfortify')
library(ggfortify)
df=as.data.frame(t(dat))
df$group=group_list 
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')


#第八步：dim函数查看数据维度，这步目前来看有点无关紧要。
dim(dat)


#第九步：SL大法。
save(dat,group_list,file = "./stage2_finished.Rdata")

#阶段总结：第二部分工作完成，准备工作结束。

###########################################
##############第三部分#####################
rm(list = ls())

#第一步：载入。
load('./stage2_finished.Rdata')


#第二步：制作设计矩阵。
#新文化翻译代码：就是把一行字符串的group_list换成矩阵形式。
#代码注释：model.matrix将一个线性模型的公式转换为一个矩阵，以便在统计分析中使用。

design_matrix <- model.matrix(~0+factor(group_list))
colnames(design) = levels(factor(group_list))
rownames(design) = colnames(new_exprSet)

