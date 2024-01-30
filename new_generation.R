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
#补充：exprs()是Bioconductor中Biobase包中的一个函数，用于获取对象为ExpressionSet属性的表达矩阵。
dat=exprs(a)


#第八步：使用pData函数来获取临床信息。
pd=pData(a)
#说明：后续需要根据临床信息表内的title来制定分组，该表结构也比较清晰易懂，可以查看。
#扩展：但其实在GEO主页看sunnmary更加易懂。

