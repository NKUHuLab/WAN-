options(java.parameters = "-Xmx80g")
library(lattice)
library(rlang)
library(rJava)
library(xlsx)
library(readxl)
#library(hydrogeo)  
library(randomForest)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr) #冲突pdp注意！
library(randomForestExplainer)
library(pdp)
library(tcltk)
library(patchwork)
library(caret)
library(ggrepel)
library(data.table)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer) 
library(remotes)
setwd('F:/unthetao/chl/11low27')
label1<-excel_sheets('F:/unthetao/chl/11low27/1.xlsx')
label<-c(label1)
colindex<-c('chl','Fe','NO3','NPPV','O2','Ph','phyc','PO4','Si','SPCO2','so',
            'uo','vo','DIC','NH4','CO3','DOC','talk','phyn','dfftdin','tob',
            'mlotst','zos','dpco2','protect','mp','human','Index')

#####模型构建#####

rf.list<-list()
for (i in 1:2){
  data<-read.xlsx('1.xlsx',i,header=TRUE)
  colnames(data)<-colindex
  # train<-read.xlsx('Ocean-data-R.xlsx',i)$X0
  # data_train<-data[train,]
  # data_test<-data[-train,]
  seed<-rep(1,2)
  set.seed(seed[i])
  rf.list[[label1[i]]]<-local({
    data=data
    randomForest(Index~.,data=data,importance=TRUE,proximity=T,ntree=700,mtry=2)    #model
  })
  rf<-rf.list[[label[i]]]
  pre_train<-predict(rf.list[[label1[i]]],data[, -ncol(data)])
  cor_train<-cor(pre_train,data$Index)
  #rmse_train<-rmse(pre_train,data$Index)
  # pre_test<-predict(rf,data_test[,1:(ncol(data_test)-1)])
  # cor_test<-cor(pre_test,data_test$Index)
  # rmse_test<-rmse(pre_test,data_test$Index)
  print(cor_train)
  # print(cor_test)
}
dir.create("Rda", showWarnings = FALSE)
save(rf.list, file = 'Rda/rflist.rda') 

#####重要性分析#####

md<-list()
mi<-list()
colindexf<-factor(colindex[-25],levels=colindex[-25])
impframe<-data.frame(index=colindexf)
for (i in 1:2){
  print(i)
  rf<-rf.list[[i]]
  imp<-importance(rf)
  incmse<-data.frame(index=names(imp[,2]),incmse=imp[,1]/max(imp[,1]))
  colnames(incmse)[2]<-paste0(i)
  impframe<-merge(impframe,incmse,by='index',all=T)
  min_depth_frame<-min_depth_distribution(rf)
  md[[label[i]]]<-min_depth_frame
  im_frame<-measure_importance(rf)
  im_frame[4]<-im_frame[4]/max(im_frame[4])
  im_frame[5]<-im_frame[5]/max(im_frame[5])
  mi[[label[i]]]<-im_frame
  }
#colnames(impframe)[2:22]<-label
impframe <- impframe[, 1:(ncol(impframe)-2)]
save(impframe,md,mi,file='Rda/multi-importance.rda')

#####重要性图#####

load(file='Rda/rflist.rda')
load(file='Rda/multi-importance.rda')
mdplot<-list()
miplot<-list()

for (i in 1:2){
  print(i)
  # 对于给定的索引 i，提取 min_depth_frame 和 im_frame
  min_depth_frame <- md[[i]]
  im_frame <- mi[[i]]
  
  # 对 min_depth_frame 进行绘图，并保存到 mdplot 列表中
  mdplot[[label[i]]] <- local({
    min_depth_frame = min_depth_frame
    plot_min_depth_distribution(min_depth_frame, k = 14) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(legend.title = element_text(size = rel(0.6)),
            legend.text = element_text(size = rel(0.5)))
  })
  ggsave(paste0('A', label[i], '.pdf'), width = 7, height = 7)  # 保存绘制的图形到 PDF 文件中
  
  # 对 im_frame 进行绘图，并保存到 miplot 列表中
  im_frame$p_value <- im_frame$p_value / 5  # 修改 p_value 的值
  miplot[[label[i]]] <- local({
    im_frame = im_frame
    plot_multi_way_importance(im_frame, x_measure = "mse_increase",
                              y_measure = "node_purity_increase",
                              size_measure = "p_value", no_of_labels = 5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.line = element_line(color = 'black'),
            axis.ticks.length = unit(0.4, 'line')) +
      coord_fixed(ratio = 1) +
      theme(legend.key.size = unit(1, 'line'), legend.position = c(0.15, 0.75))
  })
  ggsave(paste0('B', label[i], '.pdf'), width = 4, height = 4)  # 保存绘制的图形到 PDF 文件中
}

save(mdplot,miplot,file='Rda/importanceplot.rda')

#####特征交互作用分析#####

load(file='Rda/rflist.rda')
load(file='Rda/multi-importance.rda')

inter_list<-list()
for (i in 1:2){
  print(i)
  im_frame<-mi[[i]]
  rf<-rf.list[[i]]
  vars <- important_variables(im_frame, k = 5, measures = c("mean_min_depth","no_of_trees"))#计算前五个最重要的变量
  interactions_frame <- min_depth_interactions(rf, vars)
  interactions_frame <- arrange(interactions_frame,-interactions_frame[,4])
  inter_list[[label[i]]]<-interactions_frame
  #load(paste0('Rda/inter_',Label[i+2],'.rda'))
  #head(interactions_frame[order(interactions_frame$occurrences, decreasing = TRUE), ])
}
save(inter_list,file='Rda/inter1.rda')

fiplot<-list()
for (i in 1:2){
  interactions_frame<-inter_list[[i]]
  hlim<-ceiling(max(interactions_frame[1:25,3],interactions_frame[1:25,6]))
  fip<-plot_min_depth_interactions(interactions_frame,k=25)+
    scale_y_continuous(limits=c(0,hlim+1.5),expand=c(0,0))+
    scale_fill_gradient(low='#00bfc4',high='#f8766d')+                 #换色
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())   #+
  #theme(legend.position=c(0.3,0.6),legend.box="horizontal")  #(0.3,0.8),图例大小要调整，或者画在外面
  fiplot[[label[i]]]<-fip
  ggsave(paste0('interaction',label[i],'.pdf'),width=8,height=4)
}

save(fiplot,file='Rda/inter_plot.rda')

#####特征交互作用计算#####
source('min_depth_distribution.R')
source('measure_importance.R')
source('min_depth_interactions.R')
source('interaction.R')
rdlist<-list()
r_interaction<-list()
pb <- tkProgressBar("?ܽ???","?????? %", 0, 100) 
for (i in 1:length(label)){
  print(paste0(i,'/',length(label)))
  info<- sprintf("?????? %d%%", round(i*100/length(label))) 
  rf<-rf.list[[i]]
  c1<-rep(1:length(colindex),each=length(colindex))
  c11<-colindex[c1]
  c2<-rep(1:length(colindex),length(colindex))
  c22<-colindex[c2]
  rd_frame<-data.frame(c11,c22)
  colnames(rd_frame)=c('variable','root_variable')
  im_frame<-mi[[i]]
  rd_frame<-merge(rd_frame,im_frame[c(1,6)],all.x=T)
  
  pb1 <- tkProgressBar("??????","?????? %", 0, 100) 
  for (j in 1:rf$ntree){
    info1<- sprintf("?????? %d%%", round(i*100/rf$ntree)) 
    D<-calculate_max_depth(rf,1)
    interactions_frame_single<-min_depth_interactions_single(rf,j,colindex)
    rD<-calculate_rD(D,interactions_frame_single,1)
    rD<-cbind(interactions_frame_single[1:2],rD)
    rd_frame<-merge(rd_frame,rD,by=c('variable','root_variable'),all=T)
    setTkProgressBar(pb1, j*100/rf$ntree, sprintf("???? (%s)", info1),info1) 
  }
  close(pb1)
  
  rd_frame[is.na(rd_frame)]<-0
  rdlist[[label[i]]]<-rd_frame
  for (k in 1:nrow(rd_frame)){
    rd_frame[k,504]<-sum(rd_frame[k,4:503])/rd_frame[k,3]
  }
  r_frame<-rd_frame[c(1,2,504)]
  colnames(r_frame)<-c("variable" , "root_variable" ,"r")
  r_interaction[[label[i]]]<-r_frame
  #save(rd_frame, file = paste0("Rda/rd_frame_",label[i],".rda"))
  #save(r_frame,file = paste0("Rda/r_frame_",label[i],".rda"))
  setTkProgressBar(pb, i*100/length(label), sprintf("?ܽ??? (%s)", info),info) 
}
close(pb)
save(r_interaction,file='Rda/r-interaction.rda')

#----------------------------------------------
# 节点和边界文件
#----------------------------------------------
load(file='Rda/r-interaction.rda')
load(file='Rda/multi-importance.rda')

type<-data.frame(label=c(colindex[-21],label),
                 type=c(rep('M',15),rep('A',2),rep('E',3),rep('y',2)),
                 color=c(rep('#98dbef',15),rep('#a4e192',2),rep('#ffc177',3),
                         rep('#ffb6d4',2)))

for (i in 1:length(label)){
  nodes<-data.frame(id=c(1:length(colindex)),label=c(colindex[-21],label[i]))
  nodes<-merge(nodes,type,all.x=T)
  nodes<-arrange(nodes,nodes['id'])
  write.csv(nodes,paste0('network/nodes_',label[i],'.csv'),row.names=FALSE,fileEncoding='UTF-8')
  
  r_frame<-r_interaction[[i]]
  edges<-cbind(r_frame,c(rep('x-x',nrow(r_frame))))
  colnames(edges)<-c('Source','Target','Weight','Type')
  edges[is.na(edges)]<-0
  edges[3]<-edges[3]/max(edges[3])
  edges[3][edges[3]<0.5]<-0
  edges<-edges[-which(edges[3]==0),]
  edges<-edges[-which(edges[1]==edges[2]),]
  for (j in 1:nrow(edges)){
    j1<-which(edges[j,1]==edges[2])
    j2<-which(edges[j,2]==edges[1])
    j3<-intersect(j1,j2)
    if (length(j3)!=0){
      edges[j,3]<-mean(c(edges[j,3],edges[j3,3]))
      edges<-edges[-j3,]
    }
  }
  im_frame<-mi[[i]]
  if (i <13){
    x_y<-data.frame(Source=im_frame$variable,
                    Target=c(rep(label[i],20)),
                    Weight=im_frame[4],
                    Type=c(rep('x-y',20)))
  } else {
    x_y<-data.frame(Source=im_frame$variable,
                    Target=c(rep(label[i],28)),
                    Weight=im_frame[4],
                    Type=c(rep('x-y',28)))
  }
  colnames(x_y)<-c('Source','Target','Weight','Type')
  edges<-rbind(edges,x_y)
  edges[3][edges[3]<=0]<-0
  for (j in 1:nrow(nodes[1])){
    edges[edges==nodes[j,1]]<-nodes[j,2]
  }
  write.csv(edges,paste0('network/edges',label[i],'.csv'),row.names=FALSE,fileEncoding='UTF-8')
}

