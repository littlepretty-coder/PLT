library(networkD3)


data = read.delim('./metab/mediation_results.csv',sep=',')
tax = read.delim('./data/tax.txt')
tax = tax[,c(6:7)]
tax$Species <- gsub("[^a-zA-Z0-9_]", "", tax$Species) 
data = merge(data,tax,by.x = 'Bacteria',by.y = 'Species')
mapping <- list(
  fun_1 = fungi_1,
  fun_2 = fungi_2,
  fun_3 = fungi_3,
  fun_4 = fungi_4,
  bac_1 = bac_1,
  bac_2 = bac_2,
  bac_3 = bac_3,
  bac_4 = bac_4
)
data$Module <- 0
for (name in names(mapping)) {
  data$Module[data$Genus %in% mapping[[name]]] <- name
}
info = read_excel('./metab/all.xlsx')
info = info[,c(1,3)]
info$`Metabolite Name` <- gsub("[^a-zA-Z0-9_]", "", info$`Metabolite Name`)
info$`Metabolite Name`<- sapply(info$`Metabolite Name`, function(name) {
  if (grepl("^[0-9]", name)) { # 检查是否以数字开头
    paste0("X", name) # 添加前缀 "X"
  } else {
    name # 保持原列名
  }
})
colnames(info)[1] = 'Metabolite'
data0 = merge(data,info,by= 'Metabolite')
write.table(data0,'./figure/suppdata/fig4/Discovery_20metabolites_modules_mediation.csv',sep=',',row.names = F)
nodes <- data.frame(name = unique(c(data0$Module, data0$Class)))
nodes = nodes[c(7,6,1,3,2,5,8,4,9:14),,drop=F]
# 创建链接数据
links <- data.frame(
  source = match(data0$Module, nodes$name) - 1,  # 将源节点转换为索引
  target = match(data0$Class, nodes$name) - 1, # 将目标节点转换为索引
  value = 0.2 # 每个连接赋值为 1，表示权重
)

colour_mapping <- c(
  "fun_1" = "#DBB428", "fun_2" = "#D4562E",
  "fun_3" = "#FFFFB3", "fun_4" = "#FCCDE5",
  "bac_1" = "#A51C36", "bac_2" = "#7ABBDB",
  "bac_3" = "#84BA42", "bac_4" = "#682487"
)

# 将颜色映射为 JSON 格式，供 D3.js 使用
colour_scale <- sprintf(
  'd3.scaleOrdinal().domain(%s).range(%s)',
  jsonlite::toJSON(names(colour_mapping)), 
  jsonlite::toJSON(unname(colour_mapping))
)
customJS <- "
function(el, x) {
  d3.select(el).selectAll('.node text')
    .attr('x', function(d) { return d.x < x.options.width / 2 ? -10 : 10; }) // 调整 x 坐标
    .attr('text-anchor', function(d) { return d.x < x.options.width / 2 ? 'end' : 'start'; }); // 设置对齐方式
}
"

links$color <- nodes$name[links$source + 1]
sankey = sankeyNetwork(
  Links = links, 
  Nodes = nodes, 
  Source = "source", 
  Target = "target", 
  Value = "value", 
  NodeID = "name", 
  units = "Tons", 
  fontSize = 12, 
  nodeWidth = 30,
  colourScale = colour_scale, # 设置节点颜色
  LinkGroup = "color" ,
  linkType = 'bezier'
  
  # 让流线颜色与目标节点一致
)
#saveNetwork(sankey, "sankey_diagram.html")
# 显示图表
sankey$x$options$NodeJS <- htmlwidgets::JS(customJS)
sankey
install.packages("webshot")

library(webshot)
saveNetwork(sankey, "sankey_diagram.html")
webshot("sankey_diagram.html" , "sankey.pdf")
