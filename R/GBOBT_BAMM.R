##reviews for the paper for EcolLett
setwd('~/Dropbox/Work_in_progress/LDG_plants/')
#dir.create('./reviews_EcolLett/')
#dir.create('./output/trees/')
#dir.create('./output/trees/GBOTB_clades/')

library(V.PhyloMaker)
library(taxonlookup)
library(caper)
library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(picante)
library(BAMMtools)


tree<-GBOTB.extended
write.tree(GBOTB.extended,file='./raw_data/GBOTB_extended.tree')
#GBOTB is ultrametric
is.ultrametric(GBOTB.extended)
#but not binary
is.binary.phylo(tree)
#get the polytomic nodes (they have >2 children nodes)
children.nodes<-lapply(c(74532:149004),function(x)Children(tree,x))
#there are 7 nodes with polytomies
polytomic.nodes<-lapply(which(unlist(lapply(children.nodes,function(x)length(x)>2))),function(x)extract.clade(tree,node=x+74531))
#they are very small clades (terminal branches only)
# 5  3  4 55 45  7  3: number of tips in each clade
unlist(lapply(polytomic.nodes,function(x)length(x$tip.label)))
#resolve the polytomies in random order
tree.bif <- multi2di(tree, random = TRUE)
#some branches now have zero length; fix by adding 1e-05
tree.bif$edge.length[tree.bif$edge.length==0]<-1e-05
is.binary(tree.bif)
write.tree(tree.bif,file='./raw_data/GBOTB_extended_bif.tree')
#drop species that are not angiosperms and on taxonlookup
GBOTB.lookup<-lookup_table(tree.bif$tip.label,by_species = T,include_counts = T)
GBOTB.bif.angios<-drop.tip(tree.bif,tip =setdiff(tree.bif$tip.label,rownames(GBOTB.lookup)[GBOTB.lookup$group=='Angiosperms']))
GBOTB.lookup.angios<-GBOTB.lookup[GBOTB.lookup$group=='Angiosperms',]
#n=72986 angiosperms, GBOTB has been through TPL
write.tree(GBOTB.bif.angios,file='./raw_data/GBOTB_extended_bif_angios.tree')
#add species names to table
GBOTB.lookup.angios$species<-rownames(GBOTB.lookup.angios)
rownames(GBOTB.lookup.angios)<-NULL
#get families and number of species per family
PlantLookup<-plant_lookup(include_counts = TRUE)
PlantLookup<-PlantLookup[PlantLookup$group=='Angiosperms',]
PlantLookup.2<-aggregate(number.of.accepted.species~family, data=PlantLookup, sum)
colnames(PlantLookup.2)[2]<-'number.of.species.family'
#add to species level table
table.sampling<-merge(GBOTB.lookup.angios,PlantLookup.2,by='family',all.x=TRUE)
#count number of species in each family on the tree
family.count.table<-as.data.frame(table(table.sampling$family),stringsAsFactors = F)
colnames(family.count.table)<-c('family','number.of.species.tree')
family.count.table<-family.count.table[family.count.table$number.of.species.tree>0,]
table.sampling<-merge(table.sampling,family.count.table,by='family',all.x=TRUE)
#for species with number.of.species.family=0 (there are no accepted species but there are unresolved species), use number of unresolved
#they are monogeneric families  so this will do
table.sampling[table.sampling$number.of.species.family==0,]$number.of.species.family<-table.sampling[table.sampling$number.of.species.family==0,]$number.of.accepted.and.unresolved.species
table.sampling$family.sampling.fraction<-table.sampling$number.of.species.tree/table.sampling$number.of.species.family
#adding info on taxonomic discordance (when there are more tips in the tree than accepted species in a family according to taxonlookup)
#410 species in 10 families
table.sampling$taxonomic.discordance<-NA
table.sampling[table.sampling$number.of.species.family<table.sampling$number.of.species.tree,]$taxonomic.discordance<-1
table.sampling[table.sampling$number.of.species.family>=table.sampling$number.of.species.tree,]$taxonomic.discordance<-0
#remove taxonomic discordant species
#remove species with taxonomic discordance
table.sampling<-table.sampling[!table.sampling$taxonomic.discordance==1,]
tree<-drop.tip(GBOTB.bif.angios,setdiff(GBOTB.bif.angios$tip.label,table.sampling$species))
#72576 species
Family_BAMM_Sampling<-table.sampling[,c("species","family","family.sampling.fraction")]
write.table(Family_BAMM_Sampling,file='./output/tables/GBOTB_BAMM_family_sampling.txt')


####this below splits the tree into clades to run BAMM
alltips<-tree$tip.label
all.tips.length<-lapply(alltips,function(x) length(x))
all.tips.length.df<-as.data.frame(unlist(all.tips.length))
all.tips.length.df$node.number<-c(1:nrow(all.tips.length.df))
all.tips.length.df$node.number<-all.tips.length.df$node.number+length(tree$tip.label)
colnames(all.tips.length.df)[1]<-'n.descendant.tips'

node.numbers<-c(1:length(tree$node.label))
node.numbers<-node.numbers+length(tree$tip.label)
length.node.descendants<-sapply(node.numbers,function(x) length(unlist(Descendants(tree,x,type='tips'))))
length.node.descendants<-as.data.frame(length.node.descendants)
length.node.descendants$parental.node<-node.numbers
colnames(length.node.descendants)[1]<-'n.descendant.tips'

node.age(tree)->phy.age
cbind(phy.age$edge,phy.age$age, tree$edge.length)->BL.position
max(phy.age$age)-BL.position[,3]->dist.tip
cbind(BL.position,dist.tip)->BL.positions
BL.positions[,5]+BL.positions[,4]->ages
cbind(BL.positions,ages)->BL.positions
as.data.frame(BL.positions)->node.ages
names(node.ages)<-c("parental.node","daughter.node","dist.root","BL","dist.tip","mrca.age")
node.ages$node.label<-NA

node.labels.df<-as.data.frame(c(1:length(tree$node.label)),tree$node.label)
node.labels.df$node.label<-row.names(node.labels.df)
row.names(node.labels.df)<-NULL
colnames(node.labels.df)[1]<-c('node.number')
node.labels.df$node.number<-node.labels.df$node.number+length(tree$tip.label)
node.ages<-merge(node.labels.df,node.ages,by.x='node.number','parental.node')
colnames(node.ages)[1]<-'parental.node'
node.ages<-merge(node.ages,length.node.descendants,by.x='parental.node',by.y='parental.node',all.x=TRUE)

nodes.to.check<-unique(node.ages$parental.node)
#start with Lamiidae, it's 4126 species in a monophyletic clade

lamiidae.node<-node.ages[node.ages$node.label=='Lamiidae',]$parental.node[1]
chosen.nodes<-vector()
length.chosen.clades<-vector()

chosen.nodes[1]<-lamiidae.node
length.chosen.clades<-node.ages[node.ages$node.label=='Lamiidae',]$n.descendant.tips[1]

if(node.ages[node.ages$parental.node==Siblings(tree,chosen.nodes[1])]$n.descendant.tips)
  
  parent.node<-node.ages[(node.ages$daughter.node==chosen.nodes[length(chosen.nodes)]),]$parental.node[1]
sister.node<-node.ages[(node.ages$parental.node==parent.node)&(node.ages$daughter.node!=chosen.nodes[length(chosen.nodes)]),]$daughter.node

while(sum(length.chosen.clades)<length(tree$tip.label)){
  parent.node<-node.ages[(node.ages$daughter.node==chosen.nodes[length(chosen.nodes)]),]$parent.node[1]
  parent.node.df<-node.ages[node.ages$daughter.node==chosen.nodes[length(chosen.nodes)],]$parent.node[1]
}

#this function creates a set of non inclusive non monophyletic clades of <size
#starting with lamiidae.node
nodenumber<-lamiidae.node<-node.ages[node.ages$node.label=='Lamiidae',]$parental.node[1]
size<-5000

recursive_collapse_clades<-function(nodenumber,tree,size,table){
  chosen.nodes<-vector()
  length.chosen.clades<-vector()
  collapsed.species<-vector()
  while(sum(length.chosen.clades)<length(tree$tip.label)){
    #get number of descendants of node
    length.descendants<-length(unlist(Descendants(tree,nodenumber,type='tips')))
    #if number of descendants of node > size, go to left daughter
    sister<-Siblings(tree,nodenumber)
    
    length.sister.descendants<-length(unlist(Descendants(tree,sister,type='tips')))
    if (length.descendants>size){
      nodenumber<-Children(tree,nodenumber)[1]
    }
    #check if sister is already in chosen.nodes vector
    #if it is, add "descendants of"cousins"
    if(length(which(chosen.nodes==sister))>0){
      
      
      
      #get number of descendants of sister node
      
      #if sum of descendants of node + sister > size
      #keep node as one clade
      if(length.descendants+length.sister.descendants>size){
        chosen.nodes<-c(chosen.nodes,nodenumber)
        length.chosen.clades<-c(length.chosen.clades,length.descendants)
        
        node.number<-sister
      }else{
        #if sum of descendants of node + sister < size, collapse (go to ancestor)
        nodenumber<-Ancestors(tree,nodenumber,type='parent')
      }
      
      
    }else{
      
    }
    
    
  }
  
  
  #check n of descendants of sister in table
  #if sum of descendants of node and sister < size, collapse and nodenumber is the ancestor node
  #if not, store nodenumber and size in vector and use sister node as node number
  
  
}





Commelinidae<-getDescendants(tree, node=which(tree$node.label=='Commelinidae')+length(tree$tip.label))
length(Commelinidae[Commelinidae<=length(tree$tip.label)])
Commelinidae.species<-tree$tip.label[Commelinidae[Commelinidae<=length(tree$tip.label)]]
# 3367 species
Monocotyledoneae<-getDesister.nodescendants(tree, node=which(tree$node.label=='Monocotyledoneae')+length(tree$tip.label))
length(Monocotyledoneae[Monocotyledoneae<length(tree$tip.label)])
Monocotyledoneae.species<-tree$tip.label[Monocotyledoneae[Monocotyledoneae<=length(tree$tip.label)]]

Monocotyledoneae.species.noCom<-Monocotyledoneae.species[!(Monocotyledoneae.species %in% Commelinidae.species)]
#3591 species

which(tree$node.label=='Commelinidae')+length(tree$tip.label)



###select clades < 5000 species
node.ages.5000<-node.ages[node.ages$n.descendant.tips<5000,]
node.ages.5000<-node.ages.6000[,c('parental.node','n.descendant.tips')]
node.ages.5000<-unique(node.ages.5000)
node.ages.5000<-node.ages.5000[order(-node.ages.5000$n.descendant.tips),]
selected.species<-vector()
selected.nodes<-vector()
selected.lengths<-vector()
for (i in 1:nrow(node.ages.5000)){
  descendant.species<-unlist(Descendants(tree,node.ages.5000[i,'parental.node'],type='tips'))
  if(length(intersect(descendant.species,selected.species))==0){
    selected.species<-c(selected.species,descendant.species)
    selected.nodes<-c(selected.nodes,node.ages.5000[i,'parental.node'])
    selected.lengths<-c(selected.lengths,length(descendant.species))
    
  }
  cat(sum(selected.lengths),'\n')
  if(sum(selected.lengths)==length(tree$tip.label)){
    cat('done','\n')
    break
  }
}
#get one species per clade
species.selected<-tree$tip.label[sapply(selected.nodes,function(x) unlist(Descendants(tree,x,type='tips'))[1])]
tree.clades<-drop.tip(tree,setdiff(tree$tip.label,species.selected))
#add the number of species in each clade
species.selected.length<-paste(species.selected,selected.lengths,sep='_')
order.vector<-sapply(species.selected[match(tree.clades$tip.label,species.selected)],function(x) grep(x,species.selected.length))
species.selected.length<-species.selected.length[order(order)]
tree.clades$tip.label<- species.selected.length[order]
plot(tree.clades)


###select clades < 6000 species
node.ages.6000<-node.ages[node.ages$n.descendant.tips<6000,]
node.ages.6000<-node.ages.6000[,c('parental.node','n.descendant.tips')]
node.ages.6000<-unique(node.ages.6000)
node.ages.6000<-node.ages.6000[order(-node.ages.6000$n.descendant.tips),]
selected.species<-vector()
selected.nodes<-vector()
selected.lengths<-vector()
#this looks through the whole tree
for (i in 1:nrow(node.ages.6000)){
  descendant.species<-unlist(Descendants(tree,node.ages.6000[i,'parental.node'],type='tips'))
  if(length(intersect(descendant.species,selected.species))==0){
    selected.species<-c(selected.species,descendant.species)
    selected.nodes<-c(selected.nodes,node.ages.6000[i,'parental.node'])
    selected.lengths<-c(selected.lengths,length(descendant.species))
    
  }
  cat(sum(selected.lengths),'\n')
  if(sum(selected.lengths)==length(tree$tip.label)){
    cat('done','\n')
    break
  }
}
#add Amborella 
amborella.node<-which(tree$tip.label=='Amborella_trichopoda')
selected.nodes.Ambo<-c(selected.nodes,amborella.node)
selected.lengths.Ambo<-c(selected.lengths,1)
#get one species per clade
species.selected<-tree$tip.label[sapply(selected.nodes.Ambo,function(x) unlist(Descendants(tree,x,type='tips'))[1])]
tree.clades<-drop.tip(tree,setdiff(tree$tip.label,species.selected))
#add the number of species in each clade
species.selected.length<-paste(species.selected,selected.lengths.Ambo,sep='_')
order.vector<-sapply(species.selected[match(tree.clades$tip.label,species.selected)],function(x) grep(x,species.selected.length))

tree.clades$tip.label<- species.selected.length[order.vector]
plot(tree.clades)

#check monophyly of clades
species.selected.all<-sapply(selected.nodes.Ambo,function(x) unlist(Descendants(tree,x,type='tips')))
monophyly.check<-sapply(species.selected.all,function(x) is.monophyletic(tree,tree$tip.label[x]))
#they're all monophyletic
#then extract the clades > 3000 spp
big.monophyletic.clades<-lapply(selected.nodes.Ambo[selected.lengths.Ambo>2000], function(x) extract.clade(tree,node=x))
#and extract the rest of the species with the backbone analysis
big.monophyletic.clade.species<-lapply(selected.nodes.Ambo[selected.lengths.Ambo>2000],function(x) tree$tip.label[unlist(Descendants(tree,x,'tips'))])
#keep the first species on each big monophyletic clade (not dropping those to get the full backbone)
monophyletic.clade.species.to.drop<-unlist(lapply(big.monophyletic.clade.species,function(x) x[-1]))
backbone.clade<-drop.tip(tree,monophyletic.clade.species.to.drop)
#double check that there's no duplicated species in big monophyletic clades
length(unique(unlist(sapply(big.monophyletic.clades,function(x) x$tip.label))))
length(unlist(sapply(big.monophyletic.clades,function(x) x$tip.label)))
#get sampling fraction for each clade
BAMM.clades.to.run<-big.monophyletic.clades
BAMM.clades.to.run[[17]]<-backbone.clade
length(unique(unlist(sapply(BAMM.clades.to.run,function(x) x$tip.label))))
#there's 16 duplicates (one for each clade)
length(unlist(sapply(BAMM.clades.to.run,function(x) x$tip.label)))
#get the dataframes with sampling fractions
#it's easy for the big monopohyletic clades
big.monophyletic.clades.sampling<-lapply(big.monophyletic.clades,function(x) table.sampling[(table.sampling$species %in% x$tip.label),])
big.monophyletic.clades.sampling.BAMM<-lapply(big.monophyletic.clades.sampling,function(x) x[,c("species","family","family.sampling.fraction")])
big.monophyletic.clades.sampling.BAMM<-lapply(big.monophyletic.clades.sampling.BAMM,function(x) {colnames(x)[1]<-'species';return(x)})

#length of tips in each monophyletic clade
unlist(sapply(big.monophyletic.clades,function(x)length(x$tip.label)))
unlist(sapply(big.monophyletic.clades.sampling,function(x)nrow(x)))
#calculate the backbone sampling for these as well (this will be useful for the backbone clade as well)
big.monophyletic.clades.backbone.sampling<-lapply(big.monophyletic.clades.sampling,function(x) unique(x[,c('family','number.of.species.family','number.of.species.tree')]))
big.monophyletic.clades.backbone.samplingfraction<-unlist(lapply(big.monophyletic.clades.backbone.sampling,function(x) sum(x$number.of.species.tree)/sum(x$number.of.species.family)))
#now calculate the sampling fractions for the backbone clade
backbone.clades.sampling<-table.sampling[(table.sampling$species %in% backbone.clade$tip.label),]
backbone.clades.sampling.BAMM<-backbone.clades.sampling[,c("species","family","family.sampling.fraction")]
colnames(backbone.clades.sampling.BAMM)[1]<-c('species')
#substitute the sampling in the species representative of each monophyletic clade
big.monophyletic.clade.species.representative<-unlist(sapply(big.monophyletic.clade.species,function(x)x[1]))
new.sampling.fractions.clades<-1/(unlist(lapply(big.monophyletic.clades.backbone.sampling,function(x) sum(x$number.of.species.family))))
#replace in dataframe
for (i in 1:length(big.monophyletic.clade.species.representative)){
  backbone.clades.sampling.BAMM[backbone.clades.sampling.BAMM$species==big.monophyletic.clade.species.representative[i],]$family.sampling.fraction<-new.sampling.fractions.clades[i]
}
BAMM.clades.sampling.df<-big.monophyletic.clades.sampling.BAMM
BAMM.clades.sampling.df[[17]]<-backbone.clades.sampling.BAMM
#check that lengths of tree & sampling df are identical
unlist(lapply(BAMM.clades.to.run,function(x) length(x$tip.label)))
unlist(lapply(BAMM.clades.sampling.df,function(x) nrow(x)))
#calculate the sampling fraction of the backbone clade analysis
BAMM.clades.sampling.backbone.fractions<-c(big.monophyletic.clades.backbone.samplingfraction,nrow(BAMM.clades.sampling.df[[17]])/sum(PlantLookup$number.of.species.family))
#write output for BAMM
#first do checks and runsetBAMMpriors
lapply(BAMM.clades.to.run,function(x){c(is.binary.tree(x),is.ultrametric(x),min(x$edge.length),summary(duplicated(x$tip.label)))})
#they're all ok
#now run setBAMMpriors and store in files
for (i in 1:length(BAMM.clades.to.run)){
  setBAMMpriors(BAMM.clades.to.run[[i]],outfile = paste('./output/trees/GBOTB_clades/clade_',i,'_priors.txt',sep=''))
}
#write treefiles
for (i in 1:length(BAMM.clades.to.run)){
  write.tree(BAMM.clades.to.run[[i]],file = paste('./output/trees/GBOTB_clades/clade_',i,'.tree',sep=''))
}
#write sampling fractions,first backbones
for (i in 1:length(BAMM.clades.sampling.df)){
  write.table(BAMM.clades.sampling.backbone.fractions[i],file = paste('./output/trees/GBOTB_clades/clade_',i,'_sampling.txt',sep=''),quote=F,row.names=F,sep='\t',col.names = F)
}

#write sampling fractions,then table
for (i in 1:length(BAMM.clades.sampling.df)){
  write.table(BAMM.clades.sampling.df[[i]],file = paste('./output/trees/GBOTB_clades/clade_',i,'_sampling.txt',sep=''),quote=F,row.names=F,sep='\t',col.names = F,append = T)
}

#had to add extinctionProbMax = 0.99999 to the control file of clade 17 (the backbone)
#explained here: https://groups.google.com/forum/#!topic/bamm-project/X4m5BQkl7ls






