load(file = "/Users/millican/Documents/Research/Koch_Lab/Manuscripts/Turfgrass_Ecology/ITS_Forward_PhyloObject.RData")
sample_data(ps.ITS)$Location = factor(sample_data(ps.ITS)$Location, levels = c("Creeping Bent Greens", "Creeping Bent Fairway", "Tall Fescue", "Natural"))
sample_data(ps.ITS)$Section = factor(sample_data(ps.ITS)$Section, levels = c("Leaf", "Thatch", "Rhizosphere", "Bulk"))



ps.ITS.above.ground = subset_samples(ps.ITS, Section == "Leaf" | Section == "Thatch")
ps.ITS.below.ground = subset_samples(ps.ITS, Section == "Rhizosphere" | Section == "Bulk")


ps.ITS.occ.abun = Occ.Abund(ps.ITS, seed = 8675309, depth = min(sample_sums(ps.ITS)))

spp.ITS = data.frame(otu_table(rarefy_even_depth(ps.ITS, sample.size = min(sample_sums(ps.ITS)), rngseed = 8675309)))
taxon.ITS = as.vector(colnames(data.frame(otu_table(rarefy_even_depth(ps.ITS, sample.size = min(sample_sums(ps.ITS)), rngseed = 8675309)))))


obs.ITS = sncm.fit(spp.ITS, taxon.ITS, stats=FALSE, pool=NULL)
sta.ITS = sncm.fit(spp.ITS, taxon.ITS, stats=TRUE, pool=NULL)

ggplot() +
  geom_point(data = ps.ITS.occ.abun, aes(x = log10(abund), y = occ), pch = 21, fill = '#37AEC3', size = 3, alpha=.4)+
  geom_line(color = '#d60000', data = obs.ITS, size=1, aes(y = freq.pred, x = log10(p)), alpha = .75) +
  geom_line(color = '#d60000', data = obs.ITS, lty = 'twodash', size=1, aes(y = pred.upr, x = log10(p)), alpha = .75)+
  geom_line(color = '#d60000', data = obs.ITS, lty = 'twodash', size=1, aes(y = pred.lwr, x = log10(p)), alpha = .75)+
  labs(x = expression(Log[10](mean~relative~abundance)), y = "Occupancy") + 
  theme_bw(base_size = 22, base_family = "Times")
ggsave(file = "Occ_Abund_metacommunity_ITS.tiff", dpi = 900, width = 12, height = 12, units = 'in')


ps.ITS.above.ground.occ.abun = Occ.Abund(ps.ITS.above.ground, seed = 8675309, 
                                         depth = min(sample_sums(ps.ITS.above.ground)))

ps.ITS.below.ground.occ.abun = Occ.Abund(ps.ITS.below.ground, seed = 8675309, 
                                         depth = min(sample_sums(ps.ITS.below.ground)))
ggplot(ps.ITS.occ.abun) + 
  geom_point(aes(x = log10(abund), y = occ))

ggplot(ps.ITS.above.ground.occ.abun) + 
  geom_point(aes(x = log10(abund), y = occ))

ggplot(ps.ITS.below.ground.occ.abun) + 
  geom_point(aes(x = log10(abund), y = occ))


spp.ITS.above.ground = data.frame(otu_table(rarefy_even_depth(ps.ITS.above.ground, sample.size = min(sample_sums(ps.ITS.above.ground)), rngseed = 8675309)))
taxon.ITS.above.ground = as.vector(colnames(data.frame(otu_table(rarefy_even_depth(ps.ITS.above.ground, sample.size = min(sample_sums(ps.ITS.above.ground)), rngseed = 8675309)))))


obs.ITS.above.ground = sncm.fit(spp.ITS.above.ground, taxon.ITS.above.ground, stats=FALSE, pool=NULL)
sta.ITS.above.ground = sncm.fit(spp.ITS.above.ground, taxon.ITS.above.ground, stats=TRUE, pool=NULL)

ggplot() +
  geom_point(data = ps.ITS.above.ground.occ.abun, aes(x = log10(abund), y = occ), pch = 21, fill = 'white', size = 3, alpha=.6)+
  geom_line(color = 'black', data = obs.ITS.above.ground, size=1, aes(y = freq.pred, x = log10(p)), alpha = .75) +
  geom_line(color = 'black', data = obs.ITS.above.ground, lty = 'twodash', size=1, aes(y = pred.upr, x = log10(p)), alpha = .75)+
  geom_line(color = 'black', data = obs.ITS.above.ground, lty = 'twodash', size=1, aes(y = pred.lwr, x = log10(p)), alpha = .75)+
  labs(x = expression(Log[10](mean~relative~abundance)), y = "Occupancy") + 
  theme_bw(base_size = 28, base_family = "Times")


spp.ITS.below.ground = data.frame(otu_table(rarefy_even_depth(ps.ITS.below.ground, sample.size = min(sample_sums(ps.ITS.below.ground)), rngseed = 8675309)))
taxon.ITS.below.ground = as.vector(colnames(data.frame(otu_table(rarefy_even_depth(ps.ITS.below.ground, sample.size = min(sample_sums(ps.ITS.below.ground)), rngseed = 8675309)))))


obs.ITS.below.ground = sncm.fit(spp.ITS.below.ground, taxon.ITS.below.ground, stats=FALSE, pool=NULL)
sta.ITS.below.ground = sncm.fit(spp.ITS.below.ground, taxon.ITS.below.ground, stats=TRUE, pool=NULL)

ggplot() +
  geom_point(data = ps.ITS.below.ground.occ.abun, aes(x = log10(abund), y = occ), pch = 21, fill = 'white', size = 3, alpha=.6)+
  geom_line(color = 'black', data = obs.ITS.below.ground, size=1, aes(y = freq.pred, x = log10(p)), alpha = .75) +
  geom_line(color = 'black', data = obs.ITS.below.ground, lty = 'twodash', size=1, aes(y = pred.upr, x = log10(p)), alpha = .75)+
  geom_line(color = 'black', data = obs.ITS.below.ground, lty = 'twodash', size=1, aes(y = pred.lwr, x = log10(p)), alpha = .75)+
  labs(x = expression(Log[10](mean~relative~abundance)), y = "Occupancy") + 
  theme_bw(base_size = 28, base_family = "Times")

#####  Graphs for ITS Neutral model  #####

ggplot() +
  geom_point(data = ps.ITS.occ.abun, aes(x = log10(abund), y = occ), pch = 21, fill = '#37AEC3', size = 3, alpha=.4)+
  geom_line(color = '#d60000', data = obs.ITS, size=1, aes(y = freq.pred, x = log10(p)), alpha = .75) +
  geom_line(color = '#d60000', data = obs.ITS, lty = 'twodash', size=1, aes(y = pred.upr, x = log10(p)), alpha = .75)+
  geom_line(color = '#d60000', data = obs.ITS, lty = 'twodash', size=1, aes(y = pred.lwr, x = log10(p)), alpha = .75)+
  labs(x = expression(Log[10](mean~relative~abundance)), y = "Occupancy") + 
  theme_bw(base_size = 22, base_family = "Times")
ggsave(file = "Occ_Abund_metacommunity_ITS.tiff", dpi = 900, width = 12, height = 12, units = 'in')


ggplot() +
  geom_point(data = ps.ITS.above.ground.occ.abun, aes(x = log10(abund), y = occ), pch = 21, fill = '#37AEC3', size = 3, alpha=.6)+
  geom_line(color = '#d60000', data = obs.ITS.above.ground, size=1, aes(y = freq.pred, x = log10(p)), alpha = .75) +
  geom_line(color = '#d60000', data = obs.ITS.above.ground, lty = 'twodash', size=1, aes(y = pred.upr, x = log10(p)), alpha = .75)+
  geom_line(color = '#d60000', data = obs.ITS.above.ground, lty = 'twodash', size=1, aes(y = pred.lwr, x = log10(p)), alpha = .75)+
  labs(x = expression(Log[10](mean~relative~abundance)), y = "Occupancy") + 
  theme_bw(base_size = 28, base_family = "Times")
ggsave(file = "Occ_Abund_Above_Ground_ITS.tiff", dpi = 900, width = 12, height = 12, units = 'in')

ggplot() +
  geom_point(data = ps.ITS.below.ground.occ.abun, aes(x = log10(abund), y = occ), pch = 21, fill = '#37AEC3', size = 3, alpha=.6)+
  geom_line(color = '#d60000', data = obs.ITS.below.ground, size=1, aes(y = freq.pred, x = log10(p)), alpha = .75) +
  geom_line(color = '#d60000', data = obs.ITS.below.ground, lty = 'twodash', size=1, aes(y = pred.upr, x = log10(p)), alpha = .75)+
  geom_line(color = '#d60000', data = obs.ITS.below.ground, lty = 'twodash', size=1, aes(y = pred.lwr, x = log10(p)), alpha = .75)+
  labs(x = expression(Log[10](mean~relative~abundance)), y = "Occupancy") + 
  theme_bw(base_size = 28, base_family = "Times")
ggsave(file = "Occ_Abund_Below_Ground_ITS.tiff", dpi = 900, width = 12, height = 12, units = 'in')



ps.ITS.above.ground.green = subset_samples(ps.ITS.above.ground, Location == "Creeping Bent Greens")
ps.ITS.above.ground.fair = subset_samples(ps.ITS.above.ground, Location == "Creeping Bent Fairway")
ps.ITS.above.ground.tall = subset_samples(ps.ITS.above.ground, Location == "Tall Fescue")
ps.ITS.above.ground.nat = subset_samples(ps.ITS.above.ground, Location == "Natural")

ps.ITS.above.ground.green.occ.abun = Occ.Abund(ps.ITS.above.ground.green, seed = 8675309, 
                                               depth = min(sample_sums(ps.ITS.above.ground.green)))
ps.ITS.above.ground.fair.occ.abun = Occ.Abund(ps.ITS.above.ground.fair, seed = 8675309, 
                                              depth = min(sample_sums(ps.ITS.above.ground.fair)))
ps.ITS.above.ground.tall.occ.abun = Occ.Abund(ps.ITS.above.ground.tall, seed = 8675309, 
                                              depth = min(sample_sums(ps.ITS.above.ground.tall)))
ps.ITS.above.ground.nat.occ.abun = Occ.Abund(ps.ITS.above.ground.nat, seed = 8675309, 
                                             depth = min(sample_sums(ps.ITS.above.ground.nat)))


spp.ITS.above.ground.green = data.frame(otu_table(rarefy_even_depth(ps.ITS.above.ground.green, sample.size = min(sample_sums(ps.ITS.above.ground.green)), rngseed = 8675309)))
taxon.ITS.above.ground.green = as.vector(colnames(data.frame(otu_table(rarefy_even_depth(ps.ITS.above.ground.green, sample.size = min(sample_sums(ps.ITS.above.ground.green)), rngseed = 8675309)))))

sta.ITS.above.ground.green = sncm.fit(spp.ITS.above.ground.green, taxon.ITS.above.ground.green, stats=TRUE, pool=NULL)


spp.ITS.above.ground.fair = data.frame(otu_table(rarefy_even_depth(ps.ITS.above.ground.fair, sample.size = min(sample_sums(ps.ITS.above.ground.fair)), rngseed = 8675309)))
taxon.ITS.above.ground.fair = as.vector(colnames(data.frame(otu_table(rarefy_even_depth(ps.ITS.above.ground.fair, sample.size = min(sample_sums(ps.ITS.above.ground.fair)), rngseed = 8675309)))))

sta.ITS.above.ground.fair = sncm.fit(spp.ITS.above.ground.fair, taxon.ITS.above.ground.fair, stats=TRUE, pool=NULL)


spp.ITS.above.ground.tall = data.frame(otu_table(rarefy_even_depth(ps.ITS.above.ground.tall, sample.size = min(sample_sums(ps.ITS.above.ground.tall)), rngseed = 8675309)))
taxon.ITS.above.ground.tall = as.vector(colnames(data.frame(otu_table(rarefy_even_depth(ps.ITS.above.ground.tall, sample.size = min(sample_sums(ps.ITS.above.ground.tall)), rngseed = 8675309)))))

sta.ITS.above.ground.tall = sncm.fit(spp.ITS.above.ground.tall, taxon.ITS.above.ground.tall, stats=TRUE, pool=NULL)


spp.ITS.above.ground.nat = data.frame(otu_table(rarefy_even_depth(ps.ITS.above.ground.nat, sample.size = min(sample_sums(ps.ITS.above.ground.nat)), rngseed = 8675309)))
taxon.ITS.above.ground.nat = as.vector(colnames(data.frame(otu_table(rarefy_even_depth(ps.ITS.above.ground.nat, sample.size = min(sample_sums(ps.ITS.above.ground.nat)), rngseed = 8675309)))))

sta.ITS.above.ground.nat = sncm.fit(spp.ITS.above.ground.nat, taxon.ITS.above.ground.nat, stats=TRUE, pool=NULL)


sta.ITS.above.ground.green$Rsqr
sta.ITS.above.ground.fair$Rsqr
sta.ITS.above.ground.tall$Rsqr
sta.ITS.above.ground.nat$Rsqr



ps.ITS.below.ground.green = subset_samples(ps.ITS.below.ground, Location == "Creeping Bent Greens")
ps.ITS.below.ground.fair = subset_samples(ps.ITS.below.ground, Location == "Creeping Bent Fairway")
ps.ITS.below.ground.tall = subset_samples(ps.ITS.below.ground, Location == "Tall Fescue")
ps.ITS.below.ground.nat = subset_samples(ps.ITS.below.ground, Location == "Natural")

ps.ITS.below.ground.green.occ.abun = Occ.Abund(ps.ITS.below.ground.green, seed = 8675309, 
                                               depth = min(sample_sums(ps.ITS.below.ground.green)))
ps.ITS.below.ground.fair.occ.abun = Occ.Abund(ps.ITS.below.ground.fair, seed = 8675309, 
                                              depth = min(sample_sums(ps.ITS.below.ground.fair)))
ps.ITS.below.ground.tall.occ.abun = Occ.Abund(ps.ITS.below.ground.tall, seed = 8675309, 
                                              depth = min(sample_sums(ps.ITS.below.ground.tall)))
ps.ITS.below.ground.nat.occ.abun = Occ.Abund(ps.ITS.below.ground.nat, seed = 8675309, 
                                             depth = min(sample_sums(ps.ITS.below.ground.nat)))


spp.ITS.below.ground.green = data.frame(otu_table(rarefy_even_depth(ps.ITS.below.ground.green, sample.size = min(sample_sums(ps.ITS.below.ground.green)), rngseed = 8675309)))
taxon.ITS.below.ground.green = as.vector(colnames(data.frame(otu_table(rarefy_even_depth(ps.ITS.below.ground.green, sample.size = min(sample_sums(ps.ITS.below.ground.green)), rngseed = 8675309)))))

sta.ITS.below.ground.green = sncm.fit(spp.ITS.below.ground.green, taxon.ITS.below.ground.green, stats=TRUE, pool=NULL)


spp.ITS.below.ground.fair = data.frame(otu_table(rarefy_even_depth(ps.ITS.below.ground.fair, sample.size = min(sample_sums(ps.ITS.below.ground.fair)), rngseed = 8675309)))
taxon.ITS.below.ground.fair = as.vector(colnames(data.frame(otu_table(rarefy_even_depth(ps.ITS.below.ground.fair, sample.size = min(sample_sums(ps.ITS.below.ground.fair)), rngseed = 8675309)))))

sta.ITS.below.ground.fair = sncm.fit(spp.ITS.below.ground.fair, taxon.ITS.below.ground.fair, stats=TRUE, pool=NULL)


spp.ITS.below.ground.tall = data.frame(otu_table(rarefy_even_depth(ps.ITS.below.ground.tall, sample.size = min(sample_sums(ps.ITS.below.ground.tall)), rngseed = 8675309)))
taxon.ITS.below.ground.tall = as.vector(colnames(data.frame(otu_table(rarefy_even_depth(ps.ITS.below.ground.tall, sample.size = min(sample_sums(ps.ITS.below.ground.tall)), rngseed = 8675309)))))

sta.ITS.below.ground.tall = sncm.fit(spp.ITS.below.ground.tall, taxon.ITS.below.ground.tall, stats=TRUE, pool=NULL)


spp.ITS.below.ground.nat = data.frame(otu_table(rarefy_even_depth(ps.ITS.below.ground.nat, sample.size = min(sample_sums(ps.ITS.below.ground.nat)), rngseed = 8675309)))
taxon.ITS.below.ground.nat = as.vector(colnames(data.frame(otu_table(rarefy_even_depth(ps.ITS.below.ground.nat, sample.size = min(sample_sums(ps.ITS.below.ground.nat)), rngseed = 8675309)))))

sta.ITS.below.ground.nat = sncm.fit(spp.ITS.below.ground.nat, taxon.ITS.below.ground.nat, stats=TRUE, pool=NULL)


sta.ITS.below.ground.green$Rsqr
sta.ITS.below.ground.fair$Rsqr
sta.ITS.below.ground.tall$Rsqr
sta.ITS.below.ground.nat$Rsqr


df.rsqr.ITS = data.frame('Location' = c('Above','Above','Above','Above','Below','Below','Below','Below'), 
                     'Management' = c("Greens", "Fairway", 'Home Lawn', 'Unmanaged',"Greens", "Fairway", 'Home Lawn', 'Unmanaged'), 
                     'Rsqr' = c(sta.ITS.above.ground.green$Rsqr, sta.ITS.above.ground.fair$Rsqr, sta.ITS.above.ground.tall$Rsqr, sta.ITS.above.ground.nat$Rsqr, 
                                sta.ITS.below.ground.green$Rsqr, sta.ITS.below.ground.fair$Rsqr, sta.ITS.below.ground.tall$Rsqr, sta.ITS.below.ground.nat$Rsqr))
df.rsqr.ITS$type = 'ITS'
df.rsqr.ITS$Management = factor(df.rsqr.ITS$Management, levels = c("Greens", "Fairway", 'Home Lawn', 'Unmanaged'))
df.rsqr.ITS$Location = factor(df.rsqr.ITS$Location, levels = c('Above', 'Below'))

ggplot(df.rsqr.ITS) + geom_point(aes(x = Management, y = Rsqr, color = Location),position = 'jitter', size = 5, alpha = 0.7) + 
  xlab('') + ggsci::scale_color_d3() + ylim(0.25,1) +
  theme_bw(base_size = 18, base_family = "Times")


df.rsqr.ITS$type = 'ITS'
df.rsqr$type = '16s'
df.Rsq = rbind(df.rsqr, df.rsqr.ITS)
ggplot(df.Rsq) + 
  geom_point(aes(x = Management, y = Rsqr, color = Location), size = 5, alpha = 0.9) +
  facet_grid(type~.) +
  xlab('') + ggsci::scale_color_d3() + ylim(0,1) +
  theme_bw(base_size = 18, base_family = "Times") + 
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black'))
ggsave(file = "RsquaredNeutral_Fits.tiff", dpi = 900, width = 10, height = 10, units = 'in')


ggplot(df.rsqr) + geom_point(data = df.rsqr, aes(x = Management, y = Rsqr, color = Location), size = 5, alpha = 0.7)



dev.off()

