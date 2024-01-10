ps.16s.occ.abun = Occ.Abund(ps.16s, seed = 8675309, depth = min(sample_sums(ps.16s)))
ps.16s.above.ground.occ.abun = Occ.Abund(ps.16s.above.ground, seed = 8675309, depth = min(sample_sums(ps.16s.above.ground)))
ps.16s.below.ground.occ.abun = Occ.Abund(ps.16s.below.ground, seed = 8675309, depth = min(sample_sums(ps.16s.below.ground)))

ps.16s.occ.abun = sloan.neutral(ps.16s)
ps.16s.above.ground.occ.abun = sloan.neutral(ps.16s.above.ground)
ps.16s.below.ground.occ.abun = sloan.neutral(ps.16s.below.ground)
ps.16s.occ.abun$sta$Rsqr
ps.16s.above.ground.occ.abun$sta$Rsqr
ps.16s.below.ground.occ.abun$sta$Rsqr

sum(ps.16s.occ.abun$obs$freq > (ps.16s.occ.abun$obs$pred.upr), na.rm=TRUE)/ps.16s.occ.abun$sta$Richness
sum(ps.16s.occ.abun$obs$freq < (ps.16s.occ.abun$obs$pred.lwr), na.rm=TRUE)/ps.16s.occ.abun$sta$Richness

sum(ps.ITS.occ.abun$obs$freq > (ps.ITS.occ.abun$obs$pred.upr), na.rm=TRUE)/ps.ITS.occ.abun$sta$Richness
sum(ps.ITS.occ.abun$obs$freq < (ps.ITS.occ.abun$obs$pred.lwr), na.rm=TRUE)/ps.ITS.occ.abun$sta$Richness

sum(ps.16s.occ.abun$obs$freq > .50, na.rm=TRUE)/ps.16s.occ.abun$sta$Richness
sum(ps.ITS.occ.abun$obs$freq > .50, na.rm=TRUE)/ps.ITS.occ.abun$sta$Richness

ggplot() +
  geom_point(data = ps.16s.occ.abun$occ.abun, aes(x = log10(abund), y = occ), pch = 21, fill = '#37AEC3', size = 5, alpha=.4)+
  geom_line(color = '#d60000', data = ps.16s.occ.abun$obs, size=1, aes(y = freq.pred, x = log10(p)), alpha = .75) +
  geom_line(color = '#d60000', data = ps.16s.occ.abun$obs, lty = 'twodash', size=1, aes(y = pred.upr, x = log10(p)), alpha = .75)+
  geom_line(color = '#d60000', data = ps.16s.occ.abun$obs, lty = 'twodash', size=1, aes(y = pred.lwr, x = log10(p)), alpha = .75)+
  labs(x = expression(Log[10](mean~relative~abundance)), y = "Occupancy") + 
  theme_bw(base_size = 28, base_family = "Times")
ggsave(file = "Occ_Abund_metacommunity_16s.tiff", dpi = 900, width = 12, height = 12, units = 'in')


ggplot() +
  geom_point(data = ps.16s.above.ground.occ.abun$occ.abun, aes(x = log10(abund), y = occ), pch = 21, fill = '#37AEC3', size = 5, alpha=.6)+
  geom_line(color = '#d60000', data = ps.16s.above.ground.occ.abun$obs, size=1, aes(y = freq.pred, x = log10(p)), alpha = .75) +
  geom_line(color = '#d60000', data = ps.16s.above.ground.occ.abun$obs, lty = 'twodash', size=1, aes(y = pred.upr, x = log10(p)), alpha = .75)+
  geom_line(color = '#d60000', data = ps.16s.above.ground.occ.abun$obs, lty = 'twodash', size=1, aes(y = pred.lwr, x = log10(p)), alpha = .75)+
  labs(x = expression(Log[10](mean~relative~abundance)), y = "Occupancy") + 
  theme_bw(base_size = 28, base_family = "Times")
ggsave(file = "Occ_Abund_Above_Ground_16s.tiff", dpi = 900, width = 12, height = 12, units = 'in')

ggplot() +
  geom_point(data = ps.16s.below.ground.occ.abun$occ.abun, aes(x = log10(abund), y = occ), pch = 21, fill = '#37AEC3', size = 5, alpha=.6)+
  geom_line(color = '#d60000', data = ps.16s.below.ground.occ.abun$obs, size=1, aes(y = freq.pred, x = log10(p)), alpha = .75) +
  geom_line(color = '#d60000', data = ps.16s.below.ground.occ.abun$obs, lty = 'twodash', size=1, aes(y = pred.upr, x = log10(p)), alpha = .75)+
  geom_line(color = '#d60000', data = ps.16s.below.ground.occ.abun$obs, lty = 'twodash', size=1, aes(y = pred.lwr, x = log10(p)), alpha = .75)+
  labs(x = expression(Log[10](mean~relative~abundance)), y = "Occupancy") + 
  theme_bw(base_size = 28, base_family = "Times")
ggsave(file = "Occ_Abund_Below_Ground_16s.tiff", dpi = 900, width = 12, height = 12, units = 'in')



ps.
ps.ITS.occ.abun = sloan.neutral(ps.ITS)
ps.ITS.above.ground.occ.abun = sloan.neutral(ps.ITS.above.ground)
ps.ITS.below.ground.occ.abun = sloan.neutral(ps.ITS.below.ground)
ps.ITS.occ.abun$sta$Rsqr
ps.ITS.above.ground.occ.abun$sta$Rsqr
ps.ITS.below.ground.occ.abun$sta$Rsqr

ggplot() +
  geom_point(data = ps.ITS.occ.abun$occ.abun, aes(x = log10(abund), y = occ), pch = 21, fill = '#37AEC3', size = 5, alpha=.4)+
  geom_line(color = '#d60000', data = ps.ITS.occ.abun$obs, size=1, aes(y = freq.pred, x = log10(p)), alpha = .75) +
  geom_line(color = '#d60000', data = ps.ITS.occ.abun$obs, lty = 'twodash', size=1, aes(y = pred.upr, x = log10(p)), alpha = .75)+
  geom_line(color = '#d60000', data = ps.ITS.occ.abun$obs, lty = 'twodash', size=1, aes(y = pred.lwr, x = log10(p)), alpha = .75)+
  labs(x = expression(Log[10](mean~relative~abundance)), y = "Occupancy") + 
  theme_bw(base_size = 28, base_family = "Times")
ggsave(file = "Occ_Abund_metacommunity_ITS.tiff", dpi = 900, width = 12, height = 12, units = 'in')


ggplot() +
  geom_point(data = ps.ITS.above.ground.occ.abun$occ.abun, aes(x = log10(abund), y = occ), pch = 21, fill = '#37AEC3', size = 5, alpha=.6)+
  geom_line(color = '#d60000', data = ps.ITS.above.ground.occ.abun$obs, size=1, aes(y = freq.pred, x = log10(p)), alpha = .75) +
  geom_line(color = '#d60000', data = ps.ITS.above.ground.occ.abun$obs, lty = 'twodash', size=1, aes(y = pred.upr, x = log10(p)), alpha = .75)+
  geom_line(color = '#d60000', data = ps.ITS.above.ground.occ.abun$obs, lty = 'twodash', size=1, aes(y = pred.lwr, x = log10(p)), alpha = .75)+
  labs(x = expression(Log[10](mean~relative~abundance)), y = "Occupancy") + 
  theme_bw(base_size = 28, base_family = "Times")
ggsave(file = "Occ_Abund_Above_Ground_ITS.tiff", dpi = 900, width = 12, height = 12, units = 'in')

ggplot() +
  geom_point(data = ps.ITS.below.ground.occ.abun$occ.abun, aes(x = log10(abund), y = occ), pch = 21, fill = '#37AEC3', size = 5, alpha=.6)+
  geom_line(color = '#d60000', data = ps.ITS.below.ground.occ.abun$obs, size=1, aes(y = freq.pred, x = log10(p)), alpha = .75) +
  geom_line(color = '#d60000', data = ps.ITS.below.ground.occ.abun$obs, lty = 'twodash', size=1, aes(y = pred.upr, x = log10(p)), alpha = .75)+
  geom_line(color = '#d60000', data = ps.ITS.below.ground.occ.abun$obs, lty = 'twodash', size=1, aes(y = pred.lwr, x = log10(p)), alpha = .75)+
  labs(x = expression(Log[10](mean~relative~abundance)), y = "Occupancy") + 
  theme_bw(base_size = 28, base_family = "Times")
ggsave(file = "Occ_Abund_Below_Ground_ITS.tiff", dpi = 900, width = 12, height = 12, units = 'in')

