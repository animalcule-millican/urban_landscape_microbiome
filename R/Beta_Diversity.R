ps.16s.norm = transform_sample_counts(ps.16s, function(x) x / sum(x))
ps.16s.norm.nmds = ordinate(ps.16s.norm, method = 'NMDS', distance = 'bray', k = 3, maxtry = 1000)
plot_ordination(ps.16s.norm, ps.16s.norm.nmds, color = "Location", shape = 'Section') + 
  geom_point(size = 6) + ggsci::scale_color_d3() +
  theme_bw(base_size = 18, base_family = "Times")
ggsave("Beta_Div_16s.tiff", dpi = 900, width = 12, height = 7.5, units = 'in')


ps.ITS.norm = transform_sample_counts(ps.ITS, function(x) x / sum(x))
ps.ITS.norm.nmds = ordinate(ps.ITS.norm, method = 'NMDS', distance = 'bray', k = 3, maxtry = 1000)
plot_ordination(ps.ITS.norm, ps.ITS.norm.nmds, color = 'Location', shape = 'Section') + 
  geom_point(size = 6) + ggsci::scale_color_d3() +
  theme_bw(base_size = 18, base_family = "Times")
ggsave("Beta_Div_ITS.tiff", dpi = 900, width = 12, height = 7.5, units = 'in')