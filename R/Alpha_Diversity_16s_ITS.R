library('ggplot2')
library('ggsci')
library('agricolae')
library('Rmisc')
library('tidyverse')
library('lme4')
library('multcomp')
library('emmeans')
library('car')
library('janitor')
set.seed(8675309)
load(file = "/Users/millican/Documents/Research/R_functions/microbiome/Alpha_Diversity.RData")


ps.16s.rare = alpha.div(ps.16s)

ps.ITS.rare = alpha.div(ps.ITS)

sam.dat.16s = data.frame(sample_data(ps.16s.rare))
sam.dat.ITS = data.frame(sample_data(ps.ITS.rare))

shannon.16s.stats = sam.dat.16s %>% 
  nest(-Section) %>% mutate(Section = factor(Section, levels = c("Leaf", "Thatch", "Rhizosphere", "Bulk"),
                                       labels = c("Leaf", "Thatch", "Rhizosphere", "Bulk")),
                         glm = map(data, ~ glm(Shannon ~ Location, data = ., family = quasipoisson)),
                         anova = map(data, ~ Anova(glm(Shannon ~ Location, data = ., family = quasipoisson))),
                         emmeans = map(data, ~ emmeans(glm(Shannon ~ Location, data = ., family = quasipoisson),  ~  Location, type='response')),
                         stat.grp = map(data, ~ multcomp::cld(emmeans(glm(Shannon ~ Location, data = ., family = quasipoisson),  ~  Location, type='response'), Letter="ABCDEFGHIJKLMNOPQRSTUVWZYZ")), 
                         disease.CI = map(data, ~as.data.frame(summary(emmeans(glm(Shannon ~ Location, data = ., family = quasipoisson),  ~  Location, type='response'))))) %>%
  unnest(c(stat.grp))
shannon.16s.stats$.group = gsub(" ","", shannon.16s.stats$.group)

ggplot(shannon.16s.stats) + 
  geom_errorbar(aes(x=Section, ymin=asymp.LCL, ymax=asymp.UCL, 
                    color=Location), width = 0, alpha = 0.6,
                position = position_dodge(width = 0.8), 
                size = 3, show.legend = FALSE) + 
  geom_point(aes(x = Section, 
                 y = rate, 
                 color = Location), 
             position = position_dodge(width = 0.8), size = 6, alpha = 0.99) +
  geom_text(aes(label = .group, x = Section, 
                group = Location, y = asymp.UCL+(0.007*max(asymp.UCL))), 
          position=position_dodge(width=0.8), vjust=-0.5, size = 7) +
  xlab('') + ylab("Shannon Diversity Index") +
  ggsci::scale_color_jco() +
  ylim(NA, max(shannon.16s.stats$asymp.UCL+
             (0.04*max(shannon.16s.stats$asymp.UCL)))) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), lty = 'twodash', alpha = 0.4) + 
  theme_bw(base_size = 22, base_family = "Times") + 
  theme(panel.background = element_rect(fill = '#FAFAFA')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Shannon_16s.tiff", dpi = 900, width = 12, height = 7.5, units = 'in')

Observed.16s.stats = sam.dat.16s %>% 
  nest(-Section) %>% mutate(Section = factor(Section, levels = c("Leaf", "Thatch", "Rhizosphere", "Bulk"),
                                             labels = c("Leaf", "Thatch", "Rhizosphere", "Bulk")),
                            glm = map(data, ~ glm(Observed ~ Location, data = ., family = quasipoisson)),
                            anova = map(data, ~ Anova(glm(Observed ~ Location, data = ., family = quasipoisson))),
                            emmeans = map(data, ~ emmeans(glm(Observed ~ Location, data = ., family = quasipoisson),  ~  Location, type='response')),
                            stat.grp = map(data, ~ multcomp::cld(emmeans(glm(Observed ~ Location, data = ., family = quasipoisson),  ~  Location, type='response'), Letter="ABCDEFGHIJKLMNOPQRSTUVWZYZ")), 
                            disease.CI = map(data, ~as.data.frame(summary(emmeans(glm(Observed ~ Location, data = ., family = quasipoisson),  ~  Location, type='response'))))) %>%
  unnest(c(stat.grp))
Observed.16s.stats$.group = gsub(" ","", Observed.16s.stats$.group)

ggplot(Observed.16s.stats) + 
  geom_errorbar(aes(x=Section, ymin=asymp.LCL, ymax=asymp.UCL, 
                    color=Location), width = 0, alpha = 0.6,
                position = position_dodge(width = 0.8), 
                size = 3, show.legend = FALSE) + 
  geom_point(aes(x = Section, 
                 y = rate, 
                 color = Location), 
             position = position_dodge(width = 0.8), size = 6, alpha = 0.99) +
  geom_text(aes(label = .group, x = Section, 
                group = Location, y = asymp.UCL+(0.007*max(asymp.UCL))), 
            position=position_dodge(width=0.8), vjust=-0.5, size = 7) +
  xlab('') + ylab("Observed Richness") +
  ggsci::scale_color_jco() +
  ylim(NA, max(Observed.16s.stats$asymp.UCL+
                 (0.04*max(Observed.16s.stats$asymp.UCL)))) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), lty = 'twodash', alpha = 0.4) + 
  theme_bw(base_size = 22, base_family = "Times") + 
  theme(panel.background = element_rect(fill = '#FAFAFA')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Observed_16s.tiff", dpi = 900, width = 12, height = 7.5, units = 'in')

###############################################################################################################
### FUNGAL ALPHA DIVERSITY ####################################################################################
###############################################################################################################

shannon.ITS.stats = sam.dat.ITS %>% 
  nest(-Section) %>% mutate(Section = factor(Section, levels = c("Leaf", "Thatch", "Rhizosphere", "Bulk"),
                                             labels = c("Leaf", "Thatch", "Rhizosphere", "Bulk")),
                            glm = map(data, ~ glm(Shannon ~ Location, data = ., family = quasipoisson)),
                            anova = map(data, ~ Anova(glm(Shannon ~ Location, data = ., family = quasipoisson))),
                            emmeans = map(data, ~ emmeans(glm(Shannon ~ Location, data = ., family = quasipoisson),  ~  Location, type='response')),
                            stat.grp = map(data, ~ multcomp::cld(emmeans(glm(Shannon ~ Location, data = ., family = quasipoisson),  ~  Location, type='response'), Letter="ABCDEFGHIJKLMNOPQRSTUVWZYZ")), 
                            disease.CI = map(data, ~as.data.frame(summary(emmeans(glm(Shannon ~ Location, data = ., family = quasipoisson),  ~  Location, type='response'))))) %>%
  unnest(c(stat.grp))
shannon.ITS.stats$.group = gsub(" ","", shannon.ITS.stats$.group)

ggplot(shannon.ITS.stats) + 
  geom_errorbar(aes(x=Section, ymin=asymp.LCL, ymax=asymp.UCL, 
                    color=Location), width = 0, alpha = 0.6,
                position = position_dodge(width = 0.8), 
                size = 3, show.legend = FALSE) + 
  geom_point(aes(x = Section, 
                 y = rate, 
                 color = Location), 
             position = position_dodge(width = 0.8), size = 6, alpha = 0.99) +
  geom_text(aes(label = .group, x = Section, 
                group = Location, y = asymp.UCL+(0.007*max(asymp.UCL))), 
            position=position_dodge(width=0.8), vjust=-0.5, size = 7) +
  xlab('') + ylab("Shannon Diversity Index") +
  ggsci::scale_color_jco() +
  ylim(NA, max(shannon.ITS.stats$asymp.UCL+
                 (0.04*max(shannon.ITS.stats$asymp.UCL)))) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), lty = 'twodash', alpha = 0.4) + 
  theme_bw(base_size = 22, base_family = "Times") + 
  theme(panel.background = element_rect(fill = '#FAFAFA')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Shannon_ITS.tiff", dpi = 900, width = 12, height = 7.5, units = 'in')

Observed.ITS.stats = sam.dat.ITS %>% 
  nest(-Section) %>% mutate(Section = factor(Section, levels = c("Leaf", "Thatch", "Rhizosphere", "Bulk"),
                                             labels = c("Leaf", "Thatch", "Rhizosphere", "Bulk")),
                            glm = map(data, ~ glm(Observed ~ Location, data = ., family = quasipoisson)),
                            anova = map(data, ~ Anova(glm(Observed ~ Location, data = ., family = quasipoisson))),
                            emmeans = map(data, ~ emmeans(glm(Observed ~ Location, data = ., family = quasipoisson),  ~  Location, type='response')),
                            stat.grp = map(data, ~ multcomp::cld(emmeans(glm(Observed ~ Location, data = ., family = quasipoisson),  ~  Location, type='response'), Letter="ABCDEFGHIJKLMNOPQRSTUVWZYZ")), 
                            disease.CI = map(data, ~as.data.frame(summary(emmeans(glm(Observed ~ Location, data = ., family = quasipoisson),  ~  Location, type='response'))))) %>%
  unnest(c(stat.grp))
Observed.ITS.stats$.group = gsub(" ","", Observed.ITS.stats$.group)

ggplot(Observed.ITS.stats) + 
  geom_errorbar(aes(x=Section, ymin=asymp.LCL, ymax=asymp.UCL, 
                    color=Location), width = 0, alpha = 0.6,
                position = position_dodge(width = 0.8), 
                size = 3, show.legend = FALSE) + 
  geom_point(aes(x = Section, 
                 y = rate, 
                 color = Location), 
             position = position_dodge(width = 0.8), size = 6, alpha = 0.99) +
  geom_text(aes(label = .group, x = Section, 
                group = Location, y = asymp.UCL+(0.007*max(asymp.UCL))), 
            position=position_dodge(width=0.8), vjust=-0.5, size = 7) +
  xlab('') + ylab("Observed Richness") +
  ggsci::scale_color_jco() +
  ylim(NA, max(Observed.ITS.stats$asymp.UCL+
                 (0.04*max(Observed.ITS.stats$asymp.UCL)))) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), lty = 'twodash', alpha = 0.4) + 
  theme_bw(base_size = 22, base_family = "Times") + 
  theme(panel.background = element_rect(fill = '#FAFAFA')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Observed_ITS.tiff", dpi = 900, width = 12, height = 7.5, units = 'in')


