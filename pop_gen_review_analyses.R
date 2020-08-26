## Fern Population Genetics Review Analsyes 
## Pelosi, J.A. and E.B. Sessa 

library(ggplot2)
library(dplyr)
library(ggstance)
library(Hmisc)
library(wesanderson)

refs <- read.csv("updated_aug_23_references.csv")

# Number of taxa studied = 156
nrow(refs)

# Taxonomic biases 
tax_refs <- refs %>% 
  group_by(Family) %>% 
  summarise(No._taxa = n())

write.table(tax_refs, "family_biases", sep = ",")

# Overall Averages of P, H, nm for ferns
avg_p <- mean(refs$Percent_Polymorphic_Loci, na.rm = TRUE)
avg_p

avg_h <- mean(refs$Expected_Heterozygosity, na.rm = TRUE)
avg_h

avg_nm <- mean(refs$Calculated_NM, na.rm = TRUE)
avg_nm

avg_F <- mean(refs$Mean_F, na.rm = TRUE)
avg_F

avg_fst <- mean(refs$Mean_Fst, na.rm = TRUE)
avg_fst

# Calculate Average P, H, Fst, and Nm for groups

refs_by_mating_sys <- refs %>% 
  group_by(Mating_system) %>% 
  summarise(averge_P = mean(Percent_Polymorphic_Loci, na.rm = TRUE),
            standard_dev_P = sd(Percent_Polymorphic_Loci, na.rm = TRUE),
            averge_H = mean(Expected_Heterozygosity, na.rm = TRUE),
            standard_dev_H = sd(Expected_Heterozygosity, na.rm = TRUE),
            average_fst = mean(Mean_Fst, na.rm = TRUE),
            standard_dev_fst = sd(Mean_Fst, na.rm = TRUE),
            avearge_nm = mean(Calculated_NM, na.rm = TRUE),
            standard_dev_NM = sd(Calculated_NM, na.rm = TRUE),
            average_F = mean(Mean_F, na.rm = TRUE),
            standard_dev_F= sd(Mean_F, na.rm = TRUE))

refs_by_habit <- refs %>% 
  group_by(Growth_Habit) %>% 
  summarise(averge_P = mean(Percent_Polymorphic_Loci, na.rm = TRUE),
            standard_dev_P = sd(Percent_Polymorphic_Loci, na.rm = TRUE),
            averge_H = mean(Expected_Heterozygosity, na.rm = TRUE),
            standard_dev_H = sd(Expected_Heterozygosity, na.rm = TRUE),
            average_fst = mean(Mean_Fst, na.rm = TRUE),
            standard_dev_fst = sd(Mean_Fst, na.rm = TRUE),
            avearge_nm = mean(Calculated_NM, na.rm = TRUE),
            standard_dev_NM = sd(Calculated_NM, na.rm = TRUE),
            average_F = mean(Mean_F, na.rm = TRUE),
            stadard_dev_F = sd(Mean_F, na.rm= TRUE))

# Plotting P, H, Fst, and F by mating system
mating_system_ref_no_na <- refs %>% 
  filter(Mating_system != "NA")

mating_system_ref_no_na$Mating_system_fact <- factor(mating_system_ref_no_na$Mating_system,
                                                     levels = c("Outcrossing", "Mixed", "Selfing", "Apogamous"))
# Plot Mean F
# No. Outcrossing w/ F = 32
outcrossing_F <- mating_system_ref_no_na %>% 
  filter(Mating_system == "Outcrossing" & Mean_F != "NA") %>% 
  summarise(count = n())

# No. Mixed w/ F = 10
mixed_F <- mating_system_ref_no_na %>% 
  filter(Mating_system == "Mixed" & Mean_F != "NA") %>% 
  summarise(count = n())

# No. Selfing w/ F = 3  
selfing_F <- mating_system_ref_no_na %>% 
  filter(Mating_system == "Selfing" & Mean_F != "NA") %>% 
  summarise(count = n())

# No. Apogamous  w/ Fst = 0  
apogamous_F <- mating_system_ref_no_na %>% 
  filter(Mating_system == "Apogamous" & Mean_F != "NA") %>% 
  summarise(count = n())

# Statistical tests for F 
no_apog_F <- mating_system_ref_no_na %>% 
  filter(Mating_system != "Apogamous") %>% 
  filter(Mean_F != "NA")

kruskal.test(no_apog_F$Mean_Fst, no_apog_F$Mating_system)
# Kruskal-Wallis rank sum test
#
# data:  no_apog_F$Mean_Fst and no_apog_F$Mating_system
# Kruskal-Wallis chi-squared = 7.99, df = 2, p-value = 0.01841

pairwise.wilcox.test(no_apog_F$Mean_F, no_apog_F$Mating_system, p.adjust.method = 'holm')
#             mixed   outcrossing
# outcrossing 0.00223 -          
# selfing     0.02797 0.00092   


# Plot Mean Fst 
# No. Outcrossing w/ Fst = 27
outcrossing_fst <- mating_system_ref_no_na %>% 
  filter(Mating_system == "Outcrossing" & Mean_Fst != "NA") %>% 
  summarise(count = n())

# No. Mixed w/ Fst = 11
mixed_fst <- mating_system_ref_no_na %>% 
  filter(Mating_system == "Mixed" & Mean_Fst != "NA") %>% 
  summarize(count = n())

# No. Selfing w/ Fst = 7  
selfing_fst <- mating_system_ref_no_na %>% 
  filter(Mating_system == "Selfing" & Mean_Fst != "NA") %>% 
  summarise(count = n())

# No. Apogamous  w/ Fst = 0  
apogamous_fst <- mating_system_ref_no_na %>% 
  filter(Mating_system == "Apogamous" & Mean_Fst != "NA") %>% 
  summarise(count = n())

# Plot mean_fst for mating systems 

no_apog_fst <- mating_system_ref_no_na %>% 
  filter(Mating_system != "Apogamous") %>%
  filter(Mating_system != "Heterosporous- selfing") %>% 
  filter(Mean_Fst != "NA")

data_summary <- function(x) {
  m <- mean(x)
  se <- sd(x)/sqrt(length(x))
  return(c(y=m, ymin=m-se, ymax=m+se))
}

ggplot(data = no_apog_fst, mapping = aes(x = Mating_system_fact, y = Mean_Fst, fill = Mating_system)) + 
  geom_violin() +
  scale_fill_manual(values = c("#feb24c", "#ffeda0", "#e31a1c")) + theme_classic() + ylab(expression("Mean "~ F[ST]))+
  xlab("Mating System") + theme(legend.position = "none") + stat_summary(fun.data = data_summary, geom="pointrange")

# ggplot(data = mating_system_ref_no_na, mapping = aes(x = Mean_Fst, y = -0.5)) +
#  geom_boxploth(aes(fill = Mating_system), width = 0.5) +
#  geom_density(aes(x = Mean_Fst, fill = Mating_system), inherit.aes = FALSE) +
#  facet_grid(Mating_system_fact ~., labeller = 
#  labeller(Mating_system_fact = c("outcrossing" = "Outcrossing(n=27)",
#                                  "mixed" = "Mixed(n=11)", "selfing" = "Selfing (n=7)", "apogamous" = "Apogamous(n=0)")),
#  scales = "free_x") +
#  labs(y = "Density") +
#  labs(x = expression(Mean~F[ST])) +
#  scale_fill_brewer(palette = "YlOrRd")+
#  theme_light() +
#  theme(legend.position = "none")
  
ggsave("mating_system_fst_aug_23.png", dpi = 300, width = 5, height = 4.5)


# Statistical tests for Fst  
no_apog_fst <- mating_system_ref_no_na %>% 
  filter(Mating_system != "Apogamous") %>%
  filter(Mating_system != "Heterosporous- selfing") %>% 
  filter(Mean_Fst != "NA")

kruskal.test(no_apog_fst$Mean_Fst, no_apog_fst$Mating_system)
# Wallis rank sum test
#
#data:  no_apog_fst$Mean_Fst and no_apog_fst$Mating_system
#Kruskal-Wallis chi-squared = 16.712, df = 2, p-value = 0.000235

pairwise.wilcox.test(no_apog_fst$Mean_Fst, no_apog_fst$Mating_system, p.adjust.method = "holm")
#                mixed   outcrossing
#  outcrossing 0.00053      -          
#  selfing     0.96122    0.03243   

# Plot H 
# No. Outcrossing w/ H = 27
outcrossing_H <- mating_system_ref_no_na %>% 
  filter(Mating_system == "Outcrossing" & Expected_Heterozygosity != "NA") %>% 
  summarize(count = n())

# No. Mixed w/ H = 10
mixed_H <- mating_system_ref_no_na %>% 
  filter(Mating_system == "mixed" & Expected_Heterozygosity != "NA") %>% 
  summarize(count = n())

# No. Selfing w/ H = 4  
selfing_H <- mating_system_ref_no_na %>% 
  filter(Mating_system == "Selfing" & Expected_Heterozygosity != "NA") %>% 
  summarize(count = n())

# No. Apogamous  w/ H = 1  
apogamous_H <- mating_system_ref_no_na %>% 
  filter(Mating_system == "Apogamous" & Expected_Heterozygosity != "NA") %>% 
  summarize(count = n())

no_apog_he <- mating_system_ref_no_na %>% 
  filter(Mating_system != "Apogamous") %>% 
  filter(Mating_system != "Heterosporous- selfing") %>% 
  filter(Expected_Heterozygosity != "NA")

ggplot(data = no_apog_he, mapping = aes(x = Mating_system_fact, y = Expected_Heterozygosity, fill = Mating_system)) + 
  geom_violin() +
  scale_fill_manual(values = c("#feb24c", "#ffeda0", "#e31a1c")) + theme_classic() + ylab(expression("Mean Expected Heterozygosity ("~H[e]~")"))+
  xlab("Mating System") + theme(legend.position = "none") + stat_summary(fun.data = data_summary, geom="pointrange")


#ggplot(data = mating_system_ref_no_na, mapping = aes(x = Expected_Heterozygosity, y = -0.5)) +
#  geom_boxploth(aes(fill = Mating_system), width = 0.5) +
#  geom_density(aes(x = Expected_Heterozygosity, fill = Mating_system), inherit.aes = FALSE) +
#  facet_grid(Mating_system_fact ~., labeller = 
#               labeller(Mating_system_fact = c("outcrossing" = "Outcrossing(n=27)",
#                                               "mixed" = "Mixed(n=10)", "selfing" = "Selfing (n=4)", "apogamous" = "Apogamous(n=1)")),
#            scales = "free_x") +
#  labs(y = "Density") +
#  labs(x = "Expected Heterozygosity (H)") +
#  scale_fill_brewer(palette = "YlOrRd")+
#  theme_light() +
#  theme(legend.position = "none")

ggsave("mating_system_H_aug_23.png", dpi = 300, width = 5, height = 4.5)

# Statistical tests for Expected Heterozygosity Loci 
no_apog_he <- mating_system_ref_no_na %>% 
  filter(Mating_system != "Apogamous") %>% 
  filter(Mating_system != "Heterosporus- selfing") %>% 
  filter(Expected_Heterozygosity != "NA")

kruskal.test(no_apog_he$Expected_Heterozygosity, no_apog_he$Mating_system)
# Kruskal-Wallis rank sum test
#
# data:  no_apog_he$Expected_Heterozygosity and no_apog_he$Mating_system
# Kruskal-Wallis chi-squared = 2.627, df = 2, p-value = 0.2689


pairwise.wilcox.test(no_apog_he$Expected_Heterozygosity, no_apog_he$Mating_system, p.adjust.method = "holm")
#             mixed outcrossing
# outcrossing 0.57      -          
# selfing     0.48    0.40       


#Plot P

# No. Outcrossing w/ P = 33
outcrossing_P <- mating_system_ref_no_na %>% 
  filter(Mating_system == "outcrossing" & Percent_Polymorphic_Loci != "NA") %>% 
  summarize(count = n())

# No. Mixed w/ P = 11
mixed_P <- mating_system_ref_no_na %>% 
  filter(Mating_system == "mixed" & Percent_Polymorphic_Loci != "NA") %>% 
  summarize(count = n())

# No. Selfing w/ P = 4  
selfing_P <- mating_system_ref_no_na %>% 
  filter(Mating_system == "selfing" & Percent_Polymorphic_Loci != "NA") %>% 
  summarize(count = n())

# No. Apogamous  w/ P = 1  
apogamous_P <- mating_system_ref_no_na %>% 
  filter(Mating_system == "apogamous" & Percent_Polymorphic_Loci != "NA") %>% 
  summarize(count = n())

# Plotting P  w/ regard to mating system 

no_apog_p <- mating_system_ref_no_na %>% 
  filter(Mating_system != "Apogamous") %>% 
  filter(Mating_system != "Heterosporous- selfing") %>% 
  filter(Percent_Polymorphic_Loci != "NA")

ggplot(data = no_apog_p, mapping = aes(x = Mating_system_fact, y = Percent_Polymorphic_Loci, fill = Mating_system)) + 
  geom_violin() +
  scale_fill_manual(values = c("#feb24c", "#ffeda0", "#e31a1c")) + theme_classic() + ylab("Mean Percent Polymorphic Loci (%P)")+
  xlab("Mating System") + theme(legend.position = "none") + stat_summary(fun.data = data_summary, geom="pointrange")


# ggplot(data = mating_system_ref_no_na, mapping = aes(x = Percent_Polymorphic_Loci, y = -0.01)) +
#  geom_boxploth(aes(fill = Mating_system), width = 0.004) +
#  geom_density(aes(x = Percent_Polymorphic_Loci, fill = Mating_system), inherit.aes = FALSE) +
#  facet_grid(Mating_system_fact ~., labeller = 
#               labeller(Mating_system_fact = c("outcrossing" = "Outcrossing(n=33)",
#                                               "mixed" = "Mixed(n=11)", "selfing" = "Selfing (n=4)", "apogamous" = "Apogamous(n=1)")),
#             scales = "free_x") +
#  labs(y = "Density") +
#  labs(x = "Percent Polymorphic Loci (P)") +
#  scale_fill_brewer(palette = "YlOrRd")+
#  theme_light() +
#  theme(legend.position = "none")

ggsave("mating_system_P_aug_23.png", dpi = 300, width = 5, height = 4.5)

# Statistical tests for Percent Polymorphic Loci 
no_apog_p <- mating_system_ref_no_na %>% 
  filter(Mating_system != "Apogamous") %>% 
  filter(Mating_system != "Heterosporous- selfing") %>% 
  filter(Percent_Polymorphic_Loci != "NA")

kruskal.test(no_apog_p$Percent_Polymorphic_Loci, no_apog_p$Mating_system)
# Kruskal-Wallis rank sum test
#
# data:  no_apog_p$Percent_Polymorphic_Loci and no_apog_p$Mating_system
# Kruskal-Wallis chi-squared = 15.278, df = 2, p-value = 0.0004812

pairwise.wilcox.test(no_apog_p$Percent_Polymorphic_Loci, no_apog_p$Mating_system, p.adjust.method = "holm")

#             mixed   outcrossing
#   outcrossing 0.00081 -          
#  selfing     0.48938 0.03119    

# Plot Nm w/regard to mating system
no_apog_nm <- mating_system_ref_no_na %>% 
  filter(Mating_system != "Apogamous") %>% 
  filter(Mating_system != "Heterosporous- selfing")
  filter(Calculated_NM != "NA")

ggplot(data = no_apog_nm, mapping = aes(x = Mating_system_fact, y = Calculated_NM, fill = Mating_system)) + 
  geom_violin(width = 1.33) +
  scale_fill_manual(values = c("#feb24c", "#ffeda0", "#e31a1c")) + theme_classic() + ylab("Mean Number of Migrants per Generation (" ~ italic(Nm)~")")+
  xlab("Mating System") + theme(legend.position = "none") + stat_summary(fun.data = data_summary, geom="pointrange")

ggsave("mating_system_nm_aug_23.png", dpi = 300, width = 5, height = 4.5)

# Statistical test for difference in Nm based on Mating system
no_na_nm_mat <- refs %>% 
  filter(Mating_system != "NA") %>% 
  filter(Mating_system != "Heterosporous- selfing") %>% 
  filter(Calculated_NM != "NA")

summary_nm_mating <- no_na_nm_mat %>% 
  group_by(Mating_system) %>% 
  summarise(count = n())

kruskal.test(no_na_nm_mat$Calculated_NM, no_na_nm_mat$Mating_system)
# Kruskal-Wallis rank sum test
#
# data:  no_na_nm_mat$Calculated_NM and no_na_nm_mat$Mating_system
# Kruskal-Wallis chi-squared = 13.36, df = 2, p-value = 0.001256

pairwise.wilcox.test(no_na_nm_mat$Calculated_NM, no_na_nm_mat$Mating_system, p.adjust.method = "holm")
  
#               mixed  outcrossing
# outcrossing   0.0013    -          
# selfing       0.8591  0.1369   

# ggplot(data = no_na_nm_mat, mapping = aes(x = Mating_system, y = Calculated_NM)) + geom_boxplot() +
#  ylab(expression("Number of Migrants per Generation " ~italic("(Nm)"))) + 
#  xlab("Mating System") + theme_classic() + scale_x_discrete(c("outcrossing" = "Outcrossing (n = 10)",
#                                                               "mixed" = "Mixed (n = 10)", "selfing" = "Selfing (n = 6"))


# Compare genetic diversity and gene flow across growth habits


#Plot calculated NM values by growth habit
filtered_refs_nms <- refs %>% 
  filter(Growth_Habit != "NA" & Calculated_NM != "NA" & Growth_Habit != "Aquatic") 

ggplot(data = filtered_refs_nms, mapping = aes(x = Growth_Habit, y = Calculated_NM, fill = Growth_Habit)) + 
  geom_violin(width = 1.15) +
  scale_fill_manual(values = c("gray45", "green4", "tan4")) + theme_classic() + ylab("Mean Number of Migrants per Generation (" ~ italic(Nm)~")")+
  xlab("Growth Habit") + theme(legend.position = "none") + stat_summary(fun.data = data_summary, geom="pointrange")


#ggplot(data = refs, mapping = aes(x = Growth_Habit, y= Calculated_NM)) + geom_boxplot() +
#  theme_classic() + xlab("Growth Habit") + ylab(expression("Number of Migrants per Generation " ~italic("(Nm)"))) +
#  scale_x_discrete(labels = c("Aquatic" = "Aquatic (n=1)", "Epipetric" = "Epipetric (n=6)",
#                              "Epiphytic" = "Epiphytic (n=9)", "Terrestrial" = "Terrestrial (n=31)")) +
#  theme(axis.text.x = element_text(color = "black")) + theme(axis.text.y = element_text(color = "black"))

ggsave("nm_growth_habit_aug_23.png", dpi = 300, width = 5, height = 4.5)

# Statisical test for difference in Nm w/ respect to growth habit

no_na_nm <- refs %>% 
  filter(Growth_Habit != "NA" & Growth_Habit != "Aquatic") %>% 
  filter(Calculated_NM != "NA")

kruskal.test(no_na_nm$Calculated_NM, no_na_nm$Growth_Habit)
# Kruskal-Wallis rank sum test
#
# data:  no_na_nm$Calculated_NM and no_na_nm$Growth_Habit
# Kruskal-Wallis chi-squared = 10.495, df = 2, p-value = 0.005261

pairwise.wilcox.test(no_na_nm$Calculated_NM, no_na_nm$Growth_Habit, p.adjust.method = "holm")
#                  Epipetric Epiphytic
#       Epiphytic   0.014     -        
#       Terrestrial 0.031     0.045    

ggplot(data = no_na_nm, mapping = aes(x = Growth_Habit, y = Calculated_NM, fill = Growth_Habit)) + 
  geom_violin(width = 1.15) +
  scale_fill_manual(values = c("gray45", "green4", "tan4")) + theme_classic() + ylab("Mean " ~F[ST])+
  xlab("Growth Habit") + theme(legend.position = "none") + stat_summary(fun.data = data_summary, geom="pointrange")

ggsave("growth_habit_nm_aug_23.png", dpi = 300, width = 5, height = 4.5)

# Plot Mean Fst by growth habit
filtered_refs_fst <- refs %>% 
  filter(Growth_Habit != "NA" & Mean_Fst != "NA" & Growth_Habit != "Aquatic") 

ggplot(data = filtered_refs_fst, mapping = aes(x = Growth_Habit, y = Mean_Fst, fill = Growth_Habit)) + 
  geom_violin(width = 1.25) +
  scale_fill_manual(values = c("gray45", "green4", "tan4")) + theme_classic() + ylab("Mean " ~F[ST])+
  xlab("Growth Habit") + theme(legend.position = "none") + stat_summary(fun.data = data_summary, geom="pointrange")

#ggplot(data = refs, mapping = aes(x = Growth_Habit, y= Mean_Fst)) + geom_boxplot() +
#  theme_classic() + xlab("Growth Habit") + ylab(expression(Mean~ F[ST])) +
#  scale_x_discrete(labels = c("Aquatic" = "Aquatic (n=1)", "Epipetric" = "Epipetric (n=8)",
#                              "Epiphytic" = "Epiphytic (n=11)", "Terrestrial" = "Terrestrial (n=39)")) +
#  theme(axis.text.x = element_text(color = "black")) + theme(axis.text.y = element_text(color = "black"))

ggsave("growth_habit_fst_aug_23.png", dpi = 300, width = 5, height = 4.5)

kruskal.test(filtered_refs_fst$Mean_Fst, filtered_refs_fst$Growth_Habit)

# Kruskal-Wallis rank sum test

# data:  refs$Mean_Fst and refs$Growth_Habit
# Kruskal-Wallis chi-squared = 10.821, df = 2, p-value = 0.004469

pairwise.wilcox.test(filtered_refs_fst$Mean_Fst, filtered_refs_fst$Growth_Habit, p.adjust.method = "holm")
#                 Epipetric Epiphytic
#     Epiphytic   0.0076    -        
#     Terrestrial 0.0212    0.0786   

# statistical analysis of F, fixation index w re: growth habit 
filtered_refs_F <- refs %>% 
  filter(Growth_Habit != "NA" & Mean_F != "NA" & Growth_Habit != "Aquatic") 

kruskal.test(filtered_refs_F$Mean_F, filtered_refs_F$Growth_Habit)
# Kruskal-Wallis rank sum test
#
# data:  refs$Mean_F and refs$Growth_Habit
# Kruskal-Wallis chi-squared = 4.7294, df = 2, p-value = 0.09398


# Plot Mean He by growth habit
filtered_refs_he <- refs %>% 
  filter(Growth_Habit != "NA" & Expected_Heterozygosity != "NA" & Growth_Habit != "Aquatic") 

ggplot(data = filtered_refs_he, mapping = aes(x = Growth_Habit, y = Expected_Heterozygosity, fill = Growth_Habit)) + 
  geom_violin(width = 1.25) +
  scale_fill_manual(values = c("gray45", "green4", "tan4")) + theme_classic() + ylab(expression("Mean Expected Heterozygosity (" ~H[e]~")"))+
  xlab("Growth Habit") + theme(legend.position = "none") + stat_summary(fun.data = data_summary, geom="pointrange")


ggsave("growth_habit_he_aug_23.png", dpi = 300, width = 5, height = 4.5)

kruskal.test(filtered_refs_he$Expected_Heterozygosity, filtered_refs_he$Growth_Habit)

#Kruskal-Wallis rank sum test
#
#data:  refs$Expected_Heterozygosity and refs$Growth_Habit
#Kruskal-Wallis chi-squared = 0.6154, df = 2, p-value = 0.7355_

# Plot Mean P by growth habit
filtered_refs_p <- refs %>% 
  filter(Growth_Habit != "NA" & Percent_Polymorphic_Loci != "NA" & Growth_Habit != "Aquatic") 

ggplot(data = filtered_refs_p, mapping = aes(x = Growth_Habit, y = Percent_Polymorphic_Loci, fill = Growth_Habit)) + 
  geom_violin(width = 1.25) +
  scale_fill_manual(values = c("gray45", "green4", "tan4")) + theme_classic() + 
  ylab("Mean Percent Polymorphic Loic (%P)")+
  xlab("Growth Habit") + theme(legend.position = "none") + stat_summary(fun.data = data_summary, geom="pointrange")

ggsave("growth_habit_p_aug_23.png", dpi = 300, width = 5, height = 4.5)

# Statistical tests for growth habit and percent polymorphic loci
kruskal.test(filtered_refs_p$Percent_Polymorphic_Loci, filtered_refs_p$Growth_Habit)

#Kruskal-Wallis rank sum test

#data:  refs$Percent_Polymorphic_Loci and refs$Growth_Habit
#Kruskal-Wallis chi-squared = 11.94, df = 2, p-value = 0.002555

pairwise.wilcox.test(filtered_refs_p$Percent_Polymorphic_Loci, filtered_refs_p$Growth_Habit, p.adjust.method = "holm")

#Pairwise comparisons using Wilcoxon rank sum test 
#
#               Epipetric Epiphytic
#     Epiphytic   0.0090    -        
#    Terrestrial  0.2789    0.0045  


# Plot types of markers used across studies 
ref_markers <- refs %>% 
  group_by(Molecular_Marker) %>% 
  summarize(count = n())
# using output from ref_markers we make the new dataframe
Molecular_Marker <- c("AFLP", "AFLP", 
                      "Microsatellites (SSRs)", "Microsatellites (SSRs)",
                      "cp DNA", "cp DNA", 
                      "ISSRs", "ISSRs", 
                      "RADSeq", "RAPD", "RAPD",
                      "Allozymes","Allozymes","Allozymes","Allozymes")
Number_of_Taxa <- c(2, 6, 6, 13, 3, 2, 5, 1, 4, 2, 6, 18, 75, 17, 12)
Year <- c("2000-2009", "2010-2019", "2000-2009", "2010-2019", "2000-2009", "2010-2019", "2000-2009", "2010-2019", "2010-2019",
          "1990-1999", "2000-2009", "1980-1989", "1990-1999", "2000-2009", "2010-2019")
Molecular_Markers_df <- data.frame(Molecular_Marker, Number_of_Taxa, Year)

Molecular_Markers_df$Molecular_Marker <- factor(Molecular_Markers_df$Molecular_Marker, 
                                                levels = c("Allozymes","Microsatellites (SSRs)","RAPD","AFLP", "ISSRs","cp DNA","RADSeq"))

ggplot(data = Molecular_Markers_df, mapping = aes(x = Molecular_Marker, y = Number_of_Taxa, fill = Year)) +
  geom_bar(stat = "identity") + xlab("Molecular Marker") + ylab("Number of Taxa") +
  theme_classic() + theme(axis.text.x = element_text(color ="black")) + theme(axis.text.y = element_text(color = "black"))+
  scale_fill_manual(values = c("#8DD3C7", "#B3BADA", "#FB8072", "#80B1D3"))

ggsave("molecular_markers_with_year.png", dpi = 300, width = 8.5, height = 4)

ggplot(data = Molecular_Markers_df, mapping = aes(x = Year, y = Number_of_Taxa, fill = Molecular_Marker)) +
  geom_bar(stat = "identity") + xlab("Year") + ylab("Number of Taxa") +
  theme_classic() + theme(axis.text.x = element_text(color ="black")) + theme(axis.text.y = element_text(color = "black"))+
  scale_fill_brewer(palette = "Set2")
 
library(RColorBrewer)
brewer.pal(n=5, name = "Set3")
