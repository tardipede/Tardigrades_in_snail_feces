library(ggplot2)
library(patchwork)
library(BioStatR)
library(lme4)
library(boot) 
set.seed(123456789)


### 1) Analysis on laboratory experiment data
# Alive tardigrades fed on snails

# Snail.ID: id code of the snail
# alive: % of alive tardigrades defecated over tardigrades ingested
# alive_n: number of alive tardigrades recovered from feces
# tot_eaten: number of tardigrades eaten by that snail (by snail ID)
# day: day when feces were examined

data_exp = read.table("survival.txt",header=T,dec=".",sep="\t")

data_glm$day = factor(data_glm$day, levels = c("d1","d2","d3","d4"))

# Create a contrast matrix to have d1 vs. d2, d2 vs. d3, d3 vs.d4 with BACKWARD DIFFERENCE Coding
# from: https://stats.idre.ucla.edu/r/library/r-library-contrast-coding-systems-for-categorical-variables/#backward
my_contrasts = matrix(c(-3/4, 1/4, 1/4, 1/4, -1/2, -1/2, 1/2, 1/2, 
                            -1/4, -1/4, -1/4, 3/4), ncol = 3)

contrasts(data_glm$day) = my_contrasts

# GLMM with snail ID as random effect
mod1 = glmer(cbind(alive_n,tot_eaten-alive_n) ~ day + (1|Snail.ID) , family ="binomial",data = data_glm)

summary(mod1)
#Note: with the specified contrasts the slopes have the following meanings:
# day1: difference between day1 and day2
# day2: difference between day2 and day3
# day3: difference between day3 and day4

# Extract the % of tardigrades defecated alive each day
aggregate(data_exp$alive, by = list(data_exp$day), FUN = mean)


### 2a) Wild snails data plot
# Create data frame for plotting with proportion of success and 95%CI

# no T: 20/28 snails
# Only alive T: 5/28 snails
# Only dead T: 1 snail
# Alive+ dead T: 2/28 snails


data_wild = data.frame(binom.ci(c(20,5,1,2),28),
                       category = factor(c("no_T","Alive","Dead","Alive_Dead"), levels = c("no_T","Alive_Dead","Alive","Dead")))

data_wild$category = factor(data_wild$category, levels = c("no_T","Alive","Alive_Dead","Dead"))

plot_wild = ggplot(data_wild)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")+
  geom_bar(stat="identity",aes(x=category,y=PointEst*100,fill=category),alpha=0.75,color="black")+
  geom_errorbar(aes(ymin = Lower*100, ymax = Upper*100,x=category), width = 0.5)+
  xlab("")+ylab("% of wild snails ±95% C.I.")+
  ylim(c(0,100))+
  scale_fill_manual(values = c("#e66101","#fdb863","#b2abd2","#5e3c99"))+
  scale_x_discrete(labels = c("Tardigrades \n absent \n n=20/28",
                              "Tardigrades \n present \n(Alive)\n n=5/28",
                              "Tardigrades \n present \n(Alive+Dead)\n n=2/28",
                              "Tardigrades \n present \n(Dead)\n n=1/28"))+
  ggtitle("A - Tardigrades in feces\nfrom wild snails")


# Define a simple function to convert pvalues to significant code for the following plot
signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
         symbols = c("***", "**", "*", ".", " "))
}

### 2b) Experiment data plot
plot_experiment = ggplot(data_exp)+
  theme_bw()+
  ylim(c(0,100))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")+
  geom_path(aes(x=day,y=alive*100,group=Snail.ID),alpha=0.25,color="grey")+
  geom_boxplot(aes(x=day,y=alive*100),fill=NA,color="black",outlier.shape=NA)+
  geom_dotplot(aes(x=day,y=alive*100,fill=day),binaxis = "y", stackdir = "center",dotsize=0.5,color=NA,binwidth = 2,alpha=0.75)+
  scale_x_discrete(labels = c("Day 1","Day 2","Day 3","Day 4","Cumulative"))+
  scale_fill_manual(values = c("#e66101","#fdb863","#b2abd2","#5e3c99"))+
  xlab("")+ylab("% of ingested tardigrades defecated alive")+
  geom_segment(aes(y = 95, yend = 95, x = 1, xend = 2))+ # d1 vs d2
  geom_segment(aes(y = 85, yend = 85, x = 2, xend = 3))+ # d2 vs d3
  geom_segment(aes(y = 75, yend = 75, x = 3, xend = 4))+ # d3 vs d4
  annotate("text",x = 1.5, y = 97, label = signif.num(summary(mod1)$coefficients[2,4]))+
  annotate("text",x = 2.5, y = 87, label= signif.num(summary(mod1)$coefficients[3,4]))+
  annotate("text",x = 3.5, y = 77, label= signif.num(summary(mod1)$coefficients[4,4]))+
  geom_segment(data = data.frame(x = c(1,2,2,3,3,4),    # this ads the small vertical bars in the significance comparisons segments
                                 xend = c(1,2,2,3,3,4),
                                 y = c(95,95,85,85,75,75),
                                 yend = c(93,93,83,83,73,73)),
               aes(x = x, xend = xend, y = y, yend = yend))+
  ggtitle("B - Tardigrades survival\n after gut passage")




### 2c) Cultured animals plot

# 13/16 trials were successfull 

# Create data frame for plotting with proportion of success and 95%CI
data_reproduction = data.frame(binom.ci(13,16))


plot_reproduction = ggplot(data_reproduction)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")+
  geom_bar(stat="identity",aes(y = PointEst*100, x = "Successfull_trials"),alpha=0.75,color="black", fill="#fdb863")+
  geom_errorbar(aes(ymin = Lower*100, ymax = Upper*100, x = "Successfull_trials"), width = 0.5)+
  ylim(c(0,100))+
  xlab("")+ylab("% of success in reproduction trials ±95% C.I.")+
  scale_x_discrete(labels = c("Successfull \n trials (n = 13/16)"))+
  ggtitle("C - Tardigrades reproduction\n after gut passage")


## Now put the plots together
plot_final = (plot_wild|plot_experiment|plot_reproduction) + plot_layout(widths = c(2.5,2.5,1))
ggsave(plot_final,filename = "plot.jpg",width = 40,height = 20, units = "cm", dpi = 300)
ggsave(plot_final,filename = "plot.svg",width = 40,height = 20, units = "cm", dpi = 300)

