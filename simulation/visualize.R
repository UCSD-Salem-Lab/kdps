library(dplyr)
library(readxl)
library(ggplot2)
library(ggpubr)

fuzziness = read_xlsx("fuzziness_result_summary.xlsx")
phenotype = read_xlsx("phenotype_result_summary.xlsx")
relationship = read_xlsx("relationship_result_summary.xlsx")

fuzziness = filter(fuzziness, phenotypic_naive == FALSE)
fuzziness$fuzziness = as.numeric(fuzziness$fuzziness)

df = fuzziness %>%
  mutate(ratio = n_remove / total_n) %>%
  group_by(fuzziness) %>%
  summarize(average_ratio = mean(ratio),
            sd = sd(ratio), 
            n = n()) %>%
  mutate(sem = sd / sqrt(n)) 

plt1 = ggplot(df, aes(x = fuzziness, y = average_ratio)) +
  scale_x_continuous(breaks = c(0, 1, 2, 5, 10)) + 
  geom_line(linetype = "dashed") +
  geom_errorbar(aes(ymin = average_ratio - sem, ymax = average_ratio + sem), width = 0.5) +
  geom_point(size = 3, color = "black", fill = "gray", shape = 21) +
  labs(x = "Fuzziness",
       y = "Removal Ratio") + 
  theme_pubr() + 
  theme(text = element_text(size = 15)); plt1

df = fuzziness %>%
  group_by(fuzziness) %>%
  summarize(average_time = mean(comp_time),
            sd = sd(comp_time), 
            n = n()) %>%
  mutate(sem = sd / sqrt(n)) 

plt2 = ggplot(data = df, aes(x = fuzziness, y = average_time)) + 
  scale_x_continuous(breaks = c(0, 1, 2, 5, 10)) + 
  geom_line(linetype = "dashed") +
  geom_errorbar(aes(ymin = average_time - sem, ymax = average_time + sem), width = 0.5) +
  geom_point(size = 3, color = "black", fill = "gray", shape = 21) +
  labs(x = "Fuzziness",
       y = "Computational Time (s)") + 
  theme_pubr() + 
  theme(text = element_text(size = 15)); plt2

# ggplot(data = fuzziness, aes(x = fuzziness)) + 
#   geom_point(aes(y = n_remove/total_n))


df = relationship %>%
  group_by(relationship) %>%
  summarize(average_time = mean(comp_time),
            sd = sd(comp_time), 
            n = n()) %>%
  mutate(sem = sd / sqrt(n)) %>%
  mutate(relationship = as.numeric(relationship))

plt3 = ggplot(data = df, aes(x = relationship, y = average_time)) + 
  geom_line(linetype = "dashed") + 
  # geom_errorbar(aes(ymin = average_time - sem, ymax = average_time + sem), width = 0.5) +
  scale_y_log10() +
  geom_point(size = 3, color = "black", fill = "gray", shape = 21) +
  labs(x = "Number of Relatedness",
       y = "Computational Time (s)") + 
  theme_pubr() + 
  theme(text = element_text(size = 15),
        plot.margin = margin(0.1,1,0.1,0.1, "cm")); plt3

ggsave(filename = "ratio_to_fuzziness.png", device = "png", plot = plt1, dpi = 1200,
       width = 5, height = 4.5)

ggsave(filename = "time_to_fuzziness.png", device = "png", plot = plt2, dpi = 1200,
       width = 5, height = 4.5)

ggsave(filename = "ralationship_to_time.png", device = "png", plot = plt3, dpi = 1200,
       width = 5, height = 4.5)




  