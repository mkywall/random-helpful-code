

#---- run these if you don't have ggplot2 and/or reshape installed

install.packages("ggplot2")
install.packages("reshape")

#---- libraries

library(ggplot2)
library(reshape)


#---- if you are starting with a matrix do this
melt_df<- melt(your_matrix, id = c("the-columns-you-want-to-", "keep-as-columns"))


#---- edit these
df = your_dataframe # should be melt_df
row_variable = the-name-of-the-column-you-want-as-rows
column_variable = the-name-of-the-column-you-want-as-columns
size_variable = the-name-of-the-column-you-want-the-size-of-the-dot-to-correspond-to
color_variable = the-name-of-the-column-you-want-the-color-intensity-of-the-dot-to-correspond-to #can be the same as size if you want

#---- run this!
ggplot(df, aes(row_variable, column_variable)) +
  geom_point(aes(size = size_variable, colour= color_variable)) + 
  scale_colour_gradient(low = "#EBEABC", high = #585CE4) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1, fill = NA)) +
  scale_y_discrete(limits = rev(levels(df$row_variable))) 
