#Read in data file
mcmc <- read.table("proposal_bpp_mcmc.txt", header=TRUE, sep= "\t")
#Subsample 1,000 samples from the posterior
df.new = mcmc[seq(10000, 99999, 90), ]
#Write subsampled data to file
write.table(df.new, "subsampled_divtimes_1t.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#Read in data
df.new <- read.table("subsampled_divtimes_1t.txt", header = TRUE, sep = "\t")

#Calculate gdi's se editara con cada especie que se probara tener el archivo mcmc como referencia para sustituir en las thtetas y tau
gdi_S_occipitomaculata <- numeric()
for (i in 1:length(df.new$theta_1S_occipitomaculata)) {
  gdi_S_occipitomaculata <- c(gdi_S_occipitomaculata,1 - exp(((-2*df.new$tau_6S_occipitomaculataS_dekayi_MEXS_dekayi_VCCS_dekayi_USAOesteS_dekayi_USAEste[i])/df.new$theta_1S_occipitomaculata[i])))
}
df.new <- cbind(df.new,gdi_S_occipitomaculata)
gdi_SdekayiUSAeste <- numeric()
for (i in 1:length(df.new$theta_2S_dekayi_USAEste)) {
  gdi_SdekayiUSAeste <- c(gdi_SdekayiUSAeste,1 - exp(((-2*df.new$tau_7S_dekayi_MEXS_dekayi_VCCS_dekayi_USAOesteS_dekayi_USAEste[i])/df.new$theta_2S_dekayi_USAEste[i])))
}
df.new <- cbind(df.new,gdi_SdekayiUSAeste)
gdi_SdekayiMex <- numeric()
for (i in 1:length(df.new$theta_3S_dekayi_MEX)) {
  gdi_SdekayiMex <- c(gdi_SdekayiMex,1 - exp(((-2*df.new$tau_8S_dekayi_MEXS_dekayi_VCCS_dekayi_USAOeste[i])/df.new$theta_3S_dekayi_MEX[i])))
}
df.new <- cbind(df.new,gdi_SdekayiMex)
gdi_SdekayiVCC <- numeric()
for (i in 1:length(df.new$theta_4S_dekayi_VCC)) {
  gdi_SdekayiVCC <- c(gdi_SdekayiVCC,1 - exp(((-2*df.new$tau_9S_dekayi_VCCS_dekayi_USAOeste[i])/df.new$theta_4S_dekayi_VCC[i])))
}
df.new <- cbind(df.new,gdi_SdekayiVCC)
gdi_SdekayiUSAoeste <- numeric()
for (i in 1:length(df.new$theta_5S_dekayi_USAOeste)) {
  gdi_SdekayiUSAoeste <- c(gdi_SdekayiUSAoeste,1 - exp(((-2*df.new$tau_9S_dekayi_VCCS_dekayi_USAOeste[i])/df.new$theta_5S_dekayi_USAOeste[i])))
}
df.new <- cbind(df.new,gdi_SdekayiUSAoeste)

#Calculate 95% confidence intervals
divtime_table = data.frame(matrix(ncol = 2, nrow = 0))
names(divtime_table) <- c("gdi","species")
for (i in df.new$gdi_S_occipitomaculata) {
  if( (i > quantile(df.new$gdi_S_occipitomaculata,0.025)) & (i < quantile(df.new$gdi_S_occipitomaculata,0.975)) )
    divtime_table[nrow(divtime_table)+1,] <- c(i,"S_occipitomaculata")
}
for (i in df.new$gdi_SdekayiUSAeste) {
  if( (i > quantile(df.new$gdi_SdekayiUSAeste,0.025)) & (i < quantile(df.new$gdi_SdekayiUSAeste,0.975)) )
    divtime_table[nrow(divtime_table)+1,] <- c(i,"S_dekayi_USAEste")
}
for (i in df.new$gdi_SdekayiMex) {
  if( (i > quantile(df.new$gdi_SdekayiMex,0.025)) & (i < quantile(df.new$gdi_SdekayiMex,0.975)) )
    divtime_table[nrow(divtime_table)+1,] <- c(i,"S_dekayi_Mex")
}
for (i in df.new$gdi_SdekayiVCC) {
  if( (i > quantile(df.new$gdi_SdekayiVCC,0.025)) & (i < quantile(df.new$gdi_SdekayiVCC,0.975)) )
    divtime_table[nrow(divtime_table)+1,] <- c(i,"S_sp_VCC")
}
for (i in df.new$gdi_SdekayiUSAoeste) {
  if( (i > quantile(df.new$gdi_SdekayiUSAoeste,0.025)) & (i < quantile(df.new$gdi_SdekayiUSAoeste,0.975)) )
    divtime_table[nrow(divtime_table)+1,] <- c(i,"S_dekayi_USAOeste")
}

#Write gdi data to file
write.table(divtime_table, "gdi_1t_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Cargar datos
library(ggplot2)
gdi_table <- read.table("gdi_1t_all.txt", header = TRUE, sep = "\t")

# Crear gráfico solo con cajas
p <- ggplot(gdi_table, aes(x = species, y = as.numeric(gdi), fill = species)) +
  scale_y_continuous(name = "GDI", limits = c(0, 1)) +
  geom_hline(yintercept = 0.2, linetype = 'dashed', col = 'grey') +
  geom_hline(yintercept = 0.7, linetype = 'dashed', col = 'grey') +
  geom_boxplot(width = 0.1) +  # Solo cajas, ancho ajustado
  scale_x_discrete(
    name = "Species",
    limits = c("S_occipitomaculata", "S_dekayi_USAEste", "S_dekayi_Mex", "S_sp_VCC", "S_dekayi_USAOeste")
  ) + scale_fill_manual(values = c("S_occipitomaculata" = "grey", "S_dekayi_USAEste" = "purple", "S_dekayi_Mex" = "green", "S_sp_VCC" = "blue", "S_dekayi_USAOeste" = "red")) +
  theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 14), axis.title = element_text(size = 16))

# Mostrar el gráfico
print(p)

# Guardar como SVG
ggsave("gdi_plot.svg", plot = p, width = 20, height = 6, dpi = 300, device = "svg")

# Calcular estadísticos
mean(gdi_SdekayiUSAeste)  
quantile(gdi_SdekayiUSAeste, c(0.025, 0.975))  # Intervalo de confianza al 95%

mean(gdi_SdekayiMex)  
quantile(gdi_SdekayiMex, c(0.025, 0.975))  # Intervalo de confianza al 95%

mean(gdi_SdekayiVCC)  
quantile(gdi_SdekayiVCC, c(0.025, 0.975))  

mean(gdi_SdekayiUSAoeste)  
quantile(gdi_SdekayiUSAoeste, c(0.025, 0.975))  