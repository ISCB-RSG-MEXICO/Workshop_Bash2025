#------------------------------------
#Titulo:
#Fecha:
#



#cargar paqueterias 
# Instalar paquetes si es la primera vez 
install.packages("tidyverse") 


#Cargar de librerías necesarias:

library(readr) # Lectura de archivos (TSV, CSV, texto) 
library(dplyr) # Manipulación de datos (pipes, mutate, filter…) 
library(ggplot2) # Gráficos (histogramas, dispersión, boxplots)


#----Lectura de archivos generados por Bedtools----

#Archivo de coverage
#Leer archivo obtenido con bedtools coverage

coverage <- read.delim(file = "Archivos_Bash/genes_coverage.tsv", 
                       col.names = c( "chrom", "start", "end", "gene_id", "number_of_peaks", "bases_overlapped", "gene_length", "fraction_covered" ) )

#Vista general del objeto cargado
glimpse(coverage)
head(coverage)

#Archivo de intersect 

#Leer archivo de intersecciones (solo genes con al menos un pico)
intersect_genes <- read.delim("Archivos_Bash/genes_peaks_intersect.bed", 
                              col.names = c("chrom", "start", "end", "gene_id"))

#vista general del objeto cargado
glimpse(intersect_genes) 
head(intersect_genes)

##----Chequeos de consistencia----
# Verificar que la fracción cubierta sea coherente
coverage_check <- coverage %>%
  mutate(
    fraction_calc = bases_overlapped / gene_length,
    diff_fraction = abs(fraction_calc - fraction_covered)
  )

summary(coverage_check$diff_fraction)

# Si todo está bien, diff_fraction debería ser ~0 (pequeños errores de redondeo)

# Ver cuántos genes tienen al menos un pico según coverage
coverage %>%
  summarise(
    total_genes = n(),
    genes_con_picos = sum(number_of_peaks > 0),
    genes_sin_picos = sum(number_of_peaks == 0)
  )

#Ver cuántos genes hay en la lista de intersect
intersect_genes %>%
  summarise(
    genes_intersect_unicos = n_distinct(gene_id),
    filas_intersect = n()
  )



######################
#Cruzar coverage con intersect para validar
######################

# 3.1 Agregar una columna lógica: ¿el gen aparece en intersect?
coverage_validado <- coverage %>%
  mutate(
    en_intersect = gene_id %in% intersect_genes$gene_id
  )

# 3.2 Resumen de consistencia
coverage_validado %>%
  count(number_of_peaks > 0, en_intersect) %>%
  rename(tiene_picos_según_coverage = `number_of_peaks > 0`)

# Idealmente:
# - Genes con number_of_peaks > 0 deberían tener en_intersect == TRUE
# Si no es así, puede indicar diferencias en parámetros de bedtools o en archivos.

##########################
# Exploración básica de cobertura (estilo tutorial)
#########################

# Distribución de fracción cubierta por gen (plot todo churro)
ggplot(coverage, aes(x = fraction_covered)) +
  geom_histogram(bins = 50) +
  labs(
    title = "Distribución de la fracción de cobertura por gen",
    x = "Fracción del gen cubierta por picos",
    y = "Número de genes"
  )

# Comparar genes con vs sin picos 
coverage %>%
  mutate(
    tiene_picos = number_of_peaks > 0
  ) %>%
  group_by(tiene_picos) %>%
  summarise(
    n_genes = n(),
    mediana_fraccion = median(fraction_covered),
    promedio_fraccion = mean(fraction_covered),
    max_fraccion = max(fraction_covered)
  )


#boxplot final

coverage %>% 
  mutate(tiene_picos = ifelse(number_of_peaks > 0, "Con picos", "Sin picos")) %>% 
  ggplot(aes(x = tiene_picos, y = fraction_covered, fill = tiene_picos)) +
  geom_boxplot(color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("Con picos" = "#7FC97F",  # verde pastel
                               "Sin picos" = "#BEAED4")) + # morado pastel
  labs(title = "Cobertura de picos por gen",
       x = "Grupo de genes",
       y = "Fracción cubierta") +
  theme_minimal(base_size = 14)
















