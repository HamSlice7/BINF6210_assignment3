#1.load the necessary libraries----
library("tidyverse")
library("vegan")
library("fossil")
library("ggmap") 
library("ggplot2")
library("readr")

#2. Load in Aves tsv----

#read in the aves data frame from BOLD
#Changing qoutations so data can be loaded. 
df_a <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Aves&format=tsv")
#Writing Aves_BOLD_data so we dont have to download data everytime we open this project.
write_tsv(df_a, "Aves_BOLD_data.tsv")

#3.Number of samples based on latitude----

#The function serves as a tool for exploratory analysis of my data. It selects sampleid and lat columns from df_a and filters the desired latitudinal range. It then groups the data by latitude and counts the number of samples at each latitude group. It then pulls the count column and then calculates the sum of this column to get the number of samples made of the taxa at the inputed latitudinal band
num_of_samples_aves <- function(minimum_latitude, max_latitude) {
  sample_count <- df_a %>%
    select(sampleid, lat) %>%
    filter(!is.na(lat) & lat >= minimum_latitude & lat < max_latitude) %>%
    group_by(lat) %>%
    summarise(count = n()) %>%
    pull(count) %>%
    sum() %>%
    return(sample_count) }

aves_sampling_effort <- c(num_of_samples_aves(0, 15), num_of_samples_aves(16, 30), num_of_samples_aves(31, 45), num_of_samples_aves(46, 60), num_of_samples_aves(61, 75), num_of_samples_aves(76, 90))

latitude <- c('0-15', '16-30', '31-45', '46-60', '61-75', '76-90')

aves_sampling <- data.frame(latitude, aves_sampling_effort)
#4. Creating a map of sampling efforts for Aves----


#Sample data frame with latitude and longitude 
df_a_map <- df_a %>%
  select(sampleid, lat, lon) %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  filter(lat >= 0)

register_google(key = "AIzaSyBxWwRVQigeDNM_MhkZlR2fhXHN8e8Vqbk")

#create map plot of sampling efforts for Aves
map_aves <- ggmap(get_googlemap(center = c(lon = 0, lat = 40), zoom = 1, maptype = 'terrain', color = 'color')) +
  geom_point(data = df_a_map, aes(x=lon, y=lat), color = 'blue', size = 0.1) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  ) 

print(map_aves)

# Calculating summary statistics for latitude in df_a_map
summary_stats <- df_a_map %>%
  summarise(
    Mean_Latitude = mean(lat, na.rm = TRUE),
    Median_Latitude = median(lat, na.rm = TRUE),
    Min_Latitude = min(lat, na.rm = TRUE),
    Max_Latitude = max(lat, na.rm = TRUE)
  )

# Printing summary statistics
print(summary_stats)

#Introducing spatial distribution plot using the sf package and ggplot2 to offer a more detailed and insightful representation of sampling efforts across latitudinal bands. 
library(sf)
# Convert df_a_map to sf object for spatial plotting
df_a_sf <- st_as_sf(df_a_map, coords = c("lon", "lat"), crs = 4326)

# Create a spatial distribution plot
sampling_map <- ggplot() +
  geom_sf(data = df_a_sf, color = "blue", size = 0.1) +
  theme_minimal() +
  labs(
    title = 'Spatial Distribution of Aves Sampling Efforts',
    x = 'Longitude',
    y = 'Latitude'
  )
print(sampling_map)



#5.Function that estimates species richness using vegan and taxonomic identification   ----

#Here, I create a function that will select species_name and lat and filter out any na's. It will then filter the desired latitudinal range of the user. Next, it will then group by species_name and count the number of species_name groups. It will then reformat the data so that it can be used for the specnumber() function in the vegan package. The function will then return the estimated species richness.
speciesrichness_name_f <- function(min_lat, max_lat) {
  c_a <- df_a %>%
    select(species_name, lat) %>%
    filter(!is.na(species_name), !is.na(lat)) %>%
    filter(lat >= min_lat & lat < max_lat) %>%
    group_by(species_name) %>%
    count(species_name) %>%
    pivot_wider(names_from = species_name, values_from = n) %>%
    vegan::specnumber() %>%
    return(c_a) 
}
#calling the function with my latitudinal bands as input
name_vegan_015 <- speciesrichness_name_f(0, 15)
name_vegan_1630 <- speciesrichness_name_f(16, 30)
name_vegan_3145 <- speciesrichness_name_f(31, 45)
name_vegan_4660 <- speciesrichness_name_f(46, 60)
name_vegan_6175 <- speciesrichness_name_f(61, 75)
name_vegan_7690 <- speciesrichness_name_f(76, 90)



#6.Function that estimates species richness using chao1 and bin identification----

#Here, I create a function that will select bin_uri and lat and filter out any NA's. It will then filter the desired latitudinal range of the user. Next, it will then group by bin_uri and count the number of bin_uri groups. It will then reformat the data so that it can be used for the chao1 function in the fossil package. The function will then return the estimated species richness.
chao1f <- function(min_lat, max_lat) {
  c_a <- df_a %>%
    select(bin_uri, lat) %>%
    filter(!is.na(bin_uri), !is.na(lat)) %>%
    filter(lat >= min_lat & lat < max_lat) %>%
    group_by(bin_uri) %>%
    count(bin_uri) %>%
    pivot_wider(names_from = bin_uri, values_from = n) %>%
    fossil::chao1(taxa.row = TRUE) %>%
    return(c_a) 
  }

#calling the function with my latitudinal bands as input
chao1f_015 <- chao1f(0, 15)
chao1f_1630 <- chao1f(16, 30)
chao1f_3145 <- chao1f(31, 45)
chao1f_4660 <- chao1f(46, 60)
chao1f_6175 <- chao1f(61, 75)
chao1f_7690 <- chao1f(76, 90)




    


#7. Function that estimates species richness using vegan and bin identification----

#Here, I create a function that will select bin_uri and lat and filter out any NA's. It will then filter the desired latitudinal range of the user. Next, it will then group by bin_uri and count the number of bin_uri groups. It will then reformat the data so that it can be used for the specnumber() function in the vegan package. The function will then return the estimated species richness.
bin_speciesrichness_f <- function(min_lat, max_lat) {
  c_a <- df_a %>%
    select(bin_uri, lat) %>%
    filter(!is.na(bin_uri), !is.na(lat)) %>%
    filter(lat >= min_lat & lat < max_lat) %>%
    group_by(bin_uri) %>%
    count(bin_uri) %>%
    pivot_wider(names_from = bin_uri, values_from = n) %>%
    vegan::specnumber() %>%
    return(c_a) 
}
#calling the function with latitudinal bands as inputs
bin_vegan_015 <- bin_speciesrichness_f(0, 15)
bin_vegan_1630 <- bin_speciesrichness_f(16, 30)
bin_vegan_3145 <- bin_speciesrichness_f(31, 45)
bin_vegan_4660 <- bin_speciesrichness_f(46, 60)
bin_vegan_6175 <- bin_speciesrichness_f(61, 75)
bin_vegan_7690 <- bin_speciesrichness_f(76, 90)



#8. Function that estimates species richness using Choa1 and taxanomic identification----

#Here, I create a function that will select species_name and lat and filter out any na's. It will then filter the desired latitudinal range of the user. Next, it will then group by species_name and count the number of species_name groups. It will then reformat the data so that it can be used for the chao1 function in the fossil package. The function will then return the estimated species richness.
chao1_name_f <- function(min_lat, max_lat) {
  c_a <- df_a %>%
    select(species_name, lat) %>%
    filter(!is.na(species_name), !is.na(lat)) %>%
    filter(lat >= min_lat & lat < max_lat) %>%
    group_by(species_name) %>%
    count(species_name) %>%
    pivot_wider(names_from = species_name, values_from = n) %>%
    fossil::chao1(taxa.row = TRUE) %>%
    return(c_a) 
}

#calling function with latitudinal bands as inputs
chao1_namef_015 <- chao1_name_f(0, 15)
chao1_namef_1630 <- chao1_name_f(16, 30)
chao1_namef_3145 <- chao1_name_f(31, 45)
chao1_namef_4660 <- chao1_name_f(46, 60)
chao1_namef_6175 <- chao1_name_f(61, 75)
chao1_namef_7690 <- chao1_name_f(76, 90)







#9.creating a bar plot to visualize the results across the calculations.----
#creating a data frame to organize my results across the calculations
calculation_groups <- c("Chao1  TaxaID", "Vegan TaxaID", "Chao1 BinID", "Vegan BinID")

species_richness_aves <- c(chao1_namef_015, name_vegan_015,chao1f_015,bin_vegan_015, chao1_namef_1630,name_vegan_1630,chao1f_1630, bin_vegan_1630, chao1_namef_3145,name_vegan_3145,chao1f_3145, bin_vegan_3145, chao1_namef_4660, name_vegan_4660, chao1f_4660,bin_vegan_4660, chao1_namef_6175,name_vegan_6175,bin_vegan_6175, chao1f_6175, chao1_namef_7690,name_vegan_7690, chao1f_7690, bin_vegan_7690)
#creating the data frame
df_aves_species_richness_calculations <- data.frame(latitude, calculation_groups, species_richness_aves)
#creating a bar plot to visualize the results across the calculations. I am using ggplot to create a grouped bar char of my data so that I can visualize each species richness estimate at each of the latitudinal bands. I am ensuring the bar plot colours are colour blind friendly by using 'Dark2' palette.

#Correcting Variable spelling mistakes and Improving bar plot aesthetics for better visualization
aves_plot_new <- df_aves_species_richness_calculations %>%
  ggplot(aes(x = latitude, y = species_richness_aves, fill = calculation_groups)) +
  geom_bar(stat = 'identity', position = 'dodge', width = 0.7, color = 'black') +
  labs(
    title = 'Species Richness Estimates of Aves Across Latitudes',
    x = 'Latitudes',
    y = 'Species Richness Estimates',
    fill = 'Calculations'
  ) +
  
  scale_fill_brewer(palette = 'Set1') +  # Using a different color palette
  theme_minimal() +
  theme(
    legend.position = 'top',
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  coord_flip()

print(aves_plot_new)

